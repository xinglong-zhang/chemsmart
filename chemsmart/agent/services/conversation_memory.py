from __future__ import annotations

import json
import re
from collections import Counter
from typing import Any, Literal

from pydantic import BaseModel, Field

_GAUSSIAN_ROUTE_RE = re.compile(r"^\s*#.*$", re.MULTILINE)
_ORCA_ROUTE_RE = re.compile(r"^\s*!.*$", re.MULTILINE)
_MOLECULE_REPR_RE = re.compile(r"Molecule<([^,>]+)")
_DEFAULT_RECENT_TURN_LIMIT = 3
_DEFAULT_TOKEN_BUDGET = 900
_MAX_REQUEST_CHARS = 240
_MAX_RATIONALE_CHARS = 220
_MAX_RESULT_CHARS = 220


class ConversationTurn(BaseModel):
    turn_index: int
    request: str
    plan_rationale: str = ""
    intent: Literal["workflow", "advisory", "chitchat"] | None = None
    reusable_results: list[str] = Field(default_factory=list)
    status: Literal["in_progress", "completed"] = "in_progress"
    blocked: bool | None = None

    def prompt_view(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "turn": self.turn_index,
            "request": _truncate(self.request, _MAX_REQUEST_CHARS),
        }
        if self.plan_rationale:
            payload["plan_rationale"] = _truncate(
                self.plan_rationale,
                _MAX_RATIONALE_CHARS,
            )
        if self.reusable_results:
            payload["reusable_results"] = [
                _truncate(result, _MAX_RESULT_CHARS)
                for result in self.reusable_results[-3:]
            ]
        payload["status"] = self.status
        if self.blocked is not None:
            payload["blocked"] = self.blocked
        return payload

    def older_summary_line(self) -> str:
        parts = [f"Turn {self.turn_index}: {self.request}"]
        if self.reusable_results:
            parts.append(self.reusable_results[0])
        return _truncate(" — ".join(parts), _MAX_RESULT_CHARS)


class EntityMemory(BaseModel):
    last_server: str | None = None
    last_scheduler: Literal["slurm", "pbs", "sge", "lsf"] | None = None
    last_job_id: str | None = None
    last_log_path: str | None = None
    last_probe_error: str | None = None


class ConversationMemory(BaseModel):
    turns: list[ConversationTurn] = Field(default_factory=list)
    entities: EntityMemory = Field(default_factory=EntityMemory)

    @classmethod
    def from_entries(
        cls,
        entries: list[dict[str, Any]],
    ) -> "ConversationMemory":
        turns: list[ConversationTurn] = []
        entities = EntityMemory()
        current_turn: ConversationTurn | None = None
        tool_calls: dict[int, dict[str, Any]] = {}

        for entry in entries:
            kind = entry.get("kind")
            payload = entry.get("payload")
            if kind == "request":
                request = _string_value(
                    (payload or {}).get("request")
                    if isinstance(payload, dict)
                    else None
                ) or _string_value(entry.get("rationale"))
                current_turn = ConversationTurn(
                    turn_index=len(turns) + 1,
                    request=request or "",
                )
                turns.append(current_turn)
                tool_calls = {}
                continue

            if current_turn is None or not isinstance(payload, dict):
                continue

            if kind == "plan":
                current_turn.plan_rationale = (
                    _string_value(payload.get("rationale"))
                    or _string_value(entry.get("rationale"))
                    or ""
                )
                current_turn.intent = _normalize_intent(payload.get("intent"))
                continue

            if kind == "tool_call":
                step_index = _coerce_int(payload.get("step_index"))
                if step_index is not None:
                    tool_calls[step_index] = payload
                continue

            if kind == "tool_result":
                step_index = _coerce_int(payload.get("step_index"))
                req = (
                    tool_calls.get(step_index)
                    if step_index is not None
                    else None
                )
                summary = _summarize_tool_result(
                    payload,
                    req,
                )
                if summary:
                    current_turn.reusable_results.append(summary)
                _apply_mva_entity_updates(payload, req, entities)
                continue

            # run_loop path: tool_use_request carries args; tool_use_result
            # carries the handle-wrapped result. Both are keyed by step number.
            if kind == "tool_use_request":
                step = _coerce_int(payload.get("step"))
                if step is not None:
                    tool_calls[step] = payload
                continue

            if kind == "tool_use_result":
                step = _coerce_int(payload.get("step"))
                req = tool_calls.get(step) if step is not None else None
                summary = _summarize_tool_use_result(payload, req)
                if summary:
                    current_turn.reusable_results.append(summary)
                _apply_mva_entity_updates(payload, req, entities)
                continue

            if kind == "session_summary":
                current_turn.status = "completed"
                blocked = payload.get("blocked")
                if isinstance(blocked, bool):
                    current_turn.blocked = blocked

        return cls(turns=turns, entities=entities)

    def prompt_context(
        self,
        *,
        current_turn_index: int | None = None,
        recent_turn_limit: int = _DEFAULT_RECENT_TURN_LIMIT,
        token_budget: int = _DEFAULT_TOKEN_BUDGET,
    ) -> dict[str, Any]:
        entities_payload = self.entities.model_dump(exclude_none=True)
        eligible_turns = [
            turn
            for turn in self.turns
            if turn.request
            and (
                current_turn_index is None
                or turn.turn_index < current_turn_index
            )
        ]
        if not eligible_turns:
            context = {
                "recent_turns": [],
                "older_turn_summary": [],
                "approx_token_budget": token_budget,
            }
            if entities_payload:
                context["entities"] = entities_payload
            return context

        recent_turns = eligible_turns[-recent_turn_limit:]
        older_turns = eligible_turns[: -len(recent_turns)]
        context = {
            "recent_turns": [turn.prompt_view() for turn in recent_turns],
            "older_turn_summary": [
                turn.older_summary_line() for turn in older_turns
            ],
            "approx_token_budget": token_budget,
        }
        if entities_payload:
            context["entities"] = entities_payload
        _trim_context_to_budget(context, token_budget)
        return context

    def entries_for_turn(
        self,
        entries: list[dict[str, Any]],
        turn_index: int,
    ) -> list[dict[str, Any]]:
        selected: list[dict[str, Any]] = []
        current_index = 0
        for entry in entries:
            if entry.get("kind") == "request":
                current_index += 1
            if current_index == turn_index:
                selected.append(entry)
            elif current_index > turn_index and selected:
                break
        return selected


def _trim_context_to_budget(
    context: dict[str, Any],
    token_budget: int,
) -> None:
    while _estimate_tokens(context) > token_budget:
        older_summary = context.get("older_turn_summary") or []
        if older_summary:
            older_summary.pop(0)
            continue

        recent_turns = context.get("recent_turns") or []
        if len(recent_turns) > 1:
            recent_turns.pop(0)
            continue

        if not recent_turns:
            return

        recent = recent_turns[0]
        reusable_results = recent.get("reusable_results") or []
        if len(reusable_results) > 1:
            reusable_results.pop(0)
            continue

        if reusable_results:
            shortened = _truncate(reusable_results[0], 120)
            if shortened != reusable_results[0]:
                reusable_results[0] = shortened
                continue
            recent.pop("reusable_results", None)
            continue

        rationale = recent.get("plan_rationale")
        if isinstance(rationale, str) and rationale:
            shortened = _truncate(rationale, 120)
            if shortened != rationale:
                recent["plan_rationale"] = shortened
                continue
            recent.pop("plan_rationale", None)
            continue

        request = recent.get("request")
        if isinstance(request, str):
            shortened = _truncate(request, 120)
            if shortened != request:
                recent["request"] = shortened
                continue
        return


def _summarize_tool_result(
    payload: dict[str, Any],
    tool_call: dict[str, Any] | None,
) -> str | None:
    tool = _string_value(payload.get("tool"))
    if tool is None:
        return None

    result_payload = payload.get("payload")
    if not isinstance(result_payload, dict):
        return None

    args = tool_call.get("args") if isinstance(tool_call, dict) else {}
    if not isinstance(args, dict):
        args = {}

    if tool == "build_molecule":
        source = _string_value(args.get("filepath")) or _string_value(
            args.get("smiles")
        )
        formula = _molecule_formula(result_payload)
        if source and formula:
            return f"build_molecule loaded {formula} from {source}."
        if source:
            return f"build_molecule loaded the source from {source}."
        if formula:
            return f"build_molecule produced molecule {formula}."
        return "build_molecule produced a reusable molecule."

    if tool == "recommend_method":
        method = _method_summary(result_payload)
        if method:
            return f"recommend_method suggested {method}."
        return None

    if tool in {"build_gaussian_settings", "build_orca_settings"}:
        method = _settings_summary(result_payload)
        if method:
            return f"{tool} prepared {method} settings."
        return None

    if tool == "build_job":
        kind = _string_value(args.get("kind"))
        label = _string_value(result_payload.get("label")) or _string_value(
            args.get("label")
        )
        pieces = ["build_job prepared"]
        if kind:
            pieces.append(kind)
        if label:
            pieces.append(f"label={label}")
        return " ".join(pieces) + "."

    if tool == "dry_run_input":
        inputfile = _string_value(result_payload.get("inputfile"))
        route = _extract_route_line(result_payload.get("content"))
        pieces = ["dry_run_input wrote"]
        if inputfile:
            pieces.append(inputfile)
        if route:
            pieces.append(f"with route {route}")
        return " ".join(pieces) + "."

    if tool == "extract_optimized_geometry":
        formula = _molecule_formula(result_payload)
        if formula:
            return (
                "extract_optimized_geometry recovered optimized geometry "
                f"for {formula}."
            )
        return "extract_optimized_geometry recovered optimized geometry."

    if tool == "validate_runtime":
        status = _string_value(result_payload.get("ok"))
        if status and status != "ok":
            return f"validate_runtime reported {status}."
        return None

    return None


def _normalize_intent(
    value: Any,
) -> Literal["workflow", "advisory", "chitchat"] | None:
    if value in {"workflow", "advisory", "chitchat"}:
        return value
    return None


def _string_value(value: Any) -> str | None:
    if isinstance(value, str):
        stripped = value.strip()
        return stripped or None
    return None


def _coerce_int(value: Any) -> int | None:
    if isinstance(value, int):
        return value
    if isinstance(value, str) and value.isdigit():
        return int(value)
    return None


def _normalize_scheduler_literal(
    value: Any,
) -> Literal["slurm", "pbs", "sge", "lsf"] | None:
    normalized = _string_value(value)
    if normalized is None:
        return None
    lowered = normalized.lower()
    if lowered == "slurm":
        return "slurm"
    if lowered == "pbs":
        return "pbs"
    if lowered == "sge":
        return "sge"
    if lowered == "lsf":
        return "lsf"
    return None


def _summarize_entity_error(value: Any) -> str | None:
    if isinstance(value, dict):
        message = _string_value(value.get("message"))
        if message:
            return _truncate(message, 120)
        error_type = _string_value(value.get("type"))
        if error_type:
            return _truncate(error_type, 120)
        return None
    if isinstance(value, list):
        for item in value:
            summary = _summarize_entity_error(item)
            if summary:
                return summary
        return None
    summary = _string_value(value)
    if summary:
        return _truncate(summary, 120)
    return None


def _looks_like_log_path(path: str) -> bool:
    lowered = path.lower()
    return (
        lowered.endswith((".log", ".err", ".out"))
        or "stdin.o" in lowered
        or "stdin.e" in lowered
    )


def _apply_mva_entity_updates(
    payload: dict[str, Any],
    req: dict[str, Any] | None,
    entities: EntityMemory,
) -> None:
    tool = _string_value(payload.get("tool"))
    if tool is None:
        return

    args = req.get("args") if isinstance(req, dict) else {}
    if not isinstance(args, dict):
        args = {}

    result = payload.get("payload")
    if not isinstance(result, dict):
        result = {}

    server = _string_value(args.get("server"))
    if tool == "ssh_probe":
        if server:
            entities.last_server = server

        probe_name = _string_value(args.get("probe_name")) or _string_value(
            result.get("probe")
        )
        scheduler = _normalize_scheduler_literal(result.get("scheduler"))
        if probe_name and probe_name.startswith("scheduler.") and scheduler:
            entities.last_scheduler = scheduler

        if payload.get("status") != "ok":
            error_summary = _summarize_entity_error(result.get("error"))
            if error_summary is None:
                error_summary = _summarize_entity_error(
                    payload.get("reason")
                ) or _summarize_entity_error(result)
            if error_summary:
                entities.last_probe_error = error_summary
        return

    if tool == "scheduler_query":
        if server:
            entities.last_server = server

        scheduler = _normalize_scheduler_literal(result.get("scheduler"))
        if scheduler:
            entities.last_scheduler = scheduler

        job_id = _string_value(result.get("job_id")) or _string_value(
            args.get("job_id")
        )
        if job_id:
            entities.last_job_id = str(job_id)
        return

    if tool == "log_tail":
        if server:
            entities.last_server = server

        path = _string_value(args.get("path"))
        if path:
            entities.last_log_path = path

        errors = result.get("errors")
        if isinstance(errors, list) and errors:
            error_summary = _summarize_entity_error(
                errors[0].get("line")
                if isinstance(errors[0], dict)
                else errors[0]
            )
            if error_summary:
                entities.last_probe_error = error_summary
        return

    if tool == "read":
        path = _string_value(args.get("path"))
        if path and _looks_like_log_path(path):
            entities.last_log_path = path


def _truncate(text: str, limit: int) -> str:
    if len(text) <= limit:
        return text
    return text[: max(0, limit - 1)].rstrip() + "…"


def _estimate_tokens(value: Any) -> int:
    return max(1, int(round(len(json.dumps(value, sort_keys=True)) / 4)))


def _molecule_formula(payload: dict[str, Any]) -> str | None:
    symbols = payload.get("symbols")
    if not isinstance(symbols, list) or not symbols:
        return None
    counts = Counter(
        symbol
        for symbol in symbols
        if isinstance(symbol, str) and symbol.strip()
    )
    if not counts:
        return None
    formula_parts: list[str] = []
    for symbol in sorted(counts):
        count = counts[symbol]
        formula_parts.append(symbol if count == 1 else f"{symbol}{count}")
    return "".join(formula_parts)


def _method_summary(payload: dict[str, Any]) -> str | None:
    functional = _string_value(payload.get("functional"))
    basis = _string_value(payload.get("basis"))
    ab_initio = _string_value(payload.get("ab_initio"))
    solvent = _string_value(payload.get("solvent_id"))

    method = None
    if ab_initio and basis:
        method = f"{ab_initio}/{basis}"
    elif functional and basis:
        method = f"{functional}/{basis}"
    elif functional:
        method = functional

    if method and solvent:
        return f"{method} in {solvent}"
    return method


def _settings_summary(payload: dict[str, Any]) -> str | None:
    return _method_summary(payload)


def _extract_route_line(content: Any) -> str | None:
    if not isinstance(content, str):
        return None
    match = _GAUSSIAN_ROUTE_RE.search(content) or _ORCA_ROUTE_RE.search(
        content
    )
    if match is None:
        return None
    return _truncate(match.group(0).strip(), 120)


def _summarize_tool_use_result(
    payload: dict[str, Any],
    req: dict[str, Any] | None,
) -> str | None:
    """Summarize a run_loop-style tool_use_result event (display_result form)."""
    if not isinstance(payload, dict):
        return None
    tool = _string_value(payload.get("tool"))
    if tool is None:
        return None
    status = payload.get("status")
    if status not in ("ok", "partial"):
        return None

    inner = payload.get("payload")
    args = req.get("args") if isinstance(req, dict) else {}
    if not isinstance(args, dict):
        args = {}

    if tool == "build_molecule":
        source = _string_value(args.get("filepath")) or _string_value(
            args.get("smiles")
        )
        formula = None
        if isinstance(inner, dict):
            s = inner.get("summary")
            repr_str = _string_value(
                s.get("repr") if isinstance(s, dict) else None
            )
            if repr_str:
                m = _MOLECULE_REPR_RE.search(repr_str)
                if m:
                    formula = m.group(1)
        if source and formula:
            return f"build_molecule loaded {formula} from {source}."
        if source:
            return f"build_molecule loaded the source from {source}."
        if formula:
            return f"build_molecule produced molecule {formula}."
        return "build_molecule produced a reusable molecule."

    if tool == "recommend_method":
        method = _method_summary(inner) if isinstance(inner, dict) else None
        if method:
            return f"recommend_method suggested {method}."
        return None

    if tool in {"build_gaussian_settings", "build_orca_settings"}:
        method_data: dict[str, Any] = {
            "functional": args.get("functional"),
            "basis": args.get("basis"),
            "ab_initio": args.get("ab_initio"),
            "solvent_id": args.get("solvent_id"),
        }
        method = _method_summary(method_data)
        if method:
            return f"{tool} prepared {method} settings."
        return None

    if tool == "build_job":
        kind = _string_value(args.get("kind"))
        label = _string_value(args.get("label"))
        if not label and isinstance(inner, dict):
            s = inner.get("summary")
            repr_str = _string_value(
                s.get("repr") if isinstance(s, dict) else None
            )
            if repr_str:
                m = re.search(r"label=([^,>]+)", repr_str)
                if m:
                    label = m.group(1).strip()
        pieces = ["build_job prepared"]
        if kind:
            pieces.append(kind)
        if label:
            pieces.append(f"label={label}")
        return " ".join(pieces) + "."

    if tool == "dry_run_input":
        input_summary: dict[str, Any] = {}
        if isinstance(inner, dict):
            s = inner.get("summary")
            if isinstance(s, dict):
                input_summary = s
        inputfile = _string_value(input_summary.get("inputfile"))
        route = _extract_route_line(input_summary.get("content"))
        pieces = ["dry_run_input wrote"]
        if inputfile:
            pieces.append(inputfile)
        if route:
            pieces.append(f"with route {route}")
        return " ".join(pieces) + "."

    if tool == "extract_optimized_geometry":
        formula = None
        if isinstance(inner, dict):
            s = inner.get("summary")
            repr_str = _string_value(
                s.get("repr") if isinstance(s, dict) else None
            )
            if repr_str:
                m = _MOLECULE_REPR_RE.search(repr_str)
                if m:
                    formula = m.group(1)
        if formula:
            return (
                "extract_optimized_geometry recovered optimized geometry "
                f"for {formula}."
            )
        return "extract_optimized_geometry recovered optimized geometry."

    if tool == "validate_runtime":
        ok_val = None
        if isinstance(inner, dict):
            s = inner.get("summary")
            ok_val = _string_value(
                s.get("ok") if isinstance(s, dict) else inner.get("ok")
            )
        if ok_val and ok_val != "ok":
            return f"validate_runtime reported {ok_val}."
        return None

    if tool == "ssh_probe":
        inner_dict = inner if isinstance(inner, dict) else {}
        server = _string_value(args.get("server"))
        probe_name = _string_value(args.get("probe_name")) or (
            _string_value(inner_dict.get("probe"))
        )
        scheduler_or_status = None
        scheduler_or_status = _string_value(inner_dict.get("scheduler"))
        if scheduler_or_status is None:
            scheduler_or_status = _string_value(status)

        pieces = ["ssh_probe"]
        if probe_name:
            pieces.append(probe_name)
        if server:
            pieces.append(f"on {server}")
        summary_line = " ".join(pieces)
        if scheduler_or_status:
            summary_line += f" → {scheduler_or_status}"
        return summary_line

    if tool == "scheduler_query":
        server = _string_value(args.get("server"))
        inner_dict = inner if isinstance(inner, dict) else {}
        pieces = ["scheduler_query"]
        if server:
            pieces.append(f"on {server}:")
        else:
            pieces[-1] += ":"

        scheduler_details: list[str] = []
        job_id = _string_value(inner_dict.get("job_id")) or _string_value(
            args.get("job_id")
        )
        state = _string_value(inner_dict.get("state"))
        queue = _string_value(inner_dict.get("queue")) or _string_value(
            inner_dict.get("partition_or_queue")
        )
        if job_id:
            scheduler_details.append(f"job {job_id}")
        if state:
            scheduler_details.append(f"state={state}")
        if queue:
            scheduler_details.append(f"queue={queue}")
        if scheduler_details:
            return " ".join([*pieces, *scheduler_details])
        return " ".join(pieces).rstrip(":")

    if tool == "log_tail":
        inner_dict = inner if isinstance(inner, dict) else {}
        path = _string_value(args.get("path")) or (
            _string_value(inner_dict.get("path"))
        )
        lines_returned = _coerce_int(inner_dict.get("lines_returned"))
        errors = inner_dict.get("errors")
        error_count = len(errors) if isinstance(errors, list) else 0
        top_kind = None
        if isinstance(errors, list) and errors and isinstance(errors[0], dict):
            top_kind = _string_value(errors[0].get("kind"))

        log_tail_details: list[str] = []
        if lines_returned is not None:
            log_tail_details.append(f"{lines_returned}L")
        log_tail_details.append(f"{error_count} errors")
        if top_kind:
            log_tail_details[-1] += f": {top_kind}"

        summary_line = "log_tail"
        if path:
            summary_line += f" {path}"
        if log_tail_details:
            summary_line += f" ({', '.join(log_tail_details)})"
        return summary_line

    if tool == "read":
        inner_dict = inner if isinstance(inner, dict) else {}
        path = _string_value(args.get("path")) or _string_value(
            inner_dict.get("path")
        )
        start_line = _coerce_int(inner_dict.get("start_line"))
        end_line = _coerce_int(inner_dict.get("end_line"))
        total_lines = _coerce_int(inner_dict.get("total_lines"))

        summary_line = "read"
        if path:
            summary_line += f" {path}"
        if (
            start_line is not None
            and end_line is not None
            and total_lines is not None
        ):
            summary_line += f" L{start_line}-{end_line}/{total_lines}"
        return summary_line

    return None
