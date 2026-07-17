from __future__ import annotations

import json
from typing import Any, Literal

from pydantic import BaseModel, Field

from chemsmart.agent.services.memory_summaries import (
    summarize_ask_user as _summarize_ask_user,
)
from chemsmart.agent.services.memory_summaries import (
    summarize_ask_user_answer as _summarize_ask_user_answer,
)
from chemsmart.agent.services.memory_summaries import (
    summarize_tool_result as _summarize_tool_result,
)
from chemsmart.agent.services.memory_summaries import (
    summarize_tool_use_result as _summarize_tool_use_result,
)

_DEFAULT_RECENT_TURN_LIMIT = 3
_DEFAULT_TOKEN_BUDGET = 900
_MAX_REQUEST_CHARS = 240
_MAX_RATIONALE_CHARS = 220
_MAX_RESULT_CHARS = 220


class ConversationTurn(BaseModel):
    turn_index: int
    request: str
    plan_rationale: str = ""
    assistant_message: str = ""
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
    last_command: str | None = None
    last_server: str | None = None
    last_scheduler: Literal["slurm", "pbs", "sge", "lsf"] | None = None
    last_job_id: str | None = None
    last_log_path: str | None = None
    last_probe_error: str | None = None
    last_run_id: str | None = None
    last_output_path: str | None = None
    last_calculation_status: str | None = None


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

            if kind == "calculation_event" and isinstance(payload, dict):
                run = payload.get("run")
                if isinstance(run, dict):
                    run_id = _string_value(run.get("run_id"))
                    output_path = _string_value(run.get("output_path"))
                    status = _string_value(run.get("status"))
                    if run_id:
                        entities.last_run_id = run_id
                    if output_path:
                        entities.last_output_path = output_path
                        entities.last_log_path = output_path
                    if status:
                        entities.last_calculation_status = status
                    if current_turn is not None and status in {
                        "completed",
                        "chemistry_failed",
                        "process_failed",
                        "cancelled",
                        "timeout",
                    }:
                        summary = (
                            f"Calculation {run_id or 'latest'} ended with "
                            f"status={status}"
                        )
                        energy = run.get("energy")
                        if isinstance(energy, (int, float)):
                            summary += f", energy={float(energy):.12f} Eh"
                        if output_path:
                            summary += f", output={output_path}"
                        current_turn.reusable_results.append(summary + ".")
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

            if kind == "ask_user":
                summary = _summarize_ask_user(payload)
                if summary:
                    current_turn.reusable_results.append(summary)
                continue

            if kind == "ask_user_answer":
                summary = _summarize_ask_user_answer(payload)
                if summary:
                    current_turn.reusable_results.append(summary)
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
    result_summary = result.get("summary")
    if isinstance(result_summary, dict):
        result = result_summary

    if tool in {"synthesize_command", "repair_command", "dry_run_input"}:
        command = _string_value(result.get("command"))
        if command:
            entities.last_command = command

    if tool in {"execute_chemsmart_command", "inspect_calculation"}:
        calculation = result.get("calculation")
        if isinstance(calculation, dict):
            run_id = _string_value(calculation.get("run_id"))
            output_path = _string_value(calculation.get("output_path"))
            status = _string_value(calculation.get("status"))
            if run_id:
                entities.last_run_id = run_id
            if output_path:
                entities.last_output_path = output_path
                entities.last_log_path = output_path
            if status:
                entities.last_calculation_status = status
        if tool == "execute_chemsmart_command":
            command = _string_value(result.get("command"))
            if command:
                entities.last_command = command

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
