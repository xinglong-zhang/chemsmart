"""Cross-session accumulation of agent turns as SFT-ready episodes.

Every completed unified-loop turn is appended as one JSONL line to a
month-partitioned file under the repo-local ``var/agent-training/`` store.
Unlike the
per-session decision log, an episode carries the *model input* too (the
system prompt, deduplicated by hash into ``training/prompts/``) plus outcome
labels (semantic-gate verdict, execution return code, approvals/denials), so
working turns can be exported to fine-tuning data without replaying the
session. Absolute paths are masked before writing; API keys are never part
of the captured payloads.

Configured by an optional top-level ``training_log:`` block in
``~/.chemsmart/agent/agent.yaml``::

    training_log:
      enabled: true          # default: true
      dir: ~/some/other/dir  # default: <repo>/var/agent-training
      mask_paths: true       # default: true
"""

from __future__ import annotations

import hashlib
import json
import os
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

JsonDict = dict[str, Any]

EPISODE_SCHEMA_VERSION = 2
# Optional provenance attached to a new episode snapshot. Every field is a
# short opaque identifier used later to build leakage-safe train/eval splits
# (see scripts/training/build_split_manifest.py) and to attribute a row to the
# scenario/batch/teacher/fixture/CLI-schema that produced it. It is strictly
# additive: legacy rows have no such key, and the field is omitted when empty
# so the ledger stays append-only and backward compatible.
DATASET_PROVENANCE_ENV = "CHEMSMART_AGENT_DATASET_PROVENANCE"
_PROVENANCE_KEYS = (
    "scenario_id",
    "scenario_family",
    "batch_id",
    "teacher_id",
    "fixture_id",
    "schema_version",
    "collector_version",
)
_PROVENANCE_VALUE_MAX = 200
_PUBLIC_MESSAGE_KEYS = {
    "role",
    "content",
    "tool_calls",
    "tool_call_id",
    "name",
    "function_call",
}
_SECRET_PATTERNS = (
    re.compile(r"\bhf_[A-Za-z0-9]{20,}\b"),
    re.compile(r"\bsk-[A-Za-z0-9_-]{20,}\b"),
    re.compile(r"\bsk-proj-[A-Za-z0-9_-]{20,}\b"),
)


def _default_training_dir() -> Path:
    configured = os.environ.get("CHEMSMART_AGENT_TRAINING_DIR")
    if configured:
        return Path(configured).expanduser()
    # Repo-local by default: raw data stays near the code/export scripts while
    # remaining untracked because the repository already ignores `var/`.
    return Path(__file__).resolve().parents[2] / "var" / "agent-training"


@dataclass(frozen=True)
class TrainingLogConfig:
    enabled: bool = True
    dir: Path = field(default_factory=_default_training_dir)
    mask_paths: bool = True


def load_training_log_config(
    yaml_path: str | os.PathLike[str] | None = None,
) -> TrainingLogConfig:
    """Read the ``training_log:`` block from agent.yaml (missing = defaults)."""

    from chemsmart.agent.provider_config import _default_yaml_path

    if os.environ.get("PYTEST_CURRENT_TEST") and yaml_path is None:
        return TrainingLogConfig(enabled=False)

    path = Path(yaml_path) if yaml_path is not None else _default_yaml_path()
    block: Any = None
    try:
        import yaml

        document = yaml.safe_load(path.read_text(encoding="utf-8"))
        if isinstance(document, dict):
            block = document.get("training_log")
    except Exception:
        block = None
    if not isinstance(block, dict):
        return TrainingLogConfig()
    directory = block.get("dir")
    return TrainingLogConfig(
        enabled=bool(block.get("enabled", True)),
        dir=(
            Path(str(directory)).expanduser()
            if directory
            else _default_training_dir()
        ),
        mask_paths=bool(block.get("mask_paths", True)),
    )


class TrainingEpisodeWriter:
    """Append one SFT-ready episode per completed agent turn."""

    def __init__(self, config: TrainingLogConfig | None = None) -> None:
        self.config = config or load_training_log_config()

    @property
    def enabled(self) -> bool:
        return self.config.enabled

    def write_episode(
        self,
        *,
        session_id: str,
        turn: int,
        provider_name: str,
        model: str | None,
        messages: list[JsonDict],
        tool_records: list[JsonDict],
        tool_requests: list[Any] | None = None,
        approvals_count: int = 0,
        denials_count: int = 0,
        cwd: str | None = None,
        available_tools: list[str] | None = None,
        tool_names: list[str] | None = None,
        final_answer: str = "",
        workspace: JsonDict | None = None,
        terminal_state: JsonDict | None = None,
        paused: bool = False,
        dataset_provenance: JsonDict | str | None = None,
    ) -> Path | None:
        """Persist one turn; returns the episode file path or None if off.

        ``tool_records`` are ``{"tool", "status", "result"}`` dicts (results
        already JSON-safe). ``paused`` marks a turn snapshotted while waiting
        on ask_user; the resumed completion re-writes the same
        ``(session_id, turn)`` and supersedes it at export time.
        ``dataset_provenance`` (or, if unset, the ``CHEMSMART_AGENT_DATASET_
        PROVENANCE`` env var) attaches sanitized scenario/batch/teacher/fixture
        provenance for downstream split construction. Failures are swallowed:
        training capture must never break a live agent turn.
        """

        if not self.config.enabled:
            return None
        try:
            return self._write(
                session_id=session_id,
                turn=turn,
                provider_name=provider_name,
                model=model,
                messages=messages,
                tool_records=tool_records,
                tool_requests=tool_requests,
                approvals_count=approvals_count,
                denials_count=denials_count,
                cwd=cwd,
                available_tools=available_tools,
                tool_names=tool_names,
                final_answer=final_answer,
                workspace=workspace,
                terminal_state=terminal_state,
                paused=paused,
                dataset_provenance=dataset_provenance,
            )
        except Exception:
            return None

    def _write(
        self,
        *,
        session_id: str,
        turn: int,
        provider_name: str,
        model: str | None,
        messages: list[JsonDict],
        tool_records: list[JsonDict],
        tool_requests: list[Any] | None,
        approvals_count: int,
        denials_count: int,
        cwd: str | None,
        available_tools: list[str] | None,
        tool_names: list[str] | None,
        final_answer: str,
        workspace: JsonDict | None,
        terminal_state: JsonDict | None,
        paused: bool,
        dataset_provenance: JsonDict | str | None = None,
    ) -> Path:
        # Provider SDKs may attach hidden reasoning and transport metadata to
        # message objects. Keep only the public chat/tool-call contract in the
        # append-only ledger; hidden fields must never become SFT targets.
        conversation = [
            _public_message(message)
            for message in messages
            if isinstance(message, dict)
        ]
        system_sha = ""
        if conversation and conversation[0].get("role") == "system":
            system_message = conversation.pop(0)
            system_sha = self._store_prompt(
                str(system_message.get("content") or "")
            )

        masked_cwd = cwd or ""
        event_records = _tool_event_records(tool_records, tool_requests)
        invoked_tools = _invoked_tools(event_records)
        catalog_tools = sorted(available_tools or tool_names or [])
        workspace_info = workspace or _workspace_snapshot(cwd)
        if self.config.mask_paths:
            conversation = _sanitize_payload(conversation, cwd)
            tool_records = _sanitize_payload(tool_records, cwd)
            event_records = _sanitize_payload(event_records, cwd)
            final_answer = _sanitize_payload(final_answer, cwd)
            workspace_info = _sanitize_payload(workspace_info, cwd)
            terminal_state = _sanitize_payload(terminal_state or {}, cwd)
            masked_cwd = _sanitize_payload(masked_cwd, cwd)
        else:
            conversation = _mask_secrets(conversation)
            tool_records = _mask_secrets(tool_records)
            event_records = _mask_secrets(event_records)
            terminal_state = _mask_secrets(terminal_state or {})
        final_answer = _mask_secrets(final_answer)
        workspace_info = _mask_secrets(workspace_info)

        record = {
            "v": EPISODE_SCHEMA_VERSION,
            "ts": datetime.now(timezone.utc).isoformat(),
            "session_id": session_id,
            "turn": turn,
            "provider": {"name": provider_name, "model": model},
            # Back-compat alias: v1 used `tools` for the full tool catalog.
            # v2 makes the distinction explicit and keeps `tools` as invoked.
            "tools": invoked_tools,
            "available_tools": catalog_tools,
            "invoked_tools": invoked_tools,
            "tool_events": event_records,
            "system_prompt_sha": system_sha,
            "messages": conversation,
            "synthesis": _synthesis_summary(tool_records),
            "outcome": _outcome_labels(
                tool_records,
                approvals_count=approvals_count,
                denials_count=denials_count,
            ),
            "final_answer": final_answer,
            "workspace": workspace_info,
            "terminal_state": terminal_state or None,
            "cwd_masked": masked_cwd,
            "paused": paused,
        }

        provenance = _sanitize_dataset_provenance(
            dataset_provenance
            if dataset_provenance is not None
            else os.environ.get(DATASET_PROVENANCE_ENV)
        )
        if provenance:
            record["dataset_provenance"] = provenance

        episode_path = self._episode_path()
        episode_path.parent.mkdir(parents=True, exist_ok=True)
        with episode_path.open("a", encoding="utf-8") as handle:
            handle.write(
                json.dumps(record, sort_keys=True, default=str) + "\n"
            )
        return episode_path

    def _episode_path(self) -> Path:
        stamp = datetime.now(timezone.utc).strftime("%Y%m")
        return self.config.dir / "episodes" / f"{stamp}.jsonl"

    def _store_prompt(self, content: str) -> str:
        digest = hashlib.sha256(content.encode("utf-8")).hexdigest()[:12]
        prompt_dir = self.config.dir / "prompts"
        prompt_path = prompt_dir / f"{digest}.txt"
        if not prompt_path.exists():
            prompt_dir.mkdir(parents=True, exist_ok=True)
            prompt_path.write_text(content, encoding="utf-8")
        return digest


def tool_records_from_outcomes(
    outcomes: list[Any],
    requests: list[Any] | None = None,
) -> list[JsonDict]:
    """Reduce ToolOutcome objects to the JSON-safe episode form."""

    request_map = _request_lookup(requests or [])
    records: list[JsonDict] = []
    for outcome in outcomes:
        result = getattr(outcome, "raw_result", None)
        if result is None:
            result = getattr(outcome, "result", None)
        request = request_map.get(
            (
                str(getattr(outcome, "request_id", "") or ""),
                str(getattr(outcome, "provider_call_id", "") or ""),
            )
        )
        records.append(
            {
                "tool": getattr(outcome, "name", ""),
                "status": str(getattr(outcome, "status", "")),
                "request_id": getattr(outcome, "request_id", ""),
                "provider_call_id": getattr(outcome, "provider_call_id", ""),
                "args": (
                    _jsonable(getattr(request, "arguments", {}))
                    if request is not None
                    else {}
                ),
                "normalized_args": (
                    _jsonable(getattr(request, "arguments", {}))
                    if request is not None
                    else {}
                ),
                "result": _jsonable(result),
                "error_type": getattr(outcome, "error_type", None),
                "error_message": getattr(outcome, "error_message", None),
                "handle_id": getattr(outcome, "handle_id", None),
            }
        )
    return records


def _request_lookup(requests: list[Any]) -> dict[tuple[str, str], Any]:
    lookup: dict[tuple[str, str], Any] = {}
    for request in requests:
        key = (
            str(getattr(request, "request_id", "") or ""),
            str(getattr(request, "provider_call_id", "") or ""),
        )
        lookup[key] = request
    return lookup


def _jsonable(value: Any) -> Any:
    try:
        json.dumps(value)
        return value
    except (TypeError, ValueError):
        return repr(value)


def _public_message(message: JsonDict) -> JsonDict:
    return {
        key: message[key]
        for key in _PUBLIC_MESSAGE_KEYS
        if key in message
    }


def _sanitize_dataset_provenance(value: Any) -> JsonDict | None:
    """Whitelist, coerce, and secret-mask an optional provenance object.

    Accepts a dict or a JSON string (env var form). Unknown keys are dropped,
    every value is coerced to a stripped, secret-masked, length-capped string,
    and empty values are omitted. Returns ``None`` when nothing survives so the
    caller can leave the field off entirely (append-only/back-compat).
    """

    data: Any = value
    if isinstance(data, str):
        try:
            data = json.loads(data)
        except (TypeError, ValueError):
            return None
    if not isinstance(data, dict):
        return None
    clean: JsonDict = {}
    for key in _PROVENANCE_KEYS:
        raw = data.get(key)
        if raw is None:
            continue
        text = _mask_secrets(str(raw)).strip()
        if not text:
            continue
        clean[key] = text[:_PROVENANCE_VALUE_MAX]
    return clean or None


def _tool_event_records(
    tool_records: list[JsonDict],
    tool_requests: list[Any] | None,
) -> list[JsonDict]:
    """Build compact trajectory rows from tool calls and outcomes."""

    request_map = _request_lookup(tool_requests or [])
    events: list[JsonDict] = []
    for index, record in enumerate(tool_records, start=1):
        request = request_map.get(
            (
                str(record.get("request_id") or ""),
                str(record.get("provider_call_id") or ""),
            )
        )
        args = record.get("args")
        if not args and request is not None:
            args = getattr(request, "arguments", {})
        result = record.get("result")
        semantic = result.get("semantic") if isinstance(result, dict) else None
        if semantic is None and isinstance(result, dict):
            # Project YAML validation/critic tools use verdict-like fields but
            # not the command semantic shape.
            semantic = _project_like_verdict(result)
        events.append(
            {
                "index": index,
                "tool": record.get("tool") or "",
                "status": record.get("status") or "",
                "args": _jsonable(args or {}),
                "normalized_args": _jsonable(
                    record.get("normalized_args") or args or {}
                ),
                "semantic_verdict": _semantic_verdict(semantic),
                "failed_rule_ids": _semantic_failed_rule_ids(semantic),
                "generated_input_evidence": _generated_input_evidence(
                    semantic
                ),
                "terminal_state": (
                    result.get("terminal_state")
                    if isinstance(result, dict)
                    else None
                ),
                "result_summary": _compact_result_summary(result),
                "error_type": record.get("error_type"),
                "error_message": record.get("error_message"),
                "handle_id": record.get("handle_id"),
            }
        )
    return events


def _invoked_tools(events: list[JsonDict]) -> list[str]:
    seen: set[str] = set()
    names: list[str] = []
    for event in events:
        name = str(event.get("tool") or "")
        if not name or name in seen:
            continue
        seen.add(name)
        names.append(name)
    return names


def _project_like_verdict(result: JsonDict) -> JsonDict | None:
    verdict = result.get("verdict") or result.get("validation_verdict")
    if verdict is None and isinstance(result.get("ok"), bool):
        verdict = "ok" if result.get("ok") else "reject"
    if verdict is None:
        return None
    return {
        "verdict": verdict,
        "failed_rule_ids": result.get("failed_rule_ids") or [],
        "issues": result.get("issues") or [],
        "generated_inputs": result.get("generated_inputs") or [],
    }


def _semantic_verdict(semantic: Any) -> str | None:
    return (
        str(semantic.get("verdict"))
        if isinstance(semantic, dict) and semantic.get("verdict") is not None
        else None
    )


def _semantic_failed_rule_ids(semantic: Any) -> list[str]:
    if not isinstance(semantic, dict):
        return []
    raw = semantic.get("failed_rule_ids")
    if isinstance(raw, list):
        return [str(item) for item in raw]
    issues = semantic.get("issues")
    if not isinstance(issues, list):
        return []
    failed: list[str] = []
    for issue in issues:
        if isinstance(issue, dict) and issue.get("rule_id"):
            failed.append(str(issue["rule_id"]))
    return failed


def _generated_input_evidence(semantic: Any) -> list[JsonDict]:
    if not isinstance(semantic, dict):
        return []
    generated = semantic.get("generated_inputs")
    if not isinstance(generated, list):
        return []
    evidence: list[JsonDict] = []
    for item in generated:
        if not isinstance(item, dict):
            continue
        evidence.append(
            {
                "path": item.get("path"),
                "route": item.get("route"),
            }
        )
    return evidence


def _compact_result_summary(result: Any) -> Any:
    if not isinstance(result, dict):
        return result
    keys = (
        "ok",
        "status",
        "command",
        "action",
        "schema_variant",
        "project",
        "verdict",
        "validation_verdict",
        "returncode",
        "error",
        "missing_info",
        "failed_rule_ids",
        "repaired",
        "original_command",
        "terminal_state",
        "project_name",
        "program",
        "yaml_text",
    )
    summary: JsonDict = {key: result[key] for key in keys if key in result}
    semantic = result.get("semantic")
    if isinstance(semantic, dict):
        summary["semantic"] = {
            "verdict": semantic.get("verdict"),
            "failed_rule_ids": semantic.get("failed_rule_ids") or [],
            "generated_inputs": _generated_input_evidence(semantic),
        }
    return summary


def _synthesis_summary(tool_records: list[JsonDict]) -> JsonDict | None:
    """Extract the last successful synthesize_command payload, if any."""

    for record in reversed(tool_records):
        if record.get("tool") != "synthesize_command":
            continue
        result = record.get("result")
        if not isinstance(result, dict):
            continue
        semantic = result.get("semantic")
        return {
            "status": result.get("status"),
            "command": result.get("command"),
            "explanation": result.get("explanation") or "",
            "action": result.get("action"),
            "schema_variant": result.get("schema_variant"),
            "decision_trace": result.get("decision_trace"),
            "reasoning": result.get("reasoning") or "",
            "reasoning_provenance": (
                result.get("reasoning_provenance") or ""
            ),
            "semantic_verdict": (
                semantic.get("verdict") if isinstance(semantic, dict) else None
            ),
            "failed_rule_ids": _semantic_failed_rule_ids(semantic),
            "generated_input_evidence": _generated_input_evidence(semantic),
            "missing_info": result.get("missing_info") or [],
            "alternatives": result.get("alternatives") or [],
            **(
                {"raw_response": raw_response}
                if (raw_response := _safe_compact_response(
                    result.get("raw_response")
                ))
                is not None
                else {}
            ),
        }
    return None


def _safe_compact_response(value: Any) -> JsonDict | None:
    """Keep only schema-valid compact SPEC responses, never hidden reasoning."""

    candidate: Any = value
    if isinstance(value, str):
        try:
            candidate = json.loads(value)
        except (TypeError, ValueError):
            return None
    if not isinstance(candidate, dict):
        return None
    if any(
        key in candidate
        for key in ("reasoning_content", "thinking", "analysis", "<think>")
    ):
        return None
    intent = candidate.get("intent")
    if intent == "workflow" and isinstance(candidate.get("jobs"), list):
        return candidate
    if intent in {"advisory", "decline", "chitchat"} and isinstance(
        candidate.get("message"), str
    ):
        return candidate
    return None


def _outcome_labels(
    tool_records: list[JsonDict],
    *,
    approvals_count: int,
    denials_count: int,
) -> JsonDict:
    gate = "none"
    execute_rc: int | None = None
    for record in tool_records:
        result = record.get("result")
        if not isinstance(result, dict):
            continue
        if record.get("tool") == "synthesize_command":
            semantic = result.get("semantic")
            if isinstance(semantic, dict) and semantic.get("verdict"):
                gate = str(semantic["verdict"])
        elif record.get("tool") == "execute_chemsmart_command":
            returncode = result.get("returncode")
            if isinstance(returncode, int):
                execute_rc = returncode
    return {
        "gate": gate,
        "execute_rc": execute_rc,
        "approvals": approvals_count,
        "denials": denials_count,
        "denied": denials_count > 0,
    }


def _workspace_snapshot(cwd: str | None) -> JsonDict:
    try:
        from chemsmart.settings.workspace_project import (
            resolve_workspace_project,
        )

        status = resolve_workspace_project(cwd=cwd)
    except Exception as exc:
        return {
            "yaml_loaded": False,
            "project": "",
            "program": "",
            "path": None,
            "message": f"workspace project resolution failed: {exc}",
        }
    return {
        "yaml_loaded": bool(getattr(status, "loaded", False)),
        "project": str(getattr(status, "project", "") or ""),
        "program": str(getattr(status, "program", "") or ""),
        "path": (
            str(getattr(status, "path"))
            if getattr(status, "path", None) is not None
            else None
        ),
        "candidates": [
            str(path) for path in getattr(status, "candidates", ()) or ()
        ],
        "message": str(getattr(status, "message", "") or ""),
    }


def _sanitize_payload(value: Any, cwd: str | None) -> Any:
    return _mask_secrets(_mask_paths(value, cwd))


def _mask_paths(value: Any, cwd: str | None) -> Any:
    """Recursively replace the cwd and home dir in every string."""

    home = str(Path.home())
    replacements: list[tuple[str, str]] = []
    if cwd and cwd not in {home, "/"}:
        replacements.append((cwd.rstrip("/"), "."))
    replacements.append((home, "~"))

    def mask(item: Any) -> Any:
        if isinstance(item, str):
            for needle, replacement in replacements:
                item = item.replace(needle, replacement)
            return item
        if isinstance(item, dict):
            return {key: mask(child) for key, child in item.items()}
        if isinstance(item, list):
            return [mask(child) for child in item]
        return item

    return mask(value)


def _mask_secrets(value: Any) -> Any:
    def mask(item: Any) -> Any:
        if isinstance(item, str):
            for pattern in _SECRET_PATTERNS:
                item = pattern.sub("[REDACTED_SECRET]", item)
            return item
        if isinstance(item, dict):
            return {key: mask(child) for key, child in item.items()}
        if isinstance(item, list):
            return [mask(child) for child in item]
        return item

    return mask(value)
