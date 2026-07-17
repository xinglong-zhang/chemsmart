"""Persist stable session summaries, metadata, and harness evidence."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Protocol

from chemsmart.agent.harness.models import HarnessResult
from chemsmart.agent.harness.trace import write_harness_result
from chemsmart.agent.models import CriticVerdict, SessionState, utc_now_iso
from chemsmart.agent.services.runtime_metrics import (
    elapsed_ms,
    git_sha,
    schema_hash,
)

SESSION_METADATA_KEYS = tuple(
    sorted(
        {
            "advisory_only",
            "block_reason",
            "blocked",
            "critic_confidence",
            "critic_verdict",
            "ended_at",
            "executed_steps",
            "exit_status",
            "generated_commands",
            "git_sha",
            "harness_failed_rule_ids",
            "harness_issue_count",
            "harness_verdict",
            "input_file",
            "input_files",
            "intent",
            "is_chitchat",
            "plan_steps",
            "provider_name",
            "rationale",
            "request",
            "request_intent",
            "request_started_at",
            "resolved_model",
            "runtime_v2_event_count",
            "runtime_v2_mode",
            "runtime_v2_phase",
            "runtime_v2_shadow_violations",
            "schema_hash",
            "session_id",
            "started_at",
            "tools_called",
            "total_input_tokens",
            "total_output_tokens",
            "total_steps_executed",
            "total_steps_planned",
            "wall_time_ms",
        }
    )
)


class FinalizerHost(Protocol):
    state: SessionState | None
    session_dir: Path | None
    decision_log: Any | None
    registry: Any
    runtime_v2_mode: Any
    _runtime_controller: Any | None
    _run_start_time: float | None
    _last_harness_result: HarnessResult | None

    def _tools_called(self) -> list[str]: ...
    def _metadata_provider_name(self) -> str: ...
    def _metadata_resolved_model(self) -> str: ...
    def _total_llm_tokens(self, field: str) -> int: ...
    def _refresh_conversation_history(self) -> None: ...


class SessionFinalizer:
    """Write one terminal receipt for a completed or blocked session turn."""

    def __init__(self, session: FinalizerHost, *, module_path: Path) -> None:
        self.session = session
        self.module_path = module_path

    def finalize(
        self,
        *,
        verdict: CriticVerdict | None,
        blocked: bool,
        block_reason: str | None,
        dry_run_results: list[dict[str, Any]],
        run_error: Exception | None = None,
        advisory_only: bool = False,
        is_chitchat: bool = False,
        rationale: str = "",
    ) -> None:
        session = self.session
        assert session.state is not None
        assert session.session_dir is not None
        assert session.decision_log is not None

        ended_at = utc_now_iso()
        effective_blocked, effective_reason = _effective_failure(
            blocked, block_reason, run_error
        )
        wall_time_ms = elapsed_ms(
            session._run_start_time,
            started_at=session.state.request_started_at,
            ended_at=ended_at,
        )
        summary = self._summary(
            verdict=verdict,
            blocked=effective_blocked,
            block_reason=effective_reason,
            wall_time_ms=wall_time_ms,
            run_error=run_error,
            advisory_only=advisory_only,
            is_chitchat=is_chitchat,
            rationale=rationale,
        )
        session.decision_log.write(
            "session_summary", summary, rationale=rationale or ""
        )
        session._refresh_conversation_history()
        harness_result = session._last_harness_result
        if harness_result is not None:
            write_harness_result(session.session_dir, harness_result)
        metadata = self._metadata(
            summary=summary,
            verdict=verdict,
            dry_run_results=dry_run_results,
            harness_result=harness_result,
            ended_at=ended_at,
            advisory_only=advisory_only,
            is_chitchat=is_chitchat,
            rationale=rationale,
        )
        (session.session_dir / "session_metadata.json").write_text(
            json.dumps(metadata, indent=2, sort_keys=True),
            encoding="utf-8",
        )

    def _summary(
        self,
        *,
        verdict: CriticVerdict | None,
        blocked: bool,
        block_reason: str,
        wall_time_ms: int,
        run_error: Exception | None,
        advisory_only: bool,
        is_chitchat: bool,
        rationale: str,
    ) -> dict[str, Any]:
        session = self.session
        assert session.state is not None
        return {
            "total_steps_executed": session.state.current_step_index,
            "total_steps_planned": session.state.total_steps_planned,
            "blocked": blocked,
            "block_reason": block_reason,
            "wall_time_ms": wall_time_ms,
            "tools_called": session._tools_called(),
            "critic_confidence": verdict.confidence if verdict else 0.0,
            "request_intent": session.state.request_intent,
            "provider_name": session._metadata_provider_name(),
            "resolved_model": session._metadata_resolved_model(),
            "total_input_tokens": session._total_llm_tokens("input_tokens"),
            "total_output_tokens": session._total_llm_tokens("output_tokens"),
            "exit_status": _exit_status(run_error, blocked),
            "advisory_only": advisory_only,
            "is_chitchat": is_chitchat,
            "rationale": rationale or "",
        }

    def _metadata(
        self,
        *,
        summary: dict[str, Any],
        verdict: CriticVerdict | None,
        dry_run_results: list[dict[str, Any]],
        harness_result: HarnessResult | None,
        ended_at: str,
        advisory_only: bool,
        is_chitchat: bool,
        rationale: str,
    ) -> dict[str, Any]:
        metadata = self._base_metadata(
            summary=summary,
            verdict=verdict,
            dry_run_results=dry_run_results,
            ended_at=ended_at,
            advisory_only=advisory_only,
            is_chitchat=is_chitchat,
            rationale=rationale,
        )
        metadata.update(_harness_metadata(harness_result))
        metadata.update(self._runtime_metadata())
        return metadata

    def _base_metadata(
        self,
        *,
        summary: dict[str, Any],
        verdict: CriticVerdict | None,
        dry_run_results: list[dict[str, Any]],
        ended_at: str,
        advisory_only: bool,
        is_chitchat: bool,
        rationale: str,
    ) -> dict[str, Any]:
        session = self.session
        assert session.state is not None
        input_metadata = _input_metadata(dry_run_results)
        return {
            "session_id": session.state.session_id,
            "request": session.state.request or "unknown",
            "intent": session.state.request_intent,
            "request_intent": session.state.request_intent,
            "plan_steps": session.state.total_steps_planned,
            "executed_steps": session.state.current_step_index,
            "total_steps_planned": summary["total_steps_planned"],
            "total_steps_executed": summary["total_steps_executed"],
            **input_metadata,
            "critic_verdict": verdict.verdict if verdict else "unknown",
            "critic_confidence": summary["critic_confidence"],
            "blocked": summary["blocked"],
            "block_reason": summary["block_reason"],
            "wall_time_ms": summary["wall_time_ms"],
            "started_at": session.state.started_at,
            "request_started_at": session.state.request_started_at,
            "ended_at": ended_at,
            "provider_name": summary["provider_name"],
            "resolved_model": summary["resolved_model"],
            "git_sha": git_sha(self.module_path) or "unknown",
            "schema_hash": schema_hash(session.registry.openai_tool_defs()),
            "total_input_tokens": summary["total_input_tokens"],
            "total_output_tokens": summary["total_output_tokens"],
            "tools_called": summary["tools_called"],
            "exit_status": summary["exit_status"],
            "advisory_only": advisory_only,
            "is_chitchat": is_chitchat,
            "rationale": rationale or "",
        }

    def _runtime_metadata(self) -> dict[str, Any]:
        session = self.session
        controller = session._runtime_controller
        return {
            "runtime_v2_mode": session.runtime_v2_mode.value,
            "runtime_v2_phase": (
                controller.state.phase.value if controller else "not_run"
            ),
            "runtime_v2_event_count": (
                controller.state.latest_sequence if controller else 0
            ),
            "runtime_v2_shadow_violations": (
                list(controller.state.shadow_violations) if controller else []
            ),
        }


def _effective_failure(
    blocked: bool,
    block_reason: str | None,
    run_error: Exception | None,
) -> tuple[bool, str]:
    effective_blocked = blocked or run_error is not None
    if run_error is not None and block_reason is None:
        block_reason = f"exception:{run_error.__class__.__name__}"
    return effective_blocked, block_reason or "unknown"


def _exit_status(run_error: Exception | None, blocked: bool) -> str:
    if run_error is not None:
        return "error"
    return "blocked" if blocked else "ok"


def _input_metadata(
    dry_run_results: list[dict[str, Any]],
) -> dict[str, Any]:
    primary = dry_run_results[0] if dry_run_results else None
    return {
        "input_file": (
            str(primary.get("inputfile"))
            if primary and primary.get("inputfile")
            else "unknown"
        ),
        "input_files": [
            str(result["inputfile"])
            for result in dry_run_results
            if result.get("inputfile")
        ],
        "generated_commands": [
            str(result["command"])
            for result in dry_run_results
            if result.get("command")
        ],
    }


def _harness_metadata(
    result: HarnessResult | None,
) -> dict[str, Any]:
    return {
        "harness_verdict": result.verdict if result else "not_run",
        "harness_issue_count": len(result.issues) if result else 0,
        "harness_failed_rule_ids": result.failed_rule_ids if result else [],
    }


__all__ = ["SESSION_METADATA_KEYS", "SessionFinalizer"]
