"""Execute planned tool steps and persist their stable receipts."""

from __future__ import annotations

import time
from typing import Any, Protocol

from chemsmart.agent.harness.workflow_state import workflow_state_scope
from chemsmart.agent.models import Step, utc_now_iso
from chemsmart.agent.services.result_codec import preview_value, resolve_refs
from chemsmart.agent.services.runtime_metrics import elapsed_ms


class StepExecutionHost(Protocol):
    state: Any | None
    decision_log: Any | None
    handle_store: Any | None
    registry: Any

    def _write_result_artifact(self, step_index: int, result: Any) -> Any: ...
    def _store_result_handle(
        self, tool_name: str, result: Any
    ) -> str | None: ...


class StepExecutor:
    """Run registry calls while preserving decision-log and artifact contracts."""

    def __init__(self, session: StepExecutionHost) -> None:
        self.session = session

    def execute(
        self,
        step_index: int,
        step: Step,
        prior_results: list[Any],
        extra_kwargs: dict[str, Any] | None = None,
    ) -> Any:
        session = self.session
        assert session.decision_log is not None
        resolved_args = resolve_refs(
            step.args,
            prior_results,
            handle_store=session.handle_store,
        )
        if extra_kwargs:
            resolved_args.update(extra_kwargs)
        started_at = utc_now_iso()
        started = time.perf_counter()
        self._log_call(step_index, step, resolved_args, started_at)
        result = self._call(step.tool, resolved_args)
        self._raise_for_failure(step_index, step, result, started_at, started)
        artifact = session._write_result_artifact(step_index, result)
        handle_id = session._store_result_handle(step.tool, result)
        session.decision_log.write(
            "tool_result",
            {
                "step_index": step_index,
                "tool": step.tool,
                "artifact": artifact.name,
                "handle_id": handle_id,
                "payload": preview_value(result),
                "ts_end": utc_now_iso(),
                "step_wall_time_ms": elapsed_ms(started),
            },
            rationale=step.rationale,
        )
        return result

    def preview_submit(
        self,
        step_index: int,
        step: Step,
        prior_results: list[Any],
        dry_submit: bool,
    ) -> Any:
        session = self.session
        assert session.decision_log is not None
        args = resolve_refs(
            step.args,
            prior_results,
            handle_store=session.handle_store,
        )
        args["execute"] = False
        session.decision_log.write(
            "tool_preview",
            {
                "step_index": step_index,
                "tool": step.tool,
                "args": preview_value(args),
            },
            rationale=step.rationale,
        )
        result = self._call(step.tool, args)
        result = _normalize_submit_preview(result, dry_submit=dry_submit)
        session.decision_log.write(
            "tool_preview_result",
            {
                "step_index": step_index,
                "tool": step.tool,
                "payload": preview_value(result),
            },
            rationale=step.rationale,
        )
        return result

    def record_skipped(self, step_index: int, step: Step, reason: str) -> None:
        assert self.session.decision_log is not None
        self.session.decision_log.write(
            "tool_skipped",
            {
                "step_index": step_index,
                "tool": step.tool,
                "reason": reason,
            },
            rationale=step.rationale,
        )

    def _call(self, tool_name: str, args: dict[str, Any]) -> Any:
        state = self.session.state
        scope = state.session_id if state is not None else "default"
        with workflow_state_scope(scope):
            return self.session.registry.call(tool_name, args)

    def _log_call(
        self,
        step_index: int,
        step: Step,
        args: dict[str, Any],
        started_at: str,
    ) -> None:
        assert self.session.decision_log is not None
        self.session.decision_log.write(
            "tool_call",
            {
                "step_index": step_index,
                "tool": step.tool,
                "args": preview_value(args),
                "ts_start": started_at,
                "step_wall_time_ms": None,
            },
            rationale=step.rationale,
        )

    def _raise_for_failure(
        self,
        step_index: int,
        step: Step,
        result: Any,
        started_at: str,
        started: float,
    ) -> None:
        if is_tool_error(result):
            error = result["error"]
            self._log_error(
                step_index=step_index,
                step=step,
                error_type=error.get("type", "RuntimeError"),
                message=error.get("message", "Unknown tool error"),
                payload=result,
                started_at=started_at,
                started=started,
            )
            raise RuntimeError(error.get("message", "Unknown tool error"))
        if not run_local_failed(step.tool, result):
            return
        artifact = self.session._write_result_artifact(step_index, result)
        message = run_local_failure_message(result)
        self._log_error(
            step_index=step_index,
            step=step,
            error_type="RuntimeError",
            message=message,
            payload=result,
            started_at=started_at,
            started=started,
            artifact=artifact.name,
        )
        raise RuntimeError(message)

    def _log_error(
        self,
        *,
        step_index: int,
        step: Step,
        error_type: str,
        message: str,
        payload: Any,
        started_at: str,
        started: float,
        artifact: str | None = None,
    ) -> None:
        assert self.session.decision_log is not None
        receipt = {
            "step_index": step_index,
            "tool": step.tool,
            "error_type": error_type,
            "message": message,
            "payload": preview_value(payload),
            "ts_start": started_at,
            "ts_end": utc_now_iso(),
            "step_wall_time_ms": elapsed_ms(started),
        }
        if artifact is not None:
            receipt["artifact"] = artifact
        self.session.decision_log.write(
            "tool_error", receipt, rationale=step.rationale
        )


def _normalize_submit_preview(result: Any, *, dry_submit: bool) -> Any:
    if not is_tool_error(result):
        return result
    message = result["error"]["message"]
    if not (dry_submit and is_submit_server_error(message)):
        raise RuntimeError(message)
    return {
        "transport": None,
        "script_path": None,
        "script_bytes": None,
        "command_executed": None,
        "job_id": None,
        "duplicate_check": {"duplicate": False, "message": None},
        "skipped": True,
        "skip_reason": message,
    }


def is_submit_server_error(message: str) -> bool:
    return (
        "submit_hpc requires server when no configured servers are available"
        in message
        or "submit_hpc requires server when multiple configured servers are "
        in message
    )


def is_tool_error(result: Any) -> bool:
    return (
        isinstance(result, dict)
        and result.get("ok") is False
        and "error" in result
    )


def run_local_failed(tool_name: str, result: Any) -> bool:
    return (
        tool_name == "run_local"
        and isinstance(result, dict)
        and result.get("ok") is False
    )


def run_local_failure_message(result: dict[str, Any]) -> str:
    message = "run_local failed"
    returncode = result.get("returncode")
    if returncode is not None:
        message += f" with returncode {returncode}"
    stderr_path = result.get("stderr_path")
    if isinstance(stderr_path, str) and stderr_path.strip():
        message += f"; see {stderr_path}"
    return message


__all__ = [
    "StepExecutor",
    "is_submit_server_error",
    "is_tool_error",
    "run_local_failed",
    "run_local_failure_message",
]
