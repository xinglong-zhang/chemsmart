"""Agent-worker completion, resume, and usage synchronization."""

from __future__ import annotations

from textual.worker import Worker, WorkerState

from chemsmart.agent.core import Plan
from chemsmart.agent.decision_events import session_completed
from chemsmart.agent.tui.chat_helpers import _ready_command_hint
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.transcript import Transcript


class WorkerLifecycleMixin:
    def start_resume(
        self,
        session_id: str,
        *,
        cwd_override: str | None = None,
    ) -> None:
        session_dir = self.session_root / session_id
        if not session_dir.exists():
            self.post_error(
                "Resume failed", f"Unknown session id: {session_id}"
            )
            return
        self._stop_tailer()
        self._reset_request_state(clear_transcript=True, clear_session=True)
        self.query_one(FooterWidget).set_phase(Phase.PLANNING)
        self.query_one(FooterWidget).set_hint(f"Loading session {session_id}")
        self._attach_tailer(session_dir / "decision_log.jsonl")
        if session_completed(session_dir):
            self._load_agent_session(session_id, cwd_override=cwd_override)
            self.post_agent_message(
                f"Loaded completed session `{session_id}`."
            )
            self.query_one(FooterWidget).set_hint(
                "Loaded saved session transcript"
            )
            return
        self._current_worker = self.resume_agent_session(
            session_id,
            cwd_override=cwd_override,
        )

    def on_worker_state_changed(self, event: Worker.StateChanged) -> None:
        if event.worker.group == "calculation":
            self._on_calculation_worker_state(event)
            return
        if event.worker.group == "slash-tool":
            self._on_slash_tool_worker_state(event)
            return
        if event.worker.group != "agent-session":
            return
        if event.state == WorkerState.SUCCESS:
            if not self._on_agent_session_success(event):
                return
        elif event.state == WorkerState.ERROR:
            error = event.worker.error
            self.post_error(
                error.__class__.__name__ if error else "Worker error",
                str(error) if error else "Unknown worker error",
            )
            self.query_one(FooterWidget).set_phase(Phase.ERROR)
            self.query_one(FooterWidget).set_hint("Session failed")
            self._current_worker = event.worker
        elif event.state == WorkerState.CANCELLED:
            self.post_agent_message("Run cancelled.")
            self.query_one(FooterWidget).set_phase(Phase.INTERRUPTED)
            self.query_one(FooterWidget).set_hint("Run interrupted")
            self._current_worker = event.worker
        self._refresh_job_snapshot()
        if event.state in {WorkerState.SUCCESS, WorkerState.ERROR}:
            self._maybe_start_queued_prompt()

    def _on_calculation_worker_state(self, event: Worker.StateChanged) -> None:
        if event.state != WorkerState.ERROR:
            return
        error = event.worker.error
        self.post_error(
            error.__class__.__name__ if error else "Calculation worker error",
            str(error) if error else "Unknown calculation worker error",
        )
        if not self._worker_is_busy():
            self.query_one(FooterWidget).set_phase(Phase.ERROR)
            self.query_one(FooterWidget).set_hint("Calculation worker failed")

    def _on_slash_tool_worker_state(self, event: Worker.StateChanged) -> None:
        self._current_worker = event.worker
        if event.state == WorkerState.ERROR:
            error = event.worker.error
            self.post_error(
                error.__class__.__name__ if error else "Worker error",
                str(error) if error else "Unknown worker error",
            )
            self.query_one(FooterWidget).set_phase(Phase.ERROR)
            self.query_one(FooterWidget).set_hint("Slash command failed")
        if event.state in {WorkerState.SUCCESS, WorkerState.ERROR}:
            self._maybe_start_queued_prompt()

    def _on_agent_session_success(self, event: Worker.StateChanged) -> bool:
        """Handle a finished agent-session worker; False ends the event."""

        result = event.worker.result or {}
        self._current_worker = event.worker
        self._sync_footer_usage(result)
        if event.worker.name == "agent-ask":
            self._publish_synthesis_result(result)
            self._maybe_start_queued_prompt()
            return False
        if self._tailer is not None:
            self._tailer.read_available()
        self.query_one(Transcript).collapse_tool_chain(self._active_turn_id)
        self.call_after_refresh(
            self.query_one(Transcript).collapse_tool_chain,
            self._active_turn_id,
        )
        self._session_allow_tools = set(
            result.get("session_allow_tools") or []
        )
        is_execute = event.worker.name == "agent-execute"
        if not is_execute:
            self._update_legacy_dry_run_session(result)
        self._set_session_completion_footer(result, is_execute)
        return True

    def _update_legacy_dry_run_session(
        self, result: dict[str, object]
    ) -> None:
        session_id = result.get("session_id")
        advisory = result.get("advisory_only", False)
        chitchat = result.get("is_chitchat", False)
        blocked = result.get("blocked", False)
        has_legacy_dry_run = bool(
            result.get("dry_run_result") or result.get("dry_run_results")
        )
        if (
            session_id
            and has_legacy_dry_run
            and not blocked
            and not advisory
            and not chitchat
        ):
            self._last_dry_run_session_id = str(session_id)
        else:
            self._last_dry_run_session_id = None

    def _set_session_completion_footer(
        self, result: dict[str, object], is_execute: bool
    ) -> None:
        if result.get("blocked"):
            self.query_one(FooterWidget).set_hint("Execution blocked")
        elif (
            isinstance(result.get("plan"), Plan)
            and result["plan"].is_chitchat()
        ):
            self.query_one(FooterWidget).set_phase(Phase.IDLE)
            self.query_one(FooterWidget).set_hint("Ready")
        elif is_execute:
            self._last_dry_run_session_id = None
            self.query_one(FooterWidget).set_phase(Phase.FINISHED)
            self.query_one(FooterWidget).set_hint(
                "Submitted to HPC — /jobs to track"
            )
        elif self._ready_command is not None:
            self.query_one(FooterWidget).set_phase(Phase.FINISHED)
            self.query_one(FooterWidget).set_hint(
                _ready_command_hint(self._ready_command)
            )
        elif self._last_dry_run_session_id:
            self.query_one(FooterWidget).set_phase(Phase.FINISHED)
            self.query_one(FooterWidget).set_hint(
                "Dry-run complete — /execute to submit for real"
            )
        else:
            self.query_one(FooterWidget).set_phase(Phase.FINISHED)
            self.query_one(FooterWidget).set_hint("Session complete")

    def _sync_footer_usage(self, result: dict[str, object]) -> None:
        input_tokens = result.get("total_input_tokens")
        output_tokens = result.get("total_output_tokens")
        if not isinstance(input_tokens, int) and not isinstance(
            output_tokens, int
        ):
            return
        self.query_one(FooterWidget).set_usage(
            input_tokens=(
                input_tokens if isinstance(input_tokens, int) else None
            ),
            output_tokens=(
                output_tokens if isinstance(output_tokens, int) else None
            ),
        )
