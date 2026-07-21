"""Execution approval and dispatch for validated commands."""

from __future__ import annotations

from chemsmart.agent.permissions import ApprovalDecision
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.popups import TextPromptOverlay


class CommandExecutionMixin:
    """Dispatch only commands promoted by semantic and intent gates."""

    def _handle_execute_command(self, argument: str) -> None:
        if self._current_worker and not self._current_worker.is_finished:
            self.post_error(
                "Session already running",
                "Wait for the current request to finish before submitting.",
            )
            return
        if not self._last_dry_run_session_id:
            self.post_error(
                "No dry-run to execute",
                "Run a workflow dry-run first, then use /execute to submit "
                "for real.",
            )
            return
        if argument.strip().lower() in {"yes", "y"}:
            self._start_execute(self._last_dry_run_session_id)
            return
        if self.app.plain:
            self.post_error(
                "Confirmation required",
                "Use /execute yes to confirm real HPC submission.",
            )
            return
        session_id = self._last_dry_run_session_id
        self.app.push_screen(
            TextPromptOverlay(
                title="Submit to HPC",
                prompt=(
                    f"Session: {session_id}\n\n"
                    "This will submit real jobs to the HPC scheduler.\n"
                    "Type 'yes' to confirm."
                ),
            ),
            lambda value: self._handle_execute_confirmation(session_id, value),
        )

    def _handle_execute_confirmation(
        self, session_id: str, value: str | None
    ) -> None:
        if (value or "").strip().lower() not in {"y", "yes"}:
            return
        self._start_execute(session_id)

    def _start_execute(self, session_id: str) -> None:
        footer = self.query_one(FooterWidget)
        footer.set_phase(Phase.SUBMITTING)
        footer.set_hint("Submitting to HPC…")
        session_dir = self.session_root / session_id
        if session_dir.exists() and not (
            self._tailer_path and self._tailer_path.parent == session_dir
        ):
            self._stop_tailer()
            self._attach_tailer(session_dir / "decision_log.jsonl")
        self._current_worker = self.execute_agent_session(session_id)

    def _handle_ready_command_execution(
        self,
        action: str,
        argument: str,
    ) -> None:
        slash_command = "/run" if action == "run" else "/submit"
        pending = self._pending_tool_request
        if self._pending_approval and pending is not None:
            if not self._can_approve_or_execute(action):
                self.post_error(
                    "Approval target mismatch",
                    f"{slash_command} does not approve the pending "
                    f"`{pending.name}` action.",
                )
                return
            if argument:
                self._handle_inline_approval(argument)
            else:
                self._resolve_pending_approval(ApprovalDecision.ALLOW_ONCE)
            return

        if self._worker_is_busy():
            self.post_error(
                "Session already running",
                "Wait for the current request to finish before execution.",
            )
            return
        confirmation = argument.strip().lower()
        if (
            action == "run"
            and self._has_active_local_calculation()
            and confirmation not in {"y", "yes"}
        ):
            self.post_error(
                "Local calculation already running",
                "ChemSmart keeps one local calculation active by default to avoid CPU/RAM contention. Use Ctrl+B to inspect it, or use `/run yes` to explicitly start another local calculation.",
            )
            return
        if confirmation not in {"", "y", "yes"}:
            self.post_error(
                "Invalid execution confirmation",
                f"Usage: {slash_command} [yes]",
            )
            return
        problem = self._ready_command_problem(action)
        if problem is not None:
            self.post_error("Validated command unavailable", problem)
            return

        ready = self._ready_command
        if ready is None:  # narrowed by _ready_command_problem
            return
        footer = self.query_one(FooterWidget)
        footer.set_phase(
            Phase.EXECUTING if action == "run" else Phase.SUBMITTING
        )
        footer.set_hint(
            "Executing validated local command…"
            if action == "run"
            else "Submitting validated command…"
        )
        session = self.active_agent_session
        session_dir = self._tui_session_dir
        session_id = self._tui_session_dir.name
        if (
            session is not None
            and session.session_dir is not None
            and session.state is not None
        ):
            session_dir = session.session_dir
            session_id = session.state.session_id
        self.run_calculation_request(
            ready,
            session_dir=str(session_dir),
            session_id=session_id,
            turn_id=self._active_turn_id,
        )
        self._ready_command = None
