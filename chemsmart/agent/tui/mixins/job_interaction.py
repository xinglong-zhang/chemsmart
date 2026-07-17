"""Job cancellation, extraction, and resume interactions."""

from __future__ import annotations

import os
import threading

from chemsmart.agent.tui.screens.jobs_panel import JobsPanelAction
from chemsmart.agent.tui.services.job_poller import (
    JobStateReader,
    cancel_job,
    extract_run_result,
)
from chemsmart.agent.tui.widgets.cells import RunResultCell
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.popups import (
    CwdMismatchChoice,
    CwdMismatchOverlay,
    TextPromptOverlay,
)
from chemsmart.agent.tui.widgets.transcript import Transcript


class JobInteractionMixin:
    """Coordinate job actions without owning polling or rendering."""

    def _handle_jobs_panel_action(
        self, action: JobsPanelAction | None
    ) -> None:
        if action is None or not action.job_id:
            return
        if action.action == "cancel":
            self._confirm_cancel(action.job_id)
        elif action.action == "extract":
            self._extract_job_result(action.job_id)
        elif action.action == "resume":
            session_id = self._job_snapshot.get(action.job_id, {}).get(
                "session_id"
            )
            if session_id:
                self._resume_or_prompt(str(session_id))

    def _confirm_cancel(self, job_id: str) -> None:
        self.app.push_screen(
            TextPromptOverlay(
                title="Cancel job",
                prompt=f"Type yes to cancel `{job_id}`.",
            ),
            lambda value: self._handle_cancel_confirmation(job_id, value),
        )

    def _handle_cancel_confirmation(
        self, job_id: str, value: str | None
    ) -> None:
        if (value or "").strip().lower() not in {"y", "yes"}:
            return
        job = self._job_snapshot.get(job_id)
        if not job:
            self.post_error("Cancel failed", f"Unknown job id: {job_id}")
            return

        def run_cancel() -> None:
            try:
                result = cancel_job(job)
            except Exception as exc:
                self.app.call_from_thread(
                    self.post_error,
                    "Cancel failed",
                    str(exc),
                )
                return

            def publish() -> None:
                self._cancelled_job_ids.add(job_id)
                self._emit_job_update(job_id, {"status": "cancelled"})
                command = " ".join(result["command"])
                self.post_agent_message(
                    f"Cancelled `{job_id}` with `{command}`.",
                    title="Jobs",
                )

            self.app.call_from_thread(publish)

        threading.Thread(target=run_cancel, daemon=True).start()

    def _extract_job_result(self, value: str) -> None:
        try:
            result = extract_run_result(value, self._job_snapshot)
        except Exception as exc:
            self.post_error("Extract failed", str(exc))
            return
        self.query_one(Transcript).add_cell(RunResultCell(result))
        self._rendered_run_results.add(str(value))

    def _resume_or_prompt(self, session_id: str) -> None:
        try:
            state = JobStateReader.load(self.session_root / session_id)
        except RuntimeError as exc:
            self.post_error("Resume failed", str(exc))
            return
        if state is None:
            self.post_error(
                "Resume failed", f"Unknown session id: {session_id}"
            )
            return
        current_cwd = os.path.abspath(os.getcwd())
        recorded_cwd = os.path.abspath(state.cwd)
        if recorded_cwd != current_cwd:
            if self.app.plain:
                self.post_error(
                    "Resume blocked",
                    "Session cwd differs from the current cwd.",
                    {"recorded_cwd": recorded_cwd, "current_cwd": current_cwd},
                )
                return
            self.app.push_screen(
                CwdMismatchOverlay(
                    recorded_cwd=recorded_cwd,
                    current_cwd=current_cwd,
                ),
                lambda result: self._handle_cwd_choice(
                    session_id, recorded_cwd, result
                ),
            )
            return
        self.start_resume(session_id)

    def _handle_cwd_choice(
        self,
        session_id: str,
        recorded_cwd: str,
        result: CwdMismatchChoice | None,
    ) -> None:
        if result is None or result.choice == "n":
            return
        if result.choice == "c":
            os.chdir(recorded_cwd)
            self.start_resume(session_id)
            return
        if result.choice == "i":
            self.start_resume(
                session_id,
                cwd_override=os.path.abspath(os.getcwd()),
            )

    def _disarm_soft_cancel(self) -> None:
        self._quit_armed = False
        self._quit_timer = None
        self.query_one(FooterWidget).set_hint(
            "Enter to submit • /help for commands"
        )
