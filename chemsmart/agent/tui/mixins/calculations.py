"""Local calculation lifecycle and calculation-monitor interactions."""

from __future__ import annotations

from pathlib import Path

from textual import work

from chemsmart.agent.core import DecisionLog
from chemsmart.agent.runtime.calculations import (
    CalculationContext,
    CalculationEvent,
    cancel_calculation,
    inspect_calculation,
    load_calculation_runs,
)
from chemsmart.agent.tools_command import execute_chemsmart_command_observed
from chemsmart.agent.tui.chat_helpers import (
    _calculation_diagnostic_summary,
    _jobs_sort_key,
)
from chemsmart.agent.tui.chat_models import ReadyCommand
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.screens.calculations import (
    CalculationMonitor,
    CalculationMonitorAction,
)
from chemsmart.agent.tui.services.job_poller import format_jobs_table
from chemsmart.agent.tui.widgets.calculation_strip import CalculationStatusStrip
from chemsmart.agent.tui.widgets.cells import CalculationReceiptCell
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.popups import TextPromptOverlay
from chemsmart.agent.tui.widgets.transcript import Transcript


class CalculationMixin:
    @work(
        thread=True,
        exclusive=False,
        exit_on_error=False,
        group="calculation",
        name="calculation",
    )
    def run_calculation_request(
        self,
        ready: ReadyCommand,
        *,
        session_dir: str,
        session_id: str,
        turn_id: str,
    ) -> dict[str, object]:
        tool_name = "execute_chemsmart_command"
        description = "Execute validated chemsmart command"
        arguments = {"command": ready.command, "test": False}
        self.app.call_from_thread(
            self._publish_tool_call_cell,
            tool_name,
            "approved",
            description,
            arguments,
            "Approved explicitly by /run or /submit.",
        )

        def event_sink(event: CalculationEvent) -> None:
            self.app.call_from_thread(
                self._on_calculation_event,
                event.to_dict(),
            )

        result = execute_chemsmart_command_observed(
            ready.command,
            test=False,
            calculation_context=CalculationContext(
                session_dir=Path(session_dir),
                session_id=session_id,
                turn_id=turn_id,
                semantic_verdict=ready.semantic_verdict,
                intent_verdict=ready.intent_verdict,
            ),
            event_sink=event_sink,
        )
        ok = bool(result.get("ok"))
        note = (
            self._tool_success_note(tool_name, result)
            if ok
            else str(
                result.get("error")
                or result.get("status")
                or "Calculation failed."
            )
        )
        self.app.call_from_thread(
            self._publish_tool_call_cell,
            tool_name,
            "ok" if ok else "error",
            description,
            arguments,
            note,
        )
        self.app.call_from_thread(self._finish_calculation_request, result)
        return result

    def _on_calculation_event(self, payload: dict[str, object]) -> None:
        run = payload.get("run")
        if not isinstance(run, dict):
            return
        self._on_calculation_run(dict(run), persist=True, event=payload)

    def _on_calculation_run(
        self,
        run: dict[str, object],
        *,
        persist: bool,
        event: dict[str, object] | None = None,
    ) -> None:
        run_id = str(run.get("run_id") or "")
        if not run_id:
            return
        self._calculation_runs[run_id] = dict(run)
        strip = self.query_one(CalculationStatusStrip)
        strip.update_run(run)
        cell = self._calculation_cells.get(run_id)
        if cell is None and str(run.get("status") or "") != "validating":
            cell = CalculationReceiptCell(run)
            self._calculation_cells[run_id] = cell
            self.query_one(Transcript).add_cell(
                cell,
                turn_id=str(run.get("turn_id") or self._active_turn_id),
            )
        elif cell is not None:
            cell.update_run(run)
        if isinstance(self.app.screen, CalculationMonitor):
            self.app.screen.update_run(run)
        if persist and event is not None:
            self._decision_log_for_calculation(run).write(
                "calculation_event",
                event,
                rationale=str(run.get("stage") or ""),
            )
        self._update_footer_job_counts()

    def _finish_calculation_request(self, result: dict[str, object]) -> None:
        calculation = result.get("calculation")
        if isinstance(calculation, dict):
            self._on_calculation_run(dict(calculation), persist=False)
        status = (
            str(calculation.get("status"))
            if isinstance(calculation, dict)
            else str(result.get("status") or "unknown")
        )
        if not self._worker_is_busy():
            footer = self.query_one(FooterWidget)
            footer.set_phase(
                Phase.FINISHED if status == "completed" else Phase.ERROR
            )
            footer.set_hint(
                "Calculation completed · Ctrl+B for receipt"
                if status == "completed"
                else "Calculation failed · Ctrl+B for diagnostics"
            )
        self.notify(
            "Calculation completed."
            if status == "completed"
            else f"Calculation ended with status {status}.",
            severity="information" if status == "completed" else "error",
            timeout=5,
        )

    def _decision_log_for_calculation(
        self, run: dict[str, object]
    ) -> DecisionLog:
        session = self.active_agent_session
        if (
            session is not None
            and session.state is not None
            and session.decision_log is not None
            and session.state.session_id == str(run.get("session_id") or "")
        ):
            return session.decision_log
        if self._calculation_decision_log is None:
            self._calculation_decision_log = DecisionLog(
                self._tui_session_dir / "decision_log.jsonl"
            )
        return self._calculation_decision_log

    def action_open_jobs_panel(self) -> None:
        self.action_show_calculations()

    def action_show_calculations(self) -> None:
        self._refresh_job_snapshot()
        self._restore_calculation_runs()
        runs = list(self._calculation_runs.values())
        if self.app.plain:
            rows = sorted(self._job_snapshot.values(), key=_jobs_sort_key)
            calculation_lines = [
                (
                    f"{str(run.get('status') or 'unknown'):<18} "
                    f"{str(run.get('program') or 'calc'):<9} "
                    f"{str(run.get('kind') or 'job'):<10} "
                    f"{str(run.get('label') or run.get('run_id') or '')}"
                )
                for run in runs
            ]
            self.post_agent_message(
                "```\n"
                + ("\n".join(calculation_lines) or "No local calculations.")
                + "\n\n"
                + format_jobs_table(rows)
                + "\n```",
                title="Calculations",
            )
            return
        self.app.push_screen(
            CalculationMonitor(runs, self._job_snapshot),
            self._handle_calculation_monitor_action,
        )

    def _handle_calculation_monitor_action(
        self, action: CalculationMonitorAction | None
    ) -> None:
        if action is None:
            self.focus_composer()
            return
        if action.run_id.startswith("job:"):
            job_id = action.run_id.removeprefix("job:")
            if action.action == "extract":
                self._extract_job_result(job_id)
            elif action.action == "cancel":
                self._confirm_cancel(job_id)
            return
        if action.action == "extract":
            result = inspect_calculation(
                action.run_id,
                session_root=str(self.session_root),
            )
            calculation = result.get("calculation")
            if isinstance(calculation, dict):
                self._on_calculation_run(dict(calculation), persist=False)
                self.post_agent_message(
                    _calculation_diagnostic_summary(calculation),
                    title="Calculation result",
                )
            else:
                self.post_error(
                    "Result extraction failed",
                    str(result.get("error") or "No result was found."),
                )
            return
        if action.action == "cancel":
            if self.app.plain:
                self.post_error(
                    "Confirmation required",
                    f"Use /cancel {action.run_id} yes.",
                )
                return
            self.app.push_screen(
                TextPromptOverlay(
                    title="Cancel local calculation",
                    prompt=(
                        f"Type yes to terminate {action.run_id}. "
                        "The calculation process group will receive SIGTERM."
                    ),
                ),
                lambda value, run_id=action.run_id: (
                    self._confirm_local_calculation_cancel(run_id, value)
                ),
            )

    def _confirm_local_calculation_cancel(
        self, run_id: str, value: str | None
    ) -> None:
        if str(value or "").strip().lower() not in {"y", "yes"}:
            self.notify("Calculation cancellation dismissed.", timeout=2)
            return
        run = self._calculation_runs.get(run_id, {})
        pid = run.get("pid")
        cancelled = cancel_calculation(
            run_id,
            int(pid) if isinstance(pid, int) else None,
        )
        if cancelled:
            self.notify(f"Cancellation requested for {run_id}.", timeout=3)
        else:
            self.post_error(
                "Cancellation failed",
                f"No active local process was found for {run_id}.",
            )

    def _has_active_local_calculation(self) -> bool:
        return any(
            str(run.get("execution_mode") or "local") == "local"
            and str(run.get("status") or "")
            in {"validating", "starting", "running"}
            for run in self._calculation_runs.values()
        )

    def _restore_calculation_runs(self) -> None:
        current_cwd = Path.cwd().resolve()
        persisted = [
            run.to_dict()
            for run in load_calculation_runs(self.session_root)
            if Path(run.cwd).resolve() == current_cwd
        ]
        for run in persisted:
            self._calculation_runs[str(run["run_id"])] = run
        self.query_one(CalculationStatusStrip).replace_runs(
            list(self._calculation_runs.values())
        )

