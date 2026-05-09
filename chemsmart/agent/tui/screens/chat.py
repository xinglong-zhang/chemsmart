"""Main chat screen for the chemsmart agent TUI."""

from __future__ import annotations

import os
from pathlib import Path
from threading import get_ident
from typing import Iterable

from click.testing import CliRunner
from rich.table import Table
from textual import work
from textual.app import ComposeResult
from textual.binding import Binding
from textual.screen import Screen
from textual.worker import Worker, WorkerState

from chemsmart.agent.cli import agent
from chemsmart.agent.core import AgentSession, CriticVerdict, Plan
from chemsmart.agent.tui.events import (
    CriticVerdictEvent,
    DryRunInputEvent,
    ErrorEvent,
    GeometryHandoffEvent,
    IgnoredEvent,
    MethodEvent,
    PlanEvent,
    RequestEvent,
    RuntimeValidationEvent,
    SessionSummaryEvent,
    SubmissionPreviewEvent,
    ToolCallEvent,
    ToolPreviewEvent,
    parse_decision_event,
    session_completed,
)
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.screens.jobs_panel import JobsPanel, JobsPanelAction
from chemsmart.agent.tui.screens.sessions import SessionsScreen
from chemsmart.agent.tui.services.job_poller import (
    JobPollerError,
    JobPollerMixin,
    JobStateReader,
    JobStatusUpdated,
    available_server_names,
    cancel_job,
    collect_job_snapshot,
    extract_run_result,
    format_jobs_table,
    queue_snapshot,
)
from chemsmart.agent.tui.services.log_tailer import LogTailer
from chemsmart.agent.tui.services.session_runner import SessionRunnerMixin
from chemsmart.agent.tui.widgets.cells import (
    AgentMessageCell,
    CriticVerdictCell,
    DryRunInputCell,
    ErrorCell,
    GeometryHandoffCell,
    JobStatusCell,
    MethodCell,
    MoleculeCell,
    PlanCell,
    RunResultCell,
    RuntimeValidationCell,
    SubmissionPreviewCell,
    UserMessageCell,
)
from chemsmart.agent.tui.widgets.composer import Composer
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.header import ChemsmartHeader
from chemsmart.agent.tui.widgets.popups import (
    ApprovalOverlay,
    ApprovalResult,
    CwdMismatchChoice,
    CwdMismatchOverlay,
    FilePickerOverlay,
    TextPromptOverlay,
)
from chemsmart.agent.tui.widgets.transcript import Transcript
from chemsmart.io.molecules.structure import Molecule


class ChatScreen(JobPollerMixin, SessionRunnerMixin, Screen):
    BINDINGS = [Binding("j", "open_jobs_panel", "Jobs", show=False)]

    DEFAULT_CSS = """
    ChatScreen {
        layout: vertical;
    }

    #chat-body {
        height: 1fr;
    }
    """

    def __init__(
        self,
        *,
        session_root: Path,
        job_poll_interval: float = 5.0,
    ) -> None:
        super().__init__()
        self.session_root = session_root
        self.job_poll_interval = job_poll_interval
        self._tailer: LogTailer | None = None
        self._tailer_path: Path | None = None
        self._current_worker: Worker | None = None
        self._quit_armed = False
        self._quit_timer = None
        self._session_poll_timer = None
        self._user_requests: set[str] = set()
        self._current_request: str | None = None
        self._current_plan: Plan | None = None
        self._current_plan_text: str | None = None
        self._current_verdict: CriticVerdict | None = None
        self._pending_approval = False
        self._pending_risky_tool: str | None = None
        self._approval_session_granted = False
        self._latest_dry_run_content: str | None = None
        self._job_snapshot: dict[str, dict] = {}
        self._job_cells: dict[str, JobStatusCell] = {}
        self._rendered_run_results: set[str] = set()
        self._cancelled_job_ids: set[str] = set()
        self._active_server_name = self._default_active_server_name()

    def compose(self) -> ComposeResult:
        yield ChemsmartHeader(id="chat-header")
        yield Transcript(id="chat-body")
        yield Composer()
        yield FooterWidget(id="status-footer")

    def on_mount(self) -> None:
        self.post_agent_message(
            "## Welcome to the chemsmart agent\n\n"
            "Plan calculations in natural language, preview input files, "
            "and run pre-flight checks before execution.\n\n"
            "Start with: /doctor · Ask safely: ask mode · "
            "Real workflows: run mode\n\n"
            "Example: `single-point on examples/h2o.xyz at B3LYP/6-31G(d) "
            "Gaussian`"
        )
        self.focus_composer()
        self.query_one(FooterWidget).update_draft("")
        self.run_job_poller(self.job_poll_interval)
        self._refresh_job_snapshot()

    def on_composer_submitted(self, event: Composer.Submitted) -> None:
        text = event.text.strip()
        event.composer.clear_text()
        self.query_one(FooterWidget).update_draft("")
        if text.startswith("/"):
            self._handle_slash_command(text)
            return
        self.start_request(text)

    def on_text_area_changed(self, event) -> None:
        if getattr(event.text_area, "id", None) == "composer":
            self.query_one(FooterWidget).update_draft(event.text_area.text)

    def on_job_poller_error(self, message: JobPollerError) -> None:
        self.post_error(
            message.summary,
            message.summary,
            {"traceback": message.details},
        )

    def on_job_status_updated(self, message: JobStatusUpdated) -> None:
        snapshot = self._job_snapshot.setdefault(
            message.job_id, {"job_id": message.job_id}
        )
        snapshot.update(message.fields)
        if message.job_id in self._cancelled_job_ids:
            snapshot["status"] = "cancelled"
        self._update_footer_job_counts()
        active_session_dir = self.current_session_dir()
        active_session_id = (
            active_session_dir.name if active_session_dir is not None else None
        )
        track_in_transcript = snapshot.get("session_id") == active_session_id
        cell = self._job_cells.get(message.job_id)
        if cell is None:
            if not track_in_transcript:
                return
            cell = JobStatusCell(message.job_id, snapshot)
            self._job_cells[message.job_id] = cell
            self.query_one(Transcript).add_cell(cell)
        else:
            cell.apply_update(snapshot)
        if track_in_transcript and snapshot.get("status") == "done":
            self._maybe_render_run_result(message.job_id, snapshot)

    @work(
        thread=True,
        exclusive=True,
        exit_on_error=False,
        group="agent-session",
        name="agent-run",
    )
    def run_agent_session(self, request: str) -> dict[str, object]:
        self.active_resume_id = None
        if self.active_agent_session is None:
            self.active_agent_session = AgentSession(
                session_root=str(self.session_root)
            )
        return self.active_agent_session.run(
            request,
            dry_submit=True,
            pause_before_risky=True,
        )

    @work(
        thread=True,
        exclusive=True,
        exit_on_error=False,
        group="agent-session",
        name="agent-resume",
    )
    def resume_agent_session(
        self,
        session_id: str,
        *,
        cwd_override: str | None = None,
    ) -> dict[str, object]:
        self.active_resume_id = session_id
        self.active_agent_session = AgentSession.load(
            session_id,
            session_root=str(self.session_root),
            cwd_override=cwd_override,
        )
        return self.active_agent_session.continue_loaded_session(
            dry_submit=True,
            pause_before_risky=True,
            rerender_plan=False,
        )

    @work(
        thread=True,
        exclusive=True,
        exit_on_error=False,
        group="agent-session",
        name="agent-execute",
    )
    def execute_agent_session(self, session_id: str) -> dict[str, object]:
        self.active_resume_id = session_id
        self.active_agent_session = AgentSession.load(
            session_id,
            session_root=str(self.session_root),
        )
        return self.active_agent_session.continue_loaded_session(
            dry_submit=False,
            pause_before_risky=False,
            rerender_plan=False,
        )

    def _load_agent_session(
        self,
        session_id: str,
        *,
        cwd_override: str | None = None,
    ) -> AgentSession:
        self.active_resume_id = session_id
        self.active_agent_session = AgentSession.load(
            session_id,
            session_root=str(self.session_root),
            cwd_override=cwd_override,
        )
        return self.active_agent_session

    def start_request(self, text: str) -> None:
        if self._current_worker and not self._current_worker.is_finished:
            self.post_error(
                "Session already running",
                "Wait for the current request to finish before starting a new request.",
            )
            return
        keep_conversational = self.active_agent_session is not None
        if not keep_conversational:
            self._stop_tailer()
        self._reset_request_state(
            clear_transcript=True,
            keep_conversational=keep_conversational,
        )
        self._current_request = text
        self.query_one(FooterWidget).set_phase(Phase.PLANNING)
        self.query_one(FooterWidget).set_hint("Agent is planning…")
        self.query_one(Transcript).add_cell(UserMessageCell(text))
        self._user_requests.add(text)
        self._current_worker = self.run_agent_session(text)
        if self._session_poll_timer is not None:
            self._session_poll_timer.stop()
        self._session_poll_timer = self.set_interval(
            0.1, self._attach_live_tailer, pause=False
        )

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
        if event.worker.group != "agent-session":
            return
        if event.state == WorkerState.SUCCESS:
            result = event.worker.result or {}
            self._current_worker = event.worker
            if result.get("pending_approval"):
                self._pending_approval = True
                self._pending_risky_tool = result.get("next_risky_tool")
                self.query_one(FooterWidget).set_phase(Phase.DRY_RUN_READY)
                self.query_one(FooterWidget).set_hint(
                    f"Approval required for {self._pending_risky_tool}"
                )
                if self._approval_session_granted:
                    self._approve_current_request()
                return
            if result.get("blocked"):
                self.query_one(FooterWidget).set_hint("Execution blocked")
            elif (
                isinstance(result.get("plan"), Plan)
                and result["plan"].is_chitchat()
            ):
                self.query_one(FooterWidget).set_phase(Phase.IDLE)
                self.query_one(FooterWidget).set_hint("Ready")
            else:
                self.query_one(FooterWidget).set_phase(Phase.FINISHED)
                self.query_one(FooterWidget).set_hint("Session complete")
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
            self.query_one(FooterWidget).set_phase(Phase.IDLE)
            self.query_one(FooterWidget).set_hint("Ready")
            self._current_worker = event.worker
        self._refresh_job_snapshot()

    def action_open_jobs_panel(self) -> None:
        self._refresh_job_snapshot()
        if self.app.plain:
            rows = sorted(self._job_snapshot.values(), key=_jobs_sort_key)
            self.post_agent_message(
                f"```\n{format_jobs_table(rows)}\n```",
                title="Jobs",
            )
            return
        self.app.push_screen(
            JobsPanel(self._job_snapshot), self._handle_jobs_panel_action
        )

    def action_soft_cancel(self) -> None:
        if self._quit_armed:
            if self._quit_timer is not None:
                self._quit_timer.stop()
                self._quit_timer = None
            self.app.exit()
            return
        self._quit_armed = True
        self.query_one(FooterWidget).set_hint(
            "Press Ctrl+C again within 3s to quit"
        )
        if self._current_worker and not self._current_worker.is_finished:
            self.post_agent_message(
                "Soft cancel requested. The current step may still finish "
                "because the underlying agent run is not cooperatively "
                "interruptible yet."
            )
        if self._quit_timer is not None:
            self._quit_timer.stop()
        self._quit_timer = self.set_timer(3, self._disarm_soft_cancel)

    def action_quit_if_empty(self) -> None:
        composer = self.query_one(Composer)
        if composer.text.strip():
            return
        self._reset_request_state(clear_transcript=False, clear_session=True)
        self.app.exit()

    def action_refresh_screen(self) -> None:
        self.refresh(repaint=True, layout=True)

    def action_dismiss_overlay(self) -> None:
        if self.app.screen is not self:
            self.app.pop_screen()

    def focus_composer(self) -> None:
        self.query_one(Composer).focus()

    def post_agent_message(self, text: str, *, title: str = "Agent") -> None:
        self.query_one(Transcript).add_cell(
            AgentMessageCell(text, title=title)
        )

    def post_error(
        self, title: str, message: str, details: dict | None = None
    ) -> None:
        self.query_one(Transcript).add_cell(ErrorCell(title, message, details))

    def open_file_picker(self) -> None:
        self.app.push_screen(
            FilePickerOverlay(Path.cwd()), self._handle_file_pick
        )

    def edit_method_from_cell(self, recommendation: dict) -> None:
        method = recommendation.get("functional") or "manual"
        basis = recommendation.get("basis") or "manual"
        self.app.push_screen(
            TextPromptOverlay(
                title="Revise method",
                prompt=(
                    "Describe the corrected method/basis to use. "
                    f"Current: {method}/{basis}"
                ),
            ),
            self._handle_method_revision,
        )

    def _handle_file_pick(self, value: str | None) -> None:
        if not value:
            return
        self.query_one(Composer).insert_file_reference(value)

    def _handle_method_revision(self, value: str | None) -> None:
        if not value or not self._current_request:
            return
        self.start_request(self._corrected_request(value))

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
        import threading

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
        state = JobStateReader.load(self.session_root / session_id)
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

    def _attach_live_tailer(self) -> None:
        session_dir = self.current_session_dir()
        if session_dir is None:
            return
        log_path = session_dir / "decision_log.jsonl"
        if self._tailer_path == log_path:
            return
        if log_path.exists():
            self._attach_tailer(log_path)
            if self._session_poll_timer is not None:
                self._session_poll_timer.stop()
                self._session_poll_timer = None

    def _attach_tailer(self, log_path: Path) -> None:
        self._stop_tailer()
        self._tailer_path = log_path
        self._tailer = LogTailer(log_path, self._on_log_entry)
        self._tailer.start()

    def _stop_tailer(self) -> None:
        if self._tailer is not None:
            self._tailer.stop()
        self._tailer = None
        self._tailer_path = None

    def _on_log_entry(self, entry: dict) -> None:
        app_thread_id = getattr(self.app, "_thread_id", None)
        if app_thread_id == get_ident():
            self._apply_log_entry(entry)
            return
        self.app.call_from_thread(self._apply_log_entry, entry)

    def _apply_log_entry(self, entry: dict) -> None:
        event = parse_decision_event(
            entry,
            session_dir=(
                self._tailer_path.parent if self._tailer_path else None
            ),
        )
        payload = entry.get("payload") or {}
        transcript = self.query_one(Transcript)
        footer = self.query_one(FooterWidget)
        workflow_cell = getattr(self, "_workflow_cell", None)

        if isinstance(event, RequestEvent):
            self._current_request = event.request
            if not self._has_user_message(event.request):
                self._user_requests.add(event.request)
                transcript.add_cell(UserMessageCell(event.request))
            footer.set_phase(Phase.PLANNING)
            footer.set_hint("Interpreting your request…")
        elif isinstance(event, PlanEvent):
            self._current_plan = event.plan
            self._current_plan_text = event.text
            if event.plan.is_chitchat():
                self._workflow_cell = None
                transcript.add_cell(
                    AgentMessageCell(
                        event.plan.rationale or "(no reply)",
                        title="Reply",
                    )
                )
                footer.set_phase(Phase.IDLE)
                footer.set_hint("Ready")
            elif not event.plan.steps:
                self._workflow_cell = None
                transcript.add_cell(
                    AgentMessageCell(
                        event.plan.rationale or "(no rationale)",
                        title="Advisory",
                    )
                )
                footer.set_phase(Phase.FINISHED)
                footer.set_hint("Response ready")
            else:
                workflow_cell = PlanCell(event.plan)
                self._workflow_cell = workflow_cell
                transcript.add_cell(workflow_cell)
                footer.set_phase(Phase.PLANNING)
                footer.set_hint("Filling in tool steps…")
        elif isinstance(event, ToolCallEvent):
            if workflow_cell is not None:
                if event.status == "running":
                    workflow_cell.mark_started(
                        event.step_index,
                        event.tool,
                        event.args,
                    )
                else:
                    workflow_cell.mark_completed(
                        event.step_index,
                        event.tool,
                        event.payload,
                    )
            if event.tool == "run_local":
                footer.set_phase(Phase.RUNNING)
                footer.set_hint(
                    "Local run in progress…"
                    if event.status == "running"
                    else "Local run complete"
                )
            else:
                footer.set_phase(Phase.PLANNING)
                footer.set_hint(
                    f"{event.tool} in progress…"
                    if event.status == "running"
                    else f"{event.tool} complete"
                )
        elif isinstance(event, ToolPreviewEvent):
            if workflow_cell is not None and event.status == "running":
                workflow_cell.mark_started(
                    event.step_index,
                    event.tool,
                    event.args,
                )
            footer.set_phase(Phase.DRY_RUN_READY)
            footer.set_hint("Preparing submission preview…")
        elif isinstance(event, MethodEvent):
            transcript.add_cell(MethodCell(event.recommendation))
        elif isinstance(event, DryRunInputEvent):
            if workflow_cell is not None:
                workflow_cell.mark_completed(
                    event.step_index,
                    "dry_run_input",
                    {
                        "inputfile": event.inputfile,
                        "content": event.content,
                    },
                )
            transcript.add_cell(
                DryRunInputCell(
                    event.content,
                    inputfile=event.inputfile,
                    previous_content=self._latest_dry_run_content,
                )
            )
            self._latest_dry_run_content = event.content
            footer.set_phase(Phase.DRY_RUN_READY)
            footer.set_hint("Dry-run input ready")
        elif isinstance(event, RuntimeValidationEvent):
            if workflow_cell is not None:
                workflow_cell.mark_completed(
                    event.step_index,
                    "validate_runtime",
                    event.validation,
                )
            transcript.add_cell(RuntimeValidationCell(event.validation))
            footer.set_phase(Phase.DRY_RUN_READY)
            footer.set_hint("Runtime check ready")
        elif isinstance(event, SubmissionPreviewEvent):
            if workflow_cell is not None:
                workflow_cell.mark_completed(
                    event.step_index,
                    "submit_hpc",
                    event.preview,
                )
            transcript.add_cell(SubmissionPreviewCell(event.preview))
            footer.set_phase(Phase.DRY_RUN_READY)
            footer.set_hint("Submission preview ready")
        elif isinstance(event, GeometryHandoffEvent):
            transcript.add_cell(
                GeometryHandoffCell(
                    event.molecule,
                    session_dir=event.session_dir,
                )
            )
        elif isinstance(event, CriticVerdictEvent):
            if (
                self._current_plan is not None
                and self._current_plan.is_chitchat()
            ):
                footer.set_phase(Phase.IDLE)
                footer.set_hint("Ready")
                return
            self._current_verdict = event.verdict
            transcript.add_cell(CriticVerdictCell(event.verdict))
            if event.verdict.verdict == "reject":
                footer.set_phase(Phase.ERROR)
                footer.set_hint("Critic blocked execution")
            else:
                footer.set_phase(Phase.DRY_RUN_READY)
                footer.set_hint("Critic verdict ready")
        elif isinstance(event, ErrorEvent):
            details = event.details if isinstance(event.details, dict) else {}
            if workflow_cell is not None and details:
                workflow_cell.mark_failed(
                    int(details.get("step_index") or 0),
                    str(details.get("tool") or ""),
                    event.message,
                )
            transcript.add_cell(
                ErrorCell(event.title, event.message, event.details)
            )
            footer.set_phase(Phase.ERROR)
            footer.set_hint("Agent reported an error")
        elif isinstance(event, SessionSummaryEvent):
            self._pending_approval = False
            self._pending_risky_tool = None
            if payload.get("request_intent") == "chitchat":
                footer.set_phase(Phase.IDLE)
                footer.set_hint("Ready")
                return
            footer.set_phase(Phase.ERROR if event.blocked else Phase.FINISHED)
            footer.set_hint(
                "Blocked" if event.blocked else "Finished successfully"
            )
            summary = (
                f"Session finished ({event.total_steps_executed}/"
                f"{event.total_steps_planned} steps)."
            )
            if event.block_reason:
                summary += f" Block reason: {event.block_reason}."
            transcript.add_cell(AgentMessageCell(summary, title="Summary"))
        elif isinstance(event, IgnoredEvent):
            pass

        if payload.get("tool") in {"run_local", "submit_hpc"}:
            self._refresh_job_snapshot()

    def _has_user_message(self, request: str) -> bool:
        return request in self._user_requests

    def _handle_slash_command(self, raw: str) -> None:
        command, _, remainder = raw.partition(" ")
        argument = remainder.strip()

        if command == "/help":
            self._show_help()
        elif command in {"/quit", "/exit"}:
            self._reset_request_state(
                clear_transcript=False, clear_session=True
            )
            self.app.exit()
        elif command == "/clear":
            if not self._guard_phase(command, {Phase.IDLE, Phase.FINISHED}):
                return
            self._stop_tailer()
            self._reset_request_state(
                clear_transcript=True, clear_session=True
            )
            self.notify("Transcript cleared.", timeout=3)
            self.focus_composer()
        elif command == "/sessions":
            if not self._guard_phase(command, {Phase.IDLE}):
                return
            if self.app.plain:
                self._show_sessions_snapshot()
                return
            self.app.push_screen(SessionsScreen(self.session_root))
        elif command == "/resume":
            if not self._guard_phase(command, {Phase.IDLE}):
                return
            if not argument:
                if self.app.plain:
                    self._show_sessions_snapshot()
                    return
                self.app.push_screen(SessionsScreen(self.session_root))
                return
            self._resume_or_prompt(argument)
        elif command == "/jobs":
            self.action_open_jobs_panel()
        elif command == "/queue":
            self._show_queue_snapshot()
        elif command == "/server":
            self._switch_active_server(argument)
        elif command == "/cancel":
            job_id, confirmed = self._parse_cancel_argument(argument)
            if not job_id:
                self.post_error(
                    "Missing job id",
                    "Usage: /cancel <job-id> [yes]",
                )
                return
            if confirmed:
                self._handle_cancel_confirmation(job_id, "yes")
                return
            if self.app.plain:
                self.post_error(
                    "Confirmation required",
                    "Plain mode uses /cancel <job-id> yes.",
                )
                return
            self._confirm_cancel(job_id)
        elif command == "/extract":
            if not argument:
                self.post_error(
                    "Missing target",
                    "Usage: /extract <job-id|inputfile>",
                )
                return
            self._extract_job_result(argument)
        elif command == "/molecule":
            if not argument:
                self.post_error("Missing path", "Usage: /molecule <path>")
                return
            self._show_molecule(argument)
        elif command == "/tools":
            self._run_inline_cli(["tools"], title="Tools")
        elif command == "/doctor":
            self._run_inline_cli(["doctor"], title="Doctor")
        elif command == "/dryrun":
            if not self._guard_phase(command, {Phase.DRY_RUN_READY}):
                return
            if not self._current_request:
                self.post_error("No request", "There is no active request.")
                return
            self.start_request(self._current_request)
        elif command == "/critic":
            if not self._guard_phase(
                command, {Phase.PLANNING, Phase.DRY_RUN_READY}
            ):
                return
            if self._current_verdict is None:
                self.post_error(
                    "No critic verdict", "The critic has not finished yet."
                )
                return
            self.query_one(Transcript).add_cell(
                CriticVerdictCell(self._current_verdict)
            )
        elif command == "/plan":
            if not self._guard_phase(
                command, {Phase.PLANNING, Phase.DRY_RUN_READY, Phase.RUNNING}
            ):
                return
            if not self._current_plan_text:
                self.post_error("No plan", "The planner has not finished yet.")
                return
            self.query_one(Transcript).add_cell(
                PlanCell(self._current_plan_text)
            )
        elif command == "/rationale":
            rationale = (
                self._current_plan.rationale if self._current_plan else None
            )
            if not rationale:
                self.post_error(
                    "No rationale", "No planner rationale is available."
                )
                return
            self.post_agent_message(rationale, title="Rationale")
        elif command == "/submit":
            if not self._guard_phase(command, {Phase.DRY_RUN_READY}):
                return
            if (
                not self._pending_approval
                or self._pending_risky_tool != "submit_hpc"
            ):
                self._request_approval("submit_hpc")
                return
            if self._handle_inline_approval(argument):
                return
            if self.app.plain:
                self.post_error(
                    "Approval required",
                    "Plain mode uses /submit yes|session|no|revise <instruction>.",
                )
                return
            self._request_approval("submit_hpc")
        elif command == "/run":
            if not self._guard_phase(command, {Phase.DRY_RUN_READY}):
                return
            if (
                not self._pending_approval
                or self._pending_risky_tool != "run_local"
            ):
                self._request_approval("run_local")
                return
            if self._handle_inline_approval(argument):
                return
            if self.app.plain:
                self.post_error(
                    "Approval required",
                    "Plain mode uses /run yes|session|no|revise <instruction>.",
                )
                return
            self._request_approval("run_local")
        else:
            self.post_error("Unknown command", raw)

    def _show_help(self) -> None:
        table = Table(show_header=True, box=None, padding=(0, 1))
        table.add_column("Phase", style="dim", no_wrap=True)
        table.add_column("Command", style="bold")
        table.add_column("Description")
        rows = [
            ("[A]", "/help", "show this help"),
            ("[A]", "/jobs", "open the jobs panel"),
            ("[A]", "/queue", "show the current queue snapshot"),
            ("[A]", "/server <name>", "switch the active HPC server"),
            ("[A]", "/molecule <path>", "load and preview a molecule"),
            ("[A]", "/cancel <job-id> [yes]", "cancel a queued/running job"),
            ("[F]", "/extract <job-id|inputfile>", "parse a final result"),
            ("[D]", "/dryrun", "regenerate the current dry-run"),
            ("[D]", "/submit", "approve pending HPC submission"),
            ("[D]", "/run", "approve pending local execution"),
            ("[P,D]", "/critic", "show the current critic verdict"),
            ("[P,D,R]", "/plan", "show the current plan"),
            ("[A]", "/rationale", "show planner rationale"),
            ("[I,F]", "/clear", "clear the transcript"),
            ("[I]", "/sessions", "browse recent sessions"),
            ("[I]", "/resume <session-id>", "load or continue a session"),
            ("[A]", "/tools", "list registered tools"),
            ("[A]", "/doctor", "run inline diagnostics"),
            ("[A]", "/quit", "exit the TUI"),
        ]
        for row in rows:
            table.add_row(*row)
        self.query_one(Transcript).add_cell(
            AgentMessageCell(table, title="Help")
        )

    def _show_sessions_snapshot(self) -> None:
        lines = ["Sessions", "", "Use /resume <session-id> to load one."]
        if self.session_root.exists():
            session_dirs = sorted(
                [
                    path
                    for path in self.session_root.iterdir()
                    if path.is_dir()
                ],
                reverse=True,
            )
        else:
            session_dirs = []
        if not session_dirs:
            lines.extend(["", "No sessions found."])
        for session_dir in session_dirs[:10]:
            lines.append(f"- {session_dir.name}")
        body = "\n".join(lines)
        self.post_agent_message(f"```\n{body}\n```", title="Sessions")

    def _parse_cancel_argument(self, argument: str) -> tuple[str | None, bool]:
        if not argument:
            return None, False
        parts = argument.split(maxsplit=1)
        if len(parts) == 1:
            return parts[0], False
        if parts[1].strip().lower() in {"y", "yes"}:
            return parts[0], True
        self.post_error(
            "Unknown confirmation",
            "Use /cancel <job-id> yes to confirm.",
        )
        return None, False

    def _handle_inline_approval(self, argument: str) -> bool:
        if not argument:
            return False
        keyword, _, remainder = argument.partition(" ")
        keyword = keyword.strip().lower()
        corrective_text = remainder.strip() or None
        if keyword in {"y", "yes"}:
            self._handle_approval_result(ApprovalResult("y"))
            return True
        if keyword in {"n", "no"}:
            self._handle_approval_result(ApprovalResult("n"))
            return True
        if keyword in {"s", "session"}:
            self._handle_approval_result(ApprovalResult("s"))
            return True
        if keyword in {"r", "revise"}:
            if not corrective_text:
                self.post_error(
                    "Missing instruction",
                    "Usage: /run revise <instruction> or /submit revise <instruction>.",
                )
                return True
            self._handle_approval_result(
                ApprovalResult("r", corrective_text=corrective_text)
            )
            return True
        self.post_error(
            "Unknown approval response",
            "Use yes, no, session, or revise <instruction>.",
        )
        return True

    def _run_inline_cli(self, args: Iterable[str], *, title: str) -> None:
        import threading

        from chemsmart.agent.cli import sanitize_inline_cli_output

        command_args = list(args)

        def run_inline() -> None:
            try:
                result = CliRunner().invoke(
                    agent, command_args, catch_exceptions=False
                )
            except Exception as exc:
                self.app.call_from_thread(self.post_error, title, str(exc))
                return

            text = sanitize_inline_cli_output(result.output)
            text = text or f"{title} completed."

            def publish() -> None:
                if result.exit_code == 0:
                    self.post_agent_message(f"```\n{text}\n```", title=title)
                else:
                    self.post_error(title, text)

            self.app.call_from_thread(publish)

        threading.Thread(target=run_inline, daemon=True).start()

    def _guard_phase(self, command: str, allowed: set[Phase]) -> bool:
        current = self.query_one(FooterWidget).phase
        if current in allowed:
            return True
        allowed_text = ", ".join(sorted(phase.value for phase in allowed))
        self.post_error(
            "Command unavailable",
            f"{command} is only available in: {allowed_text}.",
        )
        return False

    def _request_approval(self, expected_tool: str) -> None:
        if (
            not self._pending_approval
            or self._pending_risky_tool != expected_tool
        ):
            self.post_error(
                "Approval unavailable",
                f"No pending `{expected_tool}` action is waiting for approval.",
            )
            return
        self.app.push_screen(
            ApprovalOverlay(
                action=expected_tool,
                request=self._current_request or "",
            ),
            self._handle_approval_result,
        )

    def _handle_approval_result(self, result: ApprovalResult | None) -> None:
        if result is None:
            return
        if result.choice == "n":
            self.query_one(FooterWidget).set_hint("Approval denied")
            return
        if result.choice == "s":
            self._approval_session_granted = True
            self._approve_current_request()
            return
        if result.choice == "y":
            self._approve_current_request()
            return
        if result.choice == "r" and result.corrective_text:
            self.start_request(self._corrected_request(result.corrective_text))

    def _approve_current_request(self) -> None:
        session_dir = self.current_session_dir()
        if session_dir is None:
            self.post_error("No session", "No paused session is available.")
            return
        self._pending_approval = False
        self.query_one(FooterWidget).set_phase(Phase.RUNNING)
        self.query_one(FooterWidget).set_hint(
            f"Executing {self._pending_risky_tool or 'workflow'}…"
        )
        self._current_worker = self.execute_agent_session(session_dir.name)

    def _corrected_request(self, corrective_text: str) -> str:
        original = self._current_request or ""
        return (
            f"Corrective instruction: {corrective_text}\n\n"
            f"Original request:\n{original}"
        )

    def _reset_request_state(
        self,
        *,
        clear_transcript: bool,
        clear_session: bool = False,
        keep_conversational: bool = False,
    ) -> None:
        if clear_transcript and not keep_conversational:
            self._user_requests.clear()
        if clear_session:
            self.active_agent_session = None
            self.active_resume_id = None
        self._current_request = None
        self._current_plan = None
        self._current_plan_text = None
        self._current_verdict = None
        self._workflow_cell = None
        self._pending_approval = False
        self._pending_risky_tool = None
        self._latest_dry_run_content = None
        self._job_cells.clear()
        self._rendered_run_results.clear()
        footer = self.query_one(FooterWidget)
        footer.set_phase(Phase.IDLE)
        footer.set_hint("Enter to submit • /help for commands")
        footer.reset_job_counts()
        if clear_transcript:
            transcript = self.query_one(Transcript)
            if keep_conversational:
                transcript.clear_turn_chrome()
            else:
                transcript.clear_cells()
        self._refresh_job_snapshot()

    def _emit_job_update(self, job_id: str, fields: dict) -> None:
        self.post_message(JobStatusUpdated(job_id, fields))
        if isinstance(self.app.screen, JobsPanel):
            self.app.screen.post_message(JobStatusUpdated(job_id, fields))

    def _refresh_job_snapshot(self) -> None:
        from chemsmart.agent.tui.services import job_poller as job_poller_service

        _ = collect_job_snapshot
        snapshot = (
            job_poller_service.get_cached_job_snapshot(self.session_root) or {}
        )
        job_poller_service.request_job_snapshot_refresh(self.session_root)
        if not isinstance(snapshot, dict):
            self.post_error(
                "Job snapshot failed", "Snapshot payload was not a mapping."
            )
            return
        for job_id in self._cancelled_job_ids:
            if job_id in snapshot:
                snapshot[job_id]["status"] = "cancelled"
        for job_id, current in snapshot.items():
            changed = {
                key: value
                for key, value in current.items()
                if self._job_snapshot.get(job_id, {}).get(key) != value
            }
            if changed:
                self._emit_job_update(job_id, changed)
        self._job_snapshot = snapshot
        self._update_footer_job_counts()

    def _update_footer_job_counts(self) -> None:
        counts = {"queued": 0, "running": 0, "failed": 0}
        for snapshot in self._job_snapshot.values():
            status = snapshot.get("status")
            if status in counts:
                counts[str(status)] += 1
        self.query_one(FooterWidget).set_job_counts(**counts)

    def _maybe_render_run_result(self, job_id: str, snapshot: dict) -> None:
        if job_id in self._rendered_run_results:
            return
        output_path = snapshot.get("output_path")
        if not output_path:
            return
        try:
            result = extract_run_result(str(output_path), self._job_snapshot)
        except Exception as exc:
            self.post_error("Run result parse failed", str(exc))
            return
        self.query_one(Transcript).add_cell(RunResultCell(result))
        self._rendered_run_results.add(job_id)

    def _show_queue_snapshot(self) -> None:
        self._refresh_job_snapshot()
        try:
            rows = (
                queue_snapshot(
                    self._job_snapshot or {},
                    server_name=self._active_server_name,
                )
                or []
            )
            table = format_jobs_table(rows)
        except Exception as exc:
            self.post_error("Queue snapshot failed", str(exc))
            return
        self.post_agent_message(f"```\n{table}\n```", title="Queue")

    def _switch_active_server(self, name: str) -> None:
        try:
            names = available_server_names() or []
        except Exception as exc:
            self.post_error("Server lookup failed", str(exc))
            return
        if not name:
            current = self._active_server_name or "(none)"
            self.post_agent_message(
                f"Active server: {current}\nAvailable: {', '.join(names) or '(none)'}",
                title="Server",
            )
            return
        if not names:
            self.post_error(
                "No servers configured", "No HPC servers are configured."
            )
            return
        if name not in names:
            self.post_error(
                "Unknown server",
                f"{name} is not configured. Available: {', '.join(names)}",
            )
            return
        self._active_server_name = name
        self.post_agent_message(
            f"Active server set to `{name}`.", title="Server"
        )

    def _show_molecule(self, path: str) -> None:
        source = str(Path(path).expanduser())
        try:
            molecule = Molecule.from_filepath(path)
        except Exception as exc:
            self.post_error("Molecule load failed", str(exc))
            return
        self.query_one(Transcript).add_cell(
            MoleculeCell(molecule, source=source)
        )

    def _default_active_server_name(self) -> str | None:
        try:
            names = available_server_names() or []
        except Exception:
            return None
        return next(iter(names), None) if len(names) == 1 else None


def _jobs_sort_key(job: dict) -> tuple[int, str, str]:
    order = {"running": 0, "queued": 1, "failed": 2, "cancelled": 3, "done": 4}
    return (
        order.get(str(job.get("status")), 9),
        str(job.get("raw_started") or ""),
        str(job.get("job_id") or ""),
    )
