"""Main chat screen for the chemsmart agent TUI."""

from __future__ import annotations

from pathlib import Path
from threading import get_ident
from typing import Iterable

from click.testing import CliRunner
from textual.app import ComposeResult
from textual.screen import Screen
from textual.worker import Worker, WorkerState

from chemsmart.agent.cli import agent
from chemsmart.agent.core import CriticVerdict, Plan
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
    parse_decision_event,
    session_completed,
)
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.screens.sessions import SessionsScreen
from chemsmart.agent.tui.services.log_tailer import LogTailer
from chemsmart.agent.tui.services.session_runner import SessionRunnerMixin
from chemsmart.agent.tui.widgets.cells import (
    AgentMessageCell,
    CriticVerdictCell,
    DryRunInputCell,
    ErrorCell,
    GeometryHandoffCell,
    MethodCell,
    PlanCell,
    RuntimeValidationCell,
    SubmissionPreviewCell,
    UserMessageCell,
)
from chemsmart.agent.tui.widgets.composer import Composer
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.popups import (
    ApprovalOverlay,
    ApprovalResult,
    FilePickerOverlay,
    TextPromptOverlay,
)
from chemsmart.agent.tui.widgets.transcript import Transcript


class ChatScreen(SessionRunnerMixin, Screen):
    DEFAULT_CSS = """
    ChatScreen {
        layout: vertical;
    }

    #chat-body {
        height: 1fr;
    }
    """

    def __init__(self, *, session_root: Path) -> None:
        super().__init__()
        self.session_root = session_root
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
        self._running_risky_steps: set[tuple[str, int]] = set()
        self._job_counts = {"queued": 0, "running": 0, "done": 0, "failed": 0}

    def compose(self) -> ComposeResult:
        yield Transcript(id="chat-body")
        yield Composer()
        yield FooterWidget(id="status-footer")

    def on_mount(self) -> None:
        self.post_agent_message(
            "Welcome to the chemsmart agent TUI. Type a request or use /help."
        )
        self.focus_composer()

    def on_composer_submitted(self, event: Composer.Submitted) -> None:
        text = event.text.strip()
        event.composer.clear_text()
        if text.startswith("/"):
            self._handle_slash_command(text)
            return
        self.start_request(text)

    def start_request(self, text: str) -> None:
        if self._current_worker and not self._current_worker.is_finished:
            self.post_error(
                "Session already running",
                "Wait for the current request to finish before starting a new one.",
            )
            return
        self._stop_tailer()
        self._reset_request_state(clear_transcript=True)
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

    def start_resume(self, session_id: str) -> None:
        session_dir = self.session_root / session_id
        if not session_dir.exists():
            self.post_error(
                "Resume failed", f"Unknown session id: {session_id}"
            )
            return
        self._stop_tailer()
        self._reset_request_state(clear_transcript=True)
        self.query_one(FooterWidget).set_phase(Phase.PLANNING)
        self.query_one(FooterWidget).set_hint(f"Loading session {session_id}")
        self._attach_tailer(session_dir / "decision_log.jsonl")
        if session_completed(session_dir):
            self.post_agent_message(
                f"Loaded completed session `{session_id}`."
            )
            self.query_one(FooterWidget).set_hint(
                "Loaded saved session transcript"
            )
            return
        self._current_worker = self.resume_agent_session(session_id)

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
            else:
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

    def action_soft_cancel(self) -> None:
        if self._quit_armed:
            self.app.exit()
            return
        if not self._current_worker or self._current_worker.is_finished:
            self.app.exit()
            return
        self._quit_armed = True
        self.query_one(FooterWidget).set_hint(
            "Cancellation requested • press Ctrl+C again within 3s to quit"
        )
        self.post_agent_message(
            "Soft cancel requested. The current step may still finish because "
            "the underlying agent run is not cooperatively interruptible yet."
        )
        if self._quit_timer is not None:
            self._quit_timer.stop()
        self._quit_timer = self.set_timer(3, self._disarm_soft_cancel)

    def action_quit_if_empty(self) -> None:
        composer = self.query_one(Composer)
        if composer.text.strip():
            return
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

    def _disarm_soft_cancel(self) -> None:
        self._quit_armed = False
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
        self._update_job_counts(entry)
        event = parse_decision_event(
            entry,
            session_dir=(
                self._tailer_path.parent if self._tailer_path else None
            ),
        )
        transcript = self.query_one(Transcript)
        footer = self.query_one(FooterWidget)

        if isinstance(event, RequestEvent):
            self._current_request = event.request
            if not self._has_user_message(event.request):
                self._user_requests.add(event.request)
                transcript.add_cell(UserMessageCell(event.request))
            footer.set_phase(Phase.PLANNING)
            footer.set_hint("Planning…")
        elif isinstance(event, PlanEvent):
            self._current_plan = event.plan
            self._current_plan_text = event.text
            transcript.add_cell(PlanCell(event.text))
            footer.set_phase(Phase.PLANNING)
            footer.set_hint("Generating dry-run input…")
        elif isinstance(event, MethodEvent):
            transcript.add_cell(MethodCell(event.recommendation))
        elif isinstance(event, DryRunInputEvent):
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
            transcript.add_cell(RuntimeValidationCell(event.validation))
        elif isinstance(event, SubmissionPreviewEvent):
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
            self._current_verdict = event.verdict
            transcript.add_cell(CriticVerdictCell(event.verdict))
            if event.verdict.verdict == "reject":
                footer.set_phase(Phase.ERROR)
                footer.set_hint("Critic blocked execution")
            else:
                footer.set_phase(Phase.DRY_RUN_READY)
                footer.set_hint("Critic verdict ready")
        elif isinstance(event, ErrorEvent):
            transcript.add_cell(
                ErrorCell(event.title, event.message, event.details)
            )
            footer.set_phase(Phase.ERROR)
            footer.set_hint("Agent reported an error")
        elif isinstance(event, SessionSummaryEvent):
            self._pending_approval = False
            self._pending_risky_tool = None
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
            return

    def _has_user_message(self, request: str) -> bool:
        return request in self._user_requests

    def _handle_slash_command(self, raw: str) -> None:
        command, _, remainder = raw.partition(" ")
        argument = remainder.strip()

        if command == "/help":
            self.post_agent_message(
                "\n".join(
                    [
                        "# Phase 2 slash commands",
                        "- `[A] /help` show this help",
                        "- `[D] /dryrun` regenerate the current dry-run",
                        "- `[D] /submit` approve pending HPC submission",
                        "- `[D] /run` approve pending local execution",
                        "- `[P,D] /critic` show the current critic verdict",
                        "- `[P,D,R] /plan` show the current plan",
                        "- `[A] /rationale` show planner rationale",
                        "- `[A] /quit` exit the TUI",
                        "- `[I,F] /clear` clear the transcript",
                        "- `[I] /sessions` browse recent sessions",
                        "- `[I] /resume <session-id>` load or continue a session",
                        "- `[A] /tools` list registered tools",
                        "- `[A] /doctor` run inline diagnostics",
                    ]
                )
            )
        elif command in {"/quit", "/exit"}:
            self.app.exit()
        elif command == "/clear":
            if not self._guard_phase(command, {Phase.IDLE, Phase.FINISHED}):
                return
            self._stop_tailer()
            self._reset_request_state(clear_transcript=True)
            self.post_agent_message("Transcript cleared.")
            self.focus_composer()
        elif command == "/sessions":
            if not self._guard_phase(command, {Phase.IDLE}):
                return
            self.app.push_screen(SessionsScreen(self.session_root))
        elif command == "/resume":
            if not self._guard_phase(command, {Phase.IDLE}):
                return
            if not argument:
                self.app.push_screen(SessionsScreen(self.session_root))
                return
            self.start_resume(argument)
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
            self._request_approval("submit_hpc")
        elif command == "/run":
            if not self._guard_phase(command, {Phase.DRY_RUN_READY}):
                return
            self._request_approval("run_local")
        else:
            self.post_error("Unknown command", raw)

    def _run_inline_cli(self, args: Iterable[str], *, title: str) -> None:
        result = CliRunner().invoke(agent, list(args), catch_exceptions=False)
        text = result.output.strip() or f"{title} completed."
        if result.exit_code == 0:
            self.post_agent_message(f"```\n{text}\n```", title=title)
        else:
            self.post_error(title, text)

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

    def _reset_request_state(self, *, clear_transcript: bool) -> None:
        self._user_requests.clear()
        self._current_request = None
        self._current_plan = None
        self._current_plan_text = None
        self._current_verdict = None
        self._pending_approval = False
        self._pending_risky_tool = None
        self._running_risky_steps.clear()
        self._job_counts = {"queued": 0, "running": 0, "done": 0, "failed": 0}
        footer = self.query_one(FooterWidget)
        footer.set_phase(Phase.IDLE)
        footer.set_hint("Enter to submit • /help for commands")
        footer.reset_job_counts()
        if clear_transcript:
            self.query_one(Transcript).clear_cells()

    def _update_job_counts(self, entry: dict) -> None:
        kind = entry.get("kind")
        payload = entry.get("payload") or {}
        tool = payload.get("tool")
        step_index = int(payload.get("step_index") or -1)
        key = (str(tool), step_index)

        if kind == "tool_call" and tool in {"run_local", "submit_hpc"}:
            self._running_risky_steps.add(key)
        elif kind == "tool_result" and tool == "run_local":
            self._running_risky_steps.discard(key)
            self._job_counts["done"] += 1
        elif (
            kind == "tool_result"
            and tool == "submit_hpc"
            and not payload.get("from_preview")
        ):
            self._running_risky_steps.discard(key)
            self._job_counts["queued"] += 1
        elif kind == "tool_error" and tool in {"run_local", "submit_hpc"}:
            self._running_risky_steps.discard(key)
            self._job_counts["failed"] += 1

        self._job_counts["running"] = len(self._running_risky_steps)
        self.query_one(FooterWidget).set_job_counts(**self._job_counts)
