"""Main Phase 1 chat screen."""

from __future__ import annotations

from pathlib import Path
from threading import get_ident
from typing import Iterable

from click.testing import CliRunner
from textual.app import ComposeResult
from textual.screen import Screen
from textual.worker import Worker, WorkerState

from chemsmart.agent.cli import agent
from chemsmart.agent.tui.events import (
    CriticVerdictEvent,
    DryRunInputEvent,
    ErrorEvent,
    IgnoredEvent,
    PlanEvent,
    RequestEvent,
    SessionSummaryEvent,
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
    PlanCell,
    UserMessageCell,
)
from chemsmart.agent.tui.widgets.composer import Composer
from chemsmart.agent.tui.widgets.footer import FooterWidget
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
        self._user_requests.clear()
        self.query_one(Transcript).clear_cells()
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
        self._user_requests.clear()
        self.query_one(Transcript).clear_cells()
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
            self.query_one(FooterWidget).set_hint("Session complete")
            self._current_worker = event.worker
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
        event = parse_decision_event(entry)
        transcript = self.query_one(Transcript)
        footer = self.query_one(FooterWidget)

        if isinstance(event, RequestEvent):
            if not self._has_user_message(event.request):
                self._user_requests.add(event.request)
                transcript.add_cell(UserMessageCell(event.request))
            footer.set_phase(Phase.PLANNING)
            footer.set_hint("Planning…")
        elif isinstance(event, PlanEvent):
            transcript.add_cell(PlanCell(event.text))
            footer.set_phase(Phase.PLANNING)
            footer.set_hint("Generating dry-run input…")
        elif isinstance(event, DryRunInputEvent):
            transcript.add_cell(
                DryRunInputCell(event.content, inputfile=event.inputfile)
            )
            footer.set_phase(Phase.DRY_RUN_READY)
            footer.set_hint("Dry-run input ready")
        elif isinstance(event, CriticVerdictEvent):
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
                        "# Phase 1 slash commands",
                        "- `/help` show this help",
                        "- `/quit` exit the TUI",
                        "- `/clear` clear the transcript",
                        "- `/sessions` browse recent sessions",
                        "- `/resume <session-id>` load or continue a session",
                        "- `/tools` list registered tools",
                        "- `/doctor` run inline diagnostics",
                    ]
                )
            )
        elif command in {"/quit", "/exit"}:
            self.app.exit()
        elif command == "/clear":
            if self._current_worker and not self._current_worker.is_finished:
                self.post_error(
                    "Cannot clear",
                    "Wait for the current request to finish before clearing the transcript.",
                )
                return
            self._stop_tailer()
            self.query_one(Transcript).clear_cells()
            self.query_one(FooterWidget).set_phase(Phase.IDLE)
            self.query_one(FooterWidget).set_hint(
                "Enter to submit • /help for commands"
            )
            self._user_requests.clear()
            self.post_agent_message("Transcript cleared.")
            self.focus_composer()
        elif command == "/sessions":
            self.app.push_screen(SessionsScreen(self.session_root))
        elif command == "/resume":
            if not argument:
                self.app.push_screen(SessionsScreen(self.session_root))
                return
            self.start_resume(argument)
        elif command == "/tools":
            self._run_inline_cli(["tools"], title="Tools")
        elif command == "/doctor":
            self._run_inline_cli(["doctor"], title="Doctor")
        else:
            self.post_error("Unknown command", raw)

    def _run_inline_cli(self, args: Iterable[str], *, title: str) -> None:
        result = CliRunner().invoke(agent, list(args), catch_exceptions=False)
        text = result.output.strip() or f"{title} completed."
        if result.exit_code == 0:
            self.post_agent_message(f"```\n{text}\n```", title=title)
        else:
            self.post_error(title, text)
