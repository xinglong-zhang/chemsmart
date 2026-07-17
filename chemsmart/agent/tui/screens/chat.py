"""Main chat screen for the chemsmart agent TUI."""

from __future__ import annotations

import os
import shlex
import uuid
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from threading import Event
from typing import Iterable

from click.testing import CliRunner
from rich.table import Table
from textual.app import ComposeResult
from textual.binding import Binding
from textual.screen import Screen
from textual.widgets import OptionList
from textual.worker import Worker

from chemsmart.agent.cli_commands import agent
from chemsmart.agent.core import CriticVerdict, DecisionLog, Plan
from chemsmart.agent.permissions import (
    ApprovalDecision,
    PermissionMode,
    PermissionPolicy,
)
from chemsmart.agent.model_command_parser import parse_model_command
from chemsmart.agent.provider_adapter import ToolRequest
from chemsmart.agent.harness.workflow_state import (
    current_workflow_state,
)
from chemsmart.agent.synthesis import SynthesisSession
from chemsmart.settings.workspace_project import (
    resolve_workspace_project,
    workspace_project_path,
)
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.chat_models import ReadyCommand
from chemsmart.agent.tui.slash_catalog import (
    SLASH_PALETTE_COMMANDS as _SLASH_PALETTE_COMMANDS,  # noqa: F401
)
from chemsmart.agent.tui.chat_helpers import (
    _command_details_text,
    _decision_trace_dict,
    _default_interaction_mode,
    _final_answer_text,
    _final_command_text,
    _find_project_yaml_candidate_for_write,
    _format_semantic_result,
    _latest_project_yaml_candidate as _latest_project_yaml_candidate,
    _load_tui_provider_config,
    _ready_command_hint,
)
from chemsmart.agent.tui.state import TuiState
from chemsmart.agent.tui.screens.jobs_panel import JobsPanel, JobsPanelAction
from chemsmart.agent.tui.screens.calculations import (
    CalculationMonitor,
)
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
from chemsmart.agent.tui.services.session_index import agent_session_dirs
from chemsmart.agent.tui.services.session_runner import SessionRunnerMixin
from chemsmart.agent.tui.mixins.decision_log import DecisionLogMixin
from chemsmart.agent.tui.mixins.calculations import CalculationMixin
from chemsmart.agent.tui.mixins.wizard_commands import WizardCommandsMixin
from chemsmart.agent.tui.mixins.worker_lifecycle import WorkerLifecycleMixin
from chemsmart.agent.tui.mixins.session_workers import SessionWorkersMixin
from chemsmart.agent.tui.mixins.request_flow import RequestFlowMixin
from chemsmart.agent.tui.mixins.workspace_interaction import (
    WorkspaceInteractionMixin,
)
from chemsmart.agent.tui.tool_meta import (
    tool_description,
)
from chemsmart.agent.tui.widgets.cells import (
    AgentMessageCell,
    CalculationReceiptCell,
    CommandInterpretationCell,
    CriticVerdictCell,
    DecisionTraceCell,
    FinalAnswerCell,
    JobStatusCell,
    MoleculeCell,
    PlanCell,
    RunResultCell,
    SynthesisTraceCell,
    ToolCallCell,
)
from chemsmart.agent.tui.widgets.cells.base import BaseCell
from chemsmart.agent.tui.widgets.composer import Composer
from chemsmart.agent.tui.widgets.calculation_strip import (
    CalculationStatusStrip,
)
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.header import ChemsmartHeader
from chemsmart.agent.tui.widgets.popups import (
    ApprovalResult,
    CwdMismatchChoice,
    CwdMismatchOverlay,
    FilePickerOverlay,
    PermissionModeOverlay,
    PermissionModeResult,
    ProjectWriteOverlay,
    ProjectWriteResult,
    ResponseCopyOverlay,
    TextPromptOverlay,
    build_approval_overlay,
)
from chemsmart.agent.tui.widgets.slash_palette import (
    SlashCommandPalette,
)
from chemsmart.agent.tui.widgets.transcript import Transcript
from chemsmart.io.molecules.structure import Molecule


class ChatScreen(
    CalculationMixin,
    DecisionLogMixin,
    WizardCommandsMixin,
    WorkerLifecycleMixin,
    SessionWorkersMixin,
    RequestFlowMixin,
    WorkspaceInteractionMixin,
    JobPollerMixin,
    SessionRunnerMixin,
    Screen,
):
    BINDINGS: list[Binding] = []

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
        runtime_v2: str = "active",
    ) -> None:
        super().__init__()
        self.session_root = session_root
        self.job_poll_interval = job_poll_interval
        self.runtime_v2 = runtime_v2
        self.tui_state = TuiState()
        self._tailer: LogTailer | None = None
        self._tailer_path: Path | None = None
        self._current_worker: Worker | None = None
        self._queued_prompt: str | None = None
        self._request_history: list[str] = []
        self._history_cursor = 0
        self._waiting_for_user = False
        self._quit_armed = False
        self._quit_timer = None
        self._session_poll_timer = None
        self._user_requests: set[str] = set()
        self._current_request: str | None = None
        self._current_plan: Plan | None = None
        self._current_plan_text: str | None = None
        self._current_verdict: CriticVerdict | None = None
        self._permission_mode = PermissionMode.DRIVING
        self._yolo_enabled = False
        self._session_allow_tools: set[str] = set()
        self._pending_approval = False
        self._pending_tool_request: ToolRequest | None = None
        self._pending_approval_description: str = ""
        self._pending_approval_args: dict = {}
        self._pending_approval_index: int | None = None
        self._pending_approval_total: int | None = None
        self._approval_waiter: Event | None = None
        self._approval_decision: ApprovalDecision | None = None
        self._pending_project_write_candidate: dict[str, object] | None = None
        self._latest_dry_run_content: str | None = None
        self._job_snapshot: dict[str, dict] = {}
        self._job_cells: dict[str, JobStatusCell] = {}
        self._rendered_run_results: set[str] = set()
        self._cancelled_job_ids: set[str] = set()
        self._active_server_name = self._default_active_server_name()
        self._latest_wizard_probe: dict[str, object] | None = None
        self._last_dry_run_session_id: str | None = None
        self._ready_command: ReadyCommand | None = None
        self._tool_cells: dict[str, ToolCallCell] = {}
        self._tool_order: list[str] = []
        self._direct_tool_call_ids: dict[str, str] = {}
        self._calculation_runs: dict[str, dict[str, object]] = {}
        self._calculation_cells: dict[str, CalculationReceiptCell] = {}
        self._turn_serial = 0
        self._active_turn_id = "turn-0"
        self._tui_session_dir = (
            self.session_root
            / ".runtime"
            / (
                "tui-"
                + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
                + "-"
                + uuid.uuid4().hex[:8]
            )
        )
        self._calculation_decision_log: DecisionLog | None = None
        self._active_provider_config = _load_tui_provider_config()
        self._interaction_mode = _default_interaction_mode(
            self._active_provider_config
        )
        self._build_mode = False
        selected_project = current_workflow_state().project
        self._selected_workspace_project_path = (
            Path(selected_project.path).resolve()
            if selected_project is not None
            else None
        )
        self._workspace_project_status = resolve_workspace_project(
            selected_path=self._selected_workspace_project_path
        )
        self.active_synthesis_session: SynthesisSession | None = None

    def compose(self) -> ComposeResult:
        yield ChemsmartHeader(id="chat-header")
        yield Transcript(id="chat-body")
        yield CalculationStatusStrip(id="calculation-strip")
        yield SlashCommandPalette(id="slash-palette")
        yield Composer()
        yield FooterWidget(id="status-footer", state=self.tui_state)

    def on_mount(self) -> None:
        self.post_agent_message(
            "## Welcome to the chemsmart agent\n\n"
            "Plan calculations in natural language, preview input files, "
            "and run pre-flight checks before execution. The same interface "
            "handles planning, validation, and approved execution.\n\n"
            "Start with: /doctor · Execute a validated local command: /run · "
            "Submit a validated HPC command: /submit\n\n"
            "Example: `single-point on examples/h2o.xyz at B3LYP/6-31G(d) "
            "Gaussian`"
        )
        self.focus_composer()
        footer = self.query_one(FooterWidget)
        self._sync_footer_provider()
        footer.set_server(self._active_server_name)
        footer.set_permission(
            self._permission_mode.value,
            yolo=self._yolo_enabled,
        )
        footer.update_draft("")
        config = getattr(self.app, "tui_config", None)
        if config is not None and config.issues:
            self.post_agent_message(
                "TUI configuration used safe defaults:\n\n"
                + "\n".join(f"- {issue}" for issue in config.issues),
                title="TUI configuration",
            )
        self.run_job_poller(self.job_poll_interval)
        self._refresh_job_snapshot()
        self._restore_calculation_runs()

    def on_composer_submitted(self, event: Composer.Submitted) -> None:
        text = event.text.strip()
        if text.startswith("/"):
            event.composer.clear_text()
            self.query_one(FooterWidget).update_draft("")
            self.query_one(SlashCommandPalette).hide()
            self._handle_slash_command(text)
            return
        if self._handle_plain_approval_alias(text):
            event.composer.clear_text()
            self.query_one(FooterWidget).update_draft("")
            return
        if self._worker_is_busy():
            event.composer.load_text(text)
            self.query_one(FooterWidget).set_hint(
                "Request in progress · Tab queues this draft"
            )
            return
        event.composer.clear_text()
        self.query_one(FooterWidget).update_draft("")
        self.query_one(SlashCommandPalette).hide()
        if not self._request_history or self._request_history[-1] != text:
            self._request_history.append(text)
        self._history_cursor = len(self._request_history)
        self._waiting_for_user = False
        self.start_request(text)

    def on_option_list_option_selected(
        self, event: OptionList.OptionSelected
    ) -> None:
        palette = self.query_one(SlashCommandPalette)
        if event.option_list is not palette or not event.option_id:
            return
        composer = self.query_one(Composer)
        composer.load_text(f"{event.option_id} ")
        palette.hide()
        composer.focus()

    def on_text_area_changed(self, event) -> None:
        if getattr(event.text_area, "id", None) == "composer":
            text = event.text_area.text
            self.query_one(FooterWidget).update_draft(text)
            self._update_slash_palette(text)

    def on_base_cell_copy_requested(
        self, event: BaseCell.CopyRequested
    ) -> None:
        if self.app.plain:
            self.app.copy_to_clipboard(event.text)
            self.notify("Response copied", timeout=1.5)
            return
        self.app.push_screen(
            ResponseCopyOverlay(event.text, title=event.title),
            lambda _result: self.focus_composer(),
        )

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

    def _can_approve_or_execute(self, action: str) -> bool:
        pending = self._pending_tool_request
        if self._pending_approval and pending is not None:
            if pending.name == "execute_chemsmart_command":
                pending_command = str(pending.arguments.get("command") or "")
                return parse_model_command(pending_command).action == action
            return pending.name == (
                "run_local" if action == "run" else "submit_hpc"
            )
        return bool(
            self._ready_command is not None
            and self._ready_command.action == action
        )

    def _remember_ready_command(
        self,
        *,
        command: str,
        semantic: dict[str, object] | None,
        intent: dict[str, object] | None,
        source: str,
    ) -> bool:
        self._ready_command = None
        parsed = parse_model_command(command)
        if parsed.parse_error or parsed.action not in {"run", "sub"}:
            return False
        if not isinstance(semantic, dict):
            return False
        semantic_verdict = str(semantic.get("verdict") or "").lower()
        if semantic_verdict not in {"ok", "warn"}:
            return False
        generated_inputs = semantic.get("generated_inputs")
        if not isinstance(generated_inputs, list) or not generated_inputs:
            return False
        if not isinstance(intent, dict):
            return False
        intent_verdict = str(intent.get("verdict") or "").lower()
        if intent_verdict not in {"ok", "warn"}:
            return False

        project_path: Path | None = None
        project_hash: str | None = None
        if parsed.program in {"gaussian", "orca"}:
            if not parsed.project:
                return False
            project_path = workspace_project_path(
                parsed.project,
                parsed.program,
                cwd=Path.cwd(),
            )
            if not project_path.is_file():
                return False
            project_hash = sha256(project_path.read_bytes()).hexdigest()

        self._ready_command = ReadyCommand(
            command=command,
            action=parsed.action,
            workspace=Path.cwd().resolve(),
            semantic_verdict=semantic_verdict,
            intent_verdict=intent_verdict,
            project_path=project_path,
            project_sha256=project_hash,
            source=source,
        )
        return True

    def _ready_command_problem(self, expected_action: str) -> str | None:
        ready = self._ready_command
        if ready is None:
            return "No semantic- and intent-validated command is ready."
        if ready.action != expected_action:
            expected = (
                "chemsmart run"
                if expected_action == "run"
                else "chemsmart sub"
            )
            return f"The validated command is not a `{expected}` command."
        if ready.workspace != Path.cwd().resolve():
            return (
                "The workspace changed after validation. Regenerate the command "
                "from the current workspace before execution."
            )
        if ready.project_path is not None:
            if not ready.project_path.is_file():
                return "The validated workspace project YAML no longer exists."
            current_hash = sha256(ready.project_path.read_bytes()).hexdigest()
            if current_hash != ready.project_sha256:
                return (
                    "The workspace project YAML changed after validation. "
                    "Regenerate the command so the generated-input evidence "
                    "matches the current project."
                )
        return None

    def _publish_synthesis_result(self, result: dict[str, object]) -> None:
        synthesis = result.get("synthesis")
        if not isinstance(synthesis, dict):
            self.post_error("Synthesis failed", "Provider returned no result.")
            self.query_one(FooterWidget).set_phase(Phase.ERROR)
            self.query_one(FooterWidget).set_hint("Synthesis failed")
            return

        status = str(synthesis.get("status") or "")
        self._waiting_for_user = status == "needs_clarification"
        footer = self.query_one(FooterWidget)
        provider_type = str(result.get("provider_type") or "offline")
        provider_model = str(result.get("provider_model") or "auto")
        artifact_dir = str(result.get("synthesis_artifact_dir") or "")
        semantic = result.get("semantic_result")
        semantic_dict = semantic if isinstance(semantic, dict) else None
        if status == "ready":
            command = str(synthesis.get("command") or "")
            explanation = str(synthesis.get("explanation") or "")
            confidence = str(synthesis.get("confidence") or "low")
            project = str(synthesis.get("project") or "")
            intent = synthesis.get("intent_assertion") or synthesis.get(
                "intent"
            )
            intent_dict = intent if isinstance(intent, dict) else None
            command_is_executable = self._remember_ready_command(
                command=command,
                semantic=semantic_dict,
                intent=intent_dict,
                source="local_synthesis",
            )
            transcript = self.query_one(Transcript)
            transcript.add_cell(
                SynthesisTraceCell(
                    provider_type=provider_type,
                    model=provider_model,
                    mode=self._interaction_mode,
                    status=status,
                    command=command,
                    semantic=semantic_dict,
                    decision_trace=_decision_trace_dict(synthesis),
                    artifact_dir=artifact_dir,
                )
            )
            self._publish_decision_trace(synthesis)
            transcript.add_cell(
                CommandInterpretationCell(
                    parse_model_command(command),
                    expanded=True,
                )
            )
            transcript.add_cell(
                AgentMessageCell(
                    _command_details_text(
                        explanation=(
                            explanation or "Prepared chemsmart command."
                        ),
                        confidence=confidence,
                        project=project,
                    ),
                    title="Command details",
                )
            )
            transcript.add_cell(
                FinalAnswerCell(
                    _final_command_text(command=command),
                    title="Final Command",
                ),
            )
            footer.set_phase(Phase.FINISHED)
            footer.set_hint(
                _ready_command_hint(self._ready_command)
                if command_is_executable
                else "Command shown, but execution evidence is incomplete"
            )
            transcript.collapse_tool_chain(self._active_turn_id)
            return

        if status == "informational":
            command = str(synthesis.get("command") or "")
            explanation = str(synthesis.get("explanation") or "")
            action = str(synthesis.get("action") or "explain_command")
            title = {
                "explain_command": "Command Explanation",
                "critique_command": "Command Critic",
                "repair_command": "Command Repair",
            }.get(action, "Command Analysis")
            transcript = self.query_one(Transcript)
            transcript.add_cell(
                SynthesisTraceCell(
                    provider_type=provider_type,
                    model=provider_model,
                    mode=self._interaction_mode,
                    status=status,
                    command=command,
                    semantic=semantic_dict,
                    decision_trace=_decision_trace_dict(synthesis),
                    artifact_dir=artifact_dir,
                )
            )
            self._publish_decision_trace(synthesis)
            if command:
                transcript.add_cell(
                    CommandInterpretationCell(
                        parse_model_command(command),
                        expanded=True,
                    )
                )
            transcript.add_cell(
                FinalAnswerCell(
                    _final_answer_text(
                        command=command,
                        explanation=(
                            explanation or "No explanation was generated."
                        ),
                    ),
                    title=title,
                )
            )
            footer.set_phase(Phase.FINISHED)
            footer.set_hint("Command analysis ready")
            transcript.collapse_tool_chain(self._active_turn_id)
            return

        if status == "needs_clarification":
            missing = synthesis.get("missing_info") or []
            if not isinstance(missing, list):
                missing = [str(missing)]
            lines = "\n".join(f"- {item}" for item in missing) or "- details"
            self.query_one(Transcript).add_cell(
                SynthesisTraceCell(
                    provider_type=provider_type,
                    model=provider_model,
                    mode=self._interaction_mode,
                    status=status,
                    command=str(synthesis.get("command") or ""),
                    semantic=semantic_dict,
                    decision_trace=_decision_trace_dict(synthesis),
                    artifact_dir=artifact_dir,
                )
            )
            self._publish_decision_trace(synthesis)
            self.query_one(Transcript).add_cell(
                FinalAnswerCell(
                    (
                        "I need more information before making a command:\n\n"
                        f"{lines}"
                    ),
                    title="Clarification",
                )
            )
            footer.set_phase(Phase.WAITING_USER)
            footer.set_hint("Clarification needed")
            self.query_one(Transcript).collapse_tool_chain(
                self._active_turn_id
            )
            return

        explanation = str(
            synthesis.get("explanation") or "No executable command was made."
        )
        semantic_text = _format_semantic_result(semantic_dict)
        self.query_one(Transcript).add_cell(
            SynthesisTraceCell(
                provider_type=provider_type,
                model=provider_model,
                mode=self._interaction_mode,
                status=status,
                command=str(synthesis.get("command") or ""),
                semantic=semantic_dict,
                decision_trace=_decision_trace_dict(synthesis),
                artifact_dir=artifact_dir,
            )
        )
        self._publish_decision_trace(synthesis)
        self.query_one(Transcript).add_cell(
            FinalAnswerCell(
                f"{explanation}{semantic_text}",
                title="Final Status",
            )
        )
        footer.set_phase(Phase.FINISHED)
        footer.set_hint("No command generated")
        self.query_one(Transcript).collapse_tool_chain(self._active_turn_id)

    def _publish_decision_trace(self, synthesis: dict[str, object]) -> None:
        trace = synthesis.get("decision_trace")
        if isinstance(trace, dict) and trace:
            self.query_one(Transcript).add_cell(DecisionTraceCell(trace))

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
        elif command in {"/jobs", "/runs"}:
            self.action_show_calculations()
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
        elif command == "/init":
            self._handle_init_command(argument)
        elif command == "/write-project":
            self._handle_project_write_command(argument)
        elif command == "/wizard":
            self._handle_wizard_probe_command(argument)
        elif command == "/wizard-verify":
            self._handle_wizard_verify_command(argument)
        elif command == "/wizard-refresh":
            self._handle_wizard_refresh_command(argument)
        elif command == "/wizard-write":
            self._handle_wizard_write_command(argument)
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
        elif command == "/allow":
            self._resolve_pending_approval(ApprovalDecision.ALLOW_ONCE)
        elif command == "/allow-session":
            self._resolve_pending_approval(ApprovalDecision.ALLOW_SESSION)
        elif command == "/deny":
            self._resolve_pending_approval(ApprovalDecision.DENY)
        elif command == "/permissions":
            if argument:
                if not self._apply_permission_command(argument):
                    self.post_error(
                        "Unknown permissions command",
                        "Use /permissions permission|driving.",
                    )
                return
            if self.app.plain:
                self.post_agent_message(
                    (
                        "```\n"
                        f"mode: {self._permission_mode.value}\n"
                        f"yolo: {'on' if self._yolo_enabled else 'off'}\n"
                        "Use /permissions permission|driving and /yolo on|off.\n"
                        "```"
                    ),
                    title="Permissions",
                )
                return
            self.app.push_screen(
                PermissionModeOverlay(
                    mode=self._permission_mode,
                    yolo=self._yolo_enabled,
                ),
                self._handle_permission_mode_result,
            )
        elif command == "/yolo":
            if not argument:
                self.post_error(
                    "Missing value",
                    "Usage: /yolo on|off",
                )
                return
            self._set_yolo(argument)
        elif command == "/execute":
            self._handle_execute_command(argument)
        elif command == "/submit":
            self._handle_ready_command_execution("sub", argument)
        elif command == "/run":
            self._handle_ready_command_execution("run", argument)
        else:
            self.post_error("Unknown command", raw)

    def _handle_project_write_command(self, argument: str) -> None:
        if self._current_worker and not self._current_worker.is_finished:
            self.post_error(
                "Session already running",
                "Wait for the current request to finish before writing project YAML.",
            )
            return
        try:
            parts = shlex.split(argument)
        except ValueError as exc:
            self.post_error("Invalid write-project command", str(exc))
            return

        confirmed = False
        overwrite = False
        project_name: str | None = None
        for part in parts:
            value = part.strip()
            lowered = value.lower()
            if lowered in {"yes", "y", "--yes"}:
                confirmed = True
            elif lowered in {"overwrite", "--overwrite"}:
                overwrite = True
            elif project_name is None:
                project_name = value
            else:
                self.post_error(
                    "Invalid write-project command",
                    "Usage: /write-project [name] (confirmation opens in TUI)",
                )
                return

        candidate = _find_project_yaml_candidate_for_write(
            self.session_root,
            preferred_session_dir=self.current_session_dir(),
        )
        if candidate is None:
            self.post_error(
                "Project write unavailable",
                "Run /init and build a validated project YAML first.",
            )
            return
        if project_name:
            candidate["project_name"] = project_name

        if confirmed:
            candidate["overwrite"] = overwrite
            self._start_project_yaml_write(candidate)
            return
        if self.app.plain:
            self.post_error(
                "Confirmation required",
                "Plain mode uses /write-project [name] yes [overwrite].",
            )
            return
        self._prompt_project_yaml_write(
            candidate,
            prefer_active=project_name is None,
        )

    def _prompt_project_yaml_write(
        self,
        candidate: dict[str, object],
        *,
        prefer_active: bool,
    ) -> None:
        program = str(candidate.get("program") or "gaussian").strip().lower()
        requested = Path(str(candidate.get("project_name") or "project")).stem
        status = self._resolve_workspace_project_status()
        overwrite_project = requested
        if (
            prefer_active
            and status.loaded
            and status.program == program
            and status.project
        ):
            overwrite_project = status.project
        target = workspace_project_path(overwrite_project, program)
        new_project = self._next_available_project_name(requested, program)
        self._pending_project_write_candidate = dict(candidate)
        self.app.push_screen(
            ProjectWriteOverlay(
                target=target,
                overwrite_project=overwrite_project,
                new_project=new_project,
                target_exists=target.exists(),
            ),
            self._handle_project_write_result,
        )

    def _handle_project_write_result(
        self, result: ProjectWriteResult | None
    ) -> None:
        candidate = self._pending_project_write_candidate
        self._pending_project_write_candidate = None
        if result is None or candidate is None:
            self.query_one(FooterWidget).set_hint(
                "Project YAML write cancelled"
            )
            self.focus_composer()
            return
        candidate["project_name"] = result.project_name
        candidate["overwrite"] = result.action == "overwrite"
        self._start_project_yaml_write(candidate)

    def _start_project_yaml_write(self, candidate: dict[str, object]) -> None:
        footer = self.query_one(FooterWidget)
        footer.set_phase(Phase.PLANNING)
        footer.set_hint("Writing project YAML…")
        self._current_worker = self.run_project_yaml_write(candidate)

    def _next_available_project_name(
        self, base_name: str, program: str
    ) -> str:
        base = Path(str(base_name or "project")).stem or "project"
        if not workspace_project_path(base, program).exists():
            return base
        suffix = 2
        while workspace_project_path(f"{base}-{suffix}", program).exists():
            suffix += 1
        return f"{base}-{suffix}"

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

    def _show_help(self) -> None:
        table = Table(show_header=True, box=None, padding=(0, 1))
        table.add_column("Phase", style="dim", no_wrap=True)
        table.add_column("Command", style="bold")
        table.add_column("Description")
        rows = [
            ("[A]", "/help", "show this help"),
            ("[A]", "/jobs", "open the jobs panel"),
            ("[A]", "/runs · Ctrl+B", "monitor calculations and logs"),
            ("[A]", "/queue", "show the current queue snapshot"),
            ("[A]", "/server <name>", "switch the active HPC server"),
            ("[A]", "/molecule <path>", "load and preview a molecule"),
            ("[A]", "/cancel <job-id> [yes]", "cancel a queued/running job"),
            ("[F]", "/extract <job-id|inputfile>", "parse a final result"),
            ("[D]", "/dryrun", "regenerate the current dry-run"),
            ("[A]", "/allow", "allow the focused tool request once"),
            (
                "[A]",
                "/allow-session",
                "allow the focused tool for the rest of this session",
            ),
            ("[A]", "/deny", "deny the focused tool request"),
            ("[A]", "/permissions", "toggle permission/driving mode"),
            ("[A]", "/yolo on|off", "toggle risky-tool autonomy"),
            ("[F]", "/execute [yes]", "submit the dry-run to HPC for real"),
            ("[F]", "/submit [yes]", "submit validated chemsmart sub command"),
            ("[F]", "/run [yes]", "execute validated chemsmart run command"),
            ("[P,D]", "/critic", "show the current critic verdict"),
            ("[P,D,R]", "/plan", "show the current plan"),
            ("[A]", "/rationale", "show planner rationale"),
            ("[I,F]", "/clear", "clear the transcript"),
            ("[I]", "/sessions", "browse recent sessions"),
            ("[I]", "/resume <session-id>", "load or continue a session"),
            ("[A]", "/tools", "list registered tools"),
            ("[A]", "/wizard-refresh <name> [--force]", "refresh node cache"),
            ("[A]", "/wizard-verify <name>", "verify server transport wiring"),
            ("[A]", "/doctor", "run inline diagnostics"),
            (
                "[F]",
                "/write-project [name]",
                "write or replace the latest validated project YAML",
            ),
            ("[A]", "/quit", "exit the TUI"),
        ]
        for row in rows:
            table.add_row(*row)
        self.query_one(Transcript).add_cell(
            AgentMessageCell(table, title="Help")
        )

    def _handle_init_command(self, argument: str) -> None:
        if (
            self._active_provider_config is not None
            and self._active_provider_config.type == "local"
        ):
            self.post_error(
                "Init unavailable for local provider",
                "Project YAML build mode uses the tool-loop harness, which "
                "needs a tool-calling provider (anthropic/openai).",
            )
            return
        if self._current_worker and not self._current_worker.is_finished:
            self.post_error(
                "Session already running",
                "Wait for the current request to finish before starting "
                "project YAML build mode.",
            )
            return
        argument = argument.strip()
        self._build_mode = False
        self._sync_footer_provider()
        if argument:
            self.start_request(
                "Build a workspace project YAML from this reported method. "
                "Validate and critique it before writing; ask for approval "
                f"before write_project_yaml.\n\n{argument}"
            )
            return
        self.post_agent_message(
            (
                "**Project YAML request helper**\n\n"
                "Paste your reported computational method (a paper's "
                "*Computational Details*, a method sentence, or a few facts) "
                "as a normal message, and the unified agent will use the "
                "project-YAML harness to "
                "`extract → render → validate → critique` a chemsmart project "
                "YAML, then `write` it into this workspace's "
                "`.chemsmart/<program>/` folder once you approve.\n\n"
                'Tip: name the project (e.g. "call it co2") and say the '
                "program (Gaussian or ORCA)."
            ),
            title="Init: Project YAML",
        )
        self.query_one(FooterWidget).set_hint("Paste a project-YAML request")
        self.focus_composer()

    def _show_sessions_snapshot(self) -> None:
        lines = ["Sessions", "", "Use /resume <session-id> to load one."]
        if self.session_root.exists():
            session_dirs = agent_session_dirs(self.session_root)
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
            self._resolve_pending_approval(ApprovalDecision.ALLOW_ONCE)
            return True
        if keyword in {"n", "no"}:
            self._resolve_pending_approval(ApprovalDecision.DENY)
            return True
        if keyword in {"s", "session"}:
            self._resolve_pending_approval(ApprovalDecision.ALLOW_SESSION)
            return True
        if keyword in {"r", "revise"}:
            if not corrective_text:
                self.post_error(
                    "Missing instruction",
                    "Usage: /run revise <instruction> or /submit revise <instruction>.",
                )
                return True
            self.start_request(self._corrected_request(corrective_text))
            return True
        self.post_error(
            "Unknown approval response",
            "Use yes, no, session, or revise <instruction>.",
        )
        return True

    def _run_inline_cli(self, args: Iterable[str], *, title: str) -> None:
        import threading

        from chemsmart.agent.services.cli_presenters import (
            sanitize_inline_cli_output,
        )

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

    def _handle_plain_approval_alias(self, text: str) -> bool:
        keyword = text.strip().lower()
        if keyword in {"yes", "y"}:
            self._resolve_pending_approval(ApprovalDecision.ALLOW_ONCE)
            return True
        if keyword in {"session", "s"}:
            self._resolve_pending_approval(ApprovalDecision.ALLOW_SESSION)
            return True
        if keyword in {"no", "n"}:
            self._resolve_pending_approval(ApprovalDecision.DENY)
            return True
        return False

    def _request_approval(self, expected_tool: str) -> None:
        if not self._pending_approval:
            self.post_error(
                "Approval unavailable",
                f"No pending `{expected_tool}` action is waiting for approval.",
            )
            return
        pending = self._pending_tool_request
        if pending is None or pending.name != expected_tool:
            self.post_error(
                "Approval unavailable",
                f"No pending `{expected_tool}` action is waiting for approval.",
            )
            return
        if self.app.plain:
            self.query_one(FooterWidget).set_hint(
                "PERMISSION_PENDING "
                f"({self._pending_approval_index or 1} of "
                f"{self._pending_approval_total or 1})"
            )
            return
        if pending.name == "write_project_yaml":
            program = (
                str(pending.arguments.get("program") or "gaussian")
                .strip()
                .lower()
            )
            project = Path(
                str(pending.arguments.get("project_name") or "project")
            ).stem
            target = workspace_project_path(project, program)
            if target.exists():
                self.app.push_screen(
                    ProjectWriteOverlay(
                        target=target,
                        overwrite_project=project,
                        new_project=self._next_available_project_name(
                            project, program
                        ),
                        target_exists=True,
                    ),
                    self._handle_agent_project_write_result,
                )
                return
        self.app.push_screen(
            build_approval_overlay(
                active_mode=self._permission_mode.value.upper(),
                tool_name=pending.name,
                description=(
                    self._pending_approval_description
                    or tool_description(pending.name)
                ),
                arguments=self._pending_approval_args or pending.arguments,
                session_rule_active=(
                    pending.name in self._session_allow_tools
                ),
                queue_index=self._pending_approval_index,
                queue_total=self._pending_approval_total,
            ),
            self._handle_approval_result,
        )

    def _handle_agent_project_write_result(
        self, result: ProjectWriteResult | None
    ) -> None:
        pending = self._pending_tool_request
        if (
            result is None
            or pending is None
            or pending.name != "write_project_yaml"
        ):
            self._resolve_pending_approval(ApprovalDecision.DENY)
            return
        pending.arguments["project_name"] = result.project_name
        pending.arguments["overwrite"] = result.action == "overwrite"
        self._pending_approval_args = dict(pending.arguments)
        self._resolve_pending_approval(ApprovalDecision.ALLOW_ONCE)

    def _handle_approval_result(self, result: ApprovalResult | None) -> None:
        if result is None:
            self._resolve_pending_approval(ApprovalDecision.DENY)
            return
        if result.choice == "r":
            self.app.push_screen(
                TextPromptOverlay(
                    title="Revise request",
                    prompt=(
                        "Describe the correction. The pending tool will be denied "
                        "and the corrected request will run next."
                    ),
                ),
                self._handle_approval_revision,
            )
            return
        decision = result.to_decision()
        if decision is None:
            return
        self._resolve_pending_approval(decision)

    def _handle_approval_revision(self, value: str | None) -> None:
        correction = (value or "").strip()
        if not correction:
            self._request_approval(
                self._pending_tool_request.name
                if self._pending_tool_request is not None
                else ""
            )
            return
        self._queued_prompt = self._corrected_request(correction)
        footer = self.query_one(FooterWidget)
        footer.set_queued_prompt(True)
        footer.set_hint("Correction queued; pending tool denied")
        self._resolve_pending_approval(ApprovalDecision.DENY)

    def _corrected_request(self, corrective_text: str) -> str:
        original = self._current_request or ""
        return (
            f"Corrective instruction: {corrective_text}\n\n"
            f"Original request:\n{original}"
        )

    def _await_approval(self, request: ToolRequest) -> ApprovalDecision:
        decision_ready = Event()

        def prompt() -> None:
            self._approval_waiter = decision_ready
            self._approval_decision = None
            self._pending_approval = True
            self._pending_tool_request = request
            footer = self.query_one(FooterWidget)
            footer.set_phase(Phase.APPROVAL_REQUIRED)
            footer.set_hint(
                "PERMISSION_PENDING "
                f"({self._pending_approval_index or 1} of "
                f"{self._pending_approval_total or 1})"
            )
            if not self.app.plain:
                self._request_approval(request.name)

        self.app.call_from_thread(prompt)
        decision_ready.wait()
        return self._approval_decision or ApprovalDecision.DENY

    def _resolve_pending_approval(
        self,
        decision: ApprovalDecision,
    ) -> bool:
        if not self._pending_approval or self._approval_waiter is None:
            self.post_error(
                "Approval unavailable",
                "No tool request is currently waiting for approval.",
            )
            return False
        self._approval_decision = decision
        self._pending_approval = False
        self._pending_tool_request = None
        self._pending_approval_description = ""
        self._pending_approval_args = {}
        self._pending_approval_index = None
        self._pending_approval_total = None
        waiter = self._approval_waiter
        self._approval_waiter = None
        self.query_one(FooterWidget).set_phase(Phase.TOOL_RUNNING)
        self.query_one(FooterWidget).set_hint(
            "DENIED_CONTINUING"
            if decision == ApprovalDecision.DENY
            else (
                "DRIVING_AUTO"
                if self._permission_mode == PermissionMode.DRIVING
                else "Approval recorded"
            )
        )
        waiter.set()
        return True

    def _permission_policy(
        self, *, prompt_risky: bool = False
    ) -> PermissionPolicy:
        return PermissionPolicy(
            mode=self._permission_mode,
            yolo=self._yolo_enabled,
            prompt_risky=prompt_risky,
            session_allow=set(self._session_allow_tools),
        )

    def _handle_permission_mode_result(
        self,
        result: PermissionModeResult | None,
    ) -> None:
        if result is None:
            return
        self._permission_mode = result.mode
        self._yolo_enabled = result.yolo
        self.query_one(FooterWidget).set_permission(
            self._permission_mode.value,
            yolo=self._yolo_enabled,
        )
        self.post_agent_message(
            (
                f"Permissions set to `{self._permission_mode.value}` "
                f"(yolo={'on' if self._yolo_enabled else 'off'})."
            ),
            title="Permissions",
        )

    def _apply_permission_command(self, argument: str) -> bool:
        value = argument.strip().lower()
        if value == "permission":
            self._permission_mode = PermissionMode.PERMISSION
            self.query_one(FooterWidget).set_permission(
                self._permission_mode.value,
                yolo=self._yolo_enabled,
            )
            return True
        if value == "driving":
            self._permission_mode = PermissionMode.DRIVING
            self.query_one(FooterWidget).set_permission(
                self._permission_mode.value,
                yolo=self._yolo_enabled,
            )
            return True
        return False

    def _set_yolo(self, argument: str) -> None:
        value = argument.strip().lower()
        if value not in {"on", "off"}:
            self.post_error("Unknown value", "Usage: /yolo on|off")
            return
        self._yolo_enabled = value == "on"
        self.query_one(FooterWidget).set_permission(
            self._permission_mode.value,
            yolo=self._yolo_enabled,
        )
        self.post_agent_message(
            f"YOLO {'enabled' if self._yolo_enabled else 'disabled'}.",
            title="Permissions",
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
            self.active_synthesis_session = None
            self._session_allow_tools.clear()
            self._last_dry_run_session_id = None
        self._current_request = None
        self._current_plan = None
        self._current_plan_text = None
        self._current_verdict = None
        self._workflow_cell = None
        self._ready_command = None
        self._pending_approval = False
        self._pending_tool_request = None
        self._pending_approval_description = ""
        self._pending_approval_args = {}
        self._pending_approval_index = None
        self._pending_approval_total = None
        self._approval_waiter = None
        self._approval_decision = None
        self._latest_dry_run_content = None
        self._tool_cells.clear()
        self._tool_order.clear()
        self._direct_tool_call_ids.clear()
        self._job_cells.clear()
        self._rendered_run_results.clear()
        self._latest_wizard_probe = None
        footer = self.query_one(FooterWidget)
        footer.set_phase(Phase.IDLE)
        footer.set_hint("Enter to submit • /help for commands")
        self._sync_footer_provider()
        footer.entity_status = None
        footer.reset_job_counts()
        if clear_transcript:
            transcript = self.query_one(Transcript)
            if keep_conversational:
                transcript.clear_turn_chrome()
            else:
                transcript.clear_cells()
                self._calculation_cells.clear()
        self._refresh_job_snapshot()

    def _emit_job_update(self, job_id: str, fields: dict) -> None:
        self.post_message(JobStatusUpdated(job_id, fields))
        if isinstance(self.app.screen, JobsPanel):
            self.app.screen.post_message(JobStatusUpdated(job_id, fields))
        elif isinstance(self.app.screen, CalculationMonitor):
            self.app.screen.update_job(job_id, fields)

    def _refresh_job_snapshot(self) -> None:
        from chemsmart.agent.tui.services import (
            job_poller as job_poller_service,
        )

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
        for run in self._calculation_runs.values():
            status = str(run.get("status") or "")
            if status in {"validating", "starting", "running"}:
                counts["running"] += 1
            elif status in {
                "chemistry_failed",
                "process_failed",
                "timeout",
            }:
                counts["failed"] += 1
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
        self.query_one(FooterWidget).set_server(name)
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
