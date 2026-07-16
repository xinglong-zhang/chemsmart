"""Main chat screen for the chemsmart agent TUI."""

from __future__ import annotations

import json
import os
import shlex
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from threading import Event, get_ident
from typing import Iterable

from click.testing import CliRunner
from rich.table import Table
from textual import work
from textual.app import ComposeResult
from textual.binding import Binding
from textual.screen import Screen
from textual.widgets import OptionList
from textual.worker import Worker, WorkerState

from chemsmart.agent.cli import agent
from chemsmart.agent.core import AgentSession, CriticVerdict, DecisionLog, Plan
from chemsmart.agent.permissions import (
    ApprovalDecision,
    PermissionMode,
    PermissionPolicy,
    ResolvedDecision,
)
from chemsmart.agent.provider_config import (
    AgentProviderConfig,
    AgentProviderConfigError,
    load_active_provider_config,
)
from chemsmart.agent.model_command_parser import parse_model_command
from chemsmart.agent.provider_adapter import ToolRequest
from chemsmart.agent.registry import ToolRegistry
from chemsmart.agent.runtime.calculations import (
    CalculationContext,
    CalculationEvent,
    cancel_calculation,
    inspect_calculation,
    load_calculation_runs,
)
from chemsmart.agent.harness.workflow_state import (
    current_workflow_state,
    select_workspace_project,
)
from chemsmart.agent.services.conversation_memory import ConversationMemory
from chemsmart.agent.synthesis import SynthesisSession, resolve_default_project
from chemsmart.agent.tools_command import execute_chemsmart_command_observed
from chemsmart.settings.workspace_project import (
    WorkspaceProjectStatus,
    resolve_workspace_project,
    workspace_project_path,
)
from chemsmart.agent.tui.events import (
    AssistantTurnEvent,
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
    ToolUseEvent,
    parse_decision_event,
    session_completed,
)
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.state import TuiState
from chemsmart.agent.tui.screens.jobs_panel import JobsPanel, JobsPanelAction
from chemsmart.agent.tui.screens.calculations import (
    CalculationMonitor,
    CalculationMonitorAction,
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
from chemsmart.agent.tui.tool_meta import (
    format_assumptions_banner,
    tool_description,
)
from chemsmart.agent.tui.widgets.cells import (
    AgentMessageCell,
    CalculationReceiptCell,
    CommandInterpretationCell,
    CriticVerdictCell,
    DecisionTraceCell,
    DryRunInputCell,
    ErrorCell,
    FinalAnswerCell,
    GeometryHandoffCell,
    JobStatusCell,
    MethodCell,
    MoleculeCell,
    PlanCell,
    RunResultCell,
    RuntimeValidationCell,
    SubmissionPreviewCell,
    SynthesisTraceCell,
    ToolCallCell,
    UserMessageCell,
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
    HistorySearchOverlay,
    PermissionModeOverlay,
    PermissionModeResult,
    ProjectWriteOverlay,
    ProjectWriteResult,
    ProjectYamlOverlay,
    ResponseCopyOverlay,
    ShortcutOverlay,
    TextPromptOverlay,
    ToolActivityOverlay,
    build_approval_overlay,
)
from chemsmart.agent.tui.widgets.slash_palette import (
    SlashCommandPalette,
    SlashPaletteItem,
)
from chemsmart.agent.tui.widgets.transcript import Transcript
from chemsmart.io.molecules.structure import Molecule

_SLASH_PALETTE_COMMANDS: tuple[tuple[str, str], ...] = (
    ("/help", "show available commands"),
    ("/jobs", "open the jobs panel"),
    ("/runs", "open the calculation monitor"),
    ("/queue", "show the current queue snapshot"),
    ("/server", "switch the active HPC server"),
    ("/molecule", "load and preview a molecule"),
    ("/cancel", "cancel a queued or running job"),
    ("/extract", "parse a final result"),
    ("/dryrun", "regenerate the current dry-run"),
    ("/allow", "allow the focused tool request once"),
    ("/allow-session", "allow the focused tool for this session"),
    ("/deny", "deny the focused tool request"),
    ("/permissions", "toggle permission or driving mode"),
    ("/yolo", "toggle risky-tool autonomy"),
    ("/execute", "submit the dry-run to HPC for real"),
    ("/submit", "submit the validated chemsmart sub command"),
    ("/run", "execute the validated chemsmart run command"),
    ("/critic", "show the current critic verdict"),
    ("/plan", "show the current plan"),
    ("/rationale", "show planner rationale"),
    ("/clear", "clear the transcript"),
    ("/sessions", "browse recent sessions"),
    ("/resume", "load or continue a session"),
    ("/tools", "list registered tools"),
    ("/wizard", "probe server setup and render YAML"),
    ("/wizard-refresh", "refresh node cache"),
    ("/wizard-verify", "verify server transport wiring"),
    ("/wizard-write", "write latest wizard YAML"),
    ("/doctor", "run inline diagnostics"),
    ("/mode", "show unified provider routing information"),
    ("/init", "start a project YAML request"),
    ("/write-project", "write or replace the latest validated project YAML"),
    ("/quit", "exit the TUI"),
    ("/exit", "exit the TUI"),
)


@dataclass(frozen=True)
class _ReadyCommand:
    command: str
    action: str
    workspace: Path
    semantic_verdict: str
    intent_verdict: str
    project_path: Path | None = None
    project_sha256: str | None = None
    source: str = ""


class ChatScreen(JobPollerMixin, SessionRunnerMixin, Screen):
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
        self._ready_command: _ReadyCommand | None = None
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

    @work(
        thread=True,
        exclusive=True,
        exit_on_error=False,
        group="agent-session",
        name="agent-unified",
    )
    def run_unified_session(self, request: str) -> dict[str, object]:
        self.active_resume_id = None
        if (
            self.active_agent_session is None
            or self.active_agent_session.stage_prompt != "unified_agent.md"
        ):
            self.active_agent_session = AgentSession(
                session_root=str(self.session_root),
                stage_prompt="unified_agent.md",
                runtime_v2=self.runtime_v2,
            )
        policy = self._permission_policy(prompt_risky=True)
        result = self.active_agent_session.run_loop(
            request,
            policy=policy,
            approver=lambda req: self._await_approval(req),
        )
        result["session_allow_tools"] = sorted(policy.session_allow)
        return result

    @work(
        thread=True,
        exclusive=True,
        exit_on_error=False,
        group="agent-session",
        name="agent-ask",
    )
    def run_synthesis_session(self, request: str) -> dict[str, object]:
        provider_type = _provider_type_label(self._active_provider_config)
        provider_model = _provider_model_label(self._active_provider_config)
        artifact_dir = _new_synthesis_artifact_dir(
            self.session_root,
            provider_type=provider_type,
        )
        try:
            if self.active_synthesis_session is None:
                self.active_synthesis_session = SynthesisSession()
            result = self.active_synthesis_session.prepare_command(request)
            semantic = self.active_synthesis_session._last_semantic_result
            payload = {
                "request": request,
                "provider_type": provider_type,
                "provider_model": provider_model,
                "synthesis": result,
                "raw_response": self.active_synthesis_session._last_raw_response,
                "semantic_result": (
                    semantic.to_dict() if semantic is not None else None
                ),
            }
            _write_synthesis_artifact(artifact_dir, payload)
        except Exception as exc:
            payload = {
                "synthesis": {
                    "status": "infeasible",
                    "command": "",
                    "explanation": _format_synthesis_exception(exc),
                    "confidence": "low",
                    "missing_info": [],
                    "alternatives": [],
                },
                "raw_response": "",
                "semantic_result": None,
            }
            _write_synthesis_artifact(
                artifact_dir,
                {
                    "request": request,
                    "provider_type": provider_type,
                    "provider_model": provider_model,
                    **payload,
                },
            )
            return {
                **payload,
                "provider_type": provider_type,
                "provider_model": provider_model,
                "synthesis_artifact_dir": str(artifact_dir),
            }
        return {
            "synthesis": result,
            "raw_response": self.active_synthesis_session._last_raw_response,
            "semantic_result": (
                semantic.to_dict() if semantic is not None else None
            ),
            "provider_type": provider_type,
            "provider_model": provider_model,
            "synthesis_artifact_dir": str(artifact_dir),
        }

    @work(
        thread=True,
        exclusive=True,
        exit_on_error=False,
        group="slash-tool",
        name="project-yaml-write",
    )
    def run_project_yaml_write(
        self,
        arguments: dict[str, object],
    ) -> dict[str, object]:
        tool_name = "write_project_yaml"
        registry = ToolRegistry.default()
        normalized_args = registry.normalize_args(tool_name, arguments)
        description = registry.describe_tool(tool_name)
        self.app.call_from_thread(
            self._publish_tool_call_cell,
            tool_name,
            "approved",
            description,
            normalized_args,
            "Confirmed in the workspace project write dialog.",
        )
        result = registry.call(tool_name, normalized_args)
        if isinstance(result, dict) and result.get("ok") is False:
            message = str(result.get("error") or "Project YAML write failed.")
            self.app.call_from_thread(
                self._publish_tool_call_cell,
                tool_name,
                "error",
                description,
                normalized_args,
                message,
            )
            self.app.call_from_thread(
                self._handle_slash_tool_failure,
                tool_name,
                message,
            )
            return {"tool": tool_name, "status": "error", "result": result}

        self.app.call_from_thread(
            self._publish_tool_call_cell,
            tool_name,
            "ok",
            description,
            normalized_args,
            self._tool_success_note(tool_name, result),
        )
        self.app.call_from_thread(
            self._handle_slash_tool_success,
            tool_name,
            normalized_args,
            result,
        )
        return {"tool": tool_name, "status": "ok", "result": result}

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
            runtime_v2=self.runtime_v2,
        )
        policy = self._permission_policy()
        request = self.active_agent_session.state.request or "Continue."
        result = self.active_agent_session.run_loop(
            request,
            policy=policy,
            approver=lambda req: self._await_approval(req),
        )
        result["session_allow_tools"] = sorted(policy.session_allow)
        return result

    @work(
        thread=True,
        exclusive=True,
        exit_on_error=False,
        group="agent-session",
        name="agent-execute",
    )
    def execute_agent_session(self, session_id: str) -> dict[str, object]:
        from chemsmart.agent.transport import SshQsubTransport

        self.active_resume_id = session_id
        session_dir = self.session_root / session_id
        if self.active_agent_session is None or str(
            self.active_agent_session.session_dir
        ) != str(session_dir):
            self.active_agent_session = AgentSession.load(
                session_id,
                session_root=str(self.session_root),
                runtime_v2=self.runtime_v2,
            )
        session = self.active_agent_session

        ctx = _extract_execute_context(session_dir)
        job_handle_id = ctx.get("job_handle_id")
        if job_handle_id:
            parts = [
                "The dry-run has been reviewed and approved by the user.",
                f"The job is stored as handle '{job_handle_id}'.",
            ]
            if ctx.get("inputfile"):
                parts.append(f"The input file is at '{ctx['inputfile']}'.")
            if ctx.get("server_name"):
                parts.append(f"Server name: '{ctx['server_name']}'.")
            parts.append(f"Please call submit_hpc with job='{job_handle_id}'.")
            execute_request = " ".join(parts)
        else:
            execute_request = (
                "The dry-run and all pre-flight checks have been reviewed "
                "and approved by the user. Please submit the job to the "
                "HPC scheduler now."
            )

        original_registry = session.registry
        session.registry = _ExecuteOverrideRegistry(
            original_registry, SshQsubTransport()
        )
        policy = PermissionPolicy(
            mode=PermissionMode.DRIVING,
            yolo=False,
            session_allow=set(),
            driving_denylist=set(),
        )
        try:
            result = session.run_loop(
                execute_request,
                policy=policy,
                approver=lambda req: ApprovalDecision.ALLOW_ONCE,
            )
        finally:
            session.registry = original_registry
        result["session_allow_tools"] = []
        return result

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
            runtime_v2=self.runtime_v2,
        )
        return self.active_agent_session

    def start_request(self, text: str) -> None:
        if self._worker_is_busy():
            self.post_error(
                "Session already running",
                "The current request is still running. Press Tab to queue one follow-up.",
            )
            return
        if (
            self._active_provider_config is not None
            and self._active_provider_config.type == "local"
            and _is_calculation_diagnostic_request(text)
        ):
            self._start_local_diagnostic_request(text)
            return
        if self._active_provider_config is not None and (
            self._active_provider_config.type == "local"
        ):
            self._start_synthesis_request(text)
            return
        self._start_unified_request(text)

    def _worker_is_busy(self) -> bool:
        return bool(
            self._current_worker and not self._current_worker.is_finished
        )

    def _maybe_start_queued_prompt(self) -> None:
        if self._queued_prompt is None or self._worker_is_busy():
            return
        if self._pending_approval or self._waiting_for_user:
            self.query_one(FooterWidget).set_hint(
                "Follow-up remains queued until the requested input is resolved"
            )
            return
        prompt = self._queued_prompt
        self._queued_prompt = None
        footer = self.query_one(FooterWidget)
        footer.set_queued_prompt(False)
        if not self._request_history or self._request_history[-1] != prompt:
            self._request_history.append(prompt)
        self._history_cursor = len(self._request_history)
        self.call_after_refresh(self.start_request, prompt)

    def _start_unified_request(self, text: str) -> None:
        keep_conversational = self.active_agent_session is not None
        if not keep_conversational:
            self._stop_tailer()
        self._reset_request_state(
            clear_transcript=True,
            keep_conversational=keep_conversational,
        )
        self._current_request = text
        self._begin_turn()
        self.query_one(FooterWidget).set_phase(Phase.PLANNING)
        self.query_one(FooterWidget).set_hint("Unified agent is reasoning…")
        self.query_one(Transcript).add_cell(UserMessageCell(text))
        self._user_requests.add(text)
        self._current_worker = self.run_unified_session(text)
        if self._session_poll_timer is not None:
            self._session_poll_timer.stop()
        self._session_poll_timer = self.set_interval(
            0.1, self._attach_live_tailer, pause=False
        )

    def _start_synthesis_request(self, text: str) -> None:
        keep_conversational = self.active_synthesis_session is not None
        if not keep_conversational:
            self._stop_tailer()
        self._reset_request_state(
            clear_transcript=True,
            keep_conversational=keep_conversational,
        )
        self._current_request = text
        self._begin_turn()
        self.query_one(FooterWidget).set_phase(Phase.PLANNING)
        self.query_one(FooterWidget).set_hint("Model is synthesizing…")
        self.query_one(Transcript).add_cell(UserMessageCell(text))
        self._user_requests.add(text)
        self._current_worker = self.run_synthesis_session(text)

    def _begin_turn(self) -> str:
        self._turn_serial += 1
        self._active_turn_id = f"turn-{self._turn_serial}"
        self.query_one(Transcript).start_turn(self._active_turn_id)
        return self._active_turn_id

    def _start_local_diagnostic_request(self, text: str) -> None:
        self._reset_request_state(
            clear_transcript=True,
            keep_conversational=True,
        )
        self._current_request = text
        self._begin_turn()
        transcript = self.query_one(Transcript)
        transcript.add_cell(UserMessageCell(text))
        result = inspect_calculation(
            session_root=str(self.session_root),
        )
        calculation = result.get("calculation")
        if not result.get("ok") or not isinstance(calculation, dict):
            self.post_error(
                "Calculation result unavailable",
                str(result.get("error") or "No calculation result was found."),
            )
            return
        self._on_calculation_run(dict(calculation), persist=False)
        self.post_agent_message(
            _calculation_diagnostic_summary(calculation),
            title="Deterministic calculation diagnosis",
        )
        footer = self.query_one(FooterWidget)
        footer.set_phase(Phase.FINISHED)
        footer.set_hint("Calculation diagnosis ready")

    @work(
        thread=True,
        exclusive=True,
        exit_on_error=False,
        group="slash-tool",
        name="slash-tool",
    )
    def run_slash_tool_request(
        self,
        tool_name: str,
        arguments: dict[str, object],
        *,
        explicit_approval: bool = False,
    ) -> dict[str, object]:
        registry = ToolRegistry.default()
        normalized_args = registry.normalize_args(tool_name, arguments)
        description = registry.describe_tool(tool_name)
        request = ToolRequest(
            request_id=f"slash:{tool_name}:{get_ident()}",
            provider="slash",
            provider_call_id=f"slash-{tool_name}",
            name=tool_name,
            arguments_json=json.dumps(normalized_args, sort_keys=True),
            arguments=normalized_args,
            raw={"source": "slash"},
        )
        policy = self._permission_policy()
        resolved = policy.resolve(request)

        if explicit_approval:
            if tool_name != "execute_chemsmart_command":
                raise ValueError(
                    "explicit slash approval is only valid for command execution"
                )
            self.app.call_from_thread(
                self._publish_tool_call_cell,
                tool_name,
                "approved",
                description,
                normalized_args,
                "Approved explicitly by slash command.",
            )
        elif resolved.decision == ResolvedDecision.NEEDS_USER:
            self.app.call_from_thread(
                self._set_pending_tool_context,
                description,
                normalized_args,
            )
            self.app.call_from_thread(
                self._publish_tool_call_cell,
                tool_name,
                "pending",
                description,
                normalized_args,
                "Awaiting approval.",
            )
            decision = self._await_approval(request)
            if decision == ApprovalDecision.DENY:
                self.app.call_from_thread(
                    self._publish_tool_call_cell,
                    tool_name,
                    "denied",
                    description,
                    normalized_args,
                    "Denied by user.",
                )
                return {"tool": tool_name, "status": "denied"}
            policy.record(tool_name, decision)
            self.app.call_from_thread(
                self._sync_session_allow_tools,
                policy.session_allow,
            )
            approval_note = (
                "Approved for this session."
                if decision == ApprovalDecision.ALLOW_SESSION
                else "Approved once."
            )
            self.app.call_from_thread(
                self._publish_tool_call_cell,
                tool_name,
                "approved",
                description,
                normalized_args,
                approval_note,
            )
        elif resolved.decision == ResolvedDecision.AUTO_DENY:
            self.app.call_from_thread(
                self._publish_tool_call_cell,
                tool_name,
                "denied",
                description,
                normalized_args,
                "Blocked by current permissions; enable /yolo on to proceed.",
            )
            return {"tool": tool_name, "status": "denied"}
        else:
            self.app.call_from_thread(
                self._publish_tool_call_cell,
                tool_name,
                "approved",
                description,
                normalized_args,
                "Auto-approved.",
            )

        result = registry.call(tool_name, normalized_args)
        if isinstance(result, dict) and result.get("ok") is False:
            error = result.get("error")
            if isinstance(error, dict):
                error = error.get("message")
            message = str(
                error
                or result.get("stderr_tail")
                or f"Tool returned status {result.get('status') or 'error'}."
            )
            self.app.call_from_thread(
                self._publish_tool_call_cell,
                tool_name,
                "error",
                description,
                normalized_args,
                message,
            )
            self.app.call_from_thread(
                self._handle_slash_tool_failure,
                tool_name,
                message,
                result,
            )
            return {"tool": tool_name, "status": "error", "result": result}

        self.app.call_from_thread(
            self._publish_tool_call_cell,
            tool_name,
            "ok",
            description,
            normalized_args,
            self._tool_success_note(tool_name, result),
        )
        self.app.call_from_thread(
            self._handle_slash_tool_success,
            tool_name,
            normalized_args,
            result,
        )
        return {"tool": tool_name, "status": "ok", "result": result}

    @work(
        thread=True,
        exclusive=False,
        exit_on_error=False,
        group="calculation",
        name="calculation",
    )
    def run_calculation_request(
        self,
        ready: _ReadyCommand,
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
            if event.state == WorkerState.ERROR:
                error = event.worker.error
                self.post_error(
                    error.__class__.__name__
                    if error
                    else "Calculation worker error",
                    str(error)
                    if error
                    else "Unknown calculation worker error",
                )
                if not self._worker_is_busy():
                    self.query_one(FooterWidget).set_phase(Phase.ERROR)
                    self.query_one(FooterWidget).set_hint(
                        "Calculation worker failed"
                    )
            return
        if event.worker.group == "slash-tool":
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
            return
        if event.worker.group != "agent-session":
            return
        if event.state == WorkerState.SUCCESS:
            result = event.worker.result or {}
            self._current_worker = event.worker
            self._sync_footer_usage(result)
            if event.worker.name == "agent-ask":
                self._publish_synthesis_result(result)
                self._maybe_start_queued_prompt()
                return
            if self._tailer is not None:
                self._tailer.read_available()
            self.query_one(Transcript).collapse_tool_chain(
                self._active_turn_id
            )
            self.call_after_refresh(
                self.query_one(Transcript).collapse_tool_chain,
                self._active_turn_id,
            )
            self._session_allow_tools = set(
                result.get("session_allow_tools") or []
            )
            is_execute = event.worker.name == "agent-execute"
            if not is_execute:
                session_id = result.get("session_id")
                advisory = result.get("advisory_only", False)
                chitchat = result.get("is_chitchat", False)
                blocked = result.get("blocked", False)
                has_legacy_dry_run = bool(
                    result.get("dry_run_result")
                    or result.get("dry_run_results")
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

    def _sync_footer_usage(self, result: dict[str, object]) -> None:
        input_tokens = result.get("total_input_tokens")
        output_tokens = result.get("total_output_tokens")
        if not isinstance(input_tokens, int) and not isinstance(
            output_tokens, int
        ):
            return
        self.query_one(FooterWidget).set_usage(
            input_tokens=input_tokens
            if isinstance(input_tokens, int)
            else None,
            output_tokens=(
                output_tokens if isinstance(output_tokens, int) else None
            ),
        )

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

    def action_show_project_yaml(self) -> None:
        status = self._resolve_workspace_project_status()
        self._workspace_project_status = status
        self.query_one(FooterWidget).set_yaml_status(
            loaded=status.loaded,
            label=_yaml_footer_label(status),
        )
        if not status.candidates:
            if self.app.plain:
                self.post_agent_message(
                    "YAML MISSING\n\nBuild one with /init, then save it with "
                    "/write-project.",
                    title="Workspace YAML",
                )
                return
            self.app.push_screen(
                ProjectYamlOverlay(
                    title="Workspace YAML",
                    candidates=(),
                )
            )
            return
        if self.app.plain:
            if not status.loaded or status.path is None:
                candidates = "\n".join(
                    f"- `{path.parent.name}:{path.stem}`: `{path}`"
                    for path in status.candidates
                )
                self.post_agent_message(
                    "Multiple workspace YAML files are available. Select one "
                    "with an explicit project name in the request or CLI.\n\n"
                    + candidates,
                    title="Workspace YAML selection",
                )
                return
            try:
                yaml_text = status.path.read_text(encoding="utf-8")
            except OSError as exc:
                yaml_text = f"# Failed to read {status.path}: {exc}\n"
            self.post_agent_message(
                f"path: `{status.path}`\n\n```yaml\n{yaml_text.rstrip()}\n```",
                title=f"Workspace YAML: {status.program}:{status.project}",
            )
            return
        self.app.push_screen(
            ProjectYamlOverlay(
                title="Workspace project YAML",
                candidates=status.candidates,
                selected_path=status.path,
            ),
            self._on_project_yaml_selected,
        )

    def _on_project_yaml_selected(self, path: Path | None) -> None:
        if path is not None:
            self._activate_workspace_project_path(path)
        self.focus_composer()

    def action_show_shortcuts(self) -> None:
        if self.app.plain:
            self._show_help()
            return
        config = getattr(self.app, "tui_config", None)
        bindings = config.keybindings if config is not None else {}
        self.app.push_screen(ShortcutOverlay(bindings))

    def action_toggle_transcript(self) -> None:
        transcript = self.query_one(Transcript)
        expanded = transcript.toggle_detail_mode()
        self.query_one(FooterWidget).set_hint(
            "Transcript detail expanded" if expanded else "Transcript compact"
        )

    def action_show_activity(self) -> None:
        entries = [
            self._tool_cells[call_id].activity_snapshot()
            for call_id in self._tool_order
            if call_id in self._tool_cells
        ]
        context = {
            "phase": self.tui_state.phase.value,
            "operation": self.tui_state.operation,
            "progress": self.tui_state.tool_progress or "no active tool",
        }
        if self.app.plain:
            if not entries:
                self.post_agent_message(
                    "No tool activity in this turn.", title="Tool activity"
                )
                return
            lines = [
                f"phase: {context['phase']}",
                f"operation: {context['operation']}",
                f"workflow: {context['progress']}",
                "",
            ]
            lines.extend(
                f"- {entry['status']} {entry['tool']} {entry['elapsed']}"
                for entry in entries
            )
            self.post_agent_message("\n".join(lines), title="Tool activity")
            return
        self.app.push_screen(ToolActivityOverlay(entries, context=context))

    def action_search_history(self) -> None:
        if not self._request_history:
            self.notify("No request history yet.", timeout=2)
            return
        if not self.app.plain:
            self.app.push_screen(
                HistorySearchOverlay(self._request_history),
                self._handle_history_search,
            )
            return
        self._history_cursor = max(0, self._history_cursor - 1)
        composer = self.query_one(Composer)
        composer.load_text(self._request_history[self._history_cursor])
        composer.focus()
        self.query_one(FooterWidget).set_hint(
            f"History {self._history_cursor + 1}/{len(self._request_history)}"
        )

    def _handle_history_search(self, value: str | None) -> None:
        if not value:
            self.focus_composer()
            return
        composer = self.query_one(Composer)
        composer.load_text(value)
        composer.focus()
        self.query_one(FooterWidget).set_hint("History request restored")

    def action_context_tab(self) -> None:
        palette = self.query_one(SlashCommandPalette)
        composer = self.query_one(Composer)
        if palette.is_open:
            selected = palette.selected_item()
            if selected is not None:
                composer.load_text(f"{selected.command} ")
                palette.hide()
            return
        if not self._worker_is_busy():
            return
        draft = composer.resolve_text().strip()
        footer = self.query_one(FooterWidget)
        if draft and self._queued_prompt is None:
            self._queued_prompt = draft
            composer.clear_text()
            footer.update_draft("")
            footer.set_queued_prompt(True)
            footer.set_hint(
                "Follow-up queued · Tab on an empty draft restores it"
            )
            return
        if not draft and self._queued_prompt is not None:
            composer.load_text(self._queued_prompt)
            self._queued_prompt = None
            footer.set_queued_prompt(False)
            footer.set_hint("Queued follow-up restored for editing")
            return
        if draft and self._queued_prompt is not None:
            self.notify(
                "One follow-up is already queued.",
                severity="warning",
                timeout=3,
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
            self.call_after_refresh(self.focus_composer)
            return
        self.query_one(SlashCommandPalette).hide()
        self.focus_composer()

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

    def _update_slash_palette(self, text: str) -> None:
        palette = self.query_one(SlashCommandPalette)
        raw = text.lstrip()
        if not raw.startswith("/") or " " in raw:
            palette.hide()
            return
        query = raw[1:].lower()
        matches = self._slash_palette_items(query)
        palette.show_matches(query=query, matches=matches)

    def _slash_palette_items(self, query: str) -> list[SlashPaletteItem]:
        shortcuts = {
            "/help": getattr(self.app, "tui_config", None).keybindings.get(
                "show_shortcuts", "f1"
            )
            if getattr(self.app, "tui_config", None) is not None
            else "f1",
            "/jobs": "ctrl+b",
            "/runs": "ctrl+b",
        }
        items: list[SlashPaletteItem] = []
        for command, description in _SLASH_PALETTE_COMMANDS:
            aliases = {
                "/permissions": ("perms",),
                "/quit": ("q",),
                "/resume": ("continue",),
            }.get(command, ())
            if not (
                command[1:].startswith(query)
                or any(alias.startswith(query) for alias in aliases)
            ):
                continue
            reason = self._slash_unavailable_reason(command)
            items.append(
                SlashPaletteItem(
                    command,
                    description,
                    enabled=reason is None,
                    unavailable_reason=reason or "",
                    shortcut=shortcuts.get(command, ""),
                    aliases=aliases,
                )
            )
        return items

    def _slash_unavailable_reason(self, command: str) -> str | None:
        if (
            command in {"/allow", "/allow-session", "/deny"}
            and not self._pending_approval
        ):
            return "no approval is pending"
        if command == "/execute" and not self._last_dry_run_session_id:
            return "no validated dry-run is ready"
        if command == "/run" and not self._can_approve_or_execute("run"):
            return "no validated local command is ready"
        if command == "/submit" and not self._can_approve_or_execute("sub"):
            return "no validated submission command is ready"
        if command == "/critic" and self._current_verdict is None:
            return "no critic verdict yet"
        if command in {"/plan", "/rationale"} and self._current_plan is None:
            return "no active plan yet"
        return None

    def accept_slash_palette(self) -> str | bool | None:
        palette = self.query_one(SlashCommandPalette)
        if not palette.is_open:
            return None
        selected = palette.selected_item()
        if selected is None:
            reason = next(
                (
                    item.unavailable_reason
                    for item in palette.items
                    if item.unavailable_reason
                ),
                "No available command is selected.",
            )
            self.notify(reason, severity="warning", timeout=3)
            return False
        palette.hide()
        return selected.command

    def _sync_footer_provider(self) -> None:
        config = self._active_provider_config
        self._workspace_project_status = (
            self._resolve_workspace_project_status()
        )
        if config is not None:
            role = "synthesis" if config.type == "local" else "unified"
            self.query_one(FooterWidget).set_provider_model(
                f"{role}:{config.type}",
                config.model,
                project=(
                    self._workspace_project_status.project
                    or config.project
                    or resolve_default_project()
                ),
            )
        self.query_one(FooterWidget).set_yaml_status(
            loaded=self._workspace_project_status.loaded,
            label=_yaml_footer_label(self._workspace_project_status),
        )

    def _resolve_workspace_project_status(self) -> WorkspaceProjectStatus:
        selected_path = self._selected_workspace_project_path
        state_project = current_workflow_state().project
        if selected_path is None and state_project is not None:
            selected_path = Path(state_project.path).resolve()

        status = resolve_workspace_project(selected_path=selected_path)
        if not status.loaded and status.candidates:
            configured = (
                str(self._active_provider_config.project or "").strip()
                if self._active_provider_config is not None
                else ""
            )
            matches = tuple(
                path for path in status.candidates if path.stem == configured
            )
            if len(matches) == 1:
                status = resolve_workspace_project(selected_path=matches[0])

        self._selected_workspace_project_path = (
            status.path if status.loaded else None
        )
        return status

    def _activate_workspace_project_path(self, path: Path) -> bool:
        resolved = path.expanduser().resolve()
        status = resolve_workspace_project(selected_path=resolved)
        if not status.loaded or status.path is None:
            self.post_error(
                "Project selection failed",
                f"{resolved} is not a workspace project YAML.",
            )
            return False
        state_delta = select_workspace_project(
            status.project,
            status.program,
        )
        if not state_delta.get("selected"):
            self.post_error(
                "Project selection failed",
                str(state_delta.get("rule_id") or resolved),
            )
            return False
        self._selected_workspace_project_path = status.path
        self._workspace_project_status = status
        if self.active_synthesis_session is not None:
            self.active_synthesis_session.default_project = status.project
        footer = self.query_one(FooterWidget)
        footer.set_yaml_status(
            loaded=True,
            label=_yaml_footer_label(status),
        )
        footer.set_hint(f"Active project: {status.program}:{status.project}")
        return True

    def _apply_workspace_project_tool_result(
        self,
        tool_name: str,
        result: object,
    ) -> None:
        if tool_name not in {
            "read_project_yaml",
            "write_project_yaml",
            "update_project_yaml",
        } or not isinstance(result, dict):
            return
        raw_path = result.get("written_path") or result.get("path")
        if raw_path:
            self._activate_workspace_project_path(Path(str(raw_path)))

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

        self._ready_command = _ReadyCommand(
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
            footer.entity_status = format_assumptions_banner(
                self._current_entity_snapshot(),
                None,
            )
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
        elif isinstance(event, AssistantTurnEvent):
            if event.text.strip():
                transcript.add_cell(
                    AgentMessageCell(event.text.strip(), title="Assistant")
                )
        elif isinstance(event, ToolCallEvent):
            if event.status == "done":
                self._apply_workspace_project_tool_result(
                    event.tool,
                    event.payload,
                )
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
                footer.set_phase(Phase.EXECUTING)
                footer.set_hint(
                    "Local run in progress…"
                    if event.status == "running"
                    else "Local run complete"
                )
            else:
                footer.set_phase(Phase.TOOL_RUNNING)
                footer.set_hint(
                    f"{event.tool} in progress…"
                    if event.status == "running"
                    else f"{event.tool} complete"
                )
            footer.set_tool_progress(event.tool, step=event.step_index)
        elif isinstance(event, ToolPreviewEvent):
            if workflow_cell is not None and event.status == "running":
                workflow_cell.mark_started(
                    event.step_index,
                    event.tool,
                    event.args,
                )
            footer.set_phase(Phase.VALIDATING)
            footer.set_hint("Preparing submission preview…")
            footer.set_tool_progress(event.tool, step=event.step_index)
        elif isinstance(event, ToolUseEvent):
            self._apply_tool_use_event(event)
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
                    command=event.command,
                    cli_grounded=event.cli_grounded,
                    cli_grounding_issue=event.cli_grounding_issue,
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
            self._pending_tool_request = None
            if payload.get("request_intent") == "chitchat":
                footer.set_phase(Phase.IDLE)
                footer.set_hint("Ready")
                return
            if (
                not event.blocked
                and event.total_steps_planned == 0
                and event.total_steps_executed == 0
            ):
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
            if event.blocked and event.block_reason:
                summary += f" Block reason: {event.block_reason}."
            summary_cell = AgentMessageCell(summary, title="Summary")
            if event.blocked:
                setattr(summary_cell, "_chemsmart_final_deliverable", True)
            transcript.add_cell(summary_cell)
            transcript.collapse_tool_chain(self._active_turn_id)
        elif isinstance(event, IgnoredEvent):
            pass

        if payload.get("tool") in {"run_local", "submit_hpc"}:
            self._refresh_job_snapshot()

    def _has_user_message(self, request: str) -> bool:
        return request in self._user_requests

    def _render_ask_user_prompt(self, event: ToolUseEvent) -> None:
        """Render an ask_user tool call as a plain clarification question."""
        question = str(event.args.get("question") or "").strip()
        raw_options = event.args.get("options")
        options = (
            [str(o).strip() for o in raw_options if str(o).strip()]
            if isinstance(raw_options, list)
            else []
        )
        lines = ["The agent needs a bit more information to continue.", ""]
        if question:
            lines.extend([question, ""])
        if options:
            lines.append("Options (type the text or its number):")
            lines.append("")
            for index, option in enumerate(options, start=1):
                lines.append(f"{index}. {option}")
            lines.append("")
        lines.append("_Type your answer in the prompt below._")
        self.post_agent_message("\n".join(lines), title="Clarification needed")
        footer = self.query_one(FooterWidget)
        self._waiting_for_user = True
        footer.set_phase(Phase.WAITING_USER)
        footer.set_hint("Answer the question to continue")

    def _apply_tool_use_event(self, event: ToolUseEvent) -> None:
        transcript = self.query_one(Transcript)
        footer = self.query_one(FooterWidget)
        # ask_user is a clarification, not a permission-gated action. Render it
        # as a plain question (never a "[unknown]" risky pending tool cell) and
        # leave the composer ready for the user's answer.
        if event.tool == "ask_user":
            if event.status == "pending":
                self._render_ask_user_prompt(event)
            return
        note = None
        if event.status == "approved":
            if event.scope == "session":
                note = "approved for the rest of this session"
            else:
                note = "approved"
        elif event.status == "denied":
            note = (
                "Denied by user; model will continue without this action."
                if event.reason == "user_denied"
                else f"Denied: {event.reason or 'policy blocked'}"
            )
        elif event.status not in {"pending", "approved"}:
            note = (
                event.reason
                or _tool_use_payload_summary(event.payload)
                or event.status
            )

        call_id = event.provider_call_id or (
            f"legacy:{event.step_index}:{event.tool}"
        )
        if event.status == "ok":
            self._apply_workspace_project_tool_result(
                event.tool,
                event.payload,
            )
        self._upsert_tool_cell(
            provider_call_id=call_id,
            tool=event.tool,
            status=event.status,
            description=event.description or event.tool,
            arguments=event.args,
            note=note,
            queue_index=event.queue_index,
            queue_total=event.queue_total,
            session_rule_active=(
                event.tool in self._session_allow_tools
                if event.status == "pending"
                else False
            ),
            result=_public_tool_result_payload(event.payload),
        )
        footer.set_tool_progress(
            event.tool,
            step=event.queue_index or event.step_index or None,
            total=event.queue_total,
        )
        if event.status == "ok" and event.tool in {
            "synthesize_command",
            "repair_command",
        }:
            payload = event.payload if isinstance(event.payload, dict) else {}
            self._publish_command_tool_result(event.tool, payload)
        elif (
            event.status == "ok" and event.tool == "execute_chemsmart_command"
        ):
            payload = event.payload if isinstance(event.payload, dict) else {}
            self._publish_execute_tool_result(payload)
        if event.tool == "dry_run_input" and event.status == "ok":
            summary = _tool_use_summary_payload(event.payload)
            if summary is not None:
                content = str(summary.get("content") or "")
                inputfile = (
                    str(summary.get("inputfile"))
                    if summary.get("inputfile") is not None
                    else None
                )
                command = (
                    str(summary.get("command"))
                    if summary.get("command") is not None
                    else None
                )
                transcript.add_cell(
                    DryRunInputCell(
                        content,
                        inputfile=inputfile,
                        previous_content=self._latest_dry_run_content,
                        command=command,
                        cli_grounded=bool(summary.get("cli_grounded")),
                        cli_grounding_issue=(
                            str(summary.get("cli_grounding_issue"))
                            if summary.get("cli_grounding_issue") is not None
                            else None
                        ),
                    )
                )
                self._latest_dry_run_content = content

        if event.status == "pending":
            self._pending_approval_description = (
                event.description or event.tool
            )
            self._pending_approval_args = dict(event.args)
            self._pending_approval_index = event.queue_index
            self._pending_approval_total = event.queue_total
            footer.set_phase(
                Phase.APPROVAL_REQUIRED
                if self._permission_mode == PermissionMode.PERMISSION
                else Phase.TOOL_RUNNING
            )
            if self._permission_mode == PermissionMode.PERMISSION:
                footer.set_hint(
                    "PERMISSION_PENDING "
                    f"({event.queue_index or 1} of {event.queue_total or 1})"
                )
            else:
                footer.set_hint("DRIVING_AUTO")
        elif event.status == "denied":
            footer.set_phase(Phase.TOOL_RUNNING)
            footer.set_hint("DENIED_CONTINUING")
        elif event.status in {"approved", "ok"}:
            footer.set_phase(Phase.TOOL_RUNNING)
            footer.set_hint(
                "DRIVING_AUTO"
                if self._permission_mode == PermissionMode.DRIVING
                else "Awaiting next tool decision…"
            )
        elif event.status in {"error", "skipped", "interrupted"}:
            footer.set_phase(Phase.FAILED)
            footer.set_hint("Tool call reported an error")

        if event.status in {
            "ok",
            "partial",
            "error",
            "denied",
            "skipped",
            "ask_user",
        }:
            footer.entity_status = format_assumptions_banner(
                self._current_entity_snapshot(),
                event.status,
            )

    def _publish_command_tool_result(
        self,
        tool_name: str,
        payload: dict[str, object],
    ) -> None:
        command = str(payload.get("command") or "")
        status = str(payload.get("status") or "")
        semantic = payload.get("semantic")
        semantic_dict = semantic if isinstance(semantic, dict) else None
        intent = payload.get("intent") or payload.get("intent_assertion")
        intent_dict = intent if isinstance(intent, dict) else None
        if status == "ready" and command:
            self._remember_ready_command(
                command=command,
                semantic=semantic_dict,
                intent=intent_dict,
                source=tool_name,
            )
        decision_trace = payload.get("decision_trace")
        decision_dict = (
            decision_trace if isinstance(decision_trace, dict) else None
        )
        transcript = self.query_one(Transcript)
        transcript.add_cell(
            SynthesisTraceCell(
                provider_type="api",
                model=_provider_model_label(self._active_provider_config),
                mode="unified",
                status=status,
                command=command,
                semantic=semantic_dict,
                decision_trace=decision_dict,
                artifact_dir="decision_log.jsonl",
            )
        )
        if command:
            transcript.add_cell(
                CommandInterpretationCell(
                    parse_model_command(command),
                    expanded=True,
                )
            )
        if status == "needs_clarification":
            missing = payload.get("missing_info")
            if not isinstance(missing, list):
                missing = [str(missing or "details")]
            lines = "\n".join(f"- {item}" for item in missing)
            transcript.add_cell(
                FinalAnswerCell(
                    "I need more information before making a command:\n\n"
                    f"{lines}",
                    title="Clarification",
                )
            )
            return
        title = "Final Command" if command else "Final Status"
        if tool_name == "repair_command" and command:
            title = "Repaired Command"
        if command:
            transcript.add_cell(
                AgentMessageCell(
                    _command_details_text(
                        explanation=str(payload.get("explanation") or ""),
                        confidence=str(payload.get("confidence") or "medium"),
                        project=str(payload.get("project") or ""),
                    ),
                    title="Command details",
                )
            )
        transcript.add_cell(
            FinalAnswerCell(
                (
                    _final_command_text(command=command)
                    if command
                    else str(payload.get("explanation") or payload)
                ),
                title=title,
            )
        )

    def _publish_execute_tool_result(
        self,
        payload: dict[str, object],
    ) -> None:
        calculation = payload.get("calculation")
        if isinstance(calculation, dict):
            self._on_calculation_run(dict(calculation), persist=False)
            return
        parsed = parse_model_command(str(payload.get("command") or ""))
        fallback = {
            "run_id": f"legacy-{uuid.uuid4().hex[:8]}",
            "turn_id": self._active_turn_id,
            "command": str(payload.get("command") or ""),
            "cwd": str(Path.cwd()),
            "program": parsed.program or "",
            "kind": parsed.job or "",
            "label": parsed.label or parsed.filename or "calculation",
            "project": parsed.project or "",
            "input_path": parsed.filename or "",
            "status": ("completed" if payload.get("ok") else "process_failed"),
            "stage": str(payload.get("status") or "Execution finished"),
            "returncode": payload.get("returncode"),
            "error": "\n".join(
                str(payload.get("stderr_tail") or "").splitlines()[-40:]
            ),
        }
        self._on_calculation_run(fallback, persist=False)

    def _upsert_tool_cell(
        self,
        *,
        provider_call_id: str,
        tool: str,
        status: str,
        description: str,
        arguments: dict,
        note: str | None,
        queue_index: int | None = None,
        queue_total: int | None = None,
        session_rule_active: bool = False,
        result: dict | None = None,
    ) -> ToolCallCell:
        cell = self._tool_cells.get(provider_call_id)
        if cell is None:
            config = getattr(self.app, "tui_config", None)
            expanded = bool(config and config.tool_detail == "full")
            cell = ToolCallCell(
                tool=tool,
                status=status,
                description=description,
                arguments=arguments,
                note=note,
                queue_index=queue_index,
                queue_total=queue_total,
                session_rule_active=session_rule_active,
                result=result,
                provider_call_id=provider_call_id,
                expanded=expanded,
            )
            self._tool_cells[provider_call_id] = cell
            self._tool_order.append(provider_call_id)
            self.query_one(Transcript).add_cell(cell)
            return cell
        cell.update_lifecycle(
            status=status,
            description=description,
            arguments=arguments,
            note=note,
            queue_index=queue_index,
            queue_total=queue_total,
            session_rule_active=session_rule_active,
            result=result,
        )
        return cell

    def _current_entity_snapshot(self) -> dict[str, object] | None:
        session_dir = self.current_session_dir()
        session = self.active_agent_session

        if session_dir is not None:
            log_path = session_dir / "decision_log.jsonl"
            if log_path.exists():
                try:
                    memory = ConversationMemory.from_entries(
                        DecisionLog(log_path).read_all()
                    )
                    entities = memory.entities.model_dump(exclude_none=True)
                    return entities or None
                except Exception:
                    pass

        if session is None:
            return None
        entities = session.conversation_history.entities.model_dump(
            exclude_none=True
        )
        return entities or None

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
        elif command == "/mode":
            self._handle_mode_command(argument)
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
        requested = Path(
            str(candidate.get("project_name") or "project")
        ).stem
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
            self.query_one(FooterWidget).set_hint("Project YAML write cancelled")
            self.focus_composer()
            return
        candidate["project_name"] = result.project_name
        candidate["overwrite"] = result.action == "overwrite"
        self._start_project_yaml_write(candidate)

    def _start_project_yaml_write(
        self, candidate: dict[str, object]
    ) -> None:
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

    def _handle_wizard_probe_command(self, argument: str) -> None:
        if self._current_worker and not self._current_worker.is_finished:
            self.post_error(
                "Session already running",
                "Wait for the current request to finish before starting a new request.",
            )
            return
        try:
            parts = shlex.split(argument)
        except ValueError as exc:
            self.post_error("Invalid wizard command", str(exc))
            return
        if not parts or len(parts) > 2:
            self.post_error(
                "Invalid wizard command",
                "Usage: /wizard <name> [host]",
            )
            return
        server_name = parts[0]
        ssh_host_hint = parts[1] if len(parts) == 2 else None
        footer = self.query_one(FooterWidget)
        footer.set_phase(Phase.PLANNING)
        footer.set_hint("Wizard is probing…")
        self._current_worker = self.run_slash_tool_request(
            "wizard_probe",
            {
                "server_name": server_name,
                "ssh_host_hint": ssh_host_hint,
            },
        )

    def _handle_wizard_verify_command(self, argument: str) -> None:
        if self._current_worker and not self._current_worker.is_finished:
            self.post_error(
                "Session already running",
                "Wait for the current request to finish before starting a new request.",
            )
            return
        try:
            parts = shlex.split(argument)
        except ValueError as exc:
            self.post_error("Invalid wizard-verify command", str(exc))
            return
        if len(parts) != 1:
            self.post_error(
                "Invalid wizard-verify command",
                "Usage: /wizard-verify <name>",
            )
            return
        footer = self.query_one(FooterWidget)
        footer.set_phase(Phase.PLANNING)
        footer.set_hint("Wizard verify is checking transport wiring…")
        self._current_worker = self.run_slash_tool_request(
            "wizard_verify",
            {"server_name": parts[0]},
        )

    def _handle_wizard_refresh_command(self, argument: str) -> None:
        if self._current_worker and not self._current_worker.is_finished:
            self.post_error(
                "Session already running",
                "Wait for the current request to finish before starting a new request.",
            )
            return
        try:
            parts = shlex.split(argument)
        except ValueError as exc:
            self.post_error("Invalid wizard-refresh command", str(exc))
            return
        if not parts or len(parts) > 2:
            self.post_error(
                "Invalid wizard-refresh command",
                "Usage: /wizard-refresh <name> [--force]",
            )
            return
        force = False
        if len(parts) == 2:
            if parts[1] not in {"force", "--force"}:
                self.post_error(
                    "Invalid wizard-refresh command",
                    "Usage: /wizard-refresh <name> [--force]",
                )
                return
            force = True
        footer = self.query_one(FooterWidget)
        footer.set_phase(Phase.PLANNING)
        footer.set_hint("Wizard refresh is probing cache state…")
        self._current_worker = self.run_slash_tool_request(
            "wizard_refresh",
            {"server_name": parts[0], "force": force},
        )

    def _handle_wizard_write_command(self, argument: str) -> None:
        if self._current_worker and not self._current_worker.is_finished:
            self.post_error(
                "Session already running",
                "Wait for the current request to finish before starting a new request.",
            )
            return
        try:
            parts = shlex.split(argument)
        except ValueError as exc:
            self.post_error("Invalid wizard-write command", str(exc))
            return
        overwrite = False
        if parts:
            if len(parts) == 1 and parts[0] in {"overwrite", "--overwrite"}:
                overwrite = True
            else:
                self.post_error(
                    "Invalid wizard-write command",
                    "Usage: /wizard-write [overwrite]",
                )
                return
        probe = self._latest_wizard_probe
        if not isinstance(probe, dict):
            self.post_error(
                "Wizard write unavailable",
                "Run /wizard <name> [host] first.",
            )
            return
        validation = probe.get("validation")
        if isinstance(validation, dict) and not validation.get("ok", False):
            self.post_error(
                "Wizard write unavailable",
                "The latest wizard result did not validate.",
            )
            return
        server_name = str(probe.get("server_name") or "")
        yaml_text = str(probe.get("yaml_text") or "")
        if not server_name or not yaml_text:
            self.post_error(
                "Wizard write unavailable",
                "The latest wizard result is incomplete.",
            )
            return
        footer = self.query_one(FooterWidget)
        footer.set_phase(Phase.PLANNING)
        footer.set_hint("Wizard is writing…")
        self._current_worker = self.run_slash_tool_request(
            "wizard_write",
            {
                "server_name": server_name,
                "yaml_text": yaml_text,
                "overwrite": overwrite,
            },
        )

    def _set_pending_tool_context(
        self,
        description: str,
        arguments: dict,
    ) -> None:
        self._pending_approval_description = description
        self._pending_approval_args = dict(arguments)
        self._pending_approval_index = 1
        self._pending_approval_total = 1

    def _publish_tool_call_cell(
        self,
        tool_name: str,
        status: str,
        description: str,
        arguments: dict,
        note: str,
    ) -> None:
        call_id = self._direct_tool_call_ids.get(tool_name)
        existing = self._tool_cells.get(call_id or "")
        if call_id is None or (
            status in {"pending", "approved"}
            and existing is not None
            and existing.status
            in {"ok", "error", "denied", "skipped", "interrupted"}
        ):
            call_id = f"direct:{tool_name}:{len(self._tool_order) + 1}"
            self._direct_tool_call_ids[tool_name] = call_id
        self._upsert_tool_cell(
            provider_call_id=call_id,
            tool=tool_name,
            status=status,
            description=description,
            arguments=arguments,
            note=note,
            queue_index=1,
            queue_total=1,
            session_rule_active=tool_name in self._session_allow_tools,
        )

    def _sync_session_allow_tools(self, allowed: set[str]) -> None:
        self._session_allow_tools = set(allowed)

    def _handle_slash_tool_failure(
        self,
        tool_name: str,
        message: str,
        result: dict[str, object] | None = None,
    ) -> None:
        self.query_one(FooterWidget).set_phase(Phase.ERROR)
        self.query_one(FooterWidget).set_hint("Slash command failed")
        if tool_name == "execute_chemsmart_command" and result is not None:
            self._ready_command = None
            self._publish_execute_tool_result(result)
        self.post_error(tool_name, message)

    def _handle_slash_tool_success(
        self,
        tool_name: str,
        arguments: dict[str, object],
        result: object,
    ) -> None:
        footer = self.query_one(FooterWidget)
        footer.set_phase(Phase.FINISHED)
        self._apply_workspace_project_tool_result(tool_name, result)
        if tool_name == "execute_chemsmart_command" and isinstance(
            result, dict
        ):
            self._ready_command = None
            self._publish_execute_tool_result(result)
            footer.set_hint("Command execution completed")
            return
        if tool_name == "wizard_probe" and isinstance(result, dict):
            self._latest_wizard_probe = result
            server_name = str(
                result.get("server_name")
                or arguments.get("server_name")
                or "wizard"
            )
            yaml_text = str(result.get("yaml_text") or "")
            self.post_agent_message(
                f"```yaml\n{yaml_text.rstrip()}\n```",
                title=f"Wizard: {server_name}",
            )
            validation = result.get("validation")
            if isinstance(validation, dict) and not validation.get(
                "ok", False
            ):
                errors = validation.get("errors") or []
                self.post_error(
                    "Wizard validation failed",
                    "\n".join(str(error) for error in errors)
                    or "wizard output did not validate",
                )
                footer.set_phase(Phase.ERROR)
                footer.set_hint("Wizard validation failed")
                return
            footer.set_hint("Wizard YAML ready")
            return
        if tool_name == "wizard_write" and isinstance(result, dict):
            self.post_agent_message(
                f"Wrote `{result.get('written_path')}`.",
                title="Wizard",
            )
            footer.set_hint("Wizard YAML written")
            return
        if tool_name == "write_project_yaml" and isinstance(result, dict):
            project = str(result.get("project_name") or "project")
            program = str(result.get("program") or "gaussian")
            written_path = str(result.get("written_path") or "")
            command_example = (
                f"chemsmart run {program} -p {project} "
                "-f examples/h2o.xyz -c 0 -m 1 opt"
            )
            self.post_agent_message(
                (
                    f"Wrote `{written_path}`.\n\n"
                    "Deterministic use:\n\n"
                    "```bash\n"
                    f"{command_example}\n"
                    "```\n\n"
                    f"This workspace now has `{program}:{project}` loaded. "
                    "The command-synthesis harness will attach and validate "
                    f"this project automatically, and explicit CLI commands "
                    f"can use `-p {project}` from this same workspace."
                ),
                title="Project YAML",
            )
            footer.set_hint(f"Project YAML written: {project}")
            return
        if tool_name == "wizard_refresh" and isinstance(result, dict):
            self.post_agent_message(
                self._wizard_refresh_result_table(result),
                title=f"Wizard refresh: {result.get('server_name') or 'server'}",
            )
            status = str(result.get("status") or "")
            if status == "error":
                footer.set_phase(Phase.ERROR)
                footer.set_hint("Wizard refresh failed")
                return
            if status == "stale":
                footer.set_hint("Wizard refresh preserved stale cache")
                return
            footer.set_hint("Wizard cache ready")
            return
        if tool_name == "wizard_verify" and isinstance(result, dict):
            self.post_agent_message(
                self._wizard_verify_result_table(result),
                title=f"Wizard verify: {result.get('server_name') or 'server'}",
            )
            errors = result.get("errors") or []
            warnings = result.get("warnings") or []
            if errors:
                self.post_error(
                    "Wizard verify failed",
                    "\n".join(str(error) for error in errors),
                )
                footer.set_phase(Phase.ERROR)
                footer.set_hint("Wizard verify found errors")
                return
            if warnings:
                footer.set_hint("Wizard verify completed with warnings")
                return
            footer.set_hint("Wizard verify completed")
            return
        footer.set_hint("Slash command complete")

    def _tool_success_note(self, tool_name: str, result: object) -> str:
        if tool_name == "wizard_probe" and isinstance(result, dict):
            validation = result.get("validation")
            if isinstance(validation, dict) and validation.get("ok", False):
                return "Rendered validated wizard YAML."
            return "Rendered wizard YAML with validation issues."
        if tool_name == "wizard_verify" and isinstance(result, dict):
            if result.get("errors"):
                return "Verified transport wiring and found errors."
            if result.get("warnings"):
                return "Verified transport wiring with warnings."
            return "Verified transport wiring."
        if tool_name == "wizard_refresh" and isinstance(result, dict):
            status = str(result.get("status") or "")
            if status == "error":
                return (
                    "Wizard refresh failed and recorded an error cache entry."
                )
            if status == "stale":
                return "Wizard refresh failed; preserved the last-good stale cache."
            return "Wizard refresh produced a fresh cache entry."
        if tool_name == "wizard_write" and isinstance(result, dict):
            return f"Wrote {result.get('written_path')}"
        if tool_name == "write_project_yaml" and isinstance(result, dict):
            return f"Wrote {result.get('written_path')}"
        return "Completed successfully."

    def _wizard_verify_result_table(self, result: dict[str, object]) -> Table:
        table = Table(show_header=True, box=None, padding=(0, 1))
        table.add_column("Field", style="cyan", no_wrap=True)
        table.add_column("Value")
        table.add_row("server_name", str(result.get("server_name") or "-"))
        table.add_row("host", str(result.get("host") or "-"))
        table.add_row("mode", str(result.get("mode") or "-"))
        table.add_row(
            "would_submit_via",
            str(result.get("would_submit_via") or "-"),
        )
        invocation = result.get("transport_invocation")
        table.add_row(
            "transport_invocation",
            json.dumps(invocation or []),
        )
        warnings = result.get("warnings") or []
        table.add_row(
            "warnings",
            "\n".join(str(item) for item in warnings) if warnings else "-",
        )
        errors = result.get("errors") or []
        table.add_row(
            "errors",
            "\n".join(str(item) for item in errors) if errors else "-",
        )
        return table

    def _wizard_refresh_result_table(self, result: dict[str, object]) -> Table:
        table = Table(show_header=True, box=None, padding=(0, 1))
        table.add_column("Field", style="cyan", no_wrap=True)
        table.add_column("Value")
        table.add_row("cache_path", str(result.get("cache_path") or "-"))
        table.add_row("status", str(result.get("status") or "-"))
        table.add_row("host", str(result.get("host") or "-"))
        table.add_row("mode", str(result.get("mode") or "-"))
        table.add_row("scheduler", str(result.get("scheduler") or "-"))
        table.add_row("probed_at", str(result.get("probed_at") or "-"))
        node_summary = result.get("node_summary")
        if not isinstance(node_summary, dict):
            node_summary = {}
        table.add_row(
            "selected_queue",
            str(node_summary.get("selected_queue") or "-"),
        )
        table.add_row(
            "resources",
            f"cpu={node_summary.get('cpu')} mem_gb={node_summary.get('mem_gb')} gpu={node_summary.get('gpu')}",
        )
        table.add_row(
            "project",
            str(node_summary.get("project") or "-"),
        )
        table.add_row(
            "scratch",
            f"{node_summary.get('scratch_dir') or '-'} (writable={node_summary.get('scratch_writable')})",
        )
        table.add_row("last_error", str(result.get("last_error") or "-"))
        return table

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

    def _handle_mode_command(self, argument: str) -> None:
        value = argument.strip().lower()
        if value and value not in {"ask", "run", "unified"}:
            self.post_error("Unknown mode", "Usage: /mode")
            return
        provider_type = _provider_type_label(self._active_provider_config)
        provider_role = (
            "CLI synthesis"
            if provider_type == "local"
            else "unified tool loop"
        )
        self.post_agent_message(
            (
                "```\n"
                "interface: unified\n"
                f"provider_role: {provider_role}\n"
                "planning_and_validation: natural-language request\n"
                "approved_local_execution: /run\n"
                "approved_hpc_submission: /submit\n"
                "```\n\n"
                "`/mode ask` and `/mode run` are compatibility aliases only; "
                "they no longer switch or reset the session."
            ),
            title="Unified interface",
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
            program = str(
                pending.arguments.get("program") or "gaussian"
            ).strip().lower()
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
        if result is None or pending is None or pending.name != "write_project_yaml":
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


class _ExecuteOverrideRegistry:
    """Wraps ToolRegistry to inject execute=True + transport into submit_hpc.

    Used by execute_agent_session so the LLM does not need to know about the
    internal execute flag — the framework injects it transparently.
    """

    def __init__(self, base, transport) -> None:
        self._base = base
        self._transport = transport

    def __getattr__(self, name: str):
        return getattr(self._base, name)

    def call(self, tool_name: str, args: dict) -> object:
        if tool_name == "submit_hpc":
            args = {**args, "execute": True, "transport": self._transport}
        return self._base.call(tool_name, args)


def _extract_execute_context(session_dir: Path) -> dict[str, str | None]:
    """Read the session decision log and pull out the job handle and inputfile.

    Returns a dict with 'job_handle_id', 'inputfile', and 'server_name'.
    All values may be None if not found.
    """
    log_path = session_dir / "decision_log.jsonl"
    job_handle_id: str | None = None
    inputfile: str | None = None
    server_name: str | None = None
    if not log_path.exists():
        return {
            "job_handle_id": job_handle_id,
            "inputfile": inputfile,
            "server_name": server_name,
        }
    try:
        for raw in log_path.read_text(encoding="utf-8").splitlines():
            if not raw.strip():
                continue
            entry = json.loads(raw)
            kind = entry.get("kind")
            payload = entry.get("payload") or {}
            if kind == "tool_use_result" and isinstance(payload, dict):
                tool = payload.get("tool")
                if tool == "build_job" and payload.get("status") == "ok":
                    job_handle_id = str(payload["handle_id"])
                elif tool == "dry_run_input" and payload.get("status") == "ok":
                    inner = payload.get("payload") or {}
                    summary = inner.get("summary") or {}
                    inputfile = summary.get("inputfile") or None
            elif kind == "tool_use_request" and isinstance(payload, dict):
                if payload.get("tool") == "submit_hpc":
                    args = payload.get("args") or {}
                    server_name = args.get("server_name") or server_name
    except Exception:
        pass
    return {
        "job_handle_id": job_handle_id,
        "inputfile": inputfile,
        "server_name": server_name,
    }


def _latest_project_yaml_candidate(
    session_dir: Path | None,
) -> dict[str, object] | None:
    if session_dir is None:
        return None
    log_path = session_dir / "decision_log.jsonl"
    if not log_path.exists():
        return None

    candidate: dict[str, object] | None = None
    candidate_yaml: object | None = None
    validated = False
    validation_requests: dict[str, dict[str, object]] = {}
    try:
        for raw in log_path.read_text(encoding="utf-8").splitlines():
            if not raw.strip():
                continue
            entry = json.loads(raw)
            payload = entry.get("payload") or {}
            if not isinstance(payload, dict):
                continue
            tool = payload.get("tool")
            call_id = str(payload.get("provider_call_id") or "")
            if (
                entry.get("kind") == "tool_use_request"
                and tool == "validate_project_yaml"
                and call_id
            ):
                args = payload.get("args") or {}
                yaml_text = args.get("yaml_text") if isinstance(args, dict) else None
                if isinstance(yaml_text, str) and not yaml_text.strip():
                    continue
                if not isinstance(yaml_text, (str, dict)):
                    continue
                validation_requests[call_id] = {
                    "project_name": str(
                        args.get("project_name") or "project"
                    ),
                    "program": str(args.get("program") or "gaussian"),
                    "yaml_text": yaml_text,
                }
                continue
            if payload.get("status") not in {None, "ok"}:
                continue

            result = payload.get("payload") or {}
            if isinstance(result, dict) and isinstance(
                result.get("summary"), dict
            ):
                result = result["summary"]
            if not isinstance(result, dict):
                continue
            if tool == "validate_project_yaml":
                if result.get("verdict") not in {"ok", "warn", "reject"}:
                    continue
                request_candidate = validation_requests.pop(call_id, None)
                validated = result.get("verdict") in {"ok", "warn"}
                if validated and request_candidate is not None:
                    candidate = request_candidate
                    candidate_yaml = request_candidate["yaml_text"]
                    continue
                # Older logs recorded only results. They may validate a prior
                # accepted render, but never a rejected or missing candidate.
                if candidate is None:
                    validated = False
                    continue
                if not validated:
                    candidate = None
                    candidate_yaml = None
                continue
            if tool != "render_project_yaml":
                continue
            yaml_text = result.get("yaml_text")
            if not isinstance(yaml_text, str) or not yaml_text.strip():
                continue
            validation = result.get("validation")
            render_rejected = result.get("ok") is False or (
                isinstance(validation, dict)
                and validation.get("verdict") == "reject"
            )
            if render_rejected:
                candidate = None
                candidate_yaml = None
                validated = False
                continue
            # A re-render of the identical candidate (build-mode over-iteration)
            # must NOT discard a prior successful validation — only a genuinely
            # new/changed candidate resets the validated flag and needs
            # re-validation before it can be written.
            if yaml_text == candidate_yaml:
                continue
            candidate = {
                "project_name": str(result.get("project_name") or "project"),
                "program": str(result.get("program") or "gaussian"),
                "yaml_text": yaml_text,
            }
            candidate_yaml = yaml_text
            validated = False
    except Exception:
        return None
    return candidate if validated else None


def _yaml_footer_label(status: WorkspaceProjectStatus) -> str:
    if status.loaded:
        return f"YAML OK {status.program}:{status.project}"
    if status.candidates:
        return f"YAML SELECT {len(status.candidates)}"
    return "YAML MISSING"


def _find_project_yaml_candidate_for_write(
    session_root: Path,
    *,
    preferred_session_dir: Path | None = None,
) -> dict[str, object] | None:
    candidate = _latest_project_yaml_candidate(preferred_session_dir)
    if candidate is not None:
        return candidate
    if not session_root.exists():
        return None
    try:
        session_dirs = sorted(
            [
                path
                for path in session_root.iterdir()
                if path.is_dir() and not path.name.startswith(".")
            ],
            key=lambda path: path.stat().st_mtime,
            reverse=True,
        )
    except OSError:
        return None
    for session_dir in session_dirs:
        if (
            preferred_session_dir is not None
            and session_dir == preferred_session_dir
        ):
            continue
        candidate = _latest_project_yaml_candidate(session_dir)
        if candidate is not None:
            return candidate
    return None


def _tool_use_payload_summary(payload: object) -> str | None:
    if not isinstance(payload, dict):
        return None
    if payload.get("ok") is False:
        error = payload.get("error")
        if isinstance(error, dict):
            return str(error.get("message") or "") or None
    if "handle_id" in payload:
        return f"returned {payload['handle_id']}"
    if "inputfile" in payload:
        return str(payload.get("inputfile"))
    if "ok" in payload and payload.get("ok") is True:
        return "completed successfully"
    return None


def _tool_use_summary_payload(payload: object) -> dict[str, object] | None:
    if not isinstance(payload, dict):
        return None
    summary = payload.get("summary")
    if isinstance(summary, dict):
        return summary
    if "content" in payload or "inputfile" in payload:
        return payload
    return None


def _public_tool_result_payload(payload: object) -> dict[str, object] | None:
    """Project tool output without raw model responses or private reasoning."""

    if not isinstance(payload, dict):
        return None
    summary = payload.get("summary")
    if isinstance(summary, dict):
        payload = summary
    allowed = {
        "ok",
        "status",
        "error",
        "message",
        "path",
        "handle_id",
        "inputfile",
        "command",
        "semantic",
        "verdict",
        "failed_rule_ids",
        "issues",
        "job_id",
        "server",
        "returncode",
        "cli_grounded",
        "cli_grounding_issue",
        "calculation",
    }
    projected = {
        key: value for key, value in payload.items() if key in allowed
    }
    return projected or None


def _jobs_sort_key(job: dict) -> tuple[int, str, str]:
    order = {"running": 0, "queued": 1, "failed": 2, "cancelled": 3, "done": 4}
    return (
        order.get(str(job.get("status")), 9),
        str(job.get("raw_started") or ""),
        str(job.get("job_id") or ""),
    )


def _load_tui_provider_config() -> AgentProviderConfig | None:
    # Mirror ``get_provider``: the active ``agent.yaml`` config is the source of
    # truth for which provider the harness/synthesis path actually builds. A
    # stale ``AI_PROVIDER`` env var must NOT shadow it — otherwise the TUI picks
    # a mode that does not match the provider it will instantiate (e.g. defaults
    # to run/harness while ``get_provider`` returns a local provider, crashing
    # with "Unsupported provider 'local'").
    try:
        return load_active_provider_config()
    except AgentProviderConfigError:
        return None


def _provider_type_label(config: AgentProviderConfig | None) -> str:
    if config is None:
        return "offline"
    return str(config.type or "offline").strip().lower() or "offline"


def _provider_model_label(config: AgentProviderConfig | None) -> str:
    if config is None:
        return "auto"
    return str(config.model or "auto").strip() or "auto"


def _new_synthesis_artifact_dir(
    session_root: Path,
    *,
    provider_type: str,
) -> Path:
    lane = "local" if provider_type == "local" else "api"
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = session_root / lane / "synthesis" / stamp
    path.mkdir(parents=True, exist_ok=True)
    return path


def _write_synthesis_artifact(
    artifact_dir: Path,
    payload: dict[str, object],
) -> None:
    try:
        artifact_dir.mkdir(parents=True, exist_ok=True)
        (artifact_dir / "synthesis_turn.json").write_text(
            json.dumps(payload, indent=2, sort_keys=True, default=str),
            encoding="utf-8",
        )
    except OSError:
        return


def _default_interaction_mode(config: AgentProviderConfig | None) -> str:
    return (
        "synthesis"
        if config is not None and config.type == "local"
        else "unified"
    )


def _ready_command_hint(ready: _ReadyCommand | None) -> str:
    if ready is None:
        return "Command shown, but execution evidence is incomplete"
    if ready.action == "run":
        return "Command validated — /run to execute locally"
    return "Command validated — /submit to submit for real"


def _is_calculation_diagnostic_request(request: str) -> bool:
    text = str(request or "").lower()
    return any(
        marker in text
        for marker in (
            "diagnose the result",
            "inspect the result",
            "analyze the result",
            "calculation result",
            "check the output",
            "계산 결과",
            "결과를 진단",
            "결과 분석",
            "출력 파일 확인",
            "诊断结果",
            "分析计算结果",
        )
    )


def _calculation_diagnostic_summary(calculation: dict[str, object]) -> str:
    program = str(calculation.get("program") or "calculation").upper()
    kind = str(calculation.get("kind") or "job").upper()
    status = str(calculation.get("status") or "parsed")
    lines = [f"{program} {kind} status: `{status}`"]
    energy = calculation.get("energy")
    if isinstance(energy, (int, float)):
        lines.append(f"- Final electronic energy: `{float(energy):.12f} Eh`")
    cycles = calculation.get("scf_cycles")
    if isinstance(cycles, int):
        lines.append(f"- SCF convergence: `{cycles} cycles`")
    imag = list(calculation.get("imag_freqs") or [])
    if imag:
        lines.append(
            "- Imaginary frequencies: `"
            + ", ".join(f"{float(value):.1f}" for value in imag)
            + " cm^-1`"
        )
    normal = calculation.get("normal_termination")
    if normal is not None:
        lines.append(
            "- Program termination: `normal`"
            if normal
            else "- Program termination: `not confirmed`"
        )
    output_path = str(calculation.get("output_path") or "")
    if output_path:
        lines.append(f"- Output: `{output_path}`")
    if status != "completed":
        error = str(calculation.get("error") or "").strip()
        if error:
            lines.extend(["", "Relevant diagnostic context:", "```text"])
            lines.extend(error.splitlines()[-40:])
            lines.append("```")
    return "\n".join(lines)


def _decision_trace_dict(
    synthesis: dict[str, object],
) -> dict[str, object] | None:
    trace = synthesis.get("decision_trace")
    return trace if isinstance(trace, dict) and trace else None


def _final_command_text(
    *,
    command: str,
) -> str:
    return "\n".join(("```bash", command, "```"))


def _command_details_text(
    *,
    explanation: str,
    confidence: str,
    project: str,
) -> str:
    parts = [f"confidence: `{confidence}`"]
    if project:
        parts.extend(["", f"active project: `{project}`"])
    if explanation:
        parts.extend(["", explanation])
    return "\n".join(parts)


def _final_answer_text(*, command: str, explanation: str) -> str:
    parts: list[str] = []
    if command:
        parts.extend(["Analyzed command:", "", "```bash", command, "```", ""])
    parts.append(explanation)
    return "\n".join(parts)


def _format_semantic_result(semantic: dict[str, object] | None) -> str:
    if not semantic:
        return ""
    verdict = str(semantic.get("verdict") or "unknown")
    failed = semantic.get("failed_rule_ids") or []
    if isinstance(failed, list) and failed:
        failed_text = ", ".join(str(item) for item in failed)
    else:
        failed_text = "none"
    issues = semantic.get("issues") or []
    issue_lines = []
    if isinstance(issues, list):
        for issue in issues:
            if not isinstance(issue, dict):
                continue
            rule = str(issue.get("rule_id") or "rule")
            message = str(issue.get("message") or "")
            issue_lines.append(f"- `{rule}`: {message}")
    generated = semantic.get("generated_inputs") or []
    generated_lines = []
    if isinstance(generated, list):
        for item in generated:
            if not isinstance(item, dict):
                continue
            path = str(item.get("path") or "generated input")
            route = str(item.get("route") or "").strip()
            generated_lines.append(f"- `{path}` {route}".rstrip())
    body = [
        "",
        "",
        "runtime semantic gate:",
        f"- verdict: `{verdict}`",
        f"- failed_rule_ids: `{failed_text}`",
    ]
    if issue_lines:
        body.append("- issues:")
        body.extend(issue_lines)
    if generated_lines:
        body.append("- generated input evidence:")
        body.extend(generated_lines[:3])
    return "\n".join(body)


def _format_synthesis_exception(exc: Exception) -> str:
    message = str(exc) or exc.__class__.__name__
    if "MLX runtime requires Apple Silicon/Metal and mlx-lm" in message:
        return (
            "The active local provider is configured for MLX, but this Python "
            "environment cannot load `mlx_lm` or Apple Metal. Start the TUI "
            "with the MLX-enabled interpreter, for example "
            "`/Users/hongjiseung/developer/chemsmart/.venv-mlx/bin/python -m "
            "chemsmart.cli.main agent`, or configure an API frontier provider "
            "in `agent.yaml`. Active-env install command: "
            "`python -m pip install 'mlx-lm==0.31.3'`."
        )
    if "local provider" in message.lower() or "mlx" in message.lower():
        return (
            "The local model provider could not start.\n\n"
            f"{message}\n\n"
            "Use `/doctor` to inspect the active provider, run the TUI from the "
            "MLX-enabled environment, or configure an API frontier provider "
            "in `agent.yaml`."
        )
    return message
