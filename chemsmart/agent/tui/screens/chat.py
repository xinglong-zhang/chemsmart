"""Main chat screen for the chemsmart agent TUI."""

from __future__ import annotations

import uuid
from datetime import datetime, timezone
from pathlib import Path
from threading import Event

from textual.app import ComposeResult
from textual.binding import Binding
from textual.screen import Screen
from textual.widgets import OptionList
from textual.worker import Worker

from chemsmart.agent.core import CriticVerdict, DecisionLog, Plan
from chemsmart.agent.permissions import (
    ApprovalDecision,
    PermissionMode,
)
from chemsmart.agent.provider_adapter import ToolRequest
from chemsmart.agent.harness.workflow_state import (
    current_workflow_state,
)
from chemsmart.agent.synthesis import SynthesisSession
from chemsmart.settings.workspace_project import (
    resolve_workspace_project,
)
from chemsmart.agent.tui.chat_models import ReadyCommand
from chemsmart.agent.tui.slash_catalog import (
    SLASH_PALETTE_COMMANDS as _SLASH_PALETTE_COMMANDS,  # noqa: F401
)
from chemsmart.agent.tui.chat_helpers import (
    _find_project_yaml_candidate_for_write as _find_project_yaml_candidate_for_write,  # noqa: F401,E501
    _default_interaction_mode,
    _latest_project_yaml_candidate as _latest_project_yaml_candidate,
    _load_tui_provider_config,
)
from chemsmart.agent.tui.state import TuiState
from chemsmart.agent.tui.screens.jobs_panel import JobsPanel
from chemsmart.agent.tui.screens.calculations import (
    CalculationMonitor,
)
from chemsmart.agent.tui.services.job_poller import (
    JobPollerError,
    JobPollerMixin,
    JobStatusUpdated,
    available_server_names,
    collect_job_snapshot,
    extract_run_result,
    format_jobs_table,
    queue_snapshot,
)
from chemsmart.agent.tui.services.log_tailer import LogTailer
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
from chemsmart.agent.tui.mixins.synthesis_presentation import (
    SynthesisPresentationMixin,
)
from chemsmart.agent.tui.mixins.job_interaction import JobInteractionMixin
from chemsmart.agent.tui.mixins.project_write import ProjectWriteMixin
from chemsmart.agent.tui.mixins.command_execution import CommandExecutionMixin
from chemsmart.agent.tui.mixins.slash_commands import SlashCommandsMixin
from chemsmart.agent.tui.mixins.approval_flow import ApprovalFlowMixin
from chemsmart.agent.tui.mixins.decision_log_events import DecisionEventMixin
from chemsmart.agent.tui.widgets.cells import (
    CalculationReceiptCell,
    JobStatusCell,
    MoleculeCell,
    RunResultCell,
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
    ResponseCopyOverlay,
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
    SynthesisPresentationMixin,
    JobInteractionMixin,
    ProjectWriteMixin,
    CommandExecutionMixin,
    SlashCommandsMixin,
    ApprovalFlowMixin,
    DecisionEventMixin,
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
