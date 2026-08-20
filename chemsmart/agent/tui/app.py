"""Textual application for plan, review, approval, and execution."""

from __future__ import annotations

from pathlib import Path
import re
import shlex
import threading
import time
from typing import TYPE_CHECKING, Callable

from rich.markdown import Markdown
from rich.panel import Panel
from rich.syntax import Syntax
from rich.table import Table
from rich.text import Text
from textual import events, on, work
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.widgets import Input, Static

if TYPE_CHECKING:
    from chemsmart.agent.execution import WorkflowExecutionReviewV1
    from chemsmart.agent.executor import WorkflowExecutionResultV1
    from chemsmart.agent.live_session import LiveAgentSessionResultV1

from chemsmart.agent._contracts import ContractError

from . import commands as command_registry
from .controller import AgentTuiController, AgentTuiPhase
from .monitor import (
    EventTailerV1,
    planning_feed_update,
    planning_row_key,
)
from .panels import phase_chip, wordmark
from .presentation import session_evidence_blocks
from .theme import CHEMSMART_THEME
from .transcript import ToolRow, TranscriptView


#: The one styling rule: running = text, finished = muted, failed = error.
_ROW_STYLE = {
    "running": "",
    "finished": "dim",
    "failed": "bold red",
    "note": "dim",
}


_RECEIPT_PLACEHOLDER = re.compile(
    r"<(?P<role>[a-z0-9-]+):sha256=[0-9a-f]{64}>"
)


def _human_cli_operation(argv: tuple[str, ...]) -> str:
    """Show ChemSmart input roles without exposing receipt bookkeeping."""

    return " ".join(
        _RECEIPT_PLACEHOLDER.sub(r"<\g<role>>", token) for token in argv
    )


class ScientificRequestInput(Input):
    """Single composer that preserves a pasted multi-line scientific task."""

    def _on_paste(self, event: events.Paste) -> None:
        text = " ".join(
            line.strip() for line in event.text.splitlines() if line.strip()
        )
        if text:
            self.insert_text_at_cursor(text)
        event.stop()


class ChemSmartAgentApp(App[None]):
    """Production terminal shell over the current Runtime V2 composition."""

    CSS_PATH = "styles.tcss"
    TITLE = "ChemSmart Agent"
    SUB_TITLE = "project YAML · compiled CLI · explicit approval"
    BINDINGS = [
        Binding("ctrl+c", "safe_quit", "Quit", priority=True),
        Binding("pageup", "scroll_transcript_up", "Scroll up", show=False),
        Binding(
            "pagedown", "scroll_transcript_down", "Scroll down", show=False
        ),
        Binding("ctrl+home", "scroll_transcript_top", "Top", show=False),
        Binding("ctrl+end", "scroll_transcript_bottom", "Bottom", show=False),
    ]

    def __init__(
        self, controller: AgentTuiController, *, plain: bool = False
    ) -> None:
        super().__init__()
        self.controller = controller
        self.plain = plain
        self._busy = False
        self._tail_rows: dict[str, ToolRow] = {}
        self._tail_stop = threading.Event()
        self._live_rows_seen = 0

    def compose(self) -> ComposeResult:
        with Vertical(id="agent-shell"):
            yield Static(wordmark(), id="wordmark")
            yield Static("", id="phase-bar")
            yield TranscriptView(id="transcript")
            yield ScientificRequestInput(
                placeholder=(
                    "Describe a calculation, or type /help for the approval chain"
                ),
                id="composer",
            )
            yield Static("", id="footer")

    def on_mount(self) -> None:
        self.register_theme(CHEMSMART_THEME)
        # Plain mode respects the terminal's own ANSI palette wholesale.
        self.theme = "ansi-dark" if self.plain else "chemsmart"
        self._write(
            Markdown(
                "## Ready\n\n"
                "ChemSmart turns scientific intent into canonical project YAML, "
                "compiled CLI operations, and a safe preview. A calculation "
                "runs only after the displayed scientific workflow is reviewed "
                "and the human enters /approve."
            )
        )
        self._sync_phase("Ready for a scientific request")
        self.controller.on_run_directory = self._announce_run_directory
        self.query_one("#composer", Input).focus()

    @on(Input.Submitted, "#composer")
    def submit(self, event: Input.Submitted) -> None:
        text = event.value.strip()
        if not text:
            return
        event.input.value = ""
        if self._busy:
            self.notify(
                "An operation is already in progress", severity="warning"
            )
            return
        if text.startswith("/"):
            self._dispatch_command(text)
            return
        try:
            normalized = self.controller.begin_planning(text)
        except ContractError as exc:
            self._operation_failed("Planning", exc)
            return
        self._write(Panel(Text(normalized), title="Scientific request"))
        self._busy = True
        self._sync_phase("Planning through the selected provider")
        self._run_plan(normalized)

    def _dispatch_command(self, text: str) -> None:
        try:
            arguments = shlex.split(text)
        except ValueError as exc:
            self._operation_failed("Command", exc)
            return
        command = arguments[0].lower()
        tail = arguments[1:]
        spec = command_registry.command_for(command)
        if spec is None:
            nearest = command_registry.suggest(command)
            message = f"Unknown command: {command}"
            if nearest:
                message += f". Did you mean {nearest}?"
            message += " /help lists every command."
            self._write(Panel(message, title="Command"))
            return
        handler = self._command_handlers()[spec.name]
        if tail and not spec.takes_argument:
            self._usage(f"{spec.slash} takes no arguments")
            return
        handler(tail)

    def _command_handlers(self) -> dict[str, Callable[[list[str]], None]]:
        return {
            "help": self._show_help,
            "status": self._show_status,
            "capabilities": self._show_capabilities,
            "approve": self._approve,
            "deny": self._deny,
            "revise": self._revise,
            "quit": self._quit,
        }

    def _show_help(self, tail: list[str]) -> None:
        self._write(
            Panel(command_registry.render_help(), title="Approval chain")
        )

    def _show_status(self, tail: list[str]) -> None:
        table = Table(title="Current session", show_header=False)
        table.add_column("Field", style="bold cyan")
        table.add_column("Value")
        table.add_row("phase", self.controller.phase.value)
        table.add_row("workspace", str(self.controller.config.workspace))
        table.add_row("task", self.controller.task or "not planned")
        prepared = self.controller.prepared_execution
        table.add_row(
            "workflow",
            prepared.scientific_plan.workflow_id
            if prepared is not None
            else "not awaiting approval",
        )
        table.add_row(
            "authority",
            "human /approve required"
            if prepared is not None
            else "none pending",
        )
        self._write(table)

    def _show_capabilities(self, tail: list[str]) -> None:
        from chemsmart.agent.capabilities import load_program_capabilities

        table = Table(title="Live Agent program surface")
        table.add_column("Program", style="bold cyan")
        table.add_column("Engine")
        table.add_column("Job")
        table.add_column("Preview")
        table.add_column("Approved execution")
        registry = load_program_capabilities()
        for capability in registry.programs:
            for item in capability.resolved_engine_job_capabilities:
                table.add_row(
                    capability.program,
                    item.engine,
                    item.jobtype,
                    "yes" if item.preview_supported else "no",
                    "yes" if item.execution_supported else "no",
                )
        self._write(table)
        self._write(
            "Execution support is a product path, not proof that this host is "
            "ready. Planning must still observe the selected program "
            "environment, complete safe preview, and present a ChemSmart "
            "workflow for human review."
        )

    def _approve(self, tail: list[str]) -> None:
        try:
            self.controller.begin_execution()
        except ContractError as exc:
            # The guard speaks before any approval banner can appear.
            self._operation_failed("Approval", exc)
            return
        self._write(
            Panel(
                "The human approved the displayed molecule, project YAML, "
                "CLI operations, DAG, observed execution environments and "
                "resource limits. The provider is "
                "disconnected before any engine launch.",
                title="Approve and run",
                border_style="red",
            )
        )
        self._busy = True
        self._sync_phase("Executing the approved ChemSmart DAG")
        self._run_execution()

    def _deny(self, tail: list[str]) -> None:
        self._decline_review("denied")

    def _revise(self, tail: list[str]) -> None:
        self._decline_review("revision requested")
        composer = self.query_one("#composer", Input)
        composer.value = self.controller.task
        composer.cursor_position = len(composer.value)
        composer.focus()

    def _decline_review(self, decision: str) -> None:
        self.controller.decline()
        self._write(
            Panel(
                f"Human decision: {decision}. No chemistry engine was launched.",
                title="Workflow not executed",
                border_style="yellow",
            )
        )
        self._sync_phase("Enter a revised scientific request when ready")

    def _quit(self, tail: list[str]) -> None:
        if self._busy:
            self.notify(
                "Wait for the current host operation to finish before exiting",
                severity="warning",
            )
            return
        self.exit()

    @work(thread=True, exclusive=True, group="agent-operation")
    def _run_plan(self, task: str) -> None:
        try:
            result = self.controller.run_planning(task)
        except Exception as exc:
            self.call_from_thread(self._operation_failed, "Planning", exc)
            return
        self.call_from_thread(self._plan_finished, result)

    @work(thread=True, exclusive=True, group="agent-operation")
    def _run_execution(self) -> None:
        try:
            result = self.controller.execute_begun()
        except Exception as exc:
            self.call_from_thread(self._operation_failed, "Execution", exc)
            return
        self.call_from_thread(self._execution_finished, result)

    # -- live planning feed --------------------------------------------------

    def _announce_run_directory(self, run_directory: Path) -> None:
        """Called by the live session on the worker thread, exactly once."""

        self.call_from_thread(self._start_planning_tail, Path(run_directory))

    def _start_planning_tail(self, run_directory: Path) -> None:
        self._tail_stop.clear()
        self._tail_rows.clear()
        self._run_event_tail(EventTailerV1(run_directory / "events.jsonl"))

    @work(thread=True, exclusive=True, group="event-tail")
    def _run_event_tail(self, tailer: EventTailerV1) -> None:
        while not self._tail_stop.is_set():
            events = tailer.poll()
            if events:
                self.call_from_thread(self._apply_planning_events, events)
            time.sleep(0.3)
        events = tailer.poll()
        if events:
            self.call_from_thread(self._apply_planning_events, events)

    def _apply_planning_events(self, events) -> None:
        transcript = self.query_one("#transcript", TranscriptView)
        for event in events:
            spec = planning_feed_update(event)
            if spec is None:
                continue
            text = Text(
                f"{spec.icon} {spec.text}", style=_ROW_STYLE[spec.state]
            )
            key = planning_row_key(event)
            kind = event.get("kind")
            if kind == "tool_started" and key:
                self._tail_rows[key] = transcript.add_row(
                    text, state=spec.state
                )
            elif (
                kind in {"tool_succeeded", "tool_failed"}
                and key in self._tail_rows
            ):
                transcript.settle_row(
                    self._tail_rows.pop(key), text, state=spec.state
                )
            else:
                transcript.add_row(text, state=spec.state)
            self._live_rows_seen += 1

    def _plan_finished(self, result: LiveAgentSessionResultV1) -> None:
        self._busy = False
        self._tail_stop.set()
        table = Table(title="Plan and safe-preview result")
        table.add_column("State", style="bold cyan")
        table.add_column("Value")
        table.add_row("terminal", result.terminal_state)
        table.add_row("session", result.session_id)
        table.add_row("successful tools", str(result.successful_tool_calls))
        table.add_row("failed tools", str(result.failed_tool_calls))
        self._write(table)
        for block in session_evidence_blocks(result):
            if (
                self._live_rows_seen
                and block.title != "Canonical project YAML"
            ):
                # The live feed already narrated every tool result; repeat
                # only the scientifically load-bearing project YAML. The
                # complete canonical payloads remain in the session's
                # append-only event stream.
                continue
            self._write(
                Panel(
                    Syntax(
                        block.text,
                        block.language,
                        word_wrap=True,
                        background_color="default",
                    ),
                    title=block.title,
                )
            )
        if result.final_text:
            self._write(Panel(Markdown(result.final_text), title="Agent"))
        if self.controller.phase is AgentTuiPhase.REQUEST_REVIEWED:
            prepared = self.controller.prepared_execution
            if prepared is None:  # pragma: no cover - phase owns the object
                self._operation_failed(
                    "Workflow review",
                    RuntimeError("prepared workflow is unavailable"),
                )
                return
            self._present_request(prepared)
            self._sync_phase(
                "Review shown; /approve runs once, /revise or /deny declines"
            )
        elif self.controller.phase is AgentTuiPhase.COMPLETE:
            self._sync_phase(
                "Plan or typed analysis complete; enter another scientific request"
            )
        else:
            self._sync_phase(
                "Planning stopped before approval readiness; inspect findings"
            )

    def _execution_finished(self, result: WorkflowExecutionResultV1) -> None:
        self._busy = False
        table = Table(title="ChemSmart execution result")
        table.add_column("Node")
        table.add_column("Program")
        table.add_column("Job")
        table.add_column("State")
        table.add_column("Validation")
        for node in result.nodes:
            table.add_row(
                node.node_id,
                node.program,
                node.jobtype,
                node.state,
                "validated" if node.validated else node.failure or "not validated",
            )
        self._write(table)
        self._write(
            Panel(
                f"status: {result.status}\n"
                f"provider calls: {result.provider_calls}\n"
                f"non-executable planned stages: "
                f"{', '.join(result.non_executable_node_ids) or 'none'}\n"
                f"run directory: {result.run_directory}",
                title="Execution",
                border_style=(
                    "green"
                    if result.status == "completed"
                    and not result.non_executable_node_ids
                    else "yellow"
                ),
            )
        )
        hint = (
            "Approved execution finished; non-executable planned stages remain "
            "unperformed. "
            if result.non_executable_node_ids
            else "Execution finished; "
        )
        self._sync_phase(
            hint
            + "place completed results in the task workspace and enter a new "
            "request for typed analysis"
        )

    def _present_request(self, review: WorkflowExecutionReviewV1) -> None:
        resources = review.execution_resources
        review_by_id = {item.node_id: item for item in review.node_reviews}
        deferred = set(review.non_executable_node_ids)
        overview = Table(title="ChemSmart workflow awaiting human approval")
        overview.add_column("Node", style="bold cyan")
        overview.add_column("Program / engine")
        overview.add_column("Stage")
        overview.add_column("Molecule")
        overview.add_column("Charge / multiplicity")
        overview.add_column("Execution")
        overview.add_column("Reason")
        for planned in review.scientific_plan.nodes:
            item = review_by_id.get(planned.node_id)
            if item is not None:
                identity = item.molecular_identity
                approved_names = identity.get("approved_names") or ()
                name = str(approved_names[0]) if approved_names else ""
                formula = str(identity.get("formula") or "unknown formula")
                atom_order = identity.get("atom_order") or ()
                atom_summary = "-".join(str(symbol) for symbol in atom_order)
                molecule = name or formula
                if atom_summary:
                    molecule += f" · {atom_summary}"
                program_engine = f"{item.program} / {item.engine}"
                stage = item.stage
                molecular_state = (
                    f"{identity.get('charge')} / "
                    f"{identity.get('multiplicity')}"
                )
                execution_state = "Executable"
                reason = "reviewed below"
            else:
                molecule = planned.project_role
                program_engine = f"{planned.program} / {planned.engine}"
                stage = planned.stage
                if (
                    planned.charge is not None
                    and planned.multiplicity is not None
                ):
                    molecular_state = (
                        f"{planned.charge} / {planned.multiplicity}"
                    )
                else:
                    source_ids = tuple(
                        edge.source_node_id
                        for edge in review.scientific_plan.edges
                        if edge.target_node_id == planned.node_id
                    )
                    source_review = next(
                        (
                            review_by_id[source_id]
                            for source_id in source_ids
                            if source_id in review_by_id
                        ),
                        None,
                    )
                    if source_review is not None:
                        source_identity = source_review.molecular_identity
                        molecular_state = (
                            f"{source_identity.get('charge')} / "
                            f"{source_identity.get('multiplicity')} "
                            f"(from {source_review.node_id})"
                        )
                    else:
                        molecular_state = "not separately reviewed"
                execution_state = (
                    "Deferred"
                    if planned.node_id in deferred
                    else "Not executable"
                )
                reason = planned.blocked_reason or "no reviewed command"
            overview.add_row(
                planned.node_id,
                program_engine,
                stage,
                molecule,
                molecular_state,
                execution_state,
                reason,
            )
        self._write(overview)
        environments = Table(title="Observed execution environments")
        environments.add_column("Node", style="bold cyan")
        environments.add_column("Status")
        environments.add_column("Target")
        environments.add_column("Version")
        environments.add_column("Observed by")
        for item in review.node_reviews:
            summary = item.environment_summary
            dependencies = {
                str(name): str(version)
                for name, version in summary.get("dependency_versions", ())
            }
            version = str(summary.get("observed_version") or "")
            if not version:
                version = dependencies.get(item.program, "")
            if not version and dependencies:
                version = ", ".join(
                    f"{name} {value}"
                    for name, value in sorted(dependencies.items())
                )
            environments.add_row(
                item.node_id,
                str(summary.get("status") or "unknown"),
                str(summary.get("target_kind") or "unknown"),
                version or "not reported",
                str(summary.get("observation_method") or "host probe"),
            )
        self._write(environments)
        self._write(
            Panel(
                f"cores: {resources.cores}\n"
                f"memory: {resources.memory_gb:g} GB\n"
                f"GPUs: {resources.gpu_count}\n"
                f"node timeout: {resources.node_timeout_seconds} s\n"
                f"engine-call limit: "
                f"{review.execution_envelope.get('max_engine_calls')}",
                title="Execution bounds",
            )
        )
        if review.scientific_plan.edges:
            edge_table = Table(title="Scientific data flow")
            edge_table.add_column("Producer")
            edge_table.add_column("Artifact")
            edge_table.add_column("Consumer")
            for edge in review.scientific_plan.edges:
                edge_table.add_row(
                    edge.source_node_id,
                    edge.artifact_class or edge.edge_kind,
                    edge.target_node_id,
                )
            self._write(edge_table)
        for item in review.node_reviews:
            command = _human_cli_operation(item.real_execution_argv)
            self._write(
                Panel(
                    Syntax(
                        item.project_settings_text,
                        "json",
                        word_wrap=True,
                        background_color="default",
                    ),
                    title=f"{item.node_id} · effective project settings",
                )
            )
            self._write(
                Panel(
                    Text(command),
                    title=f"{item.node_id} · ChemSmart CLI operation",
                )
            )
        executable_text = ", ".join(
            item.node_id for item in review.node_reviews
        )
        decision = (
            "Review the molecule/state, settings, CLI operations, DAG, "
            "observed execution environments and resource bounds above. "
            "Enter /approve once to execute the reviewed nodes: "
            f"{executable_text}."
        )
        if deferred:
            decision += (
                " Deferred stages remain unapproved and will not launch: "
                + ", ".join(sorted(deferred))
                + "."
            )
        decision += " A changed scientific request must be replanned."
        self._write(
            Panel(
                decision,
                title="Human decision",
                border_style="yellow",
            )
        )

    def _operation_failed(self, label: str, exc: Exception) -> None:
        self._busy = False
        self._tail_stop.set()
        self._write(
            Panel(
                str(exc) or type(exc).__name__,
                title=f"{label} stopped",
                border_style="red",
            )
        )
        self._sync_phase("Operation stopped; no unreviewed engine launch occurred")

    def _usage(self, message: str) -> None:
        self._write(
            Panel(message, title="Command usage", border_style="yellow")
        )

    def _sync_phase(self, hint: str) -> None:
        phase = self.controller.phase.value
        self.query_one("#phase-bar", Static).update(phase_chip(phase, hint))
        self.query_one("#footer", Static).update(
            command_registry.footer_hint(phase)
        )

    def _write(self, renderable) -> None:
        self.query_one("#transcript", TranscriptView).add_block(renderable)

    def action_safe_quit(self) -> None:
        self._quit([])

    def action_scroll_transcript_up(self) -> None:
        self.query_one("#transcript", TranscriptView).scroll_page_up()

    def action_scroll_transcript_down(self) -> None:
        self.query_one("#transcript", TranscriptView).scroll_page_down()

    def action_scroll_transcript_top(self) -> None:
        self.query_one("#transcript", TranscriptView).scroll_home()

    def action_scroll_transcript_bottom(self) -> None:
        self.query_one("#transcript", TranscriptView).scroll_end()


__all__ = ["ChemSmartAgentApp"]
