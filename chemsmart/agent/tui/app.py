"""Textual application for plan, review, approval binding, and execution."""

from __future__ import annotations

import re
import shlex
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
from textual.widgets import Input, RichLog, Static

if TYPE_CHECKING:
    from chemsmart.agent.execution import WorkflowExecutionReviewV1
    from chemsmart.agent.executor import WorkflowExecutionResultV1
    from chemsmart.agent.live_session import LiveAgentSessionResultV1

from .controller import AgentTuiController, AgentTuiPhase
from .presentation import session_evidence_blocks


_HELP = """\
Enter a scientific request to create a project-YAML/CLI plan and safe preview.

Commands:
  /capabilities             show the live Agent program/engine/job surface
  /approve                  approve the displayed ChemSmart workflow and run
  /deny                     decline the displayed workflow
  /revise                   decline it and enter a revised scientific request
  /status                   show the current task and evidence bindings
  /help                     show this guide
  /quit                     exit

The provider/runtime cannot grant itself authority. The single /approve action
ends planning authority and starts the displayed YAML/CLI DAG through the
provider-free ChemSmart executor. Internal receipts remain provenance; no hash
or approval-file token is required from the human.
"""


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
        Binding("ctrl+l", "refresh", "Refresh", show=False),
    ]

    def __init__(
        self, controller: AgentTuiController, *, plain: bool = False
    ) -> None:
        super().__init__()
        self.controller = controller
        self.plain = plain
        self._busy = False

    def compose(self) -> ComposeResult:
        with Vertical(id="agent-shell"):
            yield Static(self._wordmark(), id="wordmark")
            yield Static("", id="phase-bar")
            yield RichLog(
                id="transcript",
                wrap=True,
                markup=False,
                highlight=False,
                max_lines=4000,
            )
            yield ScientificRequestInput(
                placeholder=(
                    "Describe a calculation, or type /help for the approval chain"
                ),
                id="composer",
            )
            yield Static(
                "Enter plan · review YAML/CLI DAG · /approve once to run",
                id="footer",
            )

    def on_mount(self) -> None:
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
        self._write(Panel(Text(text), title="Scientific request"))
        self._busy = True
        self._sync_phase("Planning through the selected provider")
        self._run_plan(text)

    def _dispatch_command(self, text: str) -> None:
        try:
            arguments = shlex.split(text)
        except ValueError as exc:
            self._operation_failed("Command", exc)
            return
        command = arguments[0].lower()
        tail = arguments[1:]
        handlers: dict[str, Callable[[list[str]], None]] = {
            "/help": self._show_help,
            "/status": self._show_status,
            "/capabilities": self._show_capabilities,
            "/approve": self._approve,
            "/deny": self._deny,
            "/revise": self._revise,
            "/quit": self._quit,
            "/exit": self._quit,
        }
        handler = handlers.get(command)
        if handler is None:
            self._write(Panel(f"Unknown command: {command}", title="Command"))
            return
        handler(tail)

    def _show_help(self, tail: list[str]) -> None:
        if tail:
            self._usage("/help takes no arguments")
            return
        self._write(Panel(_HELP.rstrip(), title="Approval chain"))

    def _show_status(self, tail: list[str]) -> None:
        if tail:
            self._usage("/status takes no arguments")
            return
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
        if tail:
            self._usage("/capabilities takes no arguments")
            return
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
        if tail:
            self._usage("/approve takes no arguments")
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
        if tail:
            self._usage("/deny takes no arguments")
            return
        self._decline_review("denied")

    def _revise(self, tail: list[str]) -> None:
        if tail:
            self._usage("/revise takes no arguments")
            return
        self._decline_review("revision requested")

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
        if tail:
            self._usage("/quit takes no arguments")
            return
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
            result = self.controller.plan(task)
        except Exception as exc:
            self.call_from_thread(self._operation_failed, "Planning", exc)
            return
        self.call_from_thread(self._plan_finished, result)

    @work(thread=True, exclusive=True, group="agent-operation")
    def _run_execution(self) -> None:
        try:
            result = self.controller.approve_and_execute()
        except Exception as exc:
            self.call_from_thread(self._operation_failed, "Execution", exc)
            return
        self.call_from_thread(self._execution_finished, result)

    def _plan_finished(self, result: LiveAgentSessionResultV1) -> None:
        self._busy = False
        table = Table(title="Plan and safe-preview result")
        table.add_column("State", style="bold cyan")
        table.add_column("Value")
        table.add_row("terminal", result.terminal_state)
        table.add_row("session", result.session_id)
        table.add_row("successful tools", str(result.successful_tool_calls))
        table.add_row("failed tools", str(result.failed_tool_calls))
        self._write(table)
        for block in session_evidence_blocks(result):
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
                f"run directory: {result.run_directory}",
                title="Execution",
                border_style=(
                    "green" if result.status == "completed" else "yellow"
                ),
            )
        )
        self._sync_phase(
            "Execution finished; place completed results in the task "
            "workspace and enter a new request for typed analysis"
        )

    def _present_request(self, review: WorkflowExecutionReviewV1) -> None:
        resources = review.execution_resources
        overview = Table(title="ChemSmart workflow awaiting human approval")
        overview.add_column("Node", style="bold cyan")
        overview.add_column("Program / engine")
        overview.add_column("Stage")
        overview.add_column("Molecule")
        overview.add_column("Charge / multiplicity")
        for item in review.node_reviews:
            identity = item.molecular_identity
            approved_names = identity.get("approved_names") or ()
            name = str(approved_names[0]) if approved_names else ""
            formula = str(identity.get("formula") or "unknown formula")
            atom_order = identity.get("atom_order") or ()
            atom_summary = "-".join(str(symbol) for symbol in atom_order)
            molecule = name or formula
            if atom_summary:
                molecule += f" · {atom_summary}"
            overview.add_row(
                item.node_id,
                f"{item.program} / {item.engine}",
                item.stage,
                molecule,
                f"{identity.get('charge')} / {identity.get('multiplicity')}",
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
        self._write(
            Panel(
                "Review the molecule/state, settings, CLI operations, DAG, "
                "observed execution environments and resource bounds above. "
                "Enter /approve once to execute this whole workflow. A changed "
                "scientific request must be replanned.",
                title="Human decision",
                border_style="yellow",
            )
        )

    def _operation_failed(self, label: str, exc: Exception) -> None:
        self._busy = False
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
        self.query_one("#phase-bar", Static).update(
            Text.assemble(
                ("● ", "bold cyan"),
                (phase, "bold"),
                ("  ·  ", "dim"),
                (hint, "dim"),
            )
        )

    def _write(self, renderable) -> None:
        self.query_one("#transcript", RichLog).write(renderable)

    def _wordmark(self) -> Text:
        return Text.assemble(
            ("CHEMSMART", "bold cyan"),
            ("  AGENT", "bold"),
            ("  Runtime V2", "dim"),
        )

    def action_safe_quit(self) -> None:
        self._quit([])


__all__ = ["ChemSmartAgentApp"]
