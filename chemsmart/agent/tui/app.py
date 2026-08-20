"""Textual application for plan, review, approval, and execution."""

from __future__ import annotations

from pathlib import Path
import shlex
import threading
import time
from typing import TYPE_CHECKING, Callable

from rich.markdown import Markdown
from rich.panel import Panel
from rich.syntax import Syntax
from rich.table import Table
from rich.text import Text
from functools import partial

from textual import events, on, work
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.command import DiscoveryHit, Hit, Provider
from textual.containers import Vertical
from textual.widgets import Input, Static

if TYPE_CHECKING:
    from chemsmart.agent.execution import WorkflowExecutionReviewV1
    from chemsmart.agent.executor import WorkflowExecutionResultV1
    from chemsmart.agent.live_session import LiveAgentSessionResultV1

from chemsmart.agent._contracts import ContractError

from . import commands as command_registry
from .controller import AgentTuiController, AgentTuiPhase
from .mermaid import render_workflow_mermaid
from .monitor import (
    EventTailerV1,
    execution_signal,
    planning_feed_update,
    planning_row_key,
)
from .panels import (
    DagPanel,
    JobsPanel,
    SPECTRUM_COLORS,
    phase_chip,
    spectrum_rule,
    wordmark,
)
from .report import looks_like_host_report, render_report_for_humans
from . import resume as resume_module
from .presentation import human_cli_operation, session_evidence_blocks
from .review import render_review_blocks, resolve_project_yaml_texts
from .runs import list_runs
from .theme import CHEMSMART_THEME
from .transcript import ToolRow, TranscriptView
from .voice import human_state


#: The one styling rule: running = text, finished = muted, failed = error.
_ROW_STYLE = {
    "running": "",
    "finished": "dim",
    "failed": "bold red",
    "note": "dim",
}


class ScientificRequestInput(Input):
    """Single composer that preserves a pasted multi-line scientific task."""

    def _on_paste(self, event: events.Paste) -> None:
        text = " ".join(
            line.strip() for line in event.text.splitlines() if line.strip()
        )
        if text:
            self.insert_text_at_cursor(text)
        event.stop()


class ChemSmartCommandProvider(Provider):
    """The ctrl+p palette is a query over the same command registry."""

    async def search(self, query: str):
        matcher = self.matcher(query)
        for spec in command_registry.COMMANDS:
            label = f"{spec.slash}  {spec.title}"
            score = matcher.match(label)
            if score > 0:
                yield Hit(
                    score,
                    matcher.highlight(label),
                    partial(self.app._dispatch_command, spec.slash),
                    help=spec.category,
                )

    async def discover(self):
        for spec in command_registry.COMMANDS:
            yield DiscoveryHit(
                f"{spec.slash}  {spec.title}",
                partial(self.app._dispatch_command, spec.slash),
                help=spec.category,
            )


class ChemSmartAgentApp(App[None]):
    """Production terminal shell over the current Runtime V2 composition."""

    COMMANDS = {ChemSmartCommandProvider}

    CSS_PATH = "styles.tcss"
    TITLE = "ChemSmart Agent"
    SUB_TITLE = "project YAML · compiled CLI · explicit approval"
    BINDINGS = [
        Binding("ctrl+c", "safe_quit", "Quit", priority=True),
        Binding(
            "escape", "interrupt_planning", "Cancel planning", priority=True,
            show=False,
        ),
        Binding("pageup", "scroll_transcript_up", "Scroll up", show=False),
        Binding(
            "pagedown", "scroll_transcript_down", "Scroll down", show=False
        ),
        Binding("ctrl+home", "scroll_transcript_top", "Top", show=False),
        Binding("ctrl+end", "scroll_transcript_bottom", "Bottom", show=False),
    ]

    def __init__(
        self,
        controller: AgentTuiController,
        *,
        plain: bool = False,
        resume: bool = False,
        discovery_notes: tuple[str, ...] = (),
    ) -> None:
        super().__init__()
        self.controller = controller
        self.plain = plain
        self.resume_requested = resume
        self.discovery_notes = tuple(discovery_notes)
        self._busy = False
        self._tail_rows: dict[str, ToolRow] = {}
        self._tail_stop = threading.Event()
        self._live_rows_seen = 0
        self._workflow_nodes: dict[str, dict] = {}
        self._jobs_timer = None
        self._last_report_path: Path | None = None
        self._esc_armed = False
        self._esc_disarm_timer = None
        self._pulse_index = 0
        self._phase_words = ("ready", "")
        self._pending_skill = ""
        self.exit_summary = {
            "planning_sessions": 0,
            "executions": 0,
            "report_paths": [],
        }

    def compose(self) -> ComposeResult:
        with Vertical(id="agent-shell"):
            header = Text()
            header.append_text(wordmark())
            header.append("\n")
            header.append_text(spectrum_rule(44, plain=self.plain))
            yield Static(header, id="wordmark")
            yield Static("", id="phase-bar")
            yield TranscriptView(id="transcript")
            yield DagPanel("", id="dag-panel")
            yield JobsPanel("", id="jobs-panel")
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
                "ChemSmart turns scientific intent into readable project YAML, "
                "compiled CLI operations, and a safe preview. A calculation "
                "runs only after the displayed scientific workflow is reviewed "
                "and the human enters /approve."
            )
        )
        for note in self.discovery_notes:
            self._write(Text(note, style="dim"))
        self._sync_phase("Ready for a scientific request")
        self.controller.on_run_directory = self._announce_run_directory
        self.set_interval(1.0, self._pulse_phase_dot)
        self.query_one("#composer", Input).focus()
        workspace = self.controller.config.workspace
        if self.resume_requested:
            self._resume([])
        elif resume_module.has_history(workspace):
            self._write(
                Text(
                    "Previous work exists in this workspace — /resume "
                    "restores it.",
                    style="dim",
                )
            )

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
        if self._pending_skill:
            text = (
                f"The user tagged domain skill '{self._pending_skill}' -- "
                "consult it before planning.\n\n" + text
            )
            self._pending_skill = ""
        try:
            normalized = self.controller.begin_planning(text)
        except ContractError as exc:
            self._operation_failed("Planning", exc)
            return
        self._write(Panel(Text(normalized), title="Scientific request"))
        self._busy = True
        self._sync_phase("Planning with the selected model")
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
            "resume": self._resume,
            "skills": self._show_skills,
            "skill": self._tag_skill,
            "export": self._export_transcript,
            "dag": self._toggle_dag,
            "report": self._show_report,
            "runs": self._show_runs,
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
            prepared = self.controller.begin_execution()
        except ContractError as exc:
            # The guard speaks before any approval banner can appear.
            self._operation_failed("Approval", exc)
            return
        run_directory = self.controller.execution_run_directory
        if run_directory is not None:
            run_directory.mkdir(parents=True, exist_ok=True)
            (run_directory / "workflow.mmd").write_text(
                render_workflow_mermaid(prepared), encoding="utf-8"
            )
        self._seed_workflow_nodes(prepared)
        jobs = self.query_one("#jobs-panel", JobsPanel)
        jobs.display = True
        jobs.refresh_from(self._workflow_nodes)
        self._refresh_dag()
        self._jobs_timer = self.set_interval(1.0, self._refresh_live_panels)
        if run_directory is not None:
            self._tail_stop.clear()
            self._run_event_tail(
                EventTailerV1(run_directory / "events.jsonl"),
                mode="execution",
            )
        self._write(
            Panel(
                "The human approved the displayed molecule, project YAML, "
                "CLI operations, DAG, observed execution environments and "
                "resource limits. The model is "
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
    def _run_event_tail(
        self, tailer: EventTailerV1, mode: str = "planning"
    ) -> None:
        apply = (
            self._apply_planning_events
            if mode == "planning"
            else self._apply_execution_events
        )
        while not self._tail_stop.is_set():
            events = tailer.poll()
            if events:
                self.call_from_thread(apply, events)
            time.sleep(0.3)
        events = tailer.poll()
        if events:
            self.call_from_thread(apply, events)

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

    # -- execution surface ---------------------------------------------------

    def _seed_workflow_nodes(self, review) -> None:
        review_by_id = {item.node_id: item for item in review.node_reviews}
        deferred = set(review.non_executable_node_ids)
        nodes: dict[str, dict] = {}
        for planned in review.scientific_plan.nodes:
            item = review_by_id.get(planned.node_id)
            if item is not None:
                formula = str(
                    item.molecular_identity.get("formula") or ""
                )
                label = f"{item.program} {item.stage}"
                if formula:
                    label += f" · {formula}"
                nodes[planned.node_id] = {
                    "kind": "calc",
                    "label": label,
                    "argv": human_cli_operation(item.real_execution_argv),
                    "state": "queued",
                    "detail": "",
                    "started": None,
                }
            else:
                nodes[planned.node_id] = {
                    "kind": "calc",
                    "label": f"{planned.program} {planned.stage}",
                    "argv": "",
                    "state": "deferred",
                    "detail": planned.blocked_reason
                    or (
                        "deferred"
                        if planned.node_id in deferred
                        else "not executable"
                    ),
                    "started": None,
                }
        plan = review.scientific_toolchain_plan
        if plan is not None:
            for node in plan.analysis_nodes:
                blocked = node.support_state == "blocked_unsupported"
                nodes[node.node_id] = {
                    "kind": "analysis",
                    "label": node.analysis_kind,
                    "argv": "",
                    "state": "blocked_unsupported" if blocked else "queued",
                    "detail": node.blocked_reason,
                    "started": None,
                }
        self._workflow_nodes = nodes

    def _apply_execution_events(self, events) -> None:
        for event in events:
            signal = execution_signal(event)
            if signal is None:
                continue
            node = self._workflow_nodes.get(signal.node_id)
            if signal.kind == "node_launched" and node is not None:
                node["state"] = "running"
                node["started"] = time.monotonic()
            elif signal.kind == "engine_done" and node is not None:
                node["state"] = (
                    "validated"
                    if signal.detail == "validated"
                    else "engine_complete"
                )
                node["detail"] = signal.detail
            elif signal.kind == "node_state" and node is not None:
                node["state"] = signal.state
                node["detail"] = ""
            elif signal.kind == "analysis_settled" and node is not None:
                node["state"] = signal.state
                node["detail"] = signal.detail
            elif signal.kind == "report_rendered":
                self._write(
                    Text(
                        "Σ completed-analysis report rendered · /report "
                        "opens it",
                        style="dim",
                    )
                )
        self._refresh_live_panels()

    def _refresh_live_panels(self) -> None:
        jobs = self.query_one("#jobs-panel", JobsPanel)
        if jobs.display:
            jobs.refresh_from(self._workflow_nodes)
        self._refresh_dag()

    def _refresh_dag(self) -> None:
        dag = self.query_one("#dag-panel", DagPanel)
        if dag.display:
            dag.refresh_from(self._workflow_nodes)

    def _resume(self, tail: list[str]) -> None:
        if self._busy:
            self.notify(
                "Wait for the current host operation to finish first",
                severity="warning",
            )
            return
        workspace = self.controller.config.workspace
        story = resume_module.workspace_story(workspace)
        greeting = Table(
            title="Resumed workspace", show_header=False, box=None
        )
        greeting.add_column("Field", style="bold cyan")
        greeting.add_column("Value")
        greeting.add_row("workspace", str(workspace))
        if story.story is not None:
            greeting.add_row("last task", story.story.task or "(not recorded)")
            greeting.add_row("ended", story.story.ended)
        else:
            greeting.add_row("last task", "no planning session found")
        self._write(Panel(greeting, title="Resume"))
        if story.runs:
            runs_table = Table(title="Runs in this workspace (newest first)")
            runs_table.add_column("#", style="bold cyan")
            runs_table.add_column("Workflow")
            runs_table.add_column("Kind")
            runs_table.add_column("State")
            for index, row in enumerate(story.runs, start=1):
                runs_table.add_row(
                    str(index),
                    row.workflow_id or row.name,
                    row.kind,
                    row.state + (" · /report " + str(index) if row.report_path else ""),
                )
            self._write(runs_table)
        if story.story is not None and story.story.final_prose:
            prose = story.story.final_prose
            if looks_like_host_report(prose):
                self._write(
                    Panel(
                        render_report_for_humans(prose),
                        title="Analysis results",
                    )
                )
            else:
                self._write(Panel(Markdown(prose), title="Agent (previous session)"))
        pending = story.pending
        if pending is None:
            self._sync_phase(
                "Workspace restored; enter a new scientific request"
            )
            return
        if pending.status == "approved_unexecuted":
            self._write(
                Panel(
                    "A review in this workspace was approved but its run "
                    "never started. The recorded approval is one-shot and "
                    "still unconsumed; launch it exactly as approved:\n\n"
                    f"{pending.recovery}",
                    title="Approved, not yet run",
                    border_style="yellow",
                )
            )
            self._sync_phase("Workspace restored; the approved run awaits")
            return
        self._restore_pending_review(pending)

    def _restore_pending_review(self, pending) -> None:
        from chemsmart.agent.live_session import (
            inspect_workflow_execution_replay,
            load_workflow_execution_review,
        )

        workspace = self.controller.config.workspace
        try:
            report = inspect_workflow_execution_replay(
                review_file=pending.path,
                workspace=workspace,
            )
        except ContractError as exc:
            self._write(
                Panel(
                    f"A stored review exists but cannot be re-presented: "
                    f"{exc}",
                    title="Pending review",
                    border_style="yellow",
                )
            )
            return
        if report["missing_approved_artifacts"]:
            missing = ", ".join(
                str(item).split(":")[0]
                for item in report["missing_approved_artifacts"]
            )
            self._write(
                Panel(
                    "A reviewed workflow is pending, but files it was "
                    f"approved against are no longer in place ({missing}). "
                    "No decision is offered on a run that must fail; plan "
                    "the request again.",
                    title="Pending review cannot be offered",
                    border_style="yellow",
                )
            )
            return
        try:
            review = load_workflow_execution_review(pending.path)
            self.controller.restore_prepared_execution(review)
        except ContractError as exc:
            self._operation_failed("Resume", exc)
            return
        self._write(
            Panel(
                "The previous session prepared this workflow and no decision "
                "was recorded. It is re-presented below for one fresh "
                "decision.",
                title="Pending review restored",
                border_style="yellow",
            )
        )
        self._present_request(review)
        self._sync_phase(
            "Review shown; /approve runs once, /revise or /deny declines"
        )

    def _show_skills(self, tail: list[str]) -> None:
        from chemsmart.agent.skills import available_skill_ids

        table = Table(title="Consultable domain skills")
        table.add_column("Skill", style="bold cyan")
        for skill_id in available_skill_ids():
            table.add_row(skill_id)
        self._write(table)
        self._write(
            Text(
                "/skill <id> tags your next request; the session still "
                "consults it through its own tool, so the record of what "
                "was consulted is kept.",
                style="dim",
            )
        )

    def _tag_skill(self, tail: list[str]) -> None:
        from chemsmart.agent.skills import available_skill_ids

        if len(tail) != 1:
            self._usage("/skill takes exactly one skill id; /skills lists them")
            return
        skill_id = tail[0]
        known = tuple(available_skill_ids())
        if skill_id not in known:
            self._write(
                Panel(
                    f"Unknown skill '{skill_id}'. Available: "
                    + ", ".join(known),
                    title="Skill",
                    border_style="yellow",
                )
            )
            return
        self._pending_skill = skill_id
        self._write(
            Panel(
                f"The next scientific request will carry a visible tag "
                f"asking the session to consult '{skill_id}' before "
                "planning.",
                title="Skill tagged",
            )
        )

    def _export_transcript(self, tail: list[str]) -> None:
        transcript = self.query_one("#transcript", TranscriptView)
        target = (
            self.controller.config.workspace
            / time.strftime("chemsmart-transcript-%Y%m%d-%H%M%S.txt")
        )
        target.write_text(
            transcript.recorder.export_text(), encoding="utf-8"
        )
        self._write(
            Panel(f"Transcript saved to {target}", title="Export")
        )

    def _toggle_dag(self, tail: list[str]) -> None:
        dag = self.query_one("#dag-panel", DagPanel)
        if dag.display:
            dag.display = False
            return
        if not self._workflow_nodes:
            review = self.controller.prepared_execution
            if review is not None:
                self._seed_workflow_nodes(review)
        if not self._workflow_nodes:
            self._write(
                Panel(
                    "No workflow to display yet; plan one first.",
                    title="Workflow",
                )
            )
            return
        dag.display = True
        dag.refresh_from(self._workflow_nodes)

    def _show_runs(self, tail: list[str]) -> None:
        summaries = list_runs(self.controller.config.workspace)
        if not summaries:
            self._write(
                Panel(
                    "No executions or replays exist under this workspace.",
                    title="Runs",
                )
            )
            return
        table = Table(title="Workspace runs (newest first)")
        table.add_column("#", style="bold cyan")
        table.add_column("Run")
        table.add_column("Kind")
        table.add_column("Terminal state")
        table.add_column("Report")
        for index, summary in enumerate(summaries, start=1):
            table.add_row(
                str(index),
                summary.name,
                summary.kind,
                summary.terminal_state,
                "/report " + str(index) if summary.report_path else "—",
            )
        self._write(table)

    def _show_report(self, tail: list[str]) -> None:
        path: Path | None = None
        if tail:
            summaries = list_runs(self.controller.config.workspace)
            try:
                chosen = summaries[int(tail[0]) - 1]
            except (ValueError, IndexError):
                self._usage(
                    "/report takes a run number from /runs, or no argument "
                    "for the latest report"
                )
                return
            path = chosen.report_path
        else:
            path = self._last_report_path
            if path is None:
                path = next(
                    (
                        summary.report_path
                        for summary in list_runs(
                            self.controller.config.workspace
                        )
                        if summary.report_path is not None
                    ),
                    None,
                )
        if path is None or not path.exists():
            self._write(
                Panel(
                    "No completed-analysis report exists yet. A report is "
                    "rendered when an approved run's analysis chain "
                    "completes.",
                    title="Report",
                )
            )
            return
        self._write(spectrum_rule(44, plain=self.plain))
        self._write(
            Panel(
                render_report_for_humans(path.read_text(encoding="utf-8")),
                title=f"Analysis results · {path}",
            )
        )

    def _plan_finished(self, result: LiveAgentSessionResultV1) -> None:
        self._busy = False
        self._tail_stop.set()
        self.exit_summary["planning_sessions"] += 1
        table = Table(title="Plan and safe-preview result")
        table.add_column("State", style="bold cyan")
        table.add_column("Value")
        table.add_row("outcome", human_state(result.terminal_state))
        steps = f"{result.successful_tool_calls} steps"
        if result.failed_tool_calls:
            steps += f", {result.failed_tool_calls} refused"
        table.add_row("steps", steps)
        self._write(table)
        # The live feed narrated every tool result; what bears repeating is
        # the scientifically load-bearing project YAML. The complete record
        # stays in the session's append-only run evidence.
        yaml_blocks = [
            block
            for block in session_evidence_blocks(result)
            if block.language == "yaml"
        ]
        for block in yaml_blocks:
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
        if not self._live_rows_seen and not yaml_blocks:
            self._write(
                Text(
                    "The live feed captured no tool activity; the complete "
                    "record is in the session's run evidence.",
                    style="dim",
                )
            )
        if result.final_text:
            if looks_like_host_report(result.final_text):
                self._write(
                    Panel(
                        render_report_for_humans(result.final_text),
                        title="Analysis results",
                    )
                )
            else:
                self._write(
                    Panel(Markdown(result.final_text), title="Agent")
                )
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
        elif result.terminal_state == "cancelled":
            self._sync_phase(
                "Planning cancelled at a provider-turn boundary; the session "
                "remains an observed run"
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
        self._tail_stop.set()
        if self._jobs_timer is not None:
            self._jobs_timer.stop()
            self._jobs_timer = None
        self.query_one("#jobs-panel", JobsPanel).display = False
        self._refresh_dag()
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
        if result.analysis_nodes:
            analysis = Table(title="Approved analysis chain")
            analysis.add_column("Node", style="bold cyan")
            analysis.add_column("Kind")
            analysis.add_column("State")
            analysis.add_column("Reason")
            for node in result.analysis_nodes:
                analysis.add_row(
                    node.node_id,
                    node.analysis_kind,
                    node.state,
                    node.reason or "—",
                )
            self._write(analysis)
        if result.analysis_status:
            self._write(
                Text(
                    f"analysis: {result.analysis_status}",
                    style=(
                        "green"
                        if result.analysis_status == "completed"
                        else "yellow"
                    ),
                )
            )
        self.exit_summary["executions"] += 1
        if result.analysis_report_path:
            report_path = Path(result.analysis_report_path)
            self._last_report_path = report_path
            self.exit_summary["report_paths"].append(str(report_path))
            if report_path.exists():
                self._write(spectrum_rule(44, plain=self.plain))
                self._write(
                    Panel(
                        render_report_for_humans(
                            report_path.read_text(encoding="utf-8")
                        ),
                        title=f"Analysis results · {report_path}",
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
            + "interpretation and the recorded decision remain a session "
            "act; enter a new request when ready"
        )

    def _present_request(self, review: WorkflowExecutionReviewV1) -> None:
        transcript = ()
        if self.controller.plan_result is not None:
            transcript = self.controller.plan_result.public_transcript
        yaml_texts = resolve_project_yaml_texts(
            review,
            public_transcript=transcript,
            workspace=self.controller.config.workspace,
        )
        for block in render_review_blocks(review, project_yaml=yaml_texts):
            self._write(block)

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
        self._phase_words = (human_state(phase), hint)
        self.query_one("#phase-bar", Static).update(
            phase_chip(self._phase_words[0], hint)
        )
        self.query_one("#footer", Static).update(
            command_registry.footer_hint(phase)
        )

    def _pulse_phase_dot(self) -> None:
        """While busy, the phase dot breathes through the spectrum colors."""

        if not self._busy or self.plain:
            return
        self._pulse_index = (self._pulse_index + 1) % len(SPECTRUM_COLORS)
        words, hint = self._phase_words
        self.query_one("#phase-bar", Static).update(
            phase_chip(
                words,
                hint,
                dot_style=f"bold {SPECTRUM_COLORS[self._pulse_index]}",
            )
        )

    def _write(self, renderable) -> None:
        self.query_one("#transcript", TranscriptView).add_block(renderable)

    def action_safe_quit(self) -> None:
        self._quit([])

    def action_interrupt_planning(self) -> None:
        """Double-esc cancels planning at the next provider-turn boundary.

        Engines are never killable from the interface; the event is consulted
        only by the planning loop between provider turns.
        """

        if not (
            self._busy
            and self.controller.phase is AgentTuiPhase.PLANNING
        ):
            self._disarm_interrupt()
            return
        if not self._esc_armed:
            self._esc_armed = True
            self.query_one("#footer", Static).update(
                Text(
                    "esc again to cancel planning (disarms in 5 s)",
                    style="bold yellow",
                )
            )
            self._esc_disarm_timer = self.set_timer(
                5.0, self._disarm_interrupt
            )
            return
        self.controller.cancel_planning.set()
        self._esc_armed = False
        self.query_one("#footer", Static).update(
            Text(
                "cancelling at the next provider-turn boundary…",
                style="yellow",
            )
        )

    def _disarm_interrupt(self) -> None:
        if self._esc_disarm_timer is not None:
            self._esc_disarm_timer.stop()
            self._esc_disarm_timer = None
        if self._esc_armed:
            self._esc_armed = False
            self.query_one("#footer", Static).update(
                command_registry.footer_hint(self.controller.phase.value)
            )

    def action_scroll_transcript_up(self) -> None:
        self.query_one("#transcript", TranscriptView).scroll_page_up()

    def action_scroll_transcript_down(self) -> None:
        self.query_one("#transcript", TranscriptView).scroll_page_down()

    def action_scroll_transcript_top(self) -> None:
        self.query_one("#transcript", TranscriptView).scroll_home()

    def action_scroll_transcript_bottom(self) -> None:
        self.query_one("#transcript", TranscriptView).scroll_end()


__all__ = ["ChemSmartAgentApp"]
