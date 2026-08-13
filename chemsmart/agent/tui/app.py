"""Textual application for plan, review, approval binding, and execution."""

from __future__ import annotations

import shlex
from typing import TYPE_CHECKING, Callable

from rich.markdown import Markdown
from rich.panel import Panel
from rich.syntax import Syntax
from rich.table import Table
from rich.text import Text
from textual import on, work
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.widgets import Input, RichLog, Static

if TYPE_CHECKING:
    from chemsmart.agent.executor import WorkflowExecutionResultV1
    from chemsmart.agent.live_session import LiveAgentSessionResultV1

from .approval import (
    ReviewedExecutionApprovalV1,
    ReviewedExecutionReviewV1,
)
from .controller import AgentTuiController, AgentTuiPhase
from .presentation import session_evidence_blocks


_HELP = """\
Enter a scientific request to create a project-YAML/CLI plan and safe preview.

Commands:
  /capabilities             show the live Agent program/engine/job surface
  /request PATH             review an inert approval-request artifact
  /approve ACTOR OUTPUT     approve this exact review once and write a bundle
  /approval PATH            bind an approval produced by CLI or another UI
  /execute SHA256 RUN_DIR   execute only after retyping its full file digest
  /deny ACTOR               deny an unapproved review
  /revise ACTOR             request revision of an unapproved review
  /quit-review ACTOR        record quit without exiting the TUI
  /status                   show the current task and evidence bindings
  /help                     show this guide
  /quit                     exit

Only these explicit human commands can resolve a review. The provider/runtime
cannot grant itself authority. Approval and execution remain separate acts;
execution uses the deterministic executor and makes no provider call.
"""


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
            yield Input(
                placeholder=(
                    "Describe a calculation, or type /help for the approval chain"
                ),
                id="composer",
            )
            yield Static(
                "Enter plan · /request review · /approval bind · /execute confirm",
                id="footer",
            )

    def on_mount(self) -> None:
        self._write(
            Markdown(
                "## Ready\n\n"
                "ChemSmart turns scientific intent into canonical project YAML, "
                "compiled CLI operations, and a safe preview. A calculation "
                "runs only after an exact request is reviewed, an external "
                "user approval is bound, and its SHA-256 is confirmed."
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
            "/request": self._review_request,
            "/approve": self._approve,
            "/approval": self._bind_approval,
            "/execute": self._execute,
            "/deny": self._deny,
            "/revise": self._revise,
            "/quit-review": self._quit_review,
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
        result = self.controller.plan_result
        table.add_row(
            "task spec",
            result.task_spec_sha256 if result is not None else "not available",
        )
        request = self.controller.reviewed_review
        table.add_row(
            "review",
            request.review.review_sha256 if request else "not reviewed",
        )
        approval = self.controller.reviewed_approval
        table.add_row(
            "approval file",
            approval.artifact_sha256 if approval else "not bound",
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
            "environment, complete safe preview, and produce an exact review."
        )

    def _review_request(self, tail: list[str]) -> None:
        if len(tail) != 1:
            self._usage("usage: /request PATH")
            return
        try:
            reviewed = self.controller.review_request(tail[0])
        except Exception as exc:
            self._operation_failed("Approval request", exc)
            return
        self._present_request(reviewed)
        self._sync_phase("Request reviewed; create approval outside the TUI")

    def _bind_approval(self, tail: list[str]) -> None:
        if len(tail) != 1:
            self._usage("usage: /approval PATH")
            return
        try:
            reviewed = self.controller.bind_approval(tail[0])
        except Exception as exc:
            self._operation_failed("Execution approval", exc)
            return
        self._present_approval(reviewed)
        self._sync_phase(
            "Approval bound; retype full file SHA-256 with /execute to run"
        )

    def _approve(self, tail: list[str]) -> None:
        if len(tail) != 2:
            self._usage("usage: /approve ACTOR OUTPUT_PATH")
            return
        actor, output = tail
        try:
            resolution, reviewed = self.controller.resolve_review(
                decision="approve",
                actor=actor,
                output_file=output,
            )
        except Exception as exc:
            self._operation_failed("Approval", exc)
            return
        if reviewed is None:  # pragma: no cover - contract narrows this
            self._operation_failed(
                "Approval", RuntimeError("approval bundle was not created")
            )
            return
        self._write(
            Panel(
                f"actor: {resolution.actor}\n"
                f"decision: {resolution.decision}\n"
                f"resolution SHA-256: {resolution.resolution_sha256}",
                title="Human decision recorded",
                border_style="red",
            )
        )
        self._present_approval(reviewed)
        self._sync_phase(
            "One-shot approval created; retype its full file SHA-256 to execute"
        )

    def _execute(self, tail: list[str]) -> None:
        if len(tail) != 2:
            self._usage("usage: /execute APPROVAL_SHA256 RUN_DIRECTORY")
            return
        self._busy = True
        self._sync_phase("Executing exact approved DAG; provider disconnected")
        self._run_execution(tail[0], tail[1])

    def _deny(self, tail: list[str]) -> None:
        if len(tail) != 1:
            self._usage("usage: /deny ACTOR")
            return
        self._resolve_without_grant("deny", tail[0])

    def _revise(self, tail: list[str]) -> None:
        if len(tail) != 1:
            self._usage("usage: /revise ACTOR")
            return
        self._resolve_without_grant("revise", tail[0])

    def _quit_review(self, tail: list[str]) -> None:
        if len(tail) != 1:
            self._usage("usage: /quit-review ACTOR")
            return
        self._resolve_without_grant("quit", tail[0])

    def _resolve_without_grant(self, decision: str, actor: str) -> None:
        if self.controller.reviewed_approval is not None:
            self._operation_failed(
                "Review decision",
                RuntimeError(
                    "an approval bundle has already been written and cannot "
                    "be revoked by this TUI; do not execute it, or finish a "
                    "new revised review"
                ),
            )
            return
        try:
            resolution, _reviewed = self.controller.resolve_review(
                decision=decision,
                actor=actor,
            )
        except Exception as exc:
            self._operation_failed("Review decision", exc)
            return
        self._write(
            Panel(
                f"actor: {resolution.actor}\n"
                f"decision: {resolution.decision}\n"
                f"resolution SHA-256: {resolution.resolution_sha256}\n"
                "No execution approval was created.",
                title="Human decision recorded",
                border_style="yellow",
            )
        )
        self._sync_phase("Decision recorded; this review created no authority")

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
    def _run_execution(self, digest: str, run_directory: str) -> None:
        try:
            result = self.controller.execute(
                confirmation_digest=digest,
                run_directory=run_directory,
            )
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
        table.add_row("task spec", result.task_spec_sha256)
        table.add_row("successful tools", str(result.successful_tool_calls))
        table.add_row("failed tools", str(result.failed_tool_calls))
        table.add_row("event head", result.event_stream_head_sha256)
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
        if self.controller.phase is AgentTuiPhase.PREVIEW_READY:
            try:
                reviewed = self.controller.review_request(
                    self.controller.config.review_file
                )
            except Exception as exc:
                self._operation_failed("Execution review", exc)
                return
            self._present_request(reviewed)
            self._sync_phase(
                "Review shown; /approve, /deny, /revise, or /quit-review"
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
        table = Table(title="Deterministic execution result")
        table.add_column("Node")
        table.add_column("Program")
        table.add_column("Job")
        table.add_column("State")
        table.add_column("Receipt")
        for node in result.nodes:
            table.add_row(
                node.node_id,
                node.program,
                node.jobtype,
                node.state,
                node.execution_receipt_sha256 or node.failure or "-",
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

    def _present_request(self, reviewed: ReviewedExecutionReviewV1) -> None:
        review = reviewed.review
        self._write(
            Panel(
                Syntax(
                    reviewed.canonical_text,
                    "json",
                    word_wrap=True,
                    background_color="default",
                ),
                title="Unapproved execution review (inert)",
                border_style="yellow",
            )
        )
        self._write(
            f"review SHA-256: {review.review_sha256}\n"
            f"file SHA-256: {reviewed.artifact_sha256}"
        )

    def _present_approval(self, reviewed: ReviewedExecutionApprovalV1) -> None:
        approval = reviewed.bundle.workflow_approval
        self._write(
            Panel(
                Syntax(
                    reviewed.canonical_text,
                    "json",
                    word_wrap=True,
                    background_color="default",
                ),
                title="User-created execution approval",
                border_style="red",
            )
        )
        self._write(
            f"approval bundle SHA-256: {reviewed.bundle.bundle_sha256}\n"
            f"approval record SHA-256: {approval.approval_sha256}\n"
            f"approval file SHA-256: {reviewed.artifact_sha256}\n"
            f"confirmation command: /execute "
            f"{reviewed.artifact_sha256} RUN_DIRECTORY"
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
        self._sync_phase("No new execution authority was granted")

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
