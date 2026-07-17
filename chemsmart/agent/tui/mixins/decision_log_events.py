"""Decision-log event projection into transcript and footer state."""

from __future__ import annotations

from chemsmart.agent.decision_events import (
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
)
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.tool_meta import format_assumptions_banner
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
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.transcript import Transcript


class DecisionEventMixin:
    """Apply typed decision events without owning log transport."""

    def _apply_log_entry(self, entry: dict) -> None:
        event = parse_decision_event(
            entry,
            session_dir=(
                self._tailer_path.parent if self._tailer_path else None
            ),
        )
        payload = entry.get("payload") or {}
        handlers = {
            RequestEvent: self._apply_request_event,
            PlanEvent: self._apply_plan_event,
            AssistantTurnEvent: self._apply_assistant_turn_event,
            ToolCallEvent: self._apply_tool_call_event,
            ToolPreviewEvent: self._apply_tool_preview_event,
            ToolUseEvent: self._apply_tool_use_event,
            MethodEvent: self._apply_method_event,
            DryRunInputEvent: self._apply_dry_run_event,
            RuntimeValidationEvent: self._apply_runtime_validation_event,
            SubmissionPreviewEvent: self._apply_submission_preview_event,
            GeometryHandoffEvent: self._apply_geometry_handoff_event,
            CriticVerdictEvent: self._apply_critic_verdict_event,
            ErrorEvent: self._apply_error_event,
            SessionSummaryEvent: lambda summary: (
                self._apply_session_summary_event(
                    summary,
                    request_intent=str(payload.get("request_intent") or ""),
                )
            ),
            IgnoredEvent: self._ignore_decision_event,
        }
        handler = handlers.get(type(event))
        if handler is not None:
            handler(event)

        if payload.get("tool") in {"run_local", "submit_hpc"}:
            self._refresh_job_snapshot()

    def _apply_request_event(self, event: RequestEvent) -> None:
        transcript = self.query_one(Transcript)
        footer = self.query_one(FooterWidget)
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

    def _apply_plan_event(self, event: PlanEvent) -> None:
        transcript = self.query_one(Transcript)
        footer = self.query_one(FooterWidget)
        self._current_plan = event.plan
        self._current_plan_text = event.text
        if event.plan.is_chitchat():
            self._workflow_cell = None
            transcript.add_cell(
                AgentMessageCell(
                    event.plan.rationale or "(no reply)", title="Reply"
                )
            )
            footer.set_phase(Phase.IDLE)
            footer.set_hint("Ready")
            return
        if not event.plan.steps:
            self._workflow_cell = None
            transcript.add_cell(
                AgentMessageCell(
                    event.plan.rationale or "(no rationale)", title="Advisory"
                )
            )
            footer.set_phase(Phase.FINISHED)
            footer.set_hint("Response ready")
            return
        workflow_cell = PlanCell(event.plan)
        self._workflow_cell = workflow_cell
        transcript.add_cell(workflow_cell)
        footer.set_phase(Phase.PLANNING)
        footer.set_hint("Filling in tool steps…")

    def _apply_assistant_turn_event(self, event: AssistantTurnEvent) -> None:
        if event.text.strip():
            self.query_one(Transcript).add_cell(
                AgentMessageCell(event.text.strip(), title="Assistant")
            )

    def _apply_tool_call_event(self, event: ToolCallEvent) -> None:
        workflow_cell = getattr(self, "_workflow_cell", None)
        footer = self.query_one(FooterWidget)
        if event.status == "done":
            self._apply_workspace_project_tool_result(
                event.tool, event.payload
            )
        if workflow_cell is not None:
            if event.status == "running":
                workflow_cell.mark_started(
                    event.step_index, event.tool, event.args
                )
            else:
                workflow_cell.mark_completed(
                    event.step_index, event.tool, event.payload
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

    def _apply_tool_preview_event(self, event: ToolPreviewEvent) -> None:
        workflow_cell = getattr(self, "_workflow_cell", None)
        if workflow_cell is not None and event.status == "running":
            workflow_cell.mark_started(
                event.step_index, event.tool, event.args
            )
        footer = self.query_one(FooterWidget)
        footer.set_phase(Phase.VALIDATING)
        footer.set_hint("Preparing submission preview…")
        footer.set_tool_progress(event.tool, step=event.step_index)

    def _apply_method_event(self, event: MethodEvent) -> None:
        self.query_one(Transcript).add_cell(MethodCell(event.recommendation))

    def _apply_dry_run_event(self, event: DryRunInputEvent) -> None:
        workflow_cell = getattr(self, "_workflow_cell", None)
        if workflow_cell is not None:
            workflow_cell.mark_completed(
                event.step_index,
                "dry_run_input",
                {"inputfile": event.inputfile, "content": event.content},
            )
        self.query_one(Transcript).add_cell(
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
        footer = self.query_one(FooterWidget)
        footer.set_phase(Phase.DRY_RUN_READY)
        footer.set_hint("Dry-run input ready")

    def _apply_runtime_validation_event(
        self, event: RuntimeValidationEvent
    ) -> None:
        workflow_cell = getattr(self, "_workflow_cell", None)
        if workflow_cell is not None:
            workflow_cell.mark_completed(
                event.step_index, "validate_runtime", event.validation
            )
        self.query_one(Transcript).add_cell(
            RuntimeValidationCell(event.validation)
        )
        footer = self.query_one(FooterWidget)
        footer.set_phase(Phase.DRY_RUN_READY)
        footer.set_hint("Runtime check ready")

    def _apply_submission_preview_event(
        self, event: SubmissionPreviewEvent
    ) -> None:
        workflow_cell = getattr(self, "_workflow_cell", None)
        if workflow_cell is not None:
            workflow_cell.mark_completed(
                event.step_index, "submit_hpc", event.preview
            )
        self.query_one(Transcript).add_cell(
            SubmissionPreviewCell(event.preview)
        )
        footer = self.query_one(FooterWidget)
        footer.set_phase(Phase.DRY_RUN_READY)
        footer.set_hint("Submission preview ready")

    def _apply_geometry_handoff_event(
        self, event: GeometryHandoffEvent
    ) -> None:
        self.query_one(Transcript).add_cell(
            GeometryHandoffCell(event.molecule, session_dir=event.session_dir)
        )

    def _apply_critic_verdict_event(self, event: CriticVerdictEvent) -> None:
        footer = self.query_one(FooterWidget)
        if self._current_plan is not None and self._current_plan.is_chitchat():
            footer.set_phase(Phase.IDLE)
            footer.set_hint("Ready")
            return
        self._current_verdict = event.verdict
        self.query_one(Transcript).add_cell(CriticVerdictCell(event.verdict))
        if event.verdict.verdict == "reject":
            footer.set_phase(Phase.ERROR)
            footer.set_hint("Critic blocked execution")
            return
        footer.set_phase(Phase.DRY_RUN_READY)
        footer.set_hint("Critic verdict ready")

    def _apply_error_event(self, event: ErrorEvent) -> None:
        details = event.details if isinstance(event.details, dict) else {}
        workflow_cell = getattr(self, "_workflow_cell", None)
        if workflow_cell is not None and details:
            workflow_cell.mark_failed(
                int(details.get("step_index") or 0),
                str(details.get("tool") or ""),
                event.message,
            )
        self.query_one(Transcript).add_cell(
            ErrorCell(event.title, event.message, event.details)
        )
        footer = self.query_one(FooterWidget)
        footer.set_phase(Phase.ERROR)
        footer.set_hint("Agent reported an error")

    def _apply_session_summary_event(
        self,
        event: SessionSummaryEvent,
        *,
        request_intent: str = "",
    ) -> None:
        self._pending_approval = False
        self._pending_tool_request = None
        footer = self.query_one(FooterWidget)
        if request_intent == "chitchat":
            footer.set_phase(Phase.IDLE)
            footer.set_hint("Ready")
            return
        if not event.blocked and not (
            event.total_steps_planned or event.total_steps_executed
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
        transcript = self.query_one(Transcript)
        transcript.add_cell(summary_cell)
        transcript.collapse_tool_chain(self._active_turn_id)

    @staticmethod
    def _ignore_decision_event(_event: IgnoredEvent) -> None:
        return None

    def _has_user_message(self, request: str) -> bool:
        return request in self._user_requests
