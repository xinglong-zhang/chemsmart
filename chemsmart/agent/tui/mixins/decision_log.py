"""Decision-log streaming and tool-result presentation for the chat screen."""

from __future__ import annotations

import uuid
from pathlib import Path
from threading import get_ident

from chemsmart.agent.core import DecisionLog
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
from chemsmart.agent.model_command_parser import parse_model_command
from chemsmart.agent.permissions import PermissionMode
from chemsmart.agent.services.conversation_memory import ConversationMemory
from chemsmart.agent.tui.chat_helpers import (
    _command_details_text,
    _final_command_text,
    _provider_model_label,
    _public_tool_result_payload,
    _tool_use_payload_summary,
    _tool_use_summary_payload,
)
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.services.log_tailer import LogTailer
from chemsmart.agent.tui.tool_meta import format_assumptions_banner
from chemsmart.agent.tui.widgets.cells import (
    AgentMessageCell,
    CommandInterpretationCell,
    CriticVerdictCell,
    DryRunInputCell,
    ErrorCell,
    FinalAnswerCell,
    GeometryHandoffCell,
    MethodCell,
    PlanCell,
    RuntimeValidationCell,
    SubmissionPreviewCell,
    SynthesisTraceCell,
    ToolCallCell,
    UserMessageCell,
)
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.transcript import Transcript


class DecisionLogMixin:
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

