"""Tool-use lifecycle projection for the chat transcript."""

from __future__ import annotations

import uuid
from pathlib import Path

from chemsmart.agent.decision_events import ToolUseEvent
from chemsmart.agent.model_command_parser import parse_model_command
from chemsmart.agent.permissions import PermissionMode
from chemsmart.agent.tui.chat_helpers import (
    _command_details_text,
    _final_command_text,
    _provider_model_label,
    _public_tool_result_payload,
    _tool_use_payload_summary,
    _tool_use_summary_payload,
)
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.tool_meta import format_assumptions_banner
from chemsmart.agent.tui.widgets.cells import (
    AgentMessageCell,
    CommandInterpretationCell,
    DryRunInputCell,
    FinalAnswerCell,
    SynthesisTraceCell,
    ToolCallCell,
)
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.transcript import Transcript


def _tool_lifecycle_note(event: ToolUseEvent) -> str | None:
    if event.status == "approved":
        if event.scope == "session":
            return "approved for the rest of this session"
        if event.reason == "xtb_real_runs_auto":
            return "auto-approved by the xtb_real_runs policy"
        return "approved"
    if event.status == "denied":
        return (
            "Denied by user; model will continue without this action."
            if event.reason == "user_denied"
            else f"Denied: {event.reason or 'policy blocked'}"
        )
    if event.status != "pending":
        return (
            event.reason
            or _tool_use_payload_summary(event.payload)
            or event.status
        )
    return None


class ToolEventMixin:
    """Render tool calls while preserving provider call lifecycle identity."""

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
        footer = self.query_one(FooterWidget)
        # ask_user is a clarification, not a permission-gated action. Render it
        # as a plain question (never a "[unknown]" risky pending tool cell) and
        # leave the composer ready for the user's answer.
        if event.tool == "ask_user":
            if event.status == "pending":
                self._render_ask_user_prompt(event)
            return
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
            note=_tool_lifecycle_note(event),
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
            self._render_dry_run_result(event)
        self._update_footer_for_tool_status(event, footer)

    def _render_dry_run_result(self, event: ToolUseEvent) -> None:
        summary = _tool_use_summary_payload(event.payload)
        if summary is None:
            return
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
        self.query_one(Transcript).add_cell(
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

    def _update_footer_for_tool_status(
        self, event: ToolUseEvent, footer: FooterWidget
    ) -> None:
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
