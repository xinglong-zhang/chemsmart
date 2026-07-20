"""Approval state and permission policy for chat requests."""

from __future__ import annotations

from pathlib import Path
from threading import Event

from chemsmart.agent.permissions import (
    ApprovalDecision,
    PermissionMode,
    PermissionPolicy,
)
from chemsmart.agent.provider_adapter import ToolRequest
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.tool_meta import tool_description
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.popups import (
    ApprovalResult,
    PermissionModeResult,
    ProjectWriteOverlay,
    ProjectWriteResult,
    TextPromptOverlay,
    build_approval_overlay,
)
from chemsmart.agent.tui.widgets.transcript import Transcript
from chemsmart.settings.workspace_project import workspace_project_path


class ApprovalFlowMixin:
    """Own approval transitions and reset request-scoped state."""

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
