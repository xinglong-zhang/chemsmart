"""Provider routing and direct tool request flow for the chat screen."""

from __future__ import annotations

import json
from threading import get_ident

from textual import work

from chemsmart.agent.permissions import (
    ApprovalDecision,
    ResolvedDecision,
)
from chemsmart.agent.provider_adapter import ToolRequest
from chemsmart.agent.registry import ToolRegistry
from chemsmart.agent.runtime.calculations import inspect_calculation
from chemsmart.agent.tui.chat_helpers import (
    _calculation_diagnostic_summary,
    _is_calculation_diagnostic_request,
)
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.widgets.cells import UserMessageCell
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.transcript import Transcript


class RequestFlowMixin:
    """Route user requests and execute direct slash-tool calls."""

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
        if not self._resolve_slash_approval(
            request,
            policy,
            tool_name=tool_name,
            description=description,
            normalized_args=normalized_args,
            explicit_approval=explicit_approval,
        ):
            return {"tool": tool_name, "status": "denied"}
        result = registry.call(tool_name, normalized_args)
        return self._publish_slash_tool_result(
            tool_name, description, normalized_args, result
        )

    def _resolve_slash_approval(
        self,
        request: ToolRequest,
        policy,
        *,
        tool_name: str,
        description: str,
        normalized_args: dict[str, object],
        explicit_approval: bool,
    ) -> bool:
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
                return False
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
            return False
        else:
            self.app.call_from_thread(
                self._publish_tool_call_cell,
                tool_name,
                "approved",
                description,
                normalized_args,
                "Auto-approved.",
            )
        return True

    def _publish_slash_tool_result(
        self,
        tool_name: str,
        description: str,
        normalized_args: dict[str, object],
        result: object,
    ) -> dict[str, object]:
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
