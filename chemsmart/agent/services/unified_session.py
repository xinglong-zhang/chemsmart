"""Unified tool-loop orchestration behind the public AgentSession facade."""

from __future__ import annotations

import time
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Protocol

from chemsmart.agent.harness.workflow_state import (
    hydrate_workflow_state,
    workflow_state_scope,
)
from chemsmart.agent.loop import ToolLoop, ToolLoopBudgets
from chemsmart.agent.models import Plan, SessionState
from chemsmart.agent.permissions import (
    ApprovalDecision,
    PermissionMode,
    PermissionPolicy,
    RuntimePermissionMode,
)
from chemsmart.agent.prompts.identity import ensure_system_message
from chemsmart.agent.provider_adapter import ToolRequest
from chemsmart.agent.runtime.contracts import RuntimeV2Mode, TaskPhase
from chemsmart.agent.runtime.events import EventKind
from chemsmart.agent.services.plan_support import (
    is_chitchat_request,
    render_plan,
    synthetic_plan_from_tool_requests,
)

from chemsmart.agent.runtime.orchestrator import RuntimeController


class SessionLoopHost(Protocol):
    """Structural contract needed by the unified session runner."""

    state: SessionState | None
    session_dir: Path | None
    handle_store: Any | None
    decision_log: Any | None
    registry: Any
    runtime_v2_mode: RuntimeV2Mode
    _run_start_time: float | None
    _llm_stats: list[dict[str, Any]]
    _last_harness_result: Any | None

    def _start_new_session(self, request: str) -> None: ...
    def _has_pending_ask_user(self) -> bool: ...
    def _resume_pending_ask_user(self, answer: str) -> None: ...
    def _start_new_turn(self, request: str) -> None: ...
    def _save_state(self) -> None: ...
    def _refresh_conversation_history(self) -> None: ...
    def _provider_instance(self) -> Any: ...
    def _tool_defs_for_provider(self, provider_name: str) -> Any: ...
    def _ensure_runtime_controller(self) -> RuntimeController | None: ...
    def _log_loop_mode(self, policy: PermissionPolicy) -> None: ...
    def _tool_loop_system_prompt(
        self,
        policy: PermissionPolicy,
        *,
        request: str,
    ) -> str: ...
    def _write_training_episode(self, **kwargs: Any) -> None: ...
    def _finalize_session(self, **kwargs: Any) -> None: ...
    def _runtime_v2_metadata(self) -> dict[str, Any]: ...


@dataclass(frozen=True)
class LoopSetup:
    provider: Any
    provider_name: str
    policy: PermissionPolicy
    runtime_mode: RuntimePermissionMode | None
    runtime_controller: RuntimeController | None
    runtime_lifecycle: Any | None
    allowed_tools: set[str] | None
    tool_defs: list[dict[str, Any]] | None
    messages: list[dict[str, Any]]
    budgets: ToolLoopBudgets


@dataclass(frozen=True)
class LoopProjection:
    tool_requests: list[Any]
    tool_outcomes: list[Any]
    results: list[Any]
    plan: Plan
    dry_run_results: list[dict[str, Any]]
    runtime_result: dict[str, Any] | None
    intent_request: str
    is_chitchat: bool


class UnifiedSessionRunner:
    """Coordinate one API-provider tool-loop turn for an AgentSession."""

    def __init__(self, session: SessionLoopHost) -> None:
        self.session = session

    def run(
        self,
        request: str,
        *,
        budgets: ToolLoopBudgets | None = None,
        messages: list[dict[str, Any]] | None = None,
        log_raw_provider_turns: bool = False,
        policy: PermissionPolicy | None = None,
        approver: Callable[[ToolRequest], ApprovalDecision] | None = None,
    ) -> dict[str, Any]:
        continuing_ask = self._begin_turn(request)
        setup = self._prepare_setup(
            request=request,
            continuing_ask=continuing_ask,
            budgets=budgets,
            messages=messages,
            log_raw_provider_turns=log_raw_provider_turns,
            policy=policy,
        )
        loop_result = self._run_tool_loop(setup, approver)
        self._complete_runtime(setup.runtime_controller, loop_result)
        projection = self._project_result(
            request=request,
            continuing_ask=continuing_ask,
            loop_result=loop_result,
        )
        self._persist_projection(setup, projection, loop_result)

        ask_user = loop_result.get("ask_user")
        if ask_user:
            self.session._write_training_episode(
                provider_name=setup.provider_name,
                provider=setup.provider,
                loop_result=loop_result,
                paused=True,
            )
            self.session._refresh_conversation_history()
            return self._result_payload(setup, projection, loop_result)

        self._finalize_completed_turn(setup, projection, loop_result)
        return self._result_payload(setup, projection, loop_result)

    def _begin_turn(self, request: str) -> bool:
        session = self.session
        session._run_start_time = time.perf_counter()
        session._llm_stats = []
        session._last_harness_result = None
        if session.state is None or session.session_dir is None:
            session._start_new_session(request)
            continuing_ask = False
        elif session._has_pending_ask_user():
            session._resume_pending_ask_user(request)
            continuing_ask = True
        else:
            session._start_new_turn(request)
            continuing_ask = False

        assert session.state is not None
        assert session.decision_log is not None
        assert session.handle_store is not None
        session._save_state()
        if not continuing_ask:
            session.decision_log.write(
                "request", {"request": request}, rationale=request
            )
        session._refresh_conversation_history()
        return continuing_ask

    def _prepare_setup(
        self,
        *,
        request: str,
        continuing_ask: bool,
        budgets: ToolLoopBudgets | None,
        messages: list[dict[str, Any]] | None,
        log_raw_provider_turns: bool,
        policy: PermissionPolicy | None,
    ) -> LoopSetup:
        session = self.session
        provider = session._provider_instance()
        provider_name = getattr(provider, "name", None) or "openai"
        resolved_policy = policy or PermissionPolicy(
            mode=PermissionMode.DRIVING
        )
        runtime_mode = (
            resolved_policy.mode
            if isinstance(resolved_policy.mode, RuntimePermissionMode)
            else None
        )
        tool_defs = (
            None
            if runtime_mode is not None
            else session._tool_defs_for_provider(provider_name)
        )
        controller, lifecycle, allowed_tools, tool_defs = self._start_runtime(
            request, provider_name, tool_defs
        )
        session._log_loop_mode(resolved_policy)
        prepared_messages = self._prepare_messages(
            request, continuing_ask, messages, resolved_policy
        )
        prepared_budgets = self._prepare_budgets(
            budgets, log_raw_provider_turns
        )
        return LoopSetup(
            provider=provider,
            provider_name=provider_name,
            policy=resolved_policy,
            runtime_mode=runtime_mode,
            runtime_controller=controller,
            runtime_lifecycle=lifecycle,
            allowed_tools=allowed_tools,
            tool_defs=tool_defs,
            messages=prepared_messages,
            budgets=prepared_budgets,
        )

    def _start_runtime(
        self,
        request: str,
        provider_name: str,
        tool_defs: list[dict[str, Any]] | None,
    ) -> tuple[
        RuntimeController | None,
        Any | None,
        set[str] | None,
        list[dict[str, Any]] | None,
    ]:
        session = self.session
        controller = session._ensure_runtime_controller()
        if controller is None:
            return None, None, None, tool_defs
        assert session.state is not None
        with workflow_state_scope(session.state.session_id):
            durable = controller.state
            hydrated = hydrate_workflow_state(
                _durable_workflow_payload(durable, session.state.cwd),
                cwd=session.state.cwd,
                overwrite=_has_durable_workflow_state(durable),
            )
            controller.start_turn(
                request=request,
                turn_index=session.state.turn_index,
                provider_name=provider_name,
                cwd=session.state.cwd,
                workflow_state=hydrated,
            )
        lifecycle = controller.lifecycle()
        if session.runtime_v2_mode is not RuntimeV2Mode.ACTIVE:
            return controller, lifecycle, None, tool_defs
        assert controller.selection is not None
        return (
            controller,
            lifecycle,
            set(controller.selection.direct),
            controller.tool_defs(provider_name),
        )

    def _prepare_messages(
        self,
        request: str,
        continuing_ask: bool,
        messages: list[dict[str, Any]] | None,
        policy: PermissionPolicy,
    ) -> list[dict[str, Any]]:
        session = self.session
        assert session.state is not None
        if messages is None and continuing_ask:
            messages = [
                *deepcopy(session.state.pending_messages or []),
                {"role": "user", "content": request},
            ]
        elif messages is None:
            messages = [{"role": "user", "content": request}]
        return ensure_system_message(
            messages,
            session._tool_loop_system_prompt(policy, request=request),
        )

    @staticmethod
    def _prepare_budgets(
        budgets: ToolLoopBudgets | None,
        log_raw_provider_turns: bool,
    ) -> ToolLoopBudgets:
        if budgets is None:
            return ToolLoopBudgets(
                log_provider_turn_raw=log_raw_provider_turns
            )
        if not log_raw_provider_turns or budgets.log_provider_turn_raw:
            return budgets
        return ToolLoopBudgets(
            max_model_steps_per_turn=budgets.max_model_steps_per_turn,
            max_total_tool_calls_per_turn=budgets.max_total_tool_calls_per_turn,
            max_consecutive_tool_errors=budgets.max_consecutive_tool_errors,
            max_same_signature_retries=budgets.max_same_signature_retries,
            max_provider_errors_per_turn=budgets.max_provider_errors_per_turn,
            log_provider_turn_raw=True,
        )

    def _run_tool_loop(
        self,
        setup: LoopSetup,
        approver: Callable[[ToolRequest], ApprovalDecision] | None,
    ) -> dict[str, Any]:
        session = self.session
        assert session.handle_store is not None
        assert session.decision_log is not None
        assert session.state is not None
        loop = ToolLoop(
            provider=setup.provider,
            registry=session.registry,
            handle_store=session.handle_store,
            decision_log=session.decision_log,
            budgets=setup.budgets,
            policy=setup.policy,
            approver=approver,
            lifecycle=setup.runtime_lifecycle,
        )
        with workflow_state_scope(session.state.session_id):
            return loop.run_turn(
                messages=setup.messages,
                tool_defs=setup.tool_defs,
                mode=setup.runtime_mode,
                allowed_tool_names=setup.allowed_tools,
            )

    @staticmethod
    def _complete_runtime(
        controller: RuntimeController | None,
        loop_result: dict[str, Any],
    ) -> None:
        if controller is None:
            return
        ask_user = loop_result.get("ask_user")
        if ask_user and controller.state.phase is not TaskPhase.WAITING_USER:
            controller.emit(
                EventKind.CLARIFICATION_REQUESTED,
                {
                    "slots": ["clarification"],
                    "question": ask_user.get("question", ""),
                },
            )
        elif loop_result.get("limit_reason"):
            controller.block(reason=str(loop_result["limit_reason"]))
        elif not ask_user:
            controller.complete()
        notice = controller.completion_notice()
        if notice:
            text = str(loop_result.get("assistant_text") or "").rstrip()
            loop_result["assistant_text"] = (
                f"{text}\n\n{notice}" if text else notice
            )

    def _project_result(
        self,
        *,
        request: str,
        continuing_ask: bool,
        loop_result: dict[str, Any],
    ) -> LoopProjection:
        outcomes = list(loop_result["tool_outcomes"])
        requests = list(loop_result["tool_requests"])
        results = [
            outcome.raw_result
            if outcome.raw_result is not None
            else outcome.result
            for outcome in outcomes
            if outcome.status != "ask_user"
        ]
        dry_runs = [
            outcome.raw_result
            for outcome in outcomes
            if outcome.name == "dry_run_input"
            and outcome.status == "ok"
            and isinstance(outcome.raw_result, dict)
        ]
        runtime_result = next(
            (
                outcome.raw_result
                for outcome in reversed(outcomes)
                if outcome.name == "validate_runtime"
                and outcome.status == "ok"
                and isinstance(outcome.raw_result, dict)
            ),
            None,
        )
        assert self.session.state is not None
        intent_request = (
            self.session.state.request
            if continuing_ask and self.session.state.request is not None
            else request
        )
        return LoopProjection(
            tool_requests=requests,
            tool_outcomes=outcomes,
            results=results,
            plan=synthetic_plan_from_tool_requests(requests),
            dry_run_results=dry_runs,
            runtime_result=runtime_result,
            intent_request=intent_request,
            is_chitchat=(is_chitchat_request(intent_request) and not requests),
        )

    def _persist_projection(
        self,
        setup: LoopSetup,
        projection: LoopProjection,
        loop_result: dict[str, Any],
    ) -> None:
        session = self.session
        assert session.state is not None
        session.state.plan = projection.plan
        session.state.request_intent = (
            "chitchat"
            if projection.is_chitchat
            else "workflow"
            if projection.tool_requests
            else "advisory"
        )
        session.state.total_steps_planned = len(projection.tool_requests)
        session.state.current_step_index = len(projection.results)
        session.state.pending_messages = loop_result.get("messages")
        session.state.pending_ask_user = loop_result.get("ask_user")
        session._save_state()
        session._llm_stats = [
            {
                "stage": "tool_loop",
                "attempt": 1,
                "provider_name": setup.provider_name,
                "resolved_model": getattr(setup.provider, "default_model", None),
                "input_tokens": loop_result["total_input_tokens"],
                "output_tokens": loop_result["total_output_tokens"],
                "latency_ms": _elapsed_ms(session._run_start_time),
                "success": True,
            }
        ]

    def _finalize_completed_turn(
        self,
        setup: LoopSetup,
        projection: LoopProjection,
        loop_result: dict[str, Any],
    ) -> None:
        session = self.session
        session._write_training_episode(
            provider_name=setup.provider_name,
            provider=setup.provider,
            loop_result=loop_result,
        )
        assert session.state is not None
        session.state.pending_messages = None
        session.state.pending_ask_user = None
        session._save_state()
        session._finalize_session(
            verdict=None,
            blocked=False,
            block_reason=loop_result["limit_reason"],
            dry_run_results=projection.dry_run_results,
            advisory_only=not projection.tool_requests,
            is_chitchat=projection.is_chitchat,
            rationale=loop_result["assistant_text"] or "",
        )
        session._refresh_conversation_history()

    def _result_payload(
        self,
        setup: LoopSetup,
        projection: LoopProjection,
        loop_result: dict[str, Any],
    ) -> dict[str, Any]:
        session = self.session
        assert session.state is not None
        assert session.session_dir is not None
        ask_user = loop_result.get("ask_user")
        return {
            "session_id": session.state.session_id,
            "session_dir": str(session.session_dir),
            "plan": projection.plan,
            "plan_text": render_plan(projection.plan),
            "critic_verdict": None,
            "completed_steps": session.state.current_step_index,
            "blocked": False,
            "dry_run_result": _primary_dry_run_result(
                projection.dry_run_results
            ),
            "dry_run_results": projection.dry_run_results,
            "runtime_result": projection.runtime_result,
            "preview_submit": None,
            "results": projection.results,
            "assistant_output": loop_result["assistant_text"],
            "tool_requests": projection.tool_requests,
            "tool_outcomes": projection.tool_outcomes,
            "loop_state": _loop_state(loop_result, projection.tool_outcomes),
            "final_message": loop_result["assistant_text"],
            "limit_reason": loop_result["limit_reason"],
            "provider_errors": loop_result["provider_errors"],
            "advisory_only": not projection.tool_requests,
            "is_chitchat": projection.is_chitchat,
            "approval_mode": setup.policy.mode.value,
            "driving_mode": setup.policy.mode == PermissionMode.DRIVING,
            "yolo": setup.policy.yolo,
            "denials_count": loop_result["denials_count"],
            "approvals_count": loop_result["approvals_count"],
            "ask_user_question": ask_user,
            "runtime_v2": session._runtime_v2_metadata(),
        }


def _durable_workflow_payload(durable: Any, cwd: str) -> dict[str, Any]:
    return {
        "cwd": cwd,
        "project": (
            durable.active_project.model_dump(mode="json")
            if durable.active_project is not None
            else None
        ),
        "server": (
            durable.active_server.model_dump(mode="json")
            if durable.active_server is not None
            else None
        ),
        "previous_command": durable.previous_command,
        "unresolved_slots": durable.unresolved_slots,
    }


def _has_durable_workflow_state(durable: Any) -> bool:
    return any(
        (
            durable.active_project is not None,
            durable.active_server is not None,
            bool(durable.previous_command),
            bool(durable.unresolved_slots),
        )
    )


def _loop_state(
    loop_result: dict[str, Any],
    outcomes: list[Any],
) -> dict[str, Any]:
    return {
        "stop_reason": loop_result["stop_reason"],
        "model_steps": loop_result["model_steps"],
        "tool_calls": sum(
            outcome.status not in {"skipped", "ask_user"}
            for outcome in outcomes
        ),
        "limit_reason": loop_result["limit_reason"],
        "provider_errors": loop_result["provider_errors"],
    }


def _primary_dry_run_result(
    dry_run_results: list[dict[str, Any]],
) -> dict[str, Any] | None:
    return dry_run_results[0] if dry_run_results else None


def _elapsed_ms(start_time: float | None) -> int:
    if start_time is None:
        return 0
    return max(0, int(round((time.perf_counter() - start_time) * 1000)))


__all__ = ["UnifiedSessionRunner"]
