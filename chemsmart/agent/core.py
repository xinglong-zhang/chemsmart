from __future__ import annotations

import json
import os
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

from chemsmart.agent.handles import (
    HandleStore,
    store_result_handle,
)
from chemsmart.agent.harness.models import HarnessResult
from chemsmart.agent.harness.session_gates import (
    apply_deterministic_gates,
    block_reason,
    dedupe_strings,
    dry_run_job_kind,
    filter_critic_issues,
    harness_issue_messages,
    is_remote_unknown_issue,
    malformed_input_issue,
    missing_irc_keyword_issues,
    plan_requires_geometry_handoff,
    route_text,
    source_step_for_ref,
)
from chemsmart.agent.loop import (
    ToolLoopBudgets,
    registry_tool_defs_for_provider,
)
from chemsmart.agent.models import (
    CriticVerdict,
    Plan,
    SessionState,
    Step,
    utc_now_iso,
)
from chemsmart.agent.permissions import (
    ApprovalDecision,
    PermissionPolicy,
)
from chemsmart.agent.prompts import load_prompt
from chemsmart.agent.prompts.identity import (
    build_system_prompt,
)
from chemsmart.agent.provider_adapter import ToolRequest
from chemsmart.agent.providers import get_provider
from chemsmart.agent.registry import ToolRegistry
from chemsmart.agent.runtime.contracts import RuntimeV2Mode
from chemsmart.agent.runtime.orchestrator import RuntimeController
from chemsmart.agent.services.conversation_memory import ConversationMemory
from chemsmart.agent.services.plan_support import (
    classify_intent,
    is_chitchat_request,
    is_project_yaml_workflow,
    render_plan as render_plan,
    resolve_plan_intent,
    synthetic_plan_from_tool_requests,
)
from chemsmart.agent.services.planned_session import PlannedSessionRunner
from chemsmart.agent.services.provider_calls import (
    BACKOFF_SECONDS,
    MAX_ATTEMPTS,
    ProviderCallService,
    estimate_tokens,
    parse_json_response,
    resolved_model,
    stringify_response,
)
from chemsmart.agent.services.result_codec import (
    REFERENCE_RE,
    json_safe,
    preview_value,
    resolve_refs,
    restore_json_result,
)
from chemsmart.agent.services.runtime_metrics import (
    elapsed_ms,
    git_sha,
    schema_hash,
)
from chemsmart.agent.services.session_finalizer import SessionFinalizer
from chemsmart.agent.services.session_store import load_current_session_state
from chemsmart.agent.services.step_executor import (
    StepExecutor,
    is_submit_server_error,
    is_tool_error,
    run_local_failed,
    run_local_failure_message,
)
from chemsmart.agent.services.training_capture import (
    TrainingCapture,
    registry_tool_names,
)
from chemsmart.agent.services.unified_session import UnifiedSessionRunner

UTC = timezone.utc

_REFERENCE_RE = REFERENCE_RE
_LLM_MAX_ATTEMPTS = MAX_ATTEMPTS
_LLM_BACKOFF_SECONDS = BACKOFF_SECONDS


_utc_now_iso = utc_now_iso
_synthetic_plan_from_tool_requests = synthetic_plan_from_tool_requests
_is_project_yaml_workflow = is_project_yaml_workflow
_classify_intent = classify_intent
_resolve_plan_intent = resolve_plan_intent
_is_chitchat_request = is_chitchat_request


class DecisionLog:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def write(
        self,
        kind: str,
        payload: dict[str, Any],
        rationale: str = "",
    ) -> None:
        entry = {
            "ts": datetime.now(UTC).isoformat(),
            "kind": kind,
            "payload": _json_safe(payload),
            "rationale": rationale,
        }
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(entry, sort_keys=True) + "\n")

    def read_all(self) -> list[dict[str, Any]]:
        if not self.path.exists():
            return []
        with self.path.open(encoding="utf-8") as handle:
            return [json.loads(line) for line in handle if line.strip()]


class AgentSession:
    def __init__(
        self,
        provider: Any | None = None,
        registry: ToolRegistry | None = None,
        session_root: str | os.PathLike[str] | None = None,
        transport: Any | None = None,
        stage_prompt: str = "tool_loop.md",
        runtime_v2: str | bool | None = None,
    ) -> None:
        self._provider = provider
        self.registry = registry or ToolRegistry.default()
        self.transport = transport
        self._stage_prompt = stage_prompt
        self.session_root = Path(session_root or _default_session_root())
        self.session_root.mkdir(parents=True, exist_ok=True)
        self.state: SessionState | None = None
        self.session_dir: Path | None = None
        self.decision_log: DecisionLog | None = None
        self.handle_store: HandleStore | None = None
        self.conversation_history = ConversationMemory()
        self._run_start_time: float | None = None
        self._llm_stats: list[dict[str, Any]] = []
        self._loop_mode_state: tuple[str, bool] | None = None
        self._last_harness_result: HarnessResult | None = None
        self._training_writer: Any | None = None
        runtime_setting = (
            runtime_v2
            if runtime_v2 is not None
            else os.environ.get("CHEMSMART_AGENT_RUNTIME_V2")
        )
        self.runtime_v2_mode = RuntimeV2Mode.parse(runtime_setting)
        self._runtime_controller: RuntimeController | None = None

    @property
    def stage_prompt(self) -> str:
        return self._stage_prompt

    @classmethod
    def load(
        cls,
        session_id: str,
        **kwargs: Any,
    ) -> "AgentSession":
        session = cls(**_session_kwargs(kwargs))
        session._load_existing_session(session_id)
        if kwargs.get("cwd_override"):
            assert session.state is not None
            session.state.cwd = os.path.abspath(kwargs["cwd_override"])
            session._save_state()
        return session

    @classmethod
    def resume(
        cls,
        session_id: str,
        **kwargs: Any,
    ) -> dict[str, Any]:
        session = cls.load(session_id, **kwargs)
        return session.continue_loaded_session(
            dry_submit=kwargs.get("dry_submit", True),
            pause_before_risky=kwargs.get("pause_before_risky", False),
            allow_remote_unknown=kwargs.get("allow_remote_unknown", False),
            allow_critic_override=kwargs.get("allow_critic_override", False),
            rerender_plan=False,
        )

    def continue_loaded_session(
        self,
        *,
        dry_submit: bool = True,
        pause_before_risky: bool = False,
        allow_remote_unknown: bool = False,
        allow_critic_override: bool = False,
        rerender_plan: bool = False,
    ) -> dict[str, Any]:
        assert self.state is not None
        assert self.decision_log is not None
        self._run_start_time = time.perf_counter()
        self._llm_stats = []
        original_cwd = os.getcwd()
        resume_cwd = os.path.abspath(self.state.cwd)
        os.chdir(resume_cwd)
        try:
            completed_results = self._load_completed_results()
            return self._continue_run(
                completed_results=completed_results,
                dry_submit=dry_submit,
                pause_before_risky=pause_before_risky,
                allow_remote_unknown=allow_remote_unknown,
                allow_critic_override=allow_critic_override,
                rerender_plan=rerender_plan,
            )
        finally:
            os.chdir(original_cwd)

    def run(
        self,
        request: str,
        dry_submit: bool = True,
        pause_before_risky: bool = False,
        allow_remote_unknown: bool = False,
        allow_critic_override: bool = False,
    ) -> dict[str, Any]:
        self._run_start_time = time.perf_counter()
        self._llm_stats = []
        self._last_harness_result = None
        if self.state is None or self.session_dir is None:
            self._start_new_session(request)
        else:
            self._start_new_turn(request)
        assert self.state is not None
        assert self.decision_log is not None
        self._save_state()
        self.decision_log.write(
            "request", {"request": request}, rationale=request
        )
        self._refresh_conversation_history()

        plan = self._planner_call(request)
        self.state.plan = plan
        self.state.request_intent = (
            "chitchat" if plan.is_chitchat() else _classify_intent(request)
        )
        self.state.total_steps_planned = len(plan.steps)
        self._save_state()
        self.decision_log.write(
            "plan", plan.model_dump(), rationale=plan.rationale
        )
        self._refresh_conversation_history()

        return self._continue_run(
            completed_results=[],
            dry_submit=dry_submit,
            pause_before_risky=pause_before_risky,
            allow_remote_unknown=allow_remote_unknown,
            allow_critic_override=allow_critic_override,
            rerender_plan=True,
        )

    def run_loop(
        self,
        request: str,
        *,
        budgets: ToolLoopBudgets | None = None,
        messages: list[dict[str, Any]] | None = None,
        log_raw_provider_turns: bool = False,
        policy: PermissionPolicy | None = None,
        approver: Callable[[ToolRequest], ApprovalDecision] | None = None,
    ) -> dict[str, Any]:
        return UnifiedSessionRunner(self).run(
            request,
            budgets=budgets,
            messages=messages,
            log_raw_provider_turns=log_raw_provider_turns,
            policy=policy,
            approver=approver,
        )

    def _continue_run(
        self,
        completed_results: list[Any],
        dry_submit: bool,
        pause_before_risky: bool,
        allow_remote_unknown: bool,
        allow_critic_override: bool,
        rerender_plan: bool,
    ) -> dict[str, Any]:
        return PlannedSessionRunner(self).run(
            completed_results,
            dry_submit=dry_submit,
            pause_before_risky=pause_before_risky,
            allow_remote_unknown=allow_remote_unknown,
            allow_critic_override=allow_critic_override,
            rerender_plan=rerender_plan,
        )

    def _provider_instance(self) -> Any:
        if self._provider is None:
            self._provider = get_provider()
        return self._provider

    def _planner_call(self, request: str) -> Plan:
        return ProviderCallService(self).planner(request)

    def _critic_call(
        self,
        plan: Plan,
        dry_run_results: list[dict[str, Any]],
        dry_submit: bool,
    ) -> CriticVerdict:
        return ProviderCallService(self).critic(
            plan, dry_run_results, dry_submit
        )

    def _llm_json_call(
        self,
        *,
        stage: str,
        messages: list[dict[str, Any]],
        model_cls: type[Any],
    ) -> Any:
        return ProviderCallService(self).json_call(
            stage=stage, messages=messages, model_cls=model_cls
        )

    def _record_llm_stats(
        self,
        *,
        stage: str,
        attempt: int,
        provider: Any,
        messages: list[dict[str, Any]],
        response: Any,
        latency_ms: int,
    ) -> None:
        ProviderCallService(self).record_stats(
            stage=stage,
            attempt=attempt,
            provider=provider,
            messages=messages,
            response=response,
            latency_ms=latency_ms,
        )

    def _log_llm_failure(
        self,
        *,
        stage: str,
        attempt: int,
        provider: Any,
        messages: list[dict[str, Any]],
        response: Any,
        error: Exception,
        latency_ms: int,
    ) -> None:
        ProviderCallService(self).log_failure(
            stage=stage,
            attempt=attempt,
            provider=provider,
            messages=messages,
            response=response,
            error=error,
            latency_ms=latency_ms,
        )

    def _execute_step(
        self,
        step_index: int,
        step: Step,
        prior_results: list[Any],
        extra_kwargs: dict[str, Any] | None = None,
    ) -> Any:
        return StepExecutor(self).execute(
            step_index, step, prior_results, extra_kwargs
        )

    def _preview_submit_step(
        self,
        step_index: int,
        step: Step,
        prior_results: list[Any],
        dry_submit: bool,
    ) -> Any:
        return StepExecutor(self).preview_submit(
            step_index, step, prior_results, dry_submit
        )

    def _record_skipped_step(
        self,
        step_index: int,
        step: Step,
        reason: str,
    ) -> None:
        StepExecutor(self).record_skipped(step_index, step, reason)

    def _finalize_session(
        self,
        verdict: CriticVerdict | None,
        blocked: bool,
        block_reason: str | None,
        dry_run_results: list[dict[str, Any]],
        run_error: Exception | None = None,
        advisory_only: bool = False,
        is_chitchat: bool = False,
        rationale: str = "",
    ) -> None:
        SessionFinalizer(self, module_path=Path(__file__)).finalize(
            verdict=verdict,
            blocked=blocked,
            block_reason=block_reason,
            dry_run_results=dry_run_results,
            run_error=run_error,
            advisory_only=advisory_only,
            is_chitchat=is_chitchat,
            rationale=rationale,
        )

    def _tools_called(self) -> list[str]:
        tools: list[str] = []
        for entry in self._current_turn_entries():
            if entry.get("kind") not in {
                "tool_call",
                "tool_preview",
                "tool_use_request",
            }:
                continue
            tool = entry.get("payload", {}).get("tool")
            if isinstance(tool, str) and tool not in tools:
                tools.append(tool)
        return tools

    def _tool_defs_for_provider(
        self,
        provider_name: str,
    ) -> list[dict[str, Any]]:
        return registry_tool_defs_for_provider(self.registry, provider_name)

    def _ensure_runtime_controller(
        self,
    ) -> RuntimeController | None:
        if self.runtime_v2_mode is RuntimeV2Mode.OFF:
            return None
        assert self.state is not None
        assert self.session_dir is not None
        if self._runtime_controller is None:
            self._runtime_controller = RuntimeController(
                session_dir=self.session_dir,
                session_id=self.state.session_id,
                registry=self.registry,
                mode=self.runtime_v2_mode,
            )
        return self._runtime_controller

    def _runtime_v2_metadata(self) -> dict[str, Any]:
        controller = self._runtime_controller
        if controller is None:
            return {"mode": self.runtime_v2_mode.value}
        selection = controller.selection
        return {
            "mode": self.runtime_v2_mode.value,
            "phase": controller.state.phase.value,
            "exposed_tools": (
                list(selection.direct) if selection is not None else []
            ),
            "shadow_violations": list(controller.state.shadow_violations),
            "event_log": str(controller.store.path),
            "state_snapshot": str(
                controller.session_dir / "runtime_state.json"
            ),
        }

    def _has_pending_ask_user(self) -> bool:
        if self.state is None:
            return False
        if self._get_logged_summary() is not None:
            return False
        return bool(
            self.state.pending_ask_user and self.state.pending_messages
        )

    def _resume_pending_ask_user(self, answer: str) -> None:
        assert self.state is not None
        assert self.decision_log is not None
        self.state.cwd = os.path.abspath(os.getcwd())
        self.state.env_snapshot = _env_snapshot()
        pending = self.state.pending_ask_user or {}
        self.decision_log.write(
            "ask_user_answer",
            {
                "question": pending.get("question"),
                "options": pending.get("options") or [],
                "answer": answer,
            },
            rationale=answer,
        )

    def _log_loop_mode(self, policy: PermissionPolicy) -> None:
        assert self.decision_log is not None

        from_mode = None
        if self._loop_mode_state is not None:
            from_mode = self._loop_mode_state[0]

        self.decision_log.write(
            "mode_change",
            {
                "from_mode": from_mode,
                "to_mode": policy.mode.value,
                "yolo": policy.yolo,
            },
        )
        self._loop_mode_state = (policy.mode.value, policy.yolo)

    def _metadata_provider_name(self) -> str:
        if self._llm_stats:
            return self._llm_stats[-1].get("provider_name") or "unknown"
        provider = self._provider
        return getattr(provider, "name", None) or "unknown"

    def _metadata_resolved_model(self) -> str:
        if self._llm_stats:
            return self._llm_stats[-1].get("resolved_model") or "unknown"
        provider = self._provider
        return getattr(provider, "default_model", None) or "unknown"

    def _total_llm_tokens(self, field: str) -> int:
        return sum(int(stat.get(field) or 0) for stat in self._llm_stats)

    def _start_new_session(self, request: str) -> None:
        session_id = _new_session_id()
        self.session_dir = self.session_root / session_id
        self.session_dir.mkdir(parents=True, exist_ok=True)
        self.handle_store = HandleStore(self.session_dir)
        self.decision_log = DecisionLog(
            self.session_dir / "decision_log.jsonl"
        )
        self.state = SessionState(
            session_id=session_id,
            cwd=os.path.abspath(os.getcwd()),
            request_started_at=_utc_now_iso(),
            turn_index=1,
            current_step_index=0,
            request=request,
            env_snapshot=_env_snapshot(),
        )
        self.conversation_history = ConversationMemory()

    def _start_new_turn(self, request: str) -> None:
        assert self.state is not None
        assert self.decision_log is not None
        if self._get_logged_summary() is None and self.state.plan is not None:
            raise RuntimeError(
                "Cannot start a new request while the current turn is still "
                "open. Resume or finish the existing turn first."
            )
        self.state.turn_index += 1
        self.state.cwd = os.path.abspath(os.getcwd())
        self.state.request_started_at = _utc_now_iso()
        self.state.current_step_index = 0
        self.state.total_steps_planned = 0
        self.state.plan = None
        self.state.request = request
        self.state.request_intent = "unknown"
        self.state.env_snapshot = _env_snapshot()

    def _refresh_conversation_history(self) -> None:
        if self.decision_log is None:
            self.conversation_history = ConversationMemory()
            return
        self.conversation_history = ConversationMemory.from_entries(
            self.decision_log.read_all()
        )

    def _prompt_session_meta(self, **extra: Any) -> dict[str, Any]:
        meta: dict[str, Any] = dict(extra)
        if self.state is None:
            return meta

        meta.update(
            {
                "session_id": self.state.session_id,
                "turn_index": self.state.turn_index,
                "request_intent": self.state.request_intent,
            }
        )
        if self._loop_mode_state is not None:
            meta["approval_mode"] = self._loop_mode_state[0]
            meta["yolo"] = self._loop_mode_state[1]
        if self._runtime_controller is not None:
            runtime_state = self._runtime_controller.state
            if runtime_state.active_project is not None:
                meta["active_project"] = runtime_state.active_project.name
            if runtime_state.previous_command:
                meta["previous_command"] = runtime_state.previous_command
        return meta

    def _write_training_episode(
        self,
        *,
        provider_name: str,
        provider: Any,
        loop_result: dict[str, Any],
        paused: bool = False,
    ) -> None:
        TrainingCapture(self).write_episode(
            provider_name=provider_name,
            provider=provider,
            loop_result=loop_result,
            paused=paused,
        )

    def _registry_tool_names(self) -> list[str]:
        return registry_tool_names(self.registry)

    def _tool_loop_system_prompt(
        self,
        policy: PermissionPolicy,
        *,
        request: str = "",
    ) -> str:
        current_turn_index = self.state.turn_index if self.state else None
        return build_system_prompt(
            registry=self.registry,
            stage_instructions=load_prompt(self._stage_prompt),
            session_meta=self._prompt_session_meta(
                stage="tool_loop",
                approval_mode=policy.mode.value,
                yolo=policy.yolo,
            ),
            conversation_context=self.conversation_history.prompt_context(
                current_turn_index=current_turn_index
            ),
            request=request,
            max_chars=4096,
        )

    def _current_turn_entries(self) -> list[dict[str, Any]]:
        if self.decision_log is None or self.state is None:
            return []
        return self.conversation_history.entries_for_turn(
            self.decision_log.read_all(),
            self.state.turn_index,
        )

    def _artifact_paths_for_current_turn(self) -> list[Path] | None:
        assert self.state is not None
        assert self.session_dir is not None
        tool_results = [
            entry
            for entry in self._current_turn_entries()
            if entry.get("kind") == "tool_result"
        ]
        if not tool_results:
            return None
        artifact_paths: list[Path] = []
        for entry in tool_results:
            payload = entry.get("payload") or {}
            artifact_name = payload.get("artifact")
            if isinstance(artifact_name, str) and artifact_name:
                artifact_paths.append(self.session_dir / artifact_name)
        return artifact_paths or None

    def _write_result_artifact(self, step_index: int, result: Any) -> Path:
        assert self.session_dir is not None
        assert self.state is not None
        artifact_path = self.session_dir / (
            f"turn_{self.state.turn_index:02d}_step_{step_index + 1:02d}.json"
        )
        artifact_path.write_text(
            json.dumps(_json_safe(result), indent=2, sort_keys=True),
            encoding="utf-8",
        )
        return artifact_path

    def _load_existing_session(self, session_id: str) -> None:
        self.session_dir = self.session_root / session_id
        self.state = load_current_session_state(
            self.session_dir,
            required=True,
        )
        self.handle_store = HandleStore(self.session_dir)
        self.decision_log = DecisionLog(
            self.session_dir / "decision_log.jsonl"
        )
        self._refresh_conversation_history()
        expected_turn_index = max(1, len(self.conversation_history.turns))
        if self.state.turn_index != expected_turn_index:
            self.state.turn_index = expected_turn_index
            self._save_state()

    def _load_completed_results(self) -> list[Any]:
        assert self.session_dir is not None
        assert self.state is not None
        results: list[Any] = []
        artifact_paths = self._artifact_paths_for_current_turn()
        if artifact_paths is None:
            artifact_paths = [
                self.session_dir / f"step_{step_index + 1:02d}.json"
                for step_index in range(self.state.current_step_index)
            ]
        for artifact_path in artifact_paths:
            with artifact_path.open(encoding="utf-8") as handle:
                results.append(json.load(handle))
        return results

    def _store_result_handle(
        self,
        tool_name: str,
        result: Any,
    ) -> str | None:
        return store_result_handle(
            self.handle_store,
            tool_name,
            result,
            summary=_preview_value(result),
        )

    def _save_state(self) -> None:
        assert self.session_dir is not None
        assert self.state is not None
        self.state.save(self.session_dir / "session.json")
        self.state.save(self.session_dir / "state.json")

    def _collect_prior_results(
        self,
        results: list[Any],
        tool_name: str,
    ) -> list[Any]:
        assert self.state is not None
        assert self.state.plan is not None
        return [
            result
            for step, result in zip(self.state.plan.steps, results)
            if step.tool == tool_name
        ]

    def _find_prior_result(self, results: list[Any], tool_name: str) -> Any:
        assert self.state is not None
        assert self.state.plan is not None
        return next(
            (
                result
                for step, result in zip(self.state.plan.steps, results)
                if step.tool == tool_name
            ),
            None,
        )

    def _get_logged_verdict(self) -> CriticVerdict | None:
        verdict_entries = [
            entry
            for entry in self._current_turn_entries()
            if entry.get("kind") == "critic_verdict"
        ]
        if not verdict_entries:
            return None
        return CriticVerdict.model_validate(verdict_entries[-1]["payload"])

    def _get_logged_summary(self) -> dict[str, Any] | None:
        summary_entries = [
            entry
            for entry in self._current_turn_entries()
            if entry.get("kind") == "session_summary"
        ]
        if not summary_entries:
            return None
        payload = summary_entries[-1].get("payload")
        return payload if isinstance(payload, dict) else None

    def _apply_deterministic_gates(
        self,
        plan: Plan,
        verdict: CriticVerdict,
        runtime_result: dict[str, Any] | None,
        dry_run_results: list[dict[str, Any]],
        preview_submit: dict[str, Any] | None,
        dry_submit: bool,
    ) -> CriticVerdict:
        return apply_deterministic_gates(
            self,
            plan=plan,
            verdict=verdict,
            runtime_result=runtime_result,
            dry_run_results=dry_run_results,
            preview_submit=preview_submit,
            dry_submit=dry_submit,
        )

    @staticmethod
    def _block_reason(
        verdict: CriticVerdict,
        dry_submit: bool,
        allow_remote_unknown: bool,
        allow_critic_override: bool,
    ) -> str | None:
        return block_reason(
            verdict,
            dry_submit,
            allow_remote_unknown,
            allow_critic_override,
        )


def run_agent(request: str, **kwargs: Any) -> dict[str, Any]:
    return AgentSession(**_session_kwargs(kwargs)).run(
        request,
        dry_submit=kwargs.get("dry_submit", True),
        allow_remote_unknown=kwargs.get("allow_remote_unknown", False),
        allow_critic_override=kwargs.get("allow_critic_override", False),
    )


def _session_kwargs(kwargs: dict[str, Any]) -> dict[str, Any]:
    return {
        key: kwargs[key]
        for key in (
            "provider",
            "registry",
            "session_root",
            "transport",
            "runtime_v2",
        )
        if key in kwargs
    }


def _default_session_root() -> str:
    return str(Path.home() / ".chemsmart" / "agent" / "sessions")


def _new_session_id() -> str:
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    return f"{stamp}-{uuid.uuid4().hex[:8]}"


def _env_snapshot() -> dict[str, str | None]:
    return {
        "AI_PROVIDER": os.environ.get("AI_PROVIDER"),
        "PWD": os.path.abspath(os.getcwd()),
    }


_elapsed_ms = elapsed_ms


_parse_json_response = parse_json_response


_filter_critic_issues = filter_critic_issues
_plan_requires_geometry_handoff = plan_requires_geometry_handoff
_source_step_for_ref = source_step_for_ref
_is_remote_unknown_issue = is_remote_unknown_issue


_is_submit_server_error = is_submit_server_error


_stringify_response = stringify_response
_estimate_tokens = estimate_tokens
_resolved_model_for_response = resolved_model


_schema_hash = schema_hash


def _git_sha() -> str | None:
    return git_sha(Path(__file__))


_resolve_refs = resolve_refs
_json_safe = json_safe
_preview_value = preview_value
_restore_json_result = restore_json_result


_is_tool_error = is_tool_error


def _primary_dry_run_result(
    dry_run_results: list[dict[str, Any]],
) -> dict[str, Any] | None:
    return dry_run_results[0] if dry_run_results else None


_run_local_failed = run_local_failed
_run_local_failure_message = run_local_failure_message


_malformed_input_issue = malformed_input_issue
_missing_irc_keyword_issues = missing_irc_keyword_issues
_dry_run_job_kind = dry_run_job_kind
_route_text = route_text
_dedupe_strings = dedupe_strings
_harness_issue_messages = harness_issue_messages
