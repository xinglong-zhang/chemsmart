from __future__ import annotations

import json
import os
import re
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
from chemsmart.agent.harness.runner import evaluate_harness
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
from chemsmart.agent.services.unified_session import UnifiedSessionRunner

UTC = timezone.utc

_GAUSSIAN_ROUTE_RE = re.compile(r"^\s*#\s*\S+", re.MULTILINE)
_ORCA_ROUTE_RE = re.compile(r"^\s*!\s*\S+", re.MULTILINE)
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
        """Append this completed turn to the cross-session training store."""

        try:
            from chemsmart.agent.training_log import (
                TrainingEpisodeWriter,
                tool_records_from_outcomes,
            )

            if self._training_writer is None:
                self._training_writer = TrainingEpisodeWriter()
            if not self._training_writer.enabled or self.state is None:
                return
            outcomes = list(loop_result.get("tool_outcomes") or [])
            terminal_state = next(
                (
                    outcome.raw_result.get("terminal_state")
                    for outcome in reversed(outcomes)
                    if isinstance(outcome.raw_result, dict)
                    and isinstance(
                        outcome.raw_result.get("terminal_state"), dict
                    )
                ),
                None,
            )
            self._training_writer.write_episode(
                session_id=self.state.session_id,
                turn=self.state.turn_index,
                provider_name=provider_name,
                model=getattr(provider, "default_model", None),
                messages=_json_safe(loop_result.get("messages") or []),
                tool_records=tool_records_from_outcomes(
                    outcomes,
                    requests=list(loop_result.get("tool_requests") or []),
                ),
                tool_requests=list(loop_result.get("tool_requests") or []),
                approvals_count=loop_result.get("approvals_count") or 0,
                denials_count=loop_result.get("denials_count") or 0,
                cwd=self.state.cwd,
                available_tools=self._registry_tool_names(),
                final_answer=str(loop_result.get("assistant_text") or ""),
                terminal_state=terminal_state,
                paused=paused,
            )
        except Exception:
            # Training capture must never break a live turn.
            pass

    def _registry_tool_names(self) -> list[str]:
        # Custom/test registries only need `call` + tool defs; don't let a
        # missing list_tools() silently disable episode capture.
        list_tools = getattr(self.registry, "list_tools", None)
        if not callable(list_tools):
            return []
        try:
            return [tool.name for tool in list_tools()]
        except Exception:
            return []

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
        issues = _filter_critic_issues(
            plan=plan,
            issues=verdict.issues,
            dry_submit=dry_submit,
        )
        rationale_parts = [verdict.rationale] if verdict.rationale else []
        final_verdict = verdict.verdict

        if runtime_result is not None:
            runtime_ok = runtime_result.get("ok")
            if runtime_ok == "fail":
                final_verdict = "reject"
                issues.extend(runtime_result.get("local_issues", []))
                rationale_parts.append("validate_runtime returned fail")
            elif runtime_ok == "partial" and not dry_submit:
                final_verdict = (
                    "warn" if final_verdict == "ok" else final_verdict
                )
                issues.extend(runtime_result.get("remote_unknown", []))
                rationale_parts.append("validate_runtime returned partial")

        malformed_issues = [
            issue
            for result in dry_run_results
            if (issue := _malformed_input_issue(result)) is not None
        ]
        if malformed_issues:
            if final_verdict != "reject":
                final_verdict = "warn"
            issues.extend(malformed_issues)
            rationale_parts.append(
                "dry-run input route line failed basic validation"
            )

        irc_keyword_issues = _missing_irc_keyword_issues(
            plan=plan,
            dry_run_results=dry_run_results,
        )
        if irc_keyword_issues:
            final_verdict = "reject"
            issues.extend(irc_keyword_issues)
            rationale_parts.append(
                "IRC input route line missing required keyword"
            )

        harness_result = evaluate_harness(
            plan=plan,
            dry_run_results=dry_run_results,
        )
        self._last_harness_result = harness_result
        if self.decision_log is not None:
            self.decision_log.write(
                "harness_result",
                harness_result.to_dict(),
                rationale=f"runtime harness verdict: {harness_result.verdict}",
            )
        if harness_result.verdict == "reject":
            final_verdict = "reject"
            issues.extend(_harness_issue_messages(harness_result))
            rationale_parts.append(
                "software invariant harness rejected generated input"
            )
        elif harness_result.verdict == "warn":
            final_verdict = "warn" if final_verdict == "ok" else final_verdict
            issues.extend(_harness_issue_messages(harness_result))
            rationale_parts.append(
                "software invariant harness warned on generated input"
            )

        if preview_submit is not None:
            duplicate_check = preview_submit.get("duplicate_check", {})
            if duplicate_check.get("duplicate"):
                final_verdict = "reject"
                message = (
                    duplicate_check.get("message")
                    or "duplicate submission detected"
                )
                issues.append(message)
                rationale_parts.append(
                    "submit_hpc duplicate check rejected the plan"
                )

        if final_verdict == "warn" and not issues:
            final_verdict = "ok"

        return CriticVerdict(
            verdict=final_verdict,
            confidence=verdict.confidence,
            issues=_dedupe_strings(issues),
            rationale="; ".join(part for part in rationale_parts if part),
        )

    @staticmethod
    def _block_reason(
        verdict: CriticVerdict,
        dry_submit: bool,
        allow_remote_unknown: bool,
        allow_critic_override: bool,
    ) -> str | None:
        if verdict.verdict == "reject":
            return "critic_reject"
        if verdict.verdict == "ok":
            return None

        has_remote_unknown = any(
            _is_remote_unknown_issue(issue) for issue in verdict.issues
        )
        has_other_warn = any(
            not _is_remote_unknown_issue(issue) for issue in verdict.issues
        )
        if not verdict.issues:
            has_other_warn = True
        if has_remote_unknown and not dry_submit and not allow_remote_unknown:
            return "critic_warn_remote_unknown"
        if has_other_warn and not allow_critic_override:
            return "critic_warn_no_override"
        return None


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


def _filter_critic_issues(
    *,
    plan: Plan,
    issues: list[str],
    dry_submit: bool,
) -> list[str]:
    geometry_handoff_required = _plan_requires_geometry_handoff(plan)
    filtered: list[str] = []
    for issue in issues:
        normalized = issue.lower()
        if dry_submit and _is_remote_unknown_issue(issue):
            continue
        if (
            "geometry handoff missing" in normalized
            and not geometry_handoff_required
        ):
            continue
        filtered.append(issue)
    return filtered


def _plan_requires_geometry_handoff(plan: Plan) -> bool:
    for step_index, step in enumerate(plan.steps):
        if step.tool != "build_job":
            continue
        kind = step.args.get("kind")
        if not (isinstance(kind, str) and kind.endswith(".sp")):
            continue
        source_step = _source_step_for_ref(plan, step.args.get("molecule"))
        if source_step is None or source_step.tool != "build_molecule":
            continue
        if any(
            prior_step.tool == "build_job"
            and isinstance(prior_step.args.get("kind"), str)
            and prior_step.args["kind"].endswith((".opt", ".ts", ".irc"))
            for prior_step in plan.steps[:step_index]
        ):
            return True
    return False


def _source_step_for_ref(plan: Plan, value: Any) -> Step | None:
    if not isinstance(value, str):
        return None
    match = _REFERENCE_RE.match(value)
    if match is None or match.group("path"):
        return None
    ref_index = int(match.group("index")) - 1
    if ref_index < 0 or ref_index >= len(plan.steps):
        return None
    return plan.steps[ref_index]


def _is_remote_unknown_issue(issue: str) -> bool:
    return (
        issue.startswith("server.")
        or issue.endswith("on HPC")
        or issue == "ssh login reachable"
    )


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


def _malformed_input_issue(
    dry_run_result: dict[str, Any] | None,
) -> str | None:
    if not dry_run_result:
        return None
    content = dry_run_result.get("content")
    inputfile = str(dry_run_result.get("inputfile", "")).lower()
    if not isinstance(content, str):
        return None
    if inputfile.endswith(".inp"):
        if _ORCA_ROUTE_RE.search(content) is None:
            return "ORCA route line missing or malformed"
        return None
    if _GAUSSIAN_ROUTE_RE.search(content) is None:
        return "Gaussian route line missing or malformed"
    return None


def _missing_irc_keyword_issues(
    plan: Plan,
    dry_run_results: list[dict[str, Any]],
) -> list[str]:
    issues: list[str] = []
    dry_run_steps = [
        step for step in plan.steps if step.tool == "dry_run_input"
    ]
    for step, result in zip(dry_run_steps, dry_run_results):
        kind = _dry_run_job_kind(plan, step)
        if kind == "gaussian.irc":
            route_text = _route_text(result)
            if (
                route_text is None
                or re.search(r"\birc\s*=", route_text, re.IGNORECASE) is None
            ):
                issues.append("Gaussian IRC input missing irc= keyword")
        elif kind == "orca.irc":
            route_text = _route_text(result)
            if (
                route_text is None
                or re.search(r"\birc\b", route_text, re.IGNORECASE) is None
            ):
                issues.append("ORCA IRC input missing IRC keyword")
    return issues


def _dry_run_job_kind(plan: Plan, dry_run_step: Step) -> str | None:
    job_ref = dry_run_step.args.get("job")
    if not isinstance(job_ref, str):
        return None
    match = _REFERENCE_RE.match(job_ref)
    if match is None:
        return None
    job_index = int(match.group("index")) - 1
    if job_index < 0 or job_index >= len(plan.steps):
        return None
    job_step = plan.steps[job_index]
    if job_step.tool != "build_job":
        return None
    kind = job_step.args.get("kind")
    if not isinstance(kind, str):
        return None
    return kind.strip().lower()


def _route_text(dry_run_result: dict[str, Any] | None) -> str | None:
    if not dry_run_result:
        return None
    content = dry_run_result.get("content")
    inputfile = str(dry_run_result.get("inputfile", "")).lower()
    if not isinstance(content, str):
        return None
    if inputfile.endswith(".inp"):
        matches = re.findall(r"^\s*!\s*(.+)$", content, re.MULTILINE)
    else:
        matches = re.findall(r"^\s*#.*$", content, re.MULTILINE)
    if not matches:
        return None
    return " ".join(match.strip() for match in matches)


def _dedupe_strings(values: list[str]) -> list[str]:
    deduped: list[str] = []
    seen: set[str] = set()
    for value in values:
        if value not in seen:
            seen.add(value)
            deduped.append(value)
    return deduped


def _harness_issue_messages(result: HarnessResult) -> list[str]:
    return [f"{issue.rule_id}: {issue.message}" for issue in result.issues]
