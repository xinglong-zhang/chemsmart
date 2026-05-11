from __future__ import annotations

import hashlib
import importlib
import json
import os
import re
import subprocess
import time
import uuid
from datetime import datetime, timezone

UTC = timezone.utc
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable, Literal

from pydantic import BaseModel, ConfigDict, Field

from chemsmart.agent.handles import HandleStore, is_handle_id
from chemsmart.agent.loop import ToolLoop, ToolLoopBudgets
from chemsmart.agent.permissions import (
    ApprovalDecision,
    PermissionMode,
    PermissionPolicy,
)
from chemsmart.agent.provider_adapter import ToolRequest
from chemsmart.agent.providers import (
    DEFAULT_TIMEOUT_S,
    extract_response_usage,
    get_provider,
)
from chemsmart.agent.registry import ToolRegistry
from chemsmart.agent.services.conversation_memory import ConversationMemory

_RISKY_TOOLS = {"run_local", "submit_hpc"}
_GAUSSIAN_ROUTE_RE = re.compile(r"^\s*#\s*\S+", re.MULTILINE)
_ORCA_ROUTE_RE = re.compile(r"^\s*!\s*\S+", re.MULTILINE)
_REFERENCE_RE = re.compile(
    r"^\$step(?P<index>\d+)(?P<path>(?:\.[A-Za-z_][A-Za-z0-9_]*)*)$"
)
_INTENT_PATTERNS = {
    "opt": (
        r"\bopt(?:imize|imization|imisation)?\b",
        r"\bgeometry optimi[sz]ation\b",
    ),
    "ts": (r"\btransition state\b", r"\bts\b"),
    "irc": (r"\birc\b", r"\breaction path\b"),
    "sp": (
        r"\bsingle[ -]?point\b",
        r"\bsp\b",
        r"\bsingle point energy\b",
    ),
    "freq": (
        r"\bfrequenc(?:y|ies)\b",
        r"\bfreq\b",
        r"\bvibrational\b",
    ),
    "scan": (r"\bscan\b", r"\bpes\b"),
}
_CHITCHAT_EXACT_PATTERNS = (
    r"(?:hi|hello|hey)(?: there)?[!. ]*",
    r"(?:thanks|thank you|thx)[!. ]*",
    r"good (?:morning|afternoon|evening)[!. ]*",
    r"what can you do(?: for me)?\??",
    r"what do you do\??",
    r"who are you\??",
    r"help\??",
)
_CHITCHAT_IDENTITY_PATTERNS = (
    r"(?:hi|hello|hey)[!.,? ]*(?:what(?:'s| is) your name|who are you|what are you|introduce yourself)[!.,? ]*",
    r"what(?:'s| is) your name[!.,? ]*",
    r"who are you[!.,? ]*",
    r"what are you[!.,? ]*",
    r"introduce yourself[!.,? ]*",
    r"tell me about yourself[!.,? ]*",
    r"which (?:model|llm|ai) (?:are you|do you use)[!.,? ]*",
)
_CHITCHAT_TOKENS = {
    "hello",
    "hey",
    "hi",
    "thanks",
    "thank",
    "you",
    "thx",
}
_LLM_MAX_ATTEMPTS = 3
_LLM_BACKOFF_SECONDS = (1, 2, 4)


def _utc_now_iso() -> str:
    return datetime.now(UTC).isoformat()


class Step(BaseModel):
    tool: str
    args: dict[str, Any] = Field(default_factory=dict)
    rationale: str = ""


class Plan(BaseModel):
    steps: list[Step]
    rationale: str = ""
    estimated_cost: str | None = None
    intent: Literal["workflow", "advisory", "chitchat"] | None = None

    def resolved_intent(self) -> Literal["workflow", "advisory", "chitchat"]:
        if self.intent in {"workflow", "advisory", "chitchat"}:
            return self.intent
        if not self.steps:
            return "advisory"
        return "workflow"

    def is_chitchat(self) -> bool:
        return self.resolved_intent() == "chitchat"


class CriticVerdict(BaseModel):
    verdict: Literal["ok", "warn", "reject"]
    confidence: float = Field(default=0.5, ge=0.0, le=1.0)
    issues: list[str] = Field(default_factory=list)
    rationale: str = ""


class SessionState(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    session_id: str
    cwd: str
    started_at: str = Field(default_factory=_utc_now_iso)
    request_started_at: str = Field(default_factory=_utc_now_iso)
    turn_index: int = 1
    request_intent: str = "unknown"
    total_steps_planned: int = 0
    current_step_index: int = 0
    plan: Plan | None = None
    request: str | None = None
    env_snapshot: dict[str, str | None] = Field(default_factory=dict)

    def save(self, path: Path) -> None:
        path.write_text(self.model_dump_json(indent=2))

    @classmethod
    def load(cls, path: Path) -> "SessionState":
        return cls.model_validate_json(path.read_text())


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
    ) -> None:
        self._provider = provider
        self.registry = registry or ToolRegistry.default()
        self.transport = transport
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
        self._run_start_time = time.perf_counter()
        self._llm_stats = []
        if self.state is None or self.session_dir is None:
            self._start_new_session(request)
        else:
            self._start_new_turn(request)
        assert self.state is not None
        assert self.decision_log is not None
        assert self.handle_store is not None

        self._save_state()
        self.decision_log.write(
            "request", {"request": request}, rationale=request
        )
        self._refresh_conversation_history()

        provider = self._provider_instance()
        provider_name = getattr(provider, "name", None) or "openai"
        tool_defs = self._tool_defs_for_provider(provider_name)
        policy = policy or PermissionPolicy(mode=PermissionMode.DRIVING)
        self._log_loop_mode(policy)
        if messages is None:
            messages = [{"role": "user", "content": request}]
        if budgets is None:
            budgets = ToolLoopBudgets(
                log_provider_turn_raw=log_raw_provider_turns
            )
        elif log_raw_provider_turns and not budgets.log_provider_turn_raw:
            budgets = ToolLoopBudgets(
                max_model_steps_per_turn=budgets.max_model_steps_per_turn,
                max_total_tool_calls_per_turn=(
                    budgets.max_total_tool_calls_per_turn
                ),
                max_consecutive_tool_errors=(
                    budgets.max_consecutive_tool_errors
                ),
                max_same_signature_retries=(
                    budgets.max_same_signature_retries
                ),
                log_provider_turn_raw=True,
            )

        loop = ToolLoop(
            provider=provider,
            registry=self.registry,
            handle_store=self.handle_store,
            decision_log=self.decision_log,
            budgets=budgets,
            policy=policy,
            approver=approver,
        )
        loop_result = loop.run_turn(messages=messages, tool_defs=tool_defs)

        tool_requests = loop_result["tool_requests"]
        tool_outcomes = loop_result["tool_outcomes"]
        results = [
            (
                outcome.raw_result
                if outcome.raw_result is not None
                else outcome.result
            )
            for outcome in tool_outcomes
        ]
        synthetic_plan = _synthetic_plan_from_tool_requests(tool_requests)
        dry_run_results = [
            outcome.raw_result
            for outcome in tool_outcomes
            if outcome.name == "dry_run_input"
            and outcome.status == "ok"
            and isinstance(outcome.raw_result, dict)
        ]
        runtime_result = next(
            (
                outcome.raw_result
                for outcome in reversed(tool_outcomes)
                if outcome.name == "validate_runtime"
                and outcome.status == "ok"
                and isinstance(outcome.raw_result, dict)
            ),
            None,
        )

        self.state.plan = synthetic_plan
        self.state.request_intent = (
            "chitchat"
            if _is_chitchat_request(request) and not tool_requests
            else "workflow" if tool_requests else "advisory"
        )
        self.state.total_steps_planned = len(tool_requests)
        self.state.current_step_index = len(results)
        self._save_state()

        self._llm_stats = [
            {
                "stage": "tool_loop",
                "attempt": 1,
                "provider_name": provider_name,
                "resolved_model": getattr(provider, "default_model", None),
                "input_tokens": loop_result["total_input_tokens"],
                "output_tokens": loop_result["total_output_tokens"],
                "latency_ms": _elapsed_ms(self._run_start_time),
                "success": True,
            }
        ]
        self._finalize_session(
            verdict=None,
            blocked=False,
            block_reason=loop_result["limit_reason"],
            dry_run_results=dry_run_results,
            advisory_only=not tool_requests,
            is_chitchat=(_is_chitchat_request(request) and not tool_requests),
            rationale=loop_result["assistant_text"] or "",
        )
        self._refresh_conversation_history()

        return {
            "session_id": self.state.session_id,
            "session_dir": str(self.session_dir),
            "plan": synthetic_plan,
            "plan_text": render_plan(synthetic_plan),
            "critic_verdict": None,
            "completed_steps": self.state.current_step_index,
            "blocked": False,
            "dry_run_result": _primary_dry_run_result(dry_run_results),
            "dry_run_results": dry_run_results,
            "runtime_result": runtime_result,
            "preview_submit": None,
            "results": results,
            "assistant_output": loop_result["assistant_text"],
            "tool_requests": tool_requests,
            "tool_outcomes": tool_outcomes,
            "loop_state": {
                "stop_reason": loop_result["stop_reason"],
                "model_steps": loop_result["model_steps"],
                "tool_calls": len(
                    [
                        outcome
                        for outcome in tool_outcomes
                        if outcome.status != "skipped"
                    ]
                ),
                "limit_reason": loop_result["limit_reason"],
            },
            "final_message": loop_result["assistant_text"],
            "limit_reason": loop_result["limit_reason"],
            "advisory_only": not tool_requests,
            "is_chitchat": (
                _is_chitchat_request(request) and not tool_requests
            ),
            "approval_mode": policy.mode.value,
            "driving_mode": policy.mode == PermissionMode.DRIVING,
            "yolo": policy.yolo,
            "denials_count": loop_result["denials_count"],
            "approvals_count": loop_result["approvals_count"],
        }

    def _continue_run(
        self,
        completed_results: list[Any],
        dry_submit: bool,
        pause_before_risky: bool,
        allow_remote_unknown: bool,
        allow_critic_override: bool,
        rerender_plan: bool,
    ) -> dict[str, Any]:
        assert self.state is not None
        assert self.state.plan is not None
        assert self.session_dir is not None
        assert self.decision_log is not None

        plan = self.state.plan
        results = list(completed_results)
        dry_run_results = self._collect_prior_results(results, "dry_run_input")
        runtime_result = self._find_prior_result(results, "validate_runtime")
        verdict: CriticVerdict | None = self._get_logged_verdict()
        blocked = False
        block_reason: str | None = None
        run_error: Exception | None = None
        paused_for_approval = False

        if not plan.steps:
            is_chitchat = plan.is_chitchat()
            self._finalize_session(
                verdict=None,
                blocked=False,
                block_reason=None,
                dry_run_results=dry_run_results,
                advisory_only=True,
                is_chitchat=is_chitchat,
                rationale=plan.rationale,
            )
            return {
                "session_id": self.state.session_id,
                "session_dir": str(self.session_dir),
                "plan": plan,
                "plan_text": render_plan(plan) if rerender_plan else None,
                "critic_verdict": None,
                "completed_steps": self.state.current_step_index,
                "blocked": False,
                "dry_run_result": None,
                "dry_run_results": dry_run_results,
                "runtime_result": runtime_result,
                "preview_submit": None,
                "results": results,
                "advisory_only": True,
                "is_chitchat": is_chitchat,
            }
        if (
            self.state.current_step_index >= len(plan.steps)
            and self._get_logged_summary() is not None
        ):
            return {
                "session_id": self.state.session_id,
                "session_dir": str(self.session_dir),
                "plan": plan,
                "plan_text": render_plan(plan) if rerender_plan else None,
                "critic_verdict": verdict,
                "completed_steps": self.state.current_step_index,
                "blocked": False,
                "dry_run_result": _primary_dry_run_result(dry_run_results),
                "dry_run_results": dry_run_results,
                "runtime_result": runtime_result,
                "preview_submit": None,
                "results": results,
            }

        preview_submit = None
        risky_start = len(plan.steps)
        try:
            for step_index in range(
                self.state.current_step_index, len(plan.steps)
            ):
                step = plan.steps[step_index]
                if step.tool in _RISKY_TOOLS:
                    risky_start = step_index
                    break
                result = self._execute_step(step_index, step, results)
                results.append(result)
                self.state.current_step_index = step_index + 1
                self._save_state()
                if step.tool == "dry_run_input":
                    dry_run_results.append(result)
                elif step.tool == "validate_runtime":
                    runtime_result = result

            if (
                risky_start < len(plan.steps)
                and plan.steps[risky_start].tool == "submit_hpc"
                and not dry_submit
            ):
                preview_submit = self._preview_submit_step(
                    risky_start,
                    plan.steps[risky_start],
                    results,
                    dry_submit=dry_submit,
                )

            if verdict is None:
                verdict = self._critic_call(
                    plan=plan,
                    dry_run_results=dry_run_results,
                    dry_submit=dry_submit,
                )
                verdict = self._apply_deterministic_gates(
                    plan=plan,
                    verdict=verdict,
                    runtime_result=runtime_result,
                    dry_run_results=dry_run_results,
                    preview_submit=preview_submit,
                    dry_submit=dry_submit,
                )
                self.decision_log.write(
                    "critic_verdict",
                    verdict.model_dump(),
                    rationale=verdict.rationale,
                )

            block_reason = self._block_reason(
                verdict=verdict,
                dry_submit=dry_submit,
                allow_remote_unknown=allow_remote_unknown,
                allow_critic_override=allow_critic_override,
            )
            blocked = block_reason is not None
            if blocked:
                return {
                    "session_id": self.state.session_id,
                    "session_dir": str(self.session_dir),
                    "plan": plan,
                    "plan_text": render_plan(plan) if rerender_plan else None,
                    "critic_verdict": verdict,
                    "completed_steps": self.state.current_step_index,
                    "blocked": True,
                    "dry_run_result": _primary_dry_run_result(dry_run_results),
                    "dry_run_results": dry_run_results,
                    "runtime_result": runtime_result,
                    "preview_submit": preview_submit,
                }

            if (
                pause_before_risky
                and risky_start < len(plan.steps)
                and not (
                    dry_submit and plan.steps[risky_start].tool == "submit_hpc"
                )
            ):
                paused_for_approval = True
                return {
                    "session_id": self.state.session_id,
                    "session_dir": str(self.session_dir),
                    "plan": plan,
                    "plan_text": render_plan(plan) if rerender_plan else None,
                    "critic_verdict": verdict,
                    "completed_steps": self.state.current_step_index,
                    "blocked": False,
                    "dry_run_result": _primary_dry_run_result(dry_run_results),
                    "dry_run_results": dry_run_results,
                    "runtime_result": runtime_result,
                    "preview_submit": preview_submit,
                    "pending_approval": True,
                    "next_risky_tool": plan.steps[risky_start].tool,
                }

            for step_index in range(risky_start, len(plan.steps)):
                if step_index < self.state.current_step_index:
                    continue
                step = plan.steps[step_index]
                if step.tool == "submit_hpc" and dry_submit:
                    self._record_skipped_step(
                        step_index,
                        step,
                        reason="dry-submit skips remote submission",
                    )
                    preview_submit = {
                        "skipped": True,
                        "skip_reason": "dry-submit skips remote submission",
                    }
                    break
                extra_kwargs = {}
                if step.tool == "submit_hpc":
                    extra_kwargs["execute"] = not dry_submit
                    if self.transport is not None:
                        extra_kwargs["transport"] = self.transport
                result = self._execute_step(
                    step_index,
                    step,
                    results,
                    extra_kwargs=extra_kwargs,
                )
                results.append(result)
                self.state.current_step_index = step_index + 1
                self._save_state()

            return {
                "session_id": self.state.session_id,
                "session_dir": str(self.session_dir),
                "plan": plan,
                "plan_text": render_plan(plan) if rerender_plan else None,
                "critic_verdict": verdict,
                "completed_steps": self.state.current_step_index,
                "blocked": False,
                "dry_run_result": _primary_dry_run_result(dry_run_results),
                "dry_run_results": dry_run_results,
                "runtime_result": runtime_result,
                "preview_submit": preview_submit,
                "results": results,
            }
        except Exception as exc:
            run_error = exc
            raise
        finally:
            if not paused_for_approval:
                self._finalize_session(
                    verdict=verdict,
                    blocked=blocked,
                    block_reason=block_reason,
                    dry_run_results=dry_run_results,
                    run_error=run_error,
                )

    def _provider_instance(self) -> Any:
        if self._provider is None:
            self._provider = get_provider()
        return self._provider

    def _planner_call(self, request: str) -> Plan:
        current_turn_index = self.state.turn_index if self.state else None
        prompt = _load_prompt("planner.md")
        tool_defs = self.registry.openai_tool_defs()
        plan = self._llm_json_call(
            stage="planner",
            messages=[
                {"role": "system", "content": prompt},
                {
                    "role": "user",
                    "content": json.dumps(
                        {
                            "request": request,
                            "conversation_history": (
                                self.conversation_history.prompt_context(
                                    current_turn_index=current_turn_index
                                )
                            ),
                            "tools": tool_defs,
                        }
                    ),
                },
            ],
            model_cls=Plan,
        )
        for step in plan.steps:
            step.args = self.registry.normalize_args(step.tool, step.args)
        plan.intent = _resolve_plan_intent(request, plan)
        return plan

    def _critic_call(
        self,
        plan: Plan,
        dry_run_results: list[dict[str, Any]],
        dry_submit: bool,
    ) -> CriticVerdict:
        prompt = _load_prompt("critic.md")
        return self._llm_json_call(
            stage="critic",
            messages=[
                {"role": "system", "content": prompt},
                {
                    "role": "user",
                    "content": json.dumps(
                        {
                            "plan": plan.model_dump(),
                            "dry_run_inputs": dry_run_results,
                            "dry_submit": dry_submit,
                        }
                    ),
                },
            ],
            model_cls=CriticVerdict,
        )

    def _llm_json_call(
        self,
        *,
        stage: str,
        messages: list[dict[str, Any]],
        model_cls: type[BaseModel],
    ) -> Any:
        provider = self._provider_instance()
        last_error: Exception | None = None

        for attempt in range(1, _LLM_MAX_ATTEMPTS + 1):
            response = None
            started = time.perf_counter()
            try:
                response = provider.chat(
                    messages,
                    tools=None,
                    timeout_s=DEFAULT_TIMEOUT_S,
                )
                parsed = _parse_json_response(response)
                model = model_cls.model_validate(parsed)
                self._record_llm_stats(
                    stage=stage,
                    attempt=attempt,
                    provider=provider,
                    messages=messages,
                    response=response,
                    latency_ms=_elapsed_ms(started),
                )
                return model
            except Exception as exc:
                last_error = exc
                self._log_llm_failure(
                    stage=stage,
                    attempt=attempt,
                    provider=provider,
                    messages=messages,
                    response=response,
                    error=exc,
                    latency_ms=_elapsed_ms(started),
                )
                if attempt >= _LLM_MAX_ATTEMPTS:
                    break
                time.sleep(_LLM_BACKOFF_SECONDS[attempt - 1])

        assert last_error is not None
        raise last_error

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
        usage = extract_response_usage(response)
        input_tokens = usage["input_tokens"]
        output_tokens = usage["output_tokens"]
        raw_response = _stringify_response(response)
        if input_tokens is None:
            input_tokens = _estimate_tokens(messages)
        if output_tokens is None:
            output_tokens = _estimate_tokens(raw_response)

        self._llm_stats.append(
            {
                "stage": stage,
                "attempt": attempt,
                "provider_name": getattr(provider, "name", None),
                "resolved_model": _resolved_model_for_response(
                    provider, response
                ),
                "input_tokens": input_tokens,
                "output_tokens": output_tokens,
                "latency_ms": latency_ms,
                "success": True,
            }
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
        if self.decision_log is None:
            pass

        raw_response = _stringify_response(response)
        usage = extract_response_usage(response)
        input_tokens = usage["input_tokens"]
        output_tokens = usage["output_tokens"]
        if input_tokens is None:
            input_tokens = _estimate_tokens(messages)
        if output_tokens is None:
            output_tokens = _estimate_tokens(raw_response)

        self._llm_stats.append(
            {
                "stage": stage,
                "attempt": attempt,
                "provider_name": getattr(provider, "name", None),
                "resolved_model": _resolved_model_for_response(
                    provider, response
                ),
                "input_tokens": input_tokens,
                "output_tokens": output_tokens,
                "latency_ms": latency_ms,
                "success": False,
            }
        )

        if self.decision_log is None:
            return

        self.decision_log.write(
            "llm_error",
            {
                "stage": stage,
                "attempt": attempt,
                "provider_name": getattr(provider, "name", None),
                "resolved_model": _resolved_model_for_response(
                    provider, response
                ),
                "error_type": error.__class__.__name__,
                "message": str(error),
                "messages": messages,
                "raw_response": raw_response,
                "step_wall_time_ms": latency_ms,
            },
            rationale=f"{stage} attempt {attempt} failed",
        )

    def _execute_step(
        self,
        step_index: int,
        step: Step,
        prior_results: list[Any],
        extra_kwargs: dict[str, Any] | None = None,
    ) -> Any:
        assert self.session_dir is not None
        assert self.decision_log is not None
        resolved_args = _resolve_refs(
            step.args,
            prior_results,
            handle_store=self.handle_store,
        )
        if extra_kwargs:
            resolved_args.update(extra_kwargs)
        ts_start = _utc_now_iso()
        step_start_time = time.perf_counter()
        self.decision_log.write(
            "tool_call",
            {
                "step_index": step_index,
                "tool": step.tool,
                "args": _preview_value(resolved_args),
                "ts_start": ts_start,
                "step_wall_time_ms": None,
            },
            rationale=step.rationale,
        )
        result = self.registry.call(step.tool, resolved_args)
        if _is_tool_error(result):
            error = result["error"]
            self.decision_log.write(
                "tool_error",
                {
                    "step_index": step_index,
                    "tool": step.tool,
                    "error_type": error.get("type", "RuntimeError"),
                    "message": error.get("message", "Unknown tool error"),
                    "payload": _preview_value(result),
                    "ts_start": ts_start,
                    "ts_end": _utc_now_iso(),
                    "step_wall_time_ms": _elapsed_ms(step_start_time),
                },
                rationale=step.rationale,
            )
            raise RuntimeError(result["error"]["message"])
        if _run_local_failed(step.tool, result):
            artifact_path = self._write_result_artifact(step_index, result)
            message = _run_local_failure_message(result)
            self.decision_log.write(
                "tool_error",
                {
                    "step_index": step_index,
                    "tool": step.tool,
                    "artifact": artifact_path.name,
                    "error_type": "RuntimeError",
                    "message": message,
                    "payload": _preview_value(result),
                    "ts_start": ts_start,
                    "ts_end": _utc_now_iso(),
                    "step_wall_time_ms": _elapsed_ms(step_start_time),
                },
                rationale=step.rationale,
            )
            raise RuntimeError(message)
        artifact_path = self._write_result_artifact(step_index, result)
        handle_id = self._store_result_handle(step.tool, result)
        self.decision_log.write(
            "tool_result",
            {
                "step_index": step_index,
                "tool": step.tool,
                "artifact": artifact_path.name,
                "handle_id": handle_id,
                "payload": _preview_value(result),
                "ts_end": _utc_now_iso(),
                "step_wall_time_ms": _elapsed_ms(step_start_time),
            },
            rationale=step.rationale,
        )
        return result

    def _preview_submit_step(
        self,
        step_index: int,
        step: Step,
        prior_results: list[Any],
        dry_submit: bool,
    ) -> Any:
        assert self.decision_log is not None
        resolved_args = _resolve_refs(
            step.args,
            prior_results,
            handle_store=self.handle_store,
        )
        resolved_args["execute"] = False
        self.decision_log.write(
            "tool_preview",
            {
                "step_index": step_index,
                "tool": step.tool,
                "args": _preview_value(resolved_args),
            },
            rationale=step.rationale,
        )
        result = self.registry.call(step.tool, resolved_args)
        if _is_tool_error(result):
            message = result["error"]["message"]
            if dry_submit and _is_submit_server_error(message):
                result = {
                    "transport": None,
                    "script_path": None,
                    "script_bytes": None,
                    "command_executed": None,
                    "job_id": None,
                    "duplicate_check": {
                        "duplicate": False,
                        "message": None,
                    },
                    "skipped": True,
                    "skip_reason": message,
                }
            else:
                raise RuntimeError(message)
        self.decision_log.write(
            "tool_preview_result",
            {
                "step_index": step_index,
                "tool": step.tool,
                "payload": _preview_value(result),
            },
            rationale=step.rationale,
        )
        return result

    def _record_skipped_step(
        self,
        step_index: int,
        step: Step,
        reason: str,
    ) -> None:
        assert self.decision_log is not None
        self.decision_log.write(
            "tool_skipped",
            {
                "step_index": step_index,
                "tool": step.tool,
                "reason": reason,
            },
            rationale=step.rationale,
        )

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
        assert self.state is not None
        assert self.session_dir is not None
        assert self.decision_log is not None

        ended_at = _utc_now_iso()
        effective_blocked = blocked or run_error is not None
        effective_block_reason = block_reason
        if run_error is not None and effective_block_reason is None:
            effective_block_reason = (
                f"exception:{run_error.__class__.__name__}"
            )

        wall_time_ms = _elapsed_ms(
            self._run_start_time,
            started_at=self.state.request_started_at,
            ended_at=ended_at,
        )
        summary = {
            "total_steps_executed": self.state.current_step_index,
            "total_steps_planned": self.state.total_steps_planned,
            "blocked": effective_blocked,
            "block_reason": effective_block_reason or "unknown",
            "wall_time_ms": wall_time_ms,
            "tools_called": self._tools_called(),
            "critic_confidence": (
                verdict.confidence if verdict is not None else 0.0
            ),
            "request_intent": self.state.request_intent,
            "provider_name": self._metadata_provider_name(),
            "resolved_model": self._metadata_resolved_model(),
            "total_input_tokens": self._total_llm_tokens("input_tokens"),
            "total_output_tokens": self._total_llm_tokens("output_tokens"),
            "exit_status": (
                "error"
                if run_error is not None
                else "blocked" if blocked else "ok"
            ),
            "advisory_only": advisory_only,
            "is_chitchat": is_chitchat,
            "rationale": rationale or "",
        }
        self.decision_log.write(
            "session_summary",
            summary,
            rationale=rationale or "",
        )
        self._refresh_conversation_history()

        primary_dry_run_result = _primary_dry_run_result(dry_run_results)
        schema_hash = _schema_hash(self.registry.openai_tool_defs())
        metadata = {
            "session_id": self.state.session_id,
            "request": self.state.request or "unknown",
            "intent": self.state.request_intent,
            "request_intent": self.state.request_intent,
            "plan_steps": self.state.total_steps_planned,
            "executed_steps": self.state.current_step_index,
            "total_steps_planned": summary["total_steps_planned"],
            "total_steps_executed": summary["total_steps_executed"],
            "input_file": (
                str(primary_dry_run_result.get("inputfile"))
                if primary_dry_run_result
                and primary_dry_run_result.get("inputfile")
                else "unknown"
            ),
            "input_files": [
                str(result.get("inputfile"))
                for result in dry_run_results
                if result.get("inputfile")
            ],
            "critic_verdict": (
                verdict.verdict if verdict is not None else "unknown"
            ),
            "critic_confidence": summary["critic_confidence"],
            "blocked": summary["blocked"],
            "block_reason": summary["block_reason"],
            "wall_time_ms": wall_time_ms,
            "started_at": self.state.started_at,
            "request_started_at": self.state.request_started_at,
            "ended_at": ended_at,
            "provider_name": summary["provider_name"],
            "resolved_model": summary["resolved_model"],
            "git_sha": _git_sha() or "unknown",
            "schema_hash": schema_hash,
            "total_input_tokens": summary["total_input_tokens"],
            "total_output_tokens": summary["total_output_tokens"],
            "tools_called": summary["tools_called"],
            "exit_status": summary["exit_status"],
            "advisory_only": advisory_only,
            "is_chitchat": is_chitchat,
            "rationale": rationale or "",
        }
        (self.session_dir / "session_metadata.json").write_text(
            json.dumps(metadata, indent=2, sort_keys=True),
            encoding="utf-8",
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
        if hasattr(self.registry, "tool_defs_for_provider"):
            return self.registry.tool_defs_for_provider(provider_name)
        return self.registry.openai_tool_defs()

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
        state_path = self.session_dir / "session.json"
        if not state_path.exists():
            state_path = self.session_dir / "state.json"
        self.state = SessionState.load(state_path)
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
        if self.handle_store is None:
            return None
        kind = _result_handle_kind(tool_name, result)
        if kind is None:
            return None
        return self.handle_store.put(
            kind=kind,
            obj=result,
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


def render_plan(plan: Plan) -> str:
    lines = ["Plan:"]
    if plan.rationale:
        lines.append(f"Rationale: {plan.rationale}")
    if plan.estimated_cost:
        lines.append(f"Estimated cost: {plan.estimated_cost}")
    for index, step in enumerate(plan.steps, start=1):
        lines.append(
            f"{index}. {step.tool} {json.dumps(_preview_value(step.args), sort_keys=True)}"
        )
        if step.rationale:
            lines.append(f"   - {step.rationale}")
    return "\n".join(lines)


def _synthetic_plan_from_tool_requests(tool_requests: list[Any]) -> Plan:
    steps = [
        Step(
            tool=request.name,
            args=dict(request.arguments),
            rationale="",
        )
        for request in tool_requests
    ]
    return Plan(
        steps=steps,
        rationale="Synthetic plan projected from tool_use_request entries.",
        intent="workflow" if steps else "advisory",
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
        for key in ("provider", "registry", "session_root", "transport")
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


def _classify_intent(request: str) -> str:
    normalized = request.lower()
    matches = [
        intent
        for intent, patterns in _INTENT_PATTERNS.items()
        if any(re.search(pattern, normalized) for pattern in patterns)
    ]
    if not matches:
        return "unknown"
    unique_matches = list(dict.fromkeys(matches))
    if len(unique_matches) == 1:
        return unique_matches[0]

    has_composite_marker = any(
        marker in normalized for marker in ("+", " then ", " and ", " after ")
    )
    non_opt = [intent for intent in unique_matches if intent != "opt"]
    if len(non_opt) == 1 and not has_composite_marker:
        return non_opt[0]
    return "composite"


def _resolve_plan_intent(
    request: str,
    plan: Plan,
) -> Literal["workflow", "advisory", "chitchat"]:
    if plan.steps:
        return "workflow"
    if _is_chitchat_request(request):
        return "chitchat"
    if plan.intent in {"workflow", "advisory", "chitchat"}:
        return plan.intent
    return "advisory"


def _is_chitchat_request(request: str) -> bool:
    normalized = re.sub(r"\s+", " ", request).strip().lower()
    if not normalized:
        return False
    if any(
        re.fullmatch(pattern, normalized)
        for pattern in _CHITCHAT_IDENTITY_PATTERNS
    ):
        return True
    if any(
        re.fullmatch(pattern, normalized)
        for pattern in _CHITCHAT_EXACT_PATTERNS
    ):
        return True
    tokens = re.findall(r"[a-z]+", normalized)
    return (
        bool(tokens)
        and len(tokens) <= 3
        and set(tokens).issubset(_CHITCHAT_TOKENS)
    )


def _load_prompt(name: str) -> str:
    prompt_path = Path(__file__).with_name("prompts") / name
    return prompt_path.read_text(encoding="utf-8")


def _elapsed_ms(
    start_time: float | None,
    *,
    started_at: str | None = None,
    ended_at: str | None = None,
) -> int:
    if started_at is not None and ended_at is not None:
        started = datetime.fromisoformat(started_at)
        ended = datetime.fromisoformat(ended_at)
        return max(0, int(round((ended - started).total_seconds() * 1000)))
    if start_time is not None:
        return max(0, int(round((time.perf_counter() - start_time) * 1000)))
    return 0


def _parse_json_response(response: Any) -> dict[str, Any]:
    if (
        isinstance(response, dict)
        and "parsed" in response
        and isinstance(response["parsed"], dict)
    ):
        return response["parsed"]
    if (
        isinstance(response, dict)
        and "json" in response
        and isinstance(response["json"], dict)
    ):
        return response["json"]
    if (
        isinstance(response, dict)
        and "content" in response
        and isinstance(response["content"], dict)
    ):
        return response["content"]
    text = _extract_text(response)
    text = text.strip()
    if text.startswith("```"):
        text = re.sub(r"^```(?:json)?\s*", "", text)
        text = re.sub(r"\s*```$", "", text)
    try:
        return json.loads(text)
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"LLM returned invalid JSON: {exc}\nRaw: {text!r}"
        ) from exc


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


def _is_submit_server_error(message: str) -> bool:
    return (
        "submit_hpc requires server when no configured servers are available"
        in message
        or "submit_hpc requires server when multiple configured servers are "
        in message
    )


def _extract_text(response: Any) -> str:
    if isinstance(response, str):
        return response
    if isinstance(response, dict):
        content = response.get("content")
        if isinstance(content, str):
            return content
        if isinstance(content, list):
            parts: list[str] = []
            for item in content:
                if isinstance(item, dict) and item.get("type") == "text":
                    parts.append(item.get("text", ""))
            if parts:
                return "\n".join(parts)
        choices = response.get("choices") or []
        if choices:
            message = choices[0].get("message", {})
            content = message.get("content", "")
            if isinstance(content, str):
                return content
            if isinstance(content, list):
                return "\n".join(
                    item.get("text", "")
                    for item in content
                    if isinstance(item, dict)
                )
    raise ValueError("Could not extract text from provider response")


def _stringify_response(response: Any) -> str | None:
    if response is None:
        return None
    if isinstance(response, str):
        return response
    try:
        return json.dumps(_json_safe(response), ensure_ascii=False)
    except TypeError:
        return repr(response)


def _estimate_tokens(value: Any) -> int:
    if value is None:
        return 0
    if isinstance(value, str):
        text = value
    else:
        text = json.dumps(_json_safe(value), sort_keys=True)
    return max(1, int(round(len(text) / 4)))


def _resolved_model_for_response(provider: Any, response: Any) -> str | None:
    fallback = getattr(provider, "default_model", None)
    if isinstance(response, dict):
        model = response.get("model")
        if isinstance(model, str) and model.strip():
            return model
    return fallback


def _schema_hash(tool_defs: list[dict[str, Any]]) -> str:
    payload = json.dumps(tool_defs, sort_keys=True).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _git_sha() -> str | None:
    repo_root = Path(__file__).resolve().parents[2]
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_root,
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    git_sha = result.stdout.strip()
    return git_sha or None


def _resolve_refs(
    value: Any,
    prior_results: list[Any],
    handle_store: HandleStore | None = None,
) -> Any:
    if isinstance(value, str):
        return _resolve_ref_string(
            value,
            prior_results,
            handle_store=handle_store,
        )
    if isinstance(value, list):
        return [
            _resolve_refs(
                item,
                prior_results,
                handle_store=handle_store,
            )
            for item in value
        ]
    if isinstance(value, dict):
        return {
            key: _resolve_refs(
                item,
                prior_results,
                handle_store=handle_store,
            )
            for key, item in value.items()
        }
    return value


def _resolve_ref_string(
    value: str,
    prior_results: list[Any],
    handle_store: HandleStore | None = None,
) -> Any:
    match = _REFERENCE_RE.match(value)
    if match is not None:
        index = int(match.group("index")) - 1
        if index < 0 or index >= len(prior_results):
            raise IndexError(f"Step reference {value!r} is out of range")
        resolved = prior_results[index]
        for part in match.group("path").split("."):
            if not part:
                continue
            if isinstance(resolved, dict):
                resolved = resolved[part]
            else:
                resolved = getattr(resolved, part)
        return _restore_json_result(resolved)

    if handle_store is not None and is_handle_id(value):
        try:
            return handle_store.get(value)
        except KeyError:
            return _restore_json_result(handle_store.get_summary(value))

    return value


@lru_cache(maxsize=1)
def _molecule_class():
    from chemsmart.io.molecules.structure import Molecule

    return Molecule


def _is_instance_of(
    value: Any,
    *,
    module_name: str,
    class_names: tuple[str, ...],
) -> bool:
    value_type = type(value)
    return (
        value_type.__module__ == module_name
        and value_type.__name__ in class_names
    )


def _json_safe(value: Any) -> Any:
    if _is_instance_of(
        value,
        module_name="chemsmart.io.molecules.structure",
        class_names=("Molecule",),
    ):
        positions = value.positions
        if hasattr(positions, "tolist"):
            positions = positions.tolist()
        return {
            "__chemsmart_type__": "molecule",
            "symbols": _json_safe(list(value.symbols)),
            "positions": _json_safe(positions),
            "charge": value.charge,
            "multiplicity": value.multiplicity,
            "frozen_atoms": _json_safe(value.frozen_atoms),
            "pbc_conditions": _json_safe(value.pbc_conditions),
            "translation_vectors": _json_safe(value.translation_vectors),
            "energy": value.energy,
            "forces": _json_safe(value.forces),
            "velocities": _json_safe(value.velocities),
            "info": _json_safe(value.info),
        }
    if _is_instance_of(
        value,
        module_name="chemsmart.jobs.gaussian.settings",
        class_names=("GaussianJobSettings",),
    ) or _is_instance_of(
        value,
        module_name="chemsmart.jobs.orca.settings",
        class_names=("ORCAJobSettings",),
    ):
        return {
            "__chemsmart_type__": "settings",
            "module": value.__class__.__module__,
            "class": value.__class__.__name__,
            **{
                key: _json_safe(item)
                for key, item in vars(value).items()
                if not key.startswith("_")
            },
        }
    if (
        value.__class__.__module__.startswith("chemsmart.jobs.")
        and hasattr(value, "molecule")
        and hasattr(value, "settings")
    ):
        return {
            "__chemsmart_type__": "job",
            "module": value.__class__.__module__,
            "class": value.__class__.__name__,
            "molecule": _json_safe(value.molecule),
            "settings": _json_safe(value.settings),
            "label": value.label,
            "folder": value.folder,
        }
    if isinstance(value, BaseModel):
        return _json_safe(value.model_dump())
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, bytes):
        return {"type": "bytes", "length": len(value)}
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return {"type": value.__class__.__name__, "repr": repr(value)}


def _preview_value(value: Any) -> Any:
    return _json_safe(value)


def _restore_json_result(value: Any) -> Any:
    if isinstance(value, list):
        return [_restore_json_result(item) for item in value]
    if not isinstance(value, dict):
        return value

    marker = value.get("__chemsmart_type__")
    if marker == "molecule":
        molecule_cls = _molecule_class()
        return molecule_cls(
            symbols=value.get("symbols"),
            positions=value.get("positions"),
            charge=value.get("charge"),
            multiplicity=value.get("multiplicity"),
            frozen_atoms=value.get("frozen_atoms"),
            pbc_conditions=value.get("pbc_conditions"),
            translation_vectors=value.get("translation_vectors"),
            energy=value.get("energy"),
            forces=value.get("forces"),
            velocities=value.get("velocities"),
            info=value.get("info"),
        )

    if marker == "settings":
        settings_cls = _load_class(value["module"], value["class"])
        kwargs = {
            key: _restore_json_result(item)
            for key, item in value.items()
            if key not in {"__chemsmart_type__", "module", "class"}
        }
        return settings_cls(**kwargs)

    if marker == "job":
        job_cls = _load_class(value["module"], value["class"])
        job = job_cls(
            molecule=_restore_json_result(value["molecule"]),
            settings=_restore_json_result(value["settings"]),
            label=value["label"],
            jobrunner=None,
        )
        if value.get("folder"):
            job.set_folder(value["folder"])
        return job

    return {key: _restore_json_result(item) for key, item in value.items()}


def _load_class(module_name: str, class_name: str) -> type[Any]:
    module = importlib.import_module(module_name)
    return getattr(module, class_name)


def _is_tool_error(result: Any) -> bool:
    return (
        isinstance(result, dict)
        and result.get("ok") is False
        and "error" in result
    )


def _primary_dry_run_result(
    dry_run_results: list[dict[str, Any]],
) -> dict[str, Any] | None:
    return dry_run_results[0] if dry_run_results else None


def _run_local_failed(tool_name: str, result: Any) -> bool:
    return (
        tool_name == "run_local"
        and isinstance(result, dict)
        and result.get("ok") is False
    )


def _run_local_failure_message(result: dict[str, Any]) -> str:
    message = "run_local failed"
    returncode = result.get("returncode")
    if returncode is not None:
        message += f" with returncode {returncode}"
    stderr_path = result.get("stderr_path")
    if isinstance(stderr_path, str) and stderr_path.strip():
        message += f"; see {stderr_path}"
    return message


def _result_handle_kind(tool_name: str, result: Any) -> str | None:
    if tool_name == "build_molecule":
        return "mol"
    if tool_name == "build_gaussian_settings":
        return "gset"
    if tool_name == "build_orca_settings":
        return "oset"
    if tool_name == "build_job":
        return "job"
    if tool_name == "dry_run_input" and isinstance(result, dict):
        return "dryrun"
    if tool_name == "validate_runtime" and isinstance(result, dict):
        return "runtime"
    if tool_name == "run_local" and isinstance(result, dict):
        return "runresult"
    if tool_name == "extract_optimized_geometry":
        return "geom"
    if (
        tool_name == "submit_hpc"
        and isinstance(result, dict)
        and "job_id" in result
    ):
        return "submit"
    if tool_name == "recommend_method" and isinstance(result, dict):
        return "recmethod"
    return None


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
