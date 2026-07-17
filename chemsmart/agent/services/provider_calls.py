"""Provider-backed planner and critic calls for an agent session."""

from __future__ import annotations

import json
import re
import time
from typing import Any, Protocol

from pydantic import BaseModel

from chemsmart.agent.models import CriticVerdict, Plan
from chemsmart.agent.prompts import load_prompt
from chemsmart.agent.prompts.identity import build_system_prompt
from chemsmart.agent.provider_adapter import extract_response_text
from chemsmart.agent.providers import (
    DEFAULT_TIMEOUT_S,
    extract_response_usage,
)
from chemsmart.agent.services.plan_support import resolve_plan_intent
from chemsmart.agent.services.result_codec import json_safe
from chemsmart.agent.services.runtime_metrics import elapsed_ms

MAX_ATTEMPTS = 3
BACKOFF_SECONDS = (1, 2, 4)


class ProviderCallHost(Protocol):
    registry: Any
    state: Any | None
    conversation_history: Any
    decision_log: Any | None
    _llm_stats: list[dict[str, Any]]

    def _provider_instance(self) -> Any: ...
    def _prompt_session_meta(self, **extra: Any) -> dict[str, Any]: ...


class ProviderCallService:
    """Build prompts, call the configured provider, and record receipts."""

    def __init__(self, session: ProviderCallHost) -> None:
        self.session = session

    def planner(self, request: str) -> Plan:
        session = self.session
        turn_index = session.state.turn_index if session.state else None
        context = session.conversation_history.prompt_context(
            current_turn_index=turn_index
        )
        prompt = build_system_prompt(
            registry=session.registry,
            stage_instructions=load_prompt("planner.md"),
            session_meta=session._prompt_session_meta(stage="planner"),
            conversation_context=context,
        )
        plan = self.json_call(
            stage="planner",
            messages=[
                {"role": "system", "content": prompt},
                {
                    "role": "user",
                    "content": json.dumps(
                        {
                            "request": request,
                            "conversation_history": context,
                            "tools": session.registry.openai_tool_defs(),
                        }
                    ),
                },
            ],
            model_cls=Plan,
        )
        for step in plan.steps:
            step.args = session.registry.normalize_args(step.tool, step.args)
        plan.intent = resolve_plan_intent(request, plan)
        return plan

    def critic(
        self,
        plan: Plan,
        dry_run_results: list[dict[str, Any]],
        dry_submit: bool,
    ) -> CriticVerdict:
        session = self.session
        turn_index = session.state.turn_index if session.state else None
        prompt = build_system_prompt(
            registry=session.registry,
            stage_instructions=load_prompt("critic.md"),
            session_meta=session._prompt_session_meta(
                stage="critic",
                submission_mode="dry-submit" if dry_submit else "execute",
            ),
            conversation_context=session.conversation_history.prompt_context(
                current_turn_index=turn_index
            ),
        )
        return self.json_call(
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

    def json_call(
        self,
        *,
        stage: str,
        messages: list[dict[str, Any]],
        model_cls: type[BaseModel],
    ) -> Any:
        provider = self.session._provider_instance()
        last_error: Exception | None = None
        for attempt in range(1, MAX_ATTEMPTS + 1):
            response = None
            started = time.perf_counter()
            try:
                response = provider.chat(
                    messages,
                    tools=None,
                    timeout_s=DEFAULT_TIMEOUT_S,
                )
                model = model_cls.model_validate(parse_json_response(response))
                self.record_stats(
                    stage=stage,
                    attempt=attempt,
                    provider=provider,
                    messages=messages,
                    response=response,
                    latency_ms=elapsed_ms(started),
                )
                return model
            except Exception as exc:
                last_error = exc
                self.log_failure(
                    stage=stage,
                    attempt=attempt,
                    provider=provider,
                    messages=messages,
                    response=response,
                    error=exc,
                    latency_ms=elapsed_ms(started),
                )
                if attempt >= MAX_ATTEMPTS:
                    break
                time.sleep(BACKOFF_SECONDS[attempt - 1])
        assert last_error is not None
        raise last_error

    def record_stats(
        self,
        *,
        stage: str,
        attempt: int,
        provider: Any,
        messages: list[dict[str, Any]],
        response: Any,
        latency_ms: int,
    ) -> None:
        input_tokens, output_tokens = _usage(messages, response)
        self.session._llm_stats.append(
            {
                "stage": stage,
                "attempt": attempt,
                "provider_name": getattr(provider, "name", None),
                "resolved_model": resolved_model(provider, response),
                "input_tokens": input_tokens,
                "output_tokens": output_tokens,
                "latency_ms": latency_ms,
                "success": True,
            }
        )

    def log_failure(
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
        raw_response = stringify_response(response)
        input_tokens, output_tokens = _usage(messages, response)
        self.session._llm_stats.append(
            {
                "stage": stage,
                "attempt": attempt,
                "provider_name": getattr(provider, "name", None),
                "resolved_model": resolved_model(provider, response),
                "input_tokens": input_tokens,
                "output_tokens": output_tokens,
                "latency_ms": latency_ms,
                "success": False,
            }
        )
        if self.session.decision_log is None:
            return
        self.session.decision_log.write(
            "llm_error",
            {
                "stage": stage,
                "attempt": attempt,
                "provider_name": getattr(provider, "name", None),
                "resolved_model": resolved_model(provider, response),
                "error_type": error.__class__.__name__,
                "message": str(error),
                "messages": messages,
                "raw_response": raw_response,
                "step_wall_time_ms": latency_ms,
            },
            rationale=f"{stage} attempt {attempt} failed",
        )


def parse_json_response(response: Any) -> dict[str, Any]:
    for key in ("parsed", "json", "content"):
        if isinstance(response, dict) and isinstance(response.get(key), dict):
            return response[key]
    text = extract_response_text(response).strip()
    if text.startswith("```"):
        text = re.sub(r"^```(?:json)?\s*", "", text)
        text = re.sub(r"\s*```$", "", text)
    try:
        return json.loads(text)
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"LLM returned invalid JSON: {exc}\nRaw: {text!r}"
        ) from exc


def stringify_response(response: Any) -> str | None:
    if response is None:
        return None
    if isinstance(response, str):
        return response
    try:
        return json.dumps(json_safe(response), ensure_ascii=False)
    except TypeError:
        return repr(response)


def estimate_tokens(value: Any) -> int:
    if value is None:
        return 0
    text = value if isinstance(value, str) else json.dumps(
        json_safe(value), sort_keys=True
    )
    return max(1, int(round(len(text) / 4)))


def resolved_model(provider: Any, response: Any) -> str | None:
    fallback = getattr(provider, "default_model", None)
    if isinstance(response, dict):
        model = response.get("model")
        if isinstance(model, str) and model.strip():
            return model
    return fallback


def _usage(
    messages: list[dict[str, Any]], response: Any
) -> tuple[int, int]:
    usage = extract_response_usage(response)
    input_tokens = usage["input_tokens"]
    output_tokens = usage["output_tokens"]
    if input_tokens is None:
        input_tokens = estimate_tokens(messages)
    if output_tokens is None:
        output_tokens = estimate_tokens(stringify_response(response))
    return input_tokens, output_tokens


__all__ = [
    "BACKOFF_SECONDS",
    "MAX_ATTEMPTS",
    "ProviderCallService",
    "estimate_tokens",
    "parse_json_response",
    "resolved_model",
    "stringify_response",
]
