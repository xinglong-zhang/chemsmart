"""Credential-scoped provider-neutral Runtime V2 session runner."""

from __future__ import annotations

from dataclasses import replace
from typing import Any, Callable

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.request_context import (
    ProviderNetworkBudgetV1,
    RequestContextProvenanceV1,
)
from chemsmart.agent.api_access import SecretLease
from chemsmart.agent.loop import ToolLoopResultV1, ToolLoopRunner
from chemsmart.agent.runtime.contracts import TaskEnvelopeV1
from chemsmart.agent.runtime.deepseek import (
    DeepSeekHttpsTransport,
    DeepSeekV4ToolSession,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1


class UnifiedSessionRunner:
    """Run one active provider session inside a one-use secret lease."""

    def __init__(
        self,
        *,
        host: CommandCompiledToolHostV1,
        event_store: RuntimeEventStore,
        credential_lease: SecretLease,
        provider_config: Any,
    ) -> None:
        self.host = host
        self.event_store = event_store
        self.credential_lease = credential_lease
        self.provider_config = provider_config
        if credential_lease.provider != self.provider_config.provider:
            raise ContractError("credential lease belongs to another provider")

    def run(
        self,
        *,
        messages: list[dict[str, Any]],
        envelope: TaskEnvelopeV1,
        request_context: RequestContextProvenanceV1,
        provider_budget: ProviderNetworkBudgetV1,
        should_stop: Callable[[], bool] | None = None,
    ) -> ToolLoopResultV1:
        if not messages or not all(
            isinstance(item, dict) and item.get("role") in {
                "system",
                "user",
                "assistant",
            }
            for item in messages
        ):
            raise ContractError("initial provider messages are malformed")

        def _leased_run(secret: str) -> ToolLoopResultV1:
            approved_output_limit = min(
                envelope.budget.max_output_tokens_per_request,
                provider_budget.max_output_tokens_per_request,
                self.provider_config.max_output_tokens,
            )
            bound_config = replace(
                self.provider_config,
                max_output_tokens=approved_output_limit,
            )
            turn_deadlines = bound_config.turn_deadlines
            if bound_config.provider == "deepseek":
                transport = DeepSeekHttpsTransport(
                    api_key=secret,
                    endpoint=bound_config.endpoint,
                    turn_deadlines=turn_deadlines,
                )
                session = DeepSeekV4ToolSession(
                    transport=transport,
                    messages=messages,
                    config=bound_config,
                )
            elif bound_config.provider == "alibaba-token-plan":
                from chemsmart.agent.runtime.alibaba import (
                    AlibabaTokenPlanHttpsTransport,
                    AlibabaTokenPlanToolSession,
                )

                transport = AlibabaTokenPlanHttpsTransport(
                    api_key=secret,
                    endpoint=bound_config.endpoint,
                    turn_deadlines=turn_deadlines,
                )
                session = AlibabaTokenPlanToolSession(
                    transport=transport,
                    messages=messages,
                    config=bound_config,
                )
            else:
                raise ContractError("active provider has no registered runner")
            try:
                return ToolLoopRunner(
                    host=self.host,
                    event_store=self.event_store,
                ).run(
                    session=session,
                    envelope=envelope,
                    request_context=request_context,
                    provider_budget=provider_budget,
                    should_stop=should_stop,
                )
            finally:
                session.close()
                transport.close()

        return self.credential_lease.invoke(_leased_run)


__all__ = ["UnifiedSessionRunner"]
