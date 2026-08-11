"""Credential-scoped provider-neutral Runtime V2 session runner."""

from __future__ import annotations

from dataclasses import replace
from typing import Any

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.adaptive_api_campaign import (
    AdaptiveHypothesisV1,
    AdaptiveNetworkBudgetV1,
)
from chemsmart.agent.api_access import SecretLease
from chemsmart.agent.feedback import FULL_FEEDBACK_V1
from chemsmart.agent.loop import ToolLoopResultV1, ToolLoopRunner
from chemsmart.agent.runtime.contracts import TaskEnvelopeV1
from chemsmart.agent.runtime.deepseek import (
    DeepSeekHttpsTransport,
    DeepSeekV4FlashConfigV1,
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
        provider_config: Any | None = None,
    ) -> None:
        self.host = host
        self.event_store = event_store
        self.credential_lease = credential_lease
        self.provider_config = provider_config or DeepSeekV4FlashConfigV1()
        if credential_lease.provider != self.provider_config.provider:
            raise ContractError("credential lease belongs to another provider")

    def run(
        self,
        *,
        messages: list[dict[str, Any]],
        envelope: TaskEnvelopeV1,
        hypothesis: AdaptiveHypothesisV1,
        network_budget: AdaptiveNetworkBudgetV1,
        feedback_projection: str = FULL_FEEDBACK_V1,
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
                network_budget.max_output_tokens_per_request,
                self.provider_config.max_output_tokens,
            )
            bound_config = replace(
                self.provider_config,
                max_output_tokens=approved_output_limit,
            )
            if bound_config.provider == "deepseek":
                transport = DeepSeekHttpsTransport(
                    api_key=secret,
                    endpoint=bound_config.endpoint,
                    timeout_seconds=network_budget.task_wall_time_seconds,
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
                    timeout_seconds=network_budget.task_wall_time_seconds,
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
                    feedback_projection=feedback_projection,
                ).run(
                    session=session,
                    envelope=envelope,
                    hypothesis=hypothesis,
                    network_budget=network_budget,
                )
            finally:
                session.close()
                transport.close()

        return self.credential_lease.invoke(_leased_run)


__all__ = ["UnifiedSessionRunner"]
