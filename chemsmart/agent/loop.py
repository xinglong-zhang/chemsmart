"""Active Runtime V2 model/tool loop for command-compiled planning."""

from __future__ import annotations

from dataclasses import dataclass
import json
import time
from typing import Any, Callable, Mapping

from chemsmart.agent._contracts import (
    ContractError,
    canonical_data,
    canonical_sha256,
    canonical_json,
)
from chemsmart.agent.adaptive_api_campaign import (
    AdaptiveHypothesisV1,
    AdaptiveNetworkBudgetV1,
    ApiAttemptReceiptV1,
    build_api_attempt_receipt,
)
from chemsmart.agent.analysis_completion import AnalysisIncompleteError
from chemsmart.agent.runtime.contracts import TaskEnvelopeV1
from chemsmart.agent.runtime.deepseek import (
    DeepSeekProtocolError,
    DeepSeekTransportError,
    ProviderTurnReceiptV1,
    public_payload,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.runtime.events import EventKind
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.agent.feedback import (
    FEEDBACK_MODES,
    FULL_FEEDBACK_V1,
    project_tool_feedback,
)


@dataclass(frozen=True)
class ToolLoopResultV1:
    schema_version: str
    session_id: str
    turn_id: str
    terminal_state: str
    final_text: str
    public_transcript: tuple[dict[str, Any], ...]
    public_transcript_sha256: str
    public_transcript_artifact_id: str
    public_transcript_artifact_sha256: str
    provider_receipt_sha256s: tuple[str, ...]
    api_attempt_receipt_sha256s: tuple[str, ...]
    successful_tool_calls: int
    failed_tool_calls: int
    event_stream_head_sha256: str
    result_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.tool-loop-result.v1":
            raise ContractError("unsupported tool-loop result schema")
        if self.terminal_state not in {
            "complete",
            "planned",
            "failed",
            "blocked",
            "waiting_for_approval",
        }:
            raise ContractError("invalid tool-loop terminal state")
        if _contains_private_reasoning(self.public_transcript):
            raise ContractError("public transcript contains private reasoning")
        if self.public_transcript_sha256 != canonical_sha256(
            self.public_transcript
        ):
            raise ContractError("public transcript digest mismatch")
        body = {
            "schema_version": self.schema_version,
            "session_id": self.session_id,
            "turn_id": self.turn_id,
            "terminal_state": self.terminal_state,
            "final_text": self.final_text,
            "public_transcript": self.public_transcript,
            "public_transcript_sha256": self.public_transcript_sha256,
            "public_transcript_artifact_id": (
                self.public_transcript_artifact_id
            ),
            "public_transcript_artifact_sha256": (
                self.public_transcript_artifact_sha256
            ),
            "provider_receipt_sha256s": self.provider_receipt_sha256s,
            "api_attempt_receipt_sha256s": (
                self.api_attempt_receipt_sha256s
            ),
            "successful_tool_calls": self.successful_tool_calls,
            "failed_tool_calls": self.failed_tool_calls,
            "event_stream_head_sha256": self.event_stream_head_sha256,
        }
        if self.result_sha256 != canonical_sha256(body):
            raise ContractError("tool-loop result digest mismatch")


class ToolLoopRunner:
    """Drive provider-native turns through the only approved dispatcher."""

    def __init__(
        self,
        *,
        host: CommandCompiledToolHostV1,
        event_store: RuntimeEventStore,
        clock: Callable[[], float] = time.monotonic,
        feedback_projection: str = FULL_FEEDBACK_V1,
    ) -> None:
        self.host = host
        self.event_store = event_store
        self.clock = clock
        self.feedback_projection = str(feedback_projection).strip().lower()
        if self.feedback_projection not in FEEDBACK_MODES:
            raise ContractError("unsupported tool-feedback projection")

    def run(
        self,
        *,
        session: Any,
        envelope: TaskEnvelopeV1,
        hypothesis: AdaptiveHypothesisV1,
        network_budget: AdaptiveNetworkBudgetV1,
    ) -> ToolLoopResultV1:
        self._validate_run_contract(
            envelope=envelope,
            hypothesis=hypothesis,
            network_budget=network_budget,
            provider=str(session.config.provider),
        )
        approved_output_limit = min(
            envelope.budget.max_output_tokens_per_request,
            network_budget.max_output_tokens_per_request,
        )
        if session.config.max_output_tokens > approved_output_limit:
            raise ContractError("provider max_tokens exceeds approved budget")
        self.event_store.append(
            turn_id=envelope.turn_id,
            kind=EventKind.SESSION_STARTED.value,
            payload={
                "phase": envelope.phase.value,
                "task_id": envelope.task_id,
                "tool_schema_sha256": envelope.tool_schema_sha256,
            },
            idempotency_key="session-started",
        )
        self.event_store.append(
            turn_id=envelope.turn_id,
            kind=EventKind.TURN_STARTED.value,
            payload={
                "request_sha256": envelope.request_sha256,
                "phase": envelope.phase.value,
            },
            idempotency_key="turn-started:" + envelope.turn_id,
        )
        self.event_store.append(
            turn_id=envelope.turn_id,
            kind=EventKind.EXPOSURE_PLANNED.value,
            payload={
                "tools": tuple(
                    item["function"]["name"]
                    for item in self.host.surface.tool_definitions
                ),
                "tool_schema_sha256": self.host.surface.tool_schema_sha256,
            },
            idempotency_key="tool-exposure:" + envelope.turn_id,
        )
        self.host.record_seeded_evidence(envelope.turn_id)
        start = self.clock()
        attempts: list[ApiAttemptReceiptV1] = []
        provider_receipts: list[ProviderTurnReceiptV1] = []
        successful_tools = 0
        failed_tools = 0
        final_text = ""
        terminal_state = "failed"
        terminal_reason = "tool loop failed"
        completion_required: tuple[str, ...] = ()
        transport_ordinal = 0
        while True:
            elapsed = self.clock() - start
            remaining_wall_time = (
                network_budget.task_wall_time_seconds - elapsed
            )
            if remaining_wall_time <= 0:
                terminal_state = "failed"
                final_text = "provider wall-time budget exhausted"
                terminal_reason = final_text
                break
            timeout_setter = getattr(
                session, "set_turn_timeout_seconds", None
            )
            if timeout_setter is not None:
                timeout_setter(remaining_wall_time)
            transport_ordinal += 1
            request = session.request_payload(
                tools=list(self.host.surface.tool_definitions)
            )
            # Provider-private reasoning continuation must not influence any
            # persisted evidence, including request digests.
            request_sha256 = canonical_sha256(public_payload(request))
            attempt_start = self.clock()
            try:
                response, provider_receipt = session.turn(
                    tools=list(self.host.surface.tool_definitions)
                )
            except (DeepSeekTransportError, DeepSeekProtocolError) as exc:
                attempt = self._failed_attempt(
                    envelope=envelope,
                    hypothesis=hypothesis,
                    network_budget=network_budget,
                    ordinal=transport_ordinal,
                    request_sha256=request_sha256,
                    started=attempt_start,
                    error=exc,
                    provider=session.config.provider,
                )
                attempts.append(attempt)
                self._emit_attempt(envelope.turn_id, attempt)
                terminal_state = "failed"
                final_text = str(exc)
                terminal_reason = "provider transport or protocol failed"
                break
            provider_receipts.append(provider_receipt)
            attempt = build_api_attempt_receipt(
                attempt_id=f"{envelope.turn_id}.provider.{transport_ordinal}",
                provider=session.config.provider,
                endpoint_origin=network_budget.endpoint_origin,
                hypothesis_sha256=hypothesis.hypothesis_sha256,
                budget_sha256=network_budget.budget_sha256,
                request_sha256=request_sha256,
                response_sha256=canonical_sha256(response),
                status="succeeded",
                latency_ms=max(
                    0, int((self.clock() - attempt_start) * 1000)
                ),
                input_tokens=provider_receipt.input_tokens,
                output_tokens=provider_receipt.output_tokens,
                reasoning_tokens=provider_receipt.reasoning_tokens,
                reported_cost_usd="",
                retry_ordinal=0,
                nonsecret_error_class="",
            )
            attempts.append(attempt)
            self._emit_attempt(envelope.turn_id, attempt)
            self.event_store.append(
                turn_id=envelope.turn_id,
                kind=EventKind.PROVIDER_TURN_OBSERVED.value,
                payload={
                    "receipt_sha256": provider_receipt.receipt_sha256,
                    "provider": provider_receipt.provider,
                    "requested_model": provider_receipt.requested_model,
                    "observed_model": provider_receipt.observed_model,
                    "finish_reason": provider_receipt.finish_reason,
                    "tool_calls_present": (
                        provider_receipt.tool_calls_present
                    ),
                    "reasoning_continuation_present": (
                        provider_receipt.reasoning_continuation_present
                    ),
                },
                idempotency_key=(
                    "provider-turn:" + provider_receipt.receipt_sha256
                ),
            )
            budget_findings = []
            if provider_receipt.input_tokens > min(
                envelope.budget.max_input_tokens_per_request,
                network_budget.max_input_tokens_per_request,
            ):
                budget_findings.append("budget.provider_input_tokens_exceeded")
            if provider_receipt.output_tokens > approved_output_limit:
                budget_findings.append("budget.provider_output_tokens_exceeded")
            if budget_findings:
                self.event_store.append(
                    turn_id=envelope.turn_id,
                    kind=EventKind.TURN_BLOCKED.value,
                    payload={
                        "reason": "provider token budget exceeded",
                        "rule_ids": tuple(budget_findings),
                        "provider_receipt_sha256": (
                            provider_receipt.receipt_sha256
                        ),
                    },
                    idempotency_key=(
                        "token-budget:" + provider_receipt.receipt_sha256
                    ),
                )
                final_text = "provider token budget exceeded"
                terminal_state = "failed"
                terminal_reason = final_text
                break
            assistant = _public_assistant_message(response)
            tool_calls = assistant.get("tool_calls") or []
            if not tool_calls:
                final_text = str(assistant.get("content") or "")
                analysis_policy = self.host.analysis_completion_policy
                if analysis_policy is not None:
                    try:
                        completion_required = (
                            self.host.completion_receipts_for_analysis(
                                turn_id=envelope.turn_id
                            )
                        )
                        final_text = self.host.render_completed_analysis_report(
                            completion_required[0]
                        )
                        terminal_state = "complete"
                        terminal_reason = (
                            "task-owned analysis completion policy passed"
                        )
                    except AnalysisIncompleteError:
                        try:
                            completion_required = (
                                self.host.latest_workflow_draft_receipt(),
                            )
                            terminal_state = "planned"
                            terminal_reason = (
                                "workflow draft recorded; required analysis "
                                "stages remain incomplete"
                            )
                        except ContractError:
                            terminal_state = "blocked"
                            terminal_reason = (
                                "model stopped before the required analysis "
                                "completion policy passed"
                            )
                    except ContractError as exc:
                        final_text = (
                            "analysis completion evaluator failed: " + str(exc)
                        )
                        terminal_state = "failed"
                        terminal_reason = final_text
                else:
                    try:
                        completion_required = (
                            self.host.completion_receipts_for_latest_preflight()
                        )
                        terminal_state = "complete"
                        terminal_reason = "host readiness gates passed"
                    except ContractError:
                        try:
                            completion_required = (
                                self.host.latest_workflow_draft_receipt(),
                            )
                            terminal_state = "planned"
                            terminal_reason = (
                                "workflow draft recorded; action-grade gates "
                                "remain"
                            )
                        except ContractError:
                            terminal_state = "blocked"
                            terminal_reason = (
                                "model stopped before a workflow or readiness gate"
                            )
                break
            if successful_tools + failed_tools + len(tool_calls) > (
                envelope.budget.max_tool_calls
            ):
                final_text = "tool-call budget exhausted"
                terminal_state = "failed"
                terminal_reason = final_text
                break
            tool_results = []
            for tool_call in tool_calls:
                call_id, tool_name, arguments = _decode_tool_call(tool_call)
                self.event_store.append(
                    turn_id=envelope.turn_id,
                    kind=EventKind.TOOL_STARTED.value,
                    payload={
                        "request_id": call_id,
                        "tool": tool_name,
                        "arguments_sha256": canonical_sha256(arguments),
                    },
                    idempotency_key="tool-started:" + call_id,
                )
                try:
                    result = self.host.dispatch(
                        turn_id=envelope.turn_id,
                        tool_name=tool_name,
                        arguments=arguments,
                    )
                    successful_tools += 1
                    tool_event_kind = EventKind.TOOL_SUCCEEDED.value
                    tool_event_payload = {
                        "request_id": call_id,
                        "tool": tool_name,
                        "result_sha256": canonical_sha256(result),
                        "verdict": "host_observed",
                        "typed_receipt_status": "present",
                    }
                    tool_event_key = "tool-succeeded:" + call_id
                except Exception as exc:
                    failed_tools += 1
                    public_message = (
                        str(exc)
                        if isinstance(
                            exc,
                            (
                                ContractError,
                                ValueError,
                                TypeError,
                                KeyError,
                                AttributeError,
                            ),
                        )
                        else "tool host operation failed"
                    )
                    result = {
                        "schema_version": "chemsmart.tool-result.v1",
                        "tool": tool_name,
                        "status": "rejected",
                        "error_class": type(exc).__name__,
                        "message": public_message,
                    }
                    tool_event_kind = EventKind.TOOL_FAILED.value
                    tool_event_payload = {
                        "request_id": call_id,
                        "tool": tool_name,
                        "rule_ids": ("tool.dispatch.rejected",),
                        "error_class": type(exc).__name__,
                    }
                    tool_event_key = "tool-failed:" + call_id
                projected = project_tool_feedback(
                    tool=tool_name,
                    result=result,
                    mode=self.feedback_projection,
                )
                # The complete public result remains evidence even when the
                # provider receives only the compact causal projection.
                tool_event_payload.update(
                    {
                        "canonical_result": canonical_data(result),
                        "feedback_projection": self.feedback_projection,
                        "feedback_equivalence_receipt": canonical_data(
                            projected.receipt
                        ),
                    }
                )
                self.event_store.append(
                    turn_id=envelope.turn_id,
                    kind=tool_event_kind,
                    payload=tool_event_payload,
                    idempotency_key=tool_event_key,
                )
                tool_results.append(
                    {
                        "role": "tool",
                        "tool_call_id": call_id,
                        "content": canonical_json(projected.content),
                    }
                )
            session.append_tool_results(tool_results)
        public_transcript = tuple(session.public_history())
        if _contains_private_reasoning(public_transcript):
            raise ContractError("provider sanitizer left private reasoning")
        transcript_artifact = self.event_store.persist_public_transcript(
            turn_id=envelope.turn_id, transcript=public_transcript
        )
        self.event_store.append(
            turn_id=envelope.turn_id,
            kind=EventKind.ARTIFACT_RECORDED.value,
            payload={
                "artifact_id": transcript_artifact["artifact_id"],
                "kind": "public_transcript",
                "artifact_sha256": transcript_artifact["artifact_sha256"],
                "transcript_sha256": transcript_artifact["transcript_sha256"],
            },
            idempotency_key=(
                "public-transcript:" + transcript_artifact["transcript_sha256"]
            ),
        )
        if not self.event_store.state().terminal_state:
            self.event_store.terminate(
                turn_id=envelope.turn_id,
                terminal_state=terminal_state,
                reason=terminal_reason,
                required_receipt_sha256s=completion_required,
            )
        state = self.event_store.state()
        body = {
            "schema_version": "chemsmart.tool-loop-result.v1",
            "session_id": envelope.session_id,
            "turn_id": envelope.turn_id,
            "terminal_state": state.terminal_state,
            "final_text": final_text,
            "public_transcript": public_transcript,
            "public_transcript_sha256": canonical_sha256(public_transcript),
            "public_transcript_artifact_id": transcript_artifact["artifact_id"],
            "public_transcript_artifact_sha256": transcript_artifact[
                "artifact_sha256"
            ],
            "provider_receipt_sha256s": tuple(
                item.receipt_sha256 for item in provider_receipts
            ),
            "api_attempt_receipt_sha256s": tuple(
                item.receipt_sha256 for item in attempts
            ),
            "successful_tool_calls": successful_tools,
            "failed_tool_calls": failed_tools,
            "event_stream_head_sha256": state.latest_event_hash,
        }
        return ToolLoopResultV1(
            **body, result_sha256=canonical_sha256(body)
        )

    def _validate_run_contract(
        self,
        *,
        envelope: TaskEnvelopeV1,
        hypothesis: AdaptiveHypothesisV1,
        network_budget: AdaptiveNetworkBudgetV1,
        provider: str,
    ) -> None:
        if envelope.session_id != self.event_store.session_id:
            raise ContractError("task envelope belongs to another event stream")
        if envelope.tool_schema_sha256 != self.host.surface.tool_schema_sha256:
            raise ContractError("task envelope uses another tool schema")
        if hypothesis.tool_schema_sha256 != self.host.surface.tool_schema_sha256:
            raise ContractError("hypothesis uses another tool schema")
        if provider not in {"deepseek", "alibaba-token-plan"}:
            raise ContractError("active tool loop provider is not registered")
        if network_budget.allowed_provider != provider:
            raise ContractError("network budget belongs to another provider")
        if network_budget.engine_calls or network_budget.hpc_calls:
            raise ContractError("active planning loop cannot run chemistry")
        execution_profile = (
            self.host.surface.profile == "command_compiled_approved_execution"
        )
        if envelope.budget.hpc_calls:
            raise ContractError("active local tool loop cannot authorize HPC")
        if envelope.budget.chemistry_engine_calls and not execution_profile:
            raise ContractError(
                "chemistry calls require the approved-execution tool profile"
            )
        if (
            envelope.budget.max_output_tokens_per_request
            > network_budget.max_output_tokens_per_request
        ):
            raise ContractError("task output budget exceeds network budget")

    def _emit_attempt(self, turn_id: str, attempt: ApiAttemptReceiptV1) -> None:
        self.event_store.append(
            turn_id=turn_id,
            kind=EventKind.API_ATTEMPT_OBSERVED.value,
            payload={
                "receipt_sha256": attempt.receipt_sha256,
                "provider": attempt.provider,
                "endpoint_origin": attempt.endpoint_origin,
                "status": attempt.status,
                "hypothesis_sha256": attempt.hypothesis_sha256,
                "budget_sha256": attempt.budget_sha256,
                "latency_ms": attempt.latency_ms,
                "input_tokens": attempt.input_tokens,
                "output_tokens": attempt.output_tokens,
                "reasoning_tokens": attempt.reasoning_tokens,
                "nonsecret_error_class": attempt.nonsecret_error_class,
            },
            idempotency_key="api-attempt:" + attempt.attempt_id,
        )

    def _failed_attempt(
        self,
        *,
        envelope: TaskEnvelopeV1,
        hypothesis: AdaptiveHypothesisV1,
        network_budget: AdaptiveNetworkBudgetV1,
        ordinal: int,
        request_sha256: str,
        started: float,
        error: Exception,
        provider: str,
    ) -> ApiAttemptReceiptV1:
        error_class = getattr(error, "error_class", type(error).__name__)
        status = {
            "quota_exhausted": "quota_exhausted",
            "credential_invalid": "credential_invalid",
            "rate_limited": "rate_limited",
            "timeout": "timeout",
        }.get(
            str(error_class),
            "protocol_failed"
            if isinstance(error, DeepSeekProtocolError)
            else "transport_failed",
        )
        return build_api_attempt_receipt(
            attempt_id=f"{envelope.turn_id}.provider.{ordinal}",
            provider=provider,
            endpoint_origin=network_budget.endpoint_origin,
            hypothesis_sha256=hypothesis.hypothesis_sha256,
            budget_sha256=network_budget.budget_sha256,
            request_sha256=request_sha256,
            response_sha256="",
            status=status,
            latency_ms=max(0, int((self.clock() - started) * 1000)),
            retry_ordinal=0,
            nonsecret_error_class=str(error_class),
        )


def _public_assistant_message(response: Mapping[str, Any]) -> dict[str, Any]:
    choices = response.get("choices")
    if not isinstance(choices, list) or not choices:
        raise DeepSeekProtocolError("public response has no choices")
    choice = choices[0]
    message = choice.get("message") if isinstance(choice, Mapping) else None
    if not isinstance(message, Mapping):
        raise DeepSeekProtocolError("public response has no assistant message")
    return dict(message)


def _decode_tool_call(
    value: Mapping[str, Any],
) -> tuple[str, str, dict[str, Any]]:
    call_id = str(value.get("id") or "")
    function = value.get("function")
    if not call_id or not isinstance(function, Mapping):
        raise DeepSeekProtocolError("malformed tool call")
    name = str(function.get("name") or "")
    try:
        arguments = json.loads(str(function.get("arguments") or ""))
    except json.JSONDecodeError as exc:
        raise DeepSeekProtocolError("malformed tool arguments") from exc
    if not name or not isinstance(arguments, dict):
        raise DeepSeekProtocolError("tool call lacks name or object arguments")
    return call_id, name, arguments


def _contains_private_reasoning(value: Any) -> bool:
    forbidden = {"reasoning_content", "thinking", "analysis", "<think>"}
    if isinstance(value, Mapping):
        return any(
            str(key).lower() in forbidden
            or _contains_private_reasoning(item)
            for key, item in value.items()
        )
    if isinstance(value, (tuple, list)):
        return any(_contains_private_reasoning(item) for item in value)
    if isinstance(value, str):
        return "<think" in value.lower()
    return False


__all__ = ["ToolLoopResultV1", "ToolLoopRunner"]
