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
from chemsmart.agent.request_context import (
    ProviderAttemptReceiptV1,
    ProviderNetworkBudgetV1,
    RequestContextProvenanceV1,
    build_provider_attempt_receipt,
)
from chemsmart.agent.analysis_completion import AnalysisIncompleteError
from chemsmart.agent.runtime.contracts import TaskEnvelopeV1
from chemsmart.agent.runtime.deepseek import (
    DeepSeekProtocolError,
    DeepSeekTransportError,
    ProviderTurnReceiptV1,
    public_provider_request,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.runtime.events import EventKind
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1


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


#: Characters per token for the conservative pre-request size estimate.  Real
#: tokenizers vary by model and are not available to the host, so this is
#: deliberately an estimate used only to refuse a request that would clearly be
#: rejected -- never to decide anything scientific.
_CHARS_PER_TOKEN = 3.5

# Keep enough of the task envelope after any provider turn to close the stream,
# persist failure/success evidence, and render the Runtime result.  The outer
# episode timeout remains a last resort rather than the normal transport stop.
_PROVIDER_POST_TURN_RESERVE_SECONDS = 30.0


def estimate_request_input_tokens(request: Mapping[str, Any]) -> int:
    """Estimate the input size of a provider request before sending it.

    Underestimating lets a doomed request through, which is the status quo.
    Overestimating refuses a request that would have worked.  The ratio is
    therefore chosen on the low side of typical byte-per-token figures so the
    guard fires only when the conversation is clearly past the limit.
    """

    total = 0
    for message in request.get("messages") or ():
        if not isinstance(message, Mapping):
            continue
        for value in message.values():
            total += len(value) if isinstance(value, str) else len(str(value))
    for tool in request.get("tools") or ():
        total += len(str(tool))
    return int(total / _CHARS_PER_TOKEN)


class ToolLoopRunner:
    """Drive provider-native turns through the only approved dispatcher."""

    def __init__(
        self,
        *,
        host: CommandCompiledToolHostV1,
        event_store: RuntimeEventStore,
        clock: Callable[[], float] = time.monotonic,
    ) -> None:
        self.host = host
        self.event_store = event_store
        self.clock = clock

    def run(
        self,
        *,
        session: Any,
        envelope: TaskEnvelopeV1,
        request_context: RequestContextProvenanceV1,
        provider_budget: ProviderNetworkBudgetV1,
    ) -> ToolLoopResultV1:
        self._validate_run_contract(
            envelope=envelope,
            request_context=request_context,
            provider_budget=provider_budget,
            provider=str(session.config.provider),
        )
        approved_output_limit = min(
            envelope.budget.max_output_tokens_per_request,
            provider_budget.max_output_tokens_per_request,
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
        attempts: list[ProviderAttemptReceiptV1] = []
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
                provider_budget.task_wall_time_seconds - elapsed
            )
            if remaining_wall_time <= 0:
                terminal_state = "failed"
                final_text = "provider wall-time budget exhausted"
                terminal_reason = final_text
                break
            post_turn_reserve = min(
                _PROVIDER_POST_TURN_RESERVE_SECONDS,
                max(0.1, provider_budget.task_wall_time_seconds * 0.1),
            )
            provider_turn_allowance = remaining_wall_time - post_turn_reserve
            if provider_turn_allowance <= 0:
                terminal_state = "failed"
                final_text = "provider wall-time reserve reached"
                terminal_reason = final_text
                break
            timeout_setter = getattr(
                session, "set_turn_timeout_seconds", None
            )
            if timeout_setter is not None:
                timeout_setter(provider_turn_allowance)
            transport_ordinal += 1
            request = session.request_payload(
                tools=list(self.host.surface.tool_definitions)
            )
            # Provider-private reasoning continuation must not influence any
            # persisted evidence, including request digests.
            request_sha256 = canonical_sha256(
                public_provider_request(request)
            )
            context_limit = min(
                envelope.budget.max_input_tokens_per_request,
                provider_budget.max_input_tokens_per_request,
            )
            projected = estimate_request_input_tokens(request)
            if projected > context_limit:
                # The existing budget check runs on the provider receipt, which
                # means the request has already been sent and charged before
                # anyone notices it was too large.  The session history grows
                # on every turn and is never pruned, so a long workflow reaches
                # this point by construction rather than by accident: measured
                # across fourteen live sessions, input tokens rose monotonically
                # on every single turn and total cost scaled with the square of
                # the turn count.  Refusing here spends nothing and says why.
                self.event_store.append(
                    turn_id=envelope.turn_id,
                    kind=EventKind.TURN_BLOCKED.value,
                    payload={
                        "reason": (
                            "the conversation would exceed the provider "
                            f"context budget: about {projected} input tokens "
                            f"against a limit of {context_limit}. The session "
                            "history grows every turn and is never pruned, so "
                            "this is reached by accumulation, not by one large "
                            "message."
                        ),
                        "rule_ids": ("budget.context_would_be_exceeded",),
                        "projected_input_tokens": projected,
                        "context_limit": context_limit,
                    },
                    idempotency_key="context-budget:" + request_sha256,
                )
                terminal_state = "blocked"
                final_text = "conversation exceeds the provider context budget"
                terminal_reason = "context budget would be exceeded"
                break
            attempt_start = self.clock()
            try:
                response, provider_receipt = session.turn(
                    tools=list(self.host.surface.tool_definitions)
                )
            except (DeepSeekTransportError, DeepSeekProtocolError) as exc:
                attempt = self._failed_attempt(
                    envelope=envelope,
                    request_context=request_context,
                    provider_budget=provider_budget,
                    ordinal=transport_ordinal,
                    request_sha256=request_sha256,
                    started=attempt_start,
                    error=exc,
                    provider=session.config.provider,
                )
                attempts.append(attempt)
                self._emit_attempt(
                    envelope.turn_id,
                    attempt,
                    protocol_observation=(
                        exc.public_observation()
                        if isinstance(exc, DeepSeekProtocolError)
                        else None
                    ),
                    transport_observation=(
                        _public_transport_failure(session, exc)
                        if isinstance(exc, DeepSeekTransportError)
                        else None
                    ),
                )
                terminal_state = "failed"
                final_text = str(exc)
                terminal_reason = "provider transport or protocol failed"
                break
            provider_receipts.append(provider_receipt)
            attempt = build_provider_attempt_receipt(
                attempt_id=f"{envelope.turn_id}.provider.{transport_ordinal}",
                provider=session.config.provider,
                endpoint_origin=provider_budget.endpoint_origin,
                request_context_sha256=request_context.provenance_sha256,
                provider_budget_sha256=provider_budget.budget_sha256,
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
                    "requested_provider": session.config.provider,
                    "observed_provider": provider_receipt.provider,
                    "requested_model": provider_receipt.requested_model,
                    "observed_model": provider_receipt.observed_model,
                    "requested_reasoning_effort": (
                        str(session.config.reasoning_effort)
                    ),
                    # The supported endpoints echo model identity but do not
                    # report the applied reasoning effort independently.
                    "observed_reasoning_effort": "not_reported",
                    "transport_deadlines": (
                        session.public_transport_deadline_record()
                    ),
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
                provider_budget.max_input_tokens_per_request,
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
            try:
                assistant = _public_assistant_message(response)
                raw_tool_calls = assistant.get("tool_calls") or []
                if not isinstance(raw_tool_calls, list):
                    raise DeepSeekProtocolError(
                        "public assistant tool calls must be a list"
                    )
                decoded_tool_calls = tuple(
                    _decode_tool_call(tool_call)
                    for tool_call in raw_tool_calls
                )
            except DeepSeekProtocolError:
                # Decode the complete public envelope before dispatching any
                # member.  This is a host-side invariant behind the provider
                # session's raw-envelope validation, not a partial retry path.
                self.event_store.append(
                    turn_id=envelope.turn_id,
                    kind=EventKind.TURN_BLOCKED.value,
                    payload={
                        "reason": "public provider tool envelope failed atomic decoding",
                        "rule_ids": (
                            "provider.public_tool_envelope_invalid",
                        ),
                    },
                    idempotency_key=(
                        "public-tool-envelope-invalid:"
                        + canonical_sha256(response)
                    ),
                )
                final_text = "provider tool envelope failed validation"
                terminal_state = "failed"
                terminal_reason = final_text
                break
            if not decoded_tool_calls:
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
            if successful_tools + failed_tools + len(decoded_tool_calls) > (
                envelope.budget.max_tool_calls
            ):
                final_text = "tool-call budget exhausted"
                terminal_state = "failed"
                terminal_reason = final_text
                break
            tool_results = []
            for call_id, tool_name, arguments in decoded_tool_calls:
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
                wait_started = None
                wait_emitted = False
                try:
                    if tool_name == "execute_approved_program_node":
                        wait_timeout = (
                            self.host.execution_wait_timeout_seconds()
                        )
                        wait_started = self.clock()
                        self.event_store.append(
                            turn_id=envelope.turn_id,
                            kind=EventKind.TOOL_WAITING.value,
                            payload={
                                "request_id": call_id,
                                "tool": tool_name,
                                "wait_kind": "approved_program_tool_dispatch",
                                "timeout_seconds": wait_timeout,
                                "provider_calls_while_waiting": 0,
                                "continuation_state": "private_in_memory",
                            },
                            idempotency_key="tool-waiting:" + call_id,
                        )
                        wait_emitted = True
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
                if wait_emitted and wait_started is not None:
                    process_observation = _execution_process_observation(result)
                    wake_reason = _execution_wake_reason(
                        result=result,
                        process_observation=process_observation,
                    )
                    self.event_store.append(
                        turn_id=envelope.turn_id,
                        kind=EventKind.TOOL_WOKE.value,
                        payload={
                            "request_id": call_id,
                            "tool": tool_name,
                            "wake_reason": wake_reason,
                            "waited_seconds": max(
                                0.0, self.clock() - wait_started
                            ),
                            "provider_calls_while_waiting": 0,
                            "process_observation_sha256": str(
                                process_observation.get("receipt_sha256") or ""
                            ),
                        },
                        idempotency_key="tool-woke:" + call_id,
                    )
                tool_event_payload["canonical_result"] = canonical_data(result)
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
                        "content": canonical_json(result),
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
        request_context: RequestContextProvenanceV1,
        provider_budget: ProviderNetworkBudgetV1,
        provider: str,
    ) -> None:
        if envelope.session_id != self.event_store.session_id:
            raise ContractError("task envelope belongs to another event stream")
        if envelope.tool_schema_sha256 != self.host.surface.tool_schema_sha256:
            raise ContractError("task envelope uses another tool schema")
        if (
            request_context.tool_schema_sha256
            != self.host.surface.tool_schema_sha256
        ):
            raise ContractError("request context uses another tool schema")
        if request_context.prompt_sha256 != envelope.request_sha256:
            raise ContractError("request context uses another prompt")
        if (
            request_context.task_spec_sha256
            not in self.host.task_spec_sha256s
        ):
            raise ContractError("request context uses another task")
        if (
            request_context.provider_budget_sha256
            != provider_budget.budget_sha256
        ):
            raise ContractError("request context uses another provider budget")
        if provider not in {"deepseek", "alibaba-token-plan"}:
            raise ContractError("active tool loop provider is not registered")
        if provider_budget.allowed_provider != provider:
            raise ContractError("provider budget belongs to another provider")
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
            > provider_budget.max_output_tokens_per_request
        ):
            raise ContractError("task output budget exceeds provider budget")

    def _emit_attempt(
        self,
        turn_id: str,
        attempt: ProviderAttemptReceiptV1,
        *,
        protocol_observation: Mapping[str, Any] | None = None,
        transport_observation: Mapping[str, Any] | None = None,
    ) -> None:
        payload = {
            "receipt_sha256": attempt.receipt_sha256,
            "provider": attempt.provider,
            "endpoint_origin": attempt.endpoint_origin,
            "status": attempt.status,
            "request_context_sha256": attempt.request_context_sha256,
            "provider_budget_sha256": attempt.provider_budget_sha256,
            "latency_ms": attempt.latency_ms,
            "input_tokens": attempt.input_tokens,
            "output_tokens": attempt.output_tokens,
            "reasoning_tokens": attempt.reasoning_tokens,
            "nonsecret_error_class": attempt.nonsecret_error_class,
        }
        if protocol_observation is not None:
            payload["response_sha256"] = attempt.response_sha256
            payload["protocol_failure"] = canonical_data(
                protocol_observation
            )
        if transport_observation is not None:
            payload["transport_failure"] = canonical_data(
                transport_observation
            )
        self.event_store.append(
            turn_id=turn_id,
            kind=EventKind.API_ATTEMPT_OBSERVED.value,
            payload=payload,
            idempotency_key="api-attempt:" + attempt.attempt_id,
        )

    def _failed_attempt(
        self,
        *,
        envelope: TaskEnvelopeV1,
        request_context: RequestContextProvenanceV1,
        provider_budget: ProviderNetworkBudgetV1,
        ordinal: int,
        request_sha256: str,
        started: float,
        error: Exception,
        provider: str,
    ) -> ProviderAttemptReceiptV1:
        error_class = getattr(error, "error_class", type(error).__name__)
        status = {
            "quota_exhausted": "quota_exhausted",
            "credential_invalid": "credential_invalid",
            "rate_limited": "rate_limited",
            "timeout": "timeout",
            "connect_timeout": "timeout",
            "first_event_timeout": "timeout",
            "inter_event_timeout": "timeout",
            "turn_deadline_exceeded": "timeout",
        }.get(
            str(error_class),
            "protocol_failed"
            if isinstance(error, DeepSeekProtocolError)
            else "transport_failed",
        )
        response_sha256 = (
            error.response_envelope_sha256
            if isinstance(error, DeepSeekProtocolError)
            else ""
        )
        return build_provider_attempt_receipt(
            attempt_id=f"{envelope.turn_id}.provider.{ordinal}",
            provider=provider,
            endpoint_origin=provider_budget.endpoint_origin,
            request_context_sha256=request_context.provenance_sha256,
            provider_budget_sha256=provider_budget.budget_sha256,
            request_sha256=request_sha256,
            response_sha256=response_sha256,
            status=status,
            latency_ms=max(0, int((self.clock() - started) * 1000)),
            retry_ordinal=0,
            nonsecret_error_class=str(error_class),
        )


def _execution_process_observation(result: Any) -> dict[str, Any]:
    """Locate the host-owned public process receipt in a tool result."""

    if not isinstance(result, Mapping):
        return {}
    outer = result.get("result")
    outer = outer if isinstance(outer, Mapping) else result
    process = outer.get("process_observation")
    if isinstance(process, Mapping):
        return dict(process)
    execution = outer.get("execution")
    if isinstance(execution, Mapping):
        observations = execution.get("observations")
        if isinstance(observations, Mapping):
            process = observations.get("process_observation")
            return dict(process) if isinstance(process, Mapping) else {}
    return {}


def _public_transport_failure(
    session: Any, error: DeepSeekTransportError
) -> dict[str, Any]:
    """Bind a sanitized failure to requested identity and deadline policy."""

    return {
        **error.public_observation(),
        "requested_provider": str(session.config.provider),
        "observed_provider": "not_observed",
        "requested_model": str(session.config.model),
        "observed_model": "not_observed",
        "requested_reasoning_effort": str(session.config.reasoning_effort),
        "observed_reasoning_effort": "not_reported",
        "transport_deadlines": session.public_transport_deadline_record(),
    }


def _execution_wake_reason(
    *, result: Any, process_observation: Mapping[str, Any]
) -> str:
    """Classify why the synchronous approved-engine wait returned."""

    outer = result.get("result") if isinstance(result, Mapping) else None
    outer = outer if isinstance(outer, Mapping) else result
    if isinstance(outer, Mapping) and outer.get("idempotent_replay") is True:
        return "replay"
    if process_observation.get("timed_out") is True:
        return "timeout"
    if isinstance(result, Mapping) and result.get("status") == "rejected":
        return "failure"
    state = str(process_observation.get("state") or "")
    returncode = process_observation.get("returncode")
    if process_observation.get("memory_limit_exceeded") is True:
        return "failure"
    if (
        isinstance(returncode, int)
        and not isinstance(returncode, bool)
        and returncode < 0
    ):
        return "signal"
    if process_observation and (
        state != "exited" or returncode not in {0, None}
    ):
        return "failure"
    return "result"


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
    if not isinstance(value, Mapping):
        raise DeepSeekProtocolError("malformed tool call")
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
    if isinstance(value, (tuple, list)):
        return any(
            _contains_private_reasoning_at_path(
                item,
                (),
                (
                    "assistant_message"
                    if isinstance(item, Mapping)
                    and str(item.get("role") or "") == "assistant"
                    else "generic"
                ),
            )
            for item in value
        )
    return _contains_private_reasoning_at_path(value, (), "generic")


def _contains_private_reasoning_at_path(
    value: Any, path: tuple[Any, ...], context: str
) -> bool:
    forbidden = {"reasoning_content", "thinking", "analysis", "<think>"}
    if isinstance(value, Mapping):
        return any(
            str(key).lower() in forbidden
            or (
                not (
                    _is_public_tool_function_path(path, context)
                    and key == "arguments"
                )
                and _contains_private_reasoning_at_path(
                    item, (*path, key), context
                )
            )
            for key, item in value.items()
        )
    if isinstance(value, (tuple, list)):
        return any(
            _contains_private_reasoning_at_path(
                item, (*path, index), context
            )
            for index, item in enumerate(value)
        )
    if isinstance(value, str):
        return "<think" in value.lower()
    return False


def _is_public_tool_function_path(
    path: tuple[Any, ...], context: str
) -> bool:
    return (
        context == "assistant_message"
        and len(path) == 3
        and path[0] == "tool_calls"
        and isinstance(path[1], int)
        and path[2] == "function"
    )


__all__ = ["ToolLoopResultV1", "ToolLoopRunner"]
