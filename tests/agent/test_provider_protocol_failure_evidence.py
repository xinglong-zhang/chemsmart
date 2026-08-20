from __future__ import annotations

import json
from types import SimpleNamespace

import pytest

from chemsmart.agent._contracts import (
    ContractError,
    canonical_json,
    canonical_sha256,
)
from chemsmart.agent.loop import ToolLoopRunner, _execution_wake_reason
from chemsmart.agent.request_context import (
    ProviderNetworkBudgetV1,
    build_request_context_provenance,
)
from chemsmart.agent.runtime.alibaba import (
    AlibabaTokenPlanHttpsTransport,
    Qwen38MaxConfigV1,
    Qwen38MaxToolSession,
)
from chemsmart.agent.runtime.contracts import (
    ResourceBudgetV1,
    TaskEnvelopeV1,
    TaskPhase,
)
from chemsmart.agent.runtime.deepseek import (
    DeepSeekTransportError,
    public_provider_response,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.runtime.events import EventKind


class _DispatchSpyHost:
    def __init__(self) -> None:
        tool_definitions = (
            {
                "type": "function",
                "function": {
                    "name": "inspect_program_capability",
                    "description": "Inspect one program capability.",
                    "parameters": {
                        "type": "object",
                        "properties": {"program": {"type": "string"}},
                        "required": ["program"],
                    },
                },
            },
        )
        self.surface = SimpleNamespace(
            tool_definitions=tool_definitions,
            tool_schema_sha256=canonical_sha256(tool_definitions),
            profile="command_compiled_preview",
        )
        self.analysis_completion_policy = None
        self.task_spec_sha256s = frozenset(
            {canonical_sha256("synthetic task")}
        )
        self.dispatched: list[tuple[str, dict]] = []

    def record_seeded_evidence(self, _turn_id: str) -> None:
        return None

    def dispatch(self, *, tool_name: str, arguments: dict, **_kwargs):
        self.dispatched.append((tool_name, arguments))
        return {"status": "supported"}

    def completion_receipts_for_latest_preflight(self):
        raise ContractError("no workflow preflight in synthetic test")

    def latest_workflow_draft_receipt(self):
        raise ContractError("no workflow draft in synthetic test")

    def unapproved_workflow_summary(self):
        """No plan exists here, so there is nothing to say about approval."""

        return None


def _run_contracts(
    host: _DispatchSpyHost,
    config: Qwen38MaxConfigV1,
    *,
    chemistry_engine_calls: int = 0,
):
    budget = ResourceBudgetV1(
        max_input_tokens_per_request=config.context_tokens,
        max_output_tokens_per_request=config.max_output_tokens,
        max_tool_calls=8,
        wall_time_seconds=30.0,
        chemistry_engine_calls=chemistry_engine_calls,
    )
    envelope_body = {
        "schema_version": "chemsmart.task-envelope.v1",
        "task_id": "malformed-envelope-regression",
        "session_id": "protocol-session",
        "turn_id": "protocol-session.turn-1",
        "request_sha256": canonical_sha256("synthetic request"),
        "cwd_sha256": canonical_sha256("synthetic workspace"),
        "phase": TaskPhase.ROUTE,
        "budget": budget,
        "tool_schema_sha256": host.surface.tool_schema_sha256,
    }
    envelope = TaskEnvelopeV1(
        **envelope_body, envelope_sha256=canonical_sha256(envelope_body)
    )
    network_body = {
        "schema_version": "chemsmart.provider-network-budget.v1",
        "allowed_provider": config.provider,
        "endpoint_origin": config.endpoint,
        "max_concurrency": 1,
        "max_input_tokens_per_request": config.context_tokens,
        "max_output_tokens_per_request": config.max_output_tokens,
        "task_wall_time_seconds": 30.0,
    }
    network = ProviderNetworkBudgetV1(
        **network_body, budget_sha256=canonical_sha256(network_body)
    )
    request_context = build_request_context_provenance(
        task_spec_sha256=canonical_sha256("synthetic task"),
        prompt_sha256=envelope.request_sha256,
        tool_schema_sha256=host.surface.tool_schema_sha256,
        configuration_sha256=canonical_sha256("synthetic configuration"),
        provider_budget_sha256=network.budget_sha256,
    )
    return envelope, request_context, network


class _ExecutionWaitHost(_DispatchSpyHost):
    def __init__(self, timeline):
        super().__init__()
        tool_definitions = (
            {
                "type": "function",
                "function": {
                    "name": "execute_approved_program_node",
                    "description": "Execute one approved program node.",
                    "parameters": {
                        "type": "object",
                        "properties": {"node_id": {"type": "string"}},
                        "required": ["node_id"],
                    },
                },
            },
        )
        self.surface = SimpleNamespace(
            tool_definitions=tool_definitions,
            tool_schema_sha256=canonical_sha256(tool_definitions),
            profile="command_compiled_approved_execution",
        )
        self.timeline = timeline

    def execution_wait_timeout_seconds(self):
        return 120.0

    def dispatch(self, *, tool_name: str, arguments: dict, **_kwargs):
        self.timeline.append("engine-wait-start")
        self.dispatched.append((tool_name, arguments))
        self.timeline.append("engine-wait-finished")
        return {
            "execution": {"execution_state": "validated"},
            "process_observation": {
                "receipt_sha256": "c" * 64,
                "state": "exited",
                "timed_out": False,
                "termination_requested": False,
            },
        }


def test_request_context_must_bind_the_active_host_task(tmp_path):
    config = Qwen38MaxConfigV1()
    host = _DispatchSpyHost()
    envelope, _request_context, network = _run_contracts(host, config)
    mismatched = build_request_context_provenance(
        task_spec_sha256=canonical_sha256("another task"),
        prompt_sha256=envelope.request_sha256,
        tool_schema_sha256=host.surface.tool_schema_sha256,
        configuration_sha256=canonical_sha256("synthetic configuration"),
        provider_budget_sha256=network.budget_sha256,
    )
    runner = ToolLoopRunner(
        host=host,
        event_store=RuntimeEventStore(
            tmp_path / "events" / "runtime.jsonl",
            session_id="protocol-session",
        ),
    )

    with pytest.raises(ContractError, match="another task"):
        runner._validate_run_contract(
            envelope=envelope,
            request_context=mismatched,
            provider_budget=network,
            provider=config.provider,
        )


def test_alibaba_malformed_envelope_is_auditable_without_partial_dispatch(
    tmp_path,
):
    malformed_arguments = '{"api_key":"sk-secret",'
    response = {
        "id": "synthetic-malformed-response",
        "model": "qwen3.8-max",
        "choices": [
            {
                "finish_reason": "tool_calls",
                "message": {
                    "role": "assistant",
                    "content": "",
                    "reasoning_content": "PRIVATE_REASONING_MUST_NOT_PERSIST",
                    "tool_calls": [
                        {
                            "id": "valid-first-call",
                            "type": "function",
                            "function": {
                                "name": "inspect_program_capability",
                                "arguments": '{"program":"gaussian"}',
                            },
                        },
                        {
                            "id": "malformed-second-call",
                            "type": "function",
                            "function": {
                                "name": "inspect_program_capability",
                                "arguments": malformed_arguments,
                            },
                        },
                    ],
                },
            }
        ],
        "usage": {"prompt_tokens": 10, "completion_tokens": 5},
    }
    config = Qwen38MaxConfigV1()
    session = Qwen38MaxToolSession(
        transport=lambda _payload: response,
        messages=[{"role": "user", "content": "Inspect capabilities."}],
        config=config,
    )
    store = RuntimeEventStore(
        tmp_path / "events" / "runtime.jsonl",
        session_id="protocol-session",
    )
    host = _DispatchSpyHost()
    envelope, request_context, network = _run_contracts(host, config)

    result = ToolLoopRunner(host=host, event_store=store).run(
        session=session,
        envelope=envelope,
        request_context=request_context,
        provider_budget=network,
    )

    assert result.terminal_state == "failed"
    assert result.successful_tool_calls == 0
    assert result.failed_tool_calls == 0
    assert host.dispatched == []
    assert session.receipts == ()
    assert session.outstanding_tool_call_ids == ()

    attempts = [
        event
        for event in store.read_events()
        if event.kind == EventKind.API_ATTEMPT_OBSERVED.value
    ]
    # A persistently malformed provider is re-asked independently a bounded
    # number of times before the session terminates; every attempt is in the
    # trail and only the last one gives up.
    assert len(attempts) == 3
    assert [
        a.payload["protocol_failure"]["retry_decision"] for a in attempts
    ] == ["retried_independently", "retried_independently", "not_retried"]
    event = attempts[0]
    public_response = public_provider_response(response)
    expected_digest = canonical_sha256(public_response)
    failure = event.payload["protocol_failure"]
    assert event.payload["status"] == "protocol_failed"
    assert event.payload["nonsecret_error_class"] == "DeepSeekProtocolError"
    assert event.payload["response_sha256"] == expected_digest
    assert event.payload["latency_ms"] >= 0
    assert failure == {
        "schema_version": "chemsmart.provider-protocol-failure.v1",
        "error_class": "DeepSeekProtocolError",
        "failing_field": (
            "choices[0].message.tool_calls[1].function.arguments"
        ),
        "json_offset": len(malformed_arguments),
        "response_envelope_sha256": expected_digest,
        "response_envelope_bytes": len(
            canonical_json(public_response).encode("utf-8")
        ),
        "retry_decision": "retried_independently",
        "recovery_requirement": "new_independent_attempt",
    }

    persisted = "\n".join(
        path.read_text(encoding="utf-8")
        for path in sorted(tmp_path.rglob("*"))
        if path.is_file()
    )
    assert "PRIVATE_REASONING_MUST_NOT_PERSIST" not in persisted
    assert "sk-secret" not in persisted
    assert malformed_arguments not in persisted
    assert "valid-first-call" not in persisted
    assert json.loads(persisted.splitlines()[-1])["kind"] == (
        EventKind.RUNTIME_TERMINATED.value
    )


def test_valid_tool_argument_text_is_preserved_and_dispatched_atomically(
    tmp_path,
):
    responses = iter(
        (
            {
                "id": "synthetic-public-arguments",
                "model": "qwen3.8-max",
                "choices": [
                    {
                        "finish_reason": "tool_calls",
                        "message": {
                            "role": "assistant",
                            "content": "",
                            "reasoning_content": "private reasoning",
                            "metadata": {
                                "name": "not-a-tool-function",
                                "arguments": (
                                    "<think>PRIVATE_METADATA_REASONING</think>"
                                ),
                                "tool_calls": [
                                    {
                                        "function": {
                                            "arguments": (
                                                "<think>PRIVATE_NESTED_METADATA"
                                                "</think>"
                                            )
                                        }
                                    }
                                ],
                            },
                            "tool_calls": [
                                {
                                    "id": "call-one",
                                    "type": "function",
                                    "function": {
                                        "name": "inspect_program_capability",
                                        "arguments": '{"program":"gaussian"}',
                                    },
                                },
                                {
                                    "id": "call-two",
                                    "type": "function",
                                    "function": {
                                        "name": "inspect_program_capability",
                                        "arguments": '{"program":"<think>literal</think>"}',
                                    },
                                },
                            ],
                        },
                    }
                ],
                "usage": {"prompt_tokens": 10, "completion_tokens": 5},
            },
            {
                "id": "synthetic-final",
                "model": "qwen3.8-max",
                "choices": [
                    {
                        "finish_reason": "stop",
                        "message": {
                            "role": "assistant",
                            "content": "Done.",
                            "reasoning_content": "private final reasoning",
                        },
                    }
                ],
                "usage": {"prompt_tokens": 10, "completion_tokens": 5},
            },
        )
    )
    config = Qwen38MaxConfigV1()
    session = Qwen38MaxToolSession(
        transport=lambda _payload: next(responses),
        messages=[{"role": "user", "content": "Inspect capabilities."}],
        config=config,
    )
    store = RuntimeEventStore(
        tmp_path / "events" / "runtime.jsonl",
        session_id="protocol-session",
    )
    host = _DispatchSpyHost()
    envelope, request_context, network = _run_contracts(host, config)

    result = ToolLoopRunner(host=host, event_store=store).run(
        session=session,
        envelope=envelope,
        request_context=request_context,
        provider_budget=network,
    )

    assert result.terminal_state == "blocked"
    assert host.dispatched == [
        ("inspect_program_capability", {"program": "gaussian"}),
        (
            "inspect_program_capability",
            {"program": "<think>literal</think>"},
        ),
    ]
    assert result.successful_tool_calls == 2
    persisted_transcript = json.dumps(result.public_transcript)
    assert "private reasoning" not in persisted_transcript
    assert "PRIVATE_METADATA_REASONING" not in persisted_transcript
    assert "PRIVATE_NESTED_METADATA" not in persisted_transcript
    assert "<think>literal</think>" in persisted_transcript


def test_public_decode_failure_cannot_dispatch_an_earlier_call(
    tmp_path, monkeypatch
):
    import chemsmart.agent.runtime.deepseek as deepseek

    response = {
        "id": "synthetic-corrupted-public-envelope",
        "model": "qwen3.8-max",
        "choices": [
            {
                "finish_reason": "tool_calls",
                "message": {
                    "role": "assistant",
                    "content": "",
                    "tool_calls": [
                        {
                            "id": "call-one",
                            "type": "function",
                            "function": {
                                "name": "inspect_program_capability",
                                "arguments": '{"program":"gaussian"}',
                            },
                        },
                        {
                            "id": "call-two",
                            "type": "function",
                            "function": {
                                "name": "inspect_program_capability",
                                "arguments": '{"program":"orca"}',
                            },
                        },
                    ],
                },
            }
        ],
        "usage": {"prompt_tokens": 10, "completion_tokens": 5},
    }
    original_public = deepseek.public_provider_response

    def _corrupt_second_call(value):
        public = original_public(value)
        public["choices"][0]["message"]["tool_calls"][1]["function"][
            "arguments"
        ] = "{"
        return public

    monkeypatch.setattr(
        deepseek, "public_provider_response", _corrupt_second_call
    )
    config = Qwen38MaxConfigV1()
    session = Qwen38MaxToolSession(
        transport=lambda _payload: response,
        messages=[{"role": "user", "content": "Inspect capabilities."}],
        config=config,
    )
    store = RuntimeEventStore(
        tmp_path / "events" / "runtime.jsonl",
        session_id="protocol-session",
    )
    host = _DispatchSpyHost()
    envelope, request_context, network = _run_contracts(host, config)

    result = ToolLoopRunner(host=host, event_store=store).run(
        session=session,
        envelope=envelope,
        request_context=request_context,
        provider_budget=network,
    )

    assert result.terminal_state == "failed"
    assert result.successful_tool_calls == 0
    assert result.failed_tool_calls == 0
    assert host.dispatched == []
    blocked = [
        event
        for event in store.read_events()
        if event.kind == EventKind.TURN_BLOCKED.value
    ]
    assert blocked[-1].payload["rule_ids"] == [
        "provider.public_tool_envelope_invalid"
    ]


def test_approved_engine_execution_parks_provider_and_wakes_with_durable_evidence(
    tmp_path,
):
    timeline = []
    responses = iter(
        (
            {
                "model": "qwen3.8-max",
                "choices": [
                    {
                        "finish_reason": "tool_calls",
                        "message": {
                            "role": "assistant",
                            "content": "Executing the approved node.",
                            "reasoning_content": "private continuation",
                            "tool_calls": [
                                {
                                    "id": "execute-one",
                                    "type": "function",
                                    "function": {
                                        "name": "execute_approved_program_node",
                                        "arguments": '{"node_id":"node-1"}',
                                    },
                                }
                            ],
                        },
                    }
                ],
                "usage": {"prompt_tokens": 10, "completion_tokens": 5},
            },
            {
                "model": "qwen3.8-max",
                "choices": [
                    {
                        "finish_reason": "stop",
                        "message": {
                            "role": "assistant",
                            "content": "Execution observed.",
                            "reasoning_content": "private final",
                        },
                    }
                ],
                "usage": {"prompt_tokens": 20, "completion_tokens": 4},
            },
        )
    )

    def transport(_payload):
        timeline.append("provider-turn")
        return next(responses)

    config = Qwen38MaxConfigV1()
    session = Qwen38MaxToolSession(
        transport=transport,
        messages=[{"role": "user", "content": "Run the approved node."}],
        config=config,
    )
    store = RuntimeEventStore(
        tmp_path / "events" / "runtime.jsonl",
        session_id="protocol-session",
    )
    host = _ExecutionWaitHost(timeline)
    envelope, request_context, network = _run_contracts(
        host, config, chemistry_engine_calls=1
    )

    result = ToolLoopRunner(host=host, event_store=store).run(
        session=session,
        envelope=envelope,
        request_context=request_context,
        provider_budget=network,
    )

    assert timeline == [
        "provider-turn",
        "engine-wait-start",
        "engine-wait-finished",
        "provider-turn",
    ]
    assert result.successful_tool_calls == 1
    events = store.read_events()
    waiting = [event for event in events if event.kind == "tool_waiting"]
    woke = [event for event in events if event.kind == "tool_woke"]
    assert waiting[0].payload == {
        "request_id": "execute-one",
        "tool": "execute_approved_program_node",
        "wait_kind": "approved_program_tool_dispatch",
        "timeout_seconds": 120.0,
        "provider_calls_while_waiting": 0,
        "continuation_state": "private_in_memory",
    }
    assert woke[0].payload["wake_reason"] == "result"
    assert woke[0].payload["provider_calls_while_waiting"] == 0
    assert woke[0].payload["process_observation_sha256"] == "c" * 64
    replayed = store.state()
    assert replayed.active_waits == {}
    assert replayed.completed_waits[0]["wake_reason"] == "result"
    persisted = (tmp_path / "events" / "runtime.jsonl").read_text()
    assert "private continuation" not in persisted


@pytest.mark.parametrize(
    ("result", "observation", "expected"),
    (
        ({"status": "rejected"}, {}, "failure"),
        ({"result": {"idempotent_replay": True}}, {}, "replay"),
        (
            {},
            {"state": "exited", "returncode": 7, "timed_out": False},
            "failure",
        ),
        (
            {},
            {
                "state": "timed_out_terminated",
                "returncode": -15,
                "timed_out": True,
            },
            "timeout",
        ),
        (
            {},
            {"state": "exited", "returncode": 0, "timed_out": False},
            "result",
        ),
        (
            {},
            {"state": "exited", "returncode": -15, "timed_out": False},
            "signal",
        ),
    ),
)
def test_approved_program_tool_wake_reason_is_evidence_based(
    result, observation, expected
):
    assert (
        _execution_wake_reason(result=result, process_observation=observation)
        == expected
    )


def test_approved_program_dispatch_rejection_still_closes_wait_event(tmp_path):
    timeline = []

    class RejectingHost(_ExecutionWaitHost):
        def dispatch(self, **_kwargs):
            self.timeline.append("engine-wait-start")
            raise ContractError("synthetic admission rejection")

    responses = iter(
        (
            {
                "model": "qwen3.8-max",
                "choices": [
                    {
                        "finish_reason": "tool_calls",
                        "message": {
                            "role": "assistant",
                            "content": "Trying the approved node.",
                            "reasoning_content": "private continuation",
                            "tool_calls": [
                                {
                                    "id": "execute-rejected",
                                    "type": "function",
                                    "function": {
                                        "name": "execute_approved_program_node",
                                        "arguments": '{"node_id":"node-1"}',
                                    },
                                }
                            ],
                        },
                    }
                ],
                "usage": {},
            },
            {
                "model": "qwen3.8-max",
                "choices": [
                    {
                        "finish_reason": "stop",
                        "message": {
                            "role": "assistant",
                            "content": "The host rejected it.",
                        },
                    }
                ],
                "usage": {},
            },
        )
    )
    config = Qwen38MaxConfigV1()
    session = Qwen38MaxToolSession(
        transport=lambda _payload: next(responses),
        messages=[{"role": "user", "content": "Run the approved node."}],
        config=config,
    )
    store = RuntimeEventStore(
        tmp_path / "events" / "runtime.jsonl",
        session_id="protocol-session",
    )
    host = RejectingHost(timeline)
    envelope, request_context, network = _run_contracts(
        host, config, chemistry_engine_calls=1
    )

    ToolLoopRunner(host=host, event_store=store).run(
        session=session,
        envelope=envelope,
        request_context=request_context,
        provider_budget=network,
    )

    events = store.read_events()
    woke = [event for event in events if event.kind == "tool_woke"]
    assert len(woke) == 1
    assert woke[0].payload["wake_reason"] == "failure"
    assert store.state().active_waits == {}


def test_execution_wait_admission_failure_is_canonical_and_terminal(tmp_path):
    timeline = []

    class ExhaustedHost(_ExecutionWaitHost):
        def execution_wait_timeout_seconds(self):
            raise ContractError("bounded execution window is exhausted")

        def dispatch(self, **_kwargs):  # pragma: no cover - must not launch
            raise AssertionError(
                "dispatch must not follow failed wait admission"
            )

    responses = iter(
        (
            {
                "model": "qwen3.8-max",
                "choices": [
                    {
                        "finish_reason": "tool_calls",
                        "message": {
                            "role": "assistant",
                            "content": "Trying the approved node.",
                            "reasoning_content": "private continuation",
                            "tool_calls": [
                                {
                                    "id": "execute-exhausted",
                                    "type": "function",
                                    "function": {
                                        "name": "execute_approved_program_node",
                                        "arguments": '{"node_id":"node-1"}',
                                    },
                                }
                            ],
                        },
                    }
                ],
                "usage": {},
            },
            {
                "model": "qwen3.8-max",
                "choices": [
                    {
                        "finish_reason": "stop",
                        "message": {
                            "role": "assistant",
                            "content": "The execution window is exhausted.",
                        },
                    }
                ],
                "usage": {},
            },
        )
    )
    config = Qwen38MaxConfigV1()
    session = Qwen38MaxToolSession(
        transport=lambda _payload: next(responses),
        messages=[{"role": "user", "content": "Run the approved node."}],
        config=config,
    )
    store = RuntimeEventStore(
        tmp_path / "events" / "runtime.jsonl",
        session_id="protocol-session",
    )
    host = ExhaustedHost(timeline)
    envelope, request_context, network = _run_contracts(
        host, config, chemistry_engine_calls=1
    )

    result = ToolLoopRunner(host=host, event_store=store).run(
        session=session,
        envelope=envelope,
        request_context=request_context,
        provider_budget=network,
    )

    assert result.failed_tool_calls == 1
    assert result.terminal_state == "blocked"
    events = store.read_events()
    failed = [event for event in events if event.kind == "tool_failed"]
    assert len(failed) == 1
    canonical = failed[0].payload["canonical_result"]
    assert canonical["status"] == "rejected"
    assert canonical["error_class"] == "ContractError"
    assert canonical["message"] == "bounded execution window is exhausted"
    assert not any(event.kind == "tool_waiting" for event in events)
    assert not any(event.kind == "tool_woke" for event in events)
    assert events[-1].kind == EventKind.RUNTIME_TERMINATED.value
    assert store.state().active_waits == {}
    assert (
        "private continuation"
        not in (tmp_path / "events" / "runtime.jsonl").read_text()
    )


def test_provider_deadline_failure_terminalizes_inside_episode_reserve(
    tmp_path,
):
    class DeadlineTransport:
        def __init__(self):
            self.effective = 300.0

        def set_timeout_seconds(self, value):
            self.effective = min(300.0, float(value))

        def public_deadline_record(self):
            return {
                "connect_seconds": 15.0,
                "first_event_seconds": 90.0,
                "inter_event_seconds": 90.0,
                "absolute_turn_seconds": 300.0,
                "effective_absolute_turn_seconds": self.effective,
            }

        def __call__(self, _payload):
            raise DeepSeekTransportError(
                "turn_deadline_exceeded",
                timeout_phase="absolute",
                timeout_seconds=self.effective,
            )

    config = Qwen38MaxConfigV1()
    session = Qwen38MaxToolSession(
        transport=DeadlineTransport(),
        messages=[{"role": "user", "content": "Inspect capabilities."}],
        config=config,
    )
    store = RuntimeEventStore(
        tmp_path / "events" / "runtime.jsonl",
        session_id="protocol-session",
    )
    host = _DispatchSpyHost()
    envelope, request_context, network = _run_contracts(host, config)

    result = ToolLoopRunner(host=host, event_store=store).run(
        session=session,
        envelope=envelope,
        request_context=request_context,
        provider_budget=network,
    )

    assert result.terminal_state == "failed"
    assert store.read_events()[-1].kind == EventKind.RUNTIME_TERMINATED.value
    attempt = next(
        event
        for event in store.read_events()
        if event.kind == EventKind.API_ATTEMPT_OBSERVED.value
    )
    assert attempt.payload["status"] == "timeout"
    failure = attempt.payload["transport_failure"]
    assert failure["timeout_phase"] == "absolute"
    assert failure["requested_provider"] == "alibaba-token-plan"
    assert failure["requested_model"] == "qwen3.8-max"
    assert failure["requested_reasoning_effort"] == "xhigh"
    assert failure["observed_reasoning_effort"] == "not_reported"
    assert failure["transport_deadlines"][
        "effective_absolute_turn_seconds"
    ] == pytest.approx(27.0, abs=0.01)


def test_raw_reset_is_sanitized_and_terminalizes_runtime(
    tmp_path, monkeypatch
):
    import chemsmart.agent.runtime.alibaba as alibaba

    def reset_open(*_args, **_kwargs):
        raise ConnectionResetError("raw peer authorization secret")

    monkeypatch.setattr(alibaba, "open_bounded_https_response", reset_open)
    config = Qwen38MaxConfigV1()
    session = Qwen38MaxToolSession(
        transport=AlibabaTokenPlanHttpsTransport(api_key="sk-sp-test"),
        messages=[{"role": "user", "content": "Inspect capabilities."}],
        config=config,
    )
    store = RuntimeEventStore(
        tmp_path / "events" / "runtime.jsonl",
        session_id="protocol-session",
    )
    host = _DispatchSpyHost()
    envelope, request_context, network = _run_contracts(host, config)

    result = ToolLoopRunner(host=host, event_store=store).run(
        session=session,
        envelope=envelope,
        request_context=request_context,
        provider_budget=network,
    )

    assert result.terminal_state == "failed"
    events = store.read_events()
    assert events[-1].kind == EventKind.RUNTIME_TERMINATED.value
    attempt = next(
        event
        for event in events
        if event.kind == EventKind.API_ATTEMPT_OBSERVED.value
    )
    assert attempt.payload["transport_failure"]["error_class"] == "transport"
    persisted = (tmp_path / "events" / "runtime.jsonl").read_text()
    assert "authorization secret" not in persisted
