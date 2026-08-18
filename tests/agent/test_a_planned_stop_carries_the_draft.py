"""A planned stop must satisfy the event store's own termination contract.

The stop-point repair downgrades a session to `planned` when its workflow is
not approvable, so the model's "ready for approval" prose can no longer stand
uncorrected. Its first live run crashed the whole session instead: the branch
had a green preflight in hand, passed those *preflight* receipts to
`terminate`, and the event store refused -- "planned termination requires the
latest workflow draft". Observation 4 of the paper task died on it.

That is the check-the-production-shape lesson again. The unit checks for the
repair exercised the reason helper and the host method; nothing exercised the
loop against the terminate contract, which is where the two meet.

This drives the real `ToolLoopRunner` with a synthetic provider to the exact
stop that crashed: a session that planned a workflow, holds a green preflight,
is not approvable, and then sends a message with no tool calls.
"""

from __future__ import annotations

from types import SimpleNamespace

from chemsmart.agent._contracts import canonical_sha256
from chemsmart.agent.loop import ToolLoopRunner
from chemsmart.agent.request_context import (
    ProviderNetworkBudgetV1,
    build_request_context_provenance,
)
from chemsmart.agent.runtime.alibaba import (
    Qwen38MaxConfigV1,
    Qwen38MaxToolSession,
)
from chemsmart.agent.runtime.contracts import (
    ResourceBudgetV1,
    TaskEnvelopeV1,
    TaskPhase,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.runtime.events import EventKind

_DRAFT = "d" * 64
_PREFLIGHT = "f" * 64


class _UnapprovableHost:
    """The host state observation 4 actually stopped in."""

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
        self.task_spec_sha256s = frozenset({canonical_sha256("task")})

    def record_seeded_evidence(self, _turn_id: str) -> None:
        return None

    def dispatch(self, **_kwargs):  # pragma: no cover - loop never gets here
        raise AssertionError("this session makes no tool calls")

    def completion_receipts_for_latest_preflight(self):
        return (_PREFLIGHT,)

    def latest_workflow_draft_receipt(self):
        return _DRAFT

    def unapproved_workflow_summary(self):
        return {
            "workflow_id": "bal-conformers",
            "blocking_node_ids": ("cal-conf-psimin",),
            "deferred_node_ids": (),
            "workflow_blocked_reason": "",
        }


def _run(tmp_path, host):
    config = Qwen38MaxConfigV1()
    budget = ResourceBudgetV1(
        max_input_tokens_per_request=config.context_tokens,
        max_output_tokens_per_request=config.max_output_tokens,
        max_tool_calls=8,
        wall_time_seconds=30.0,
        chemistry_engine_calls=0,
    )
    envelope_body = {
        "schema_version": "chemsmart.task-envelope.v1",
        "task_id": "planned-stop-regression",
        "session_id": "planned-stop",
        "turn_id": "planned-stop.turn-1",
        "request_sha256": canonical_sha256("request"),
        "cwd_sha256": canonical_sha256("workspace"),
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
        task_spec_sha256=canonical_sha256("task"),
        prompt_sha256=envelope.request_sha256,
        tool_schema_sha256=host.surface.tool_schema_sha256,
        configuration_sha256=canonical_sha256("configuration"),
        provider_budget_sha256=network.budget_sha256,
    )
    store = RuntimeEventStore(
        tmp_path / "events" / "runtime.jsonl", session_id="planned-stop"
    )
    # The draft this session planned, as the reducer learns it in production.
    store.append(
        turn_id=envelope.turn_id,
        kind=EventKind.WORKFLOW_PLANNED.value,
        payload={"receipt_sha256": _DRAFT, "status": "planned"},
    )
    responses = iter(
        (
            {
                "id": "final",
                "model": "qwen3.8-max",
                "choices": [
                    {
                        "finish_reason": "stop",
                        "message": {
                            "role": "assistant",
                            "content": (
                                "The workflow is ready for human approval."
                            ),
                        },
                    }
                ],
                "usage": {"prompt_tokens": 5, "completion_tokens": 5},
            },
        )
    )
    session = Qwen38MaxToolSession(
        transport=lambda _payload: next(responses),
        messages=[{"role": "user", "content": "Plan the workflow."}],
        config=config,
    )
    result = ToolLoopRunner(host=host, event_store=store).run(
        session=session,
        envelope=envelope,
        request_context=request_context,
        provider_budget=network,
    )
    return result, store


def test_the_stop_is_planned_and_the_blocking_node_is_named(tmp_path):
    import json

    result, store = _run(tmp_path, _UnapprovableHost())

    assert result.terminal_state == "planned"
    terminate = next(
        json.loads(line)
        for line in reversed(store.path.read_text().splitlines())
        if '"runtime-terminal"' in line or '"reason"' in line
    )
    reason = str(terminate.get("payload", {}).get("reason", ""))
    assert "not approvable" in reason
    assert "cal-conf-psimin" in reason


def test_the_terminate_record_carries_the_draft_not_the_preflight(tmp_path):
    """The exact contract violation that killed observation 4."""

    _result, store = _run(tmp_path, _UnapprovableHost())

    state = store.state()
    assert state.terminal_state == "planned"
    # The event store accepted the termination at all, which it refused when
    # the branch offered the preflight receipt instead of the draft.
    assert state.workflow_receipts[-1] == _DRAFT
