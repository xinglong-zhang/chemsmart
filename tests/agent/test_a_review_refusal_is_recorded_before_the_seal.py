"""A refusal is an output; a session that produced one settles on it.

REACH-1 ino3 (2026-09-06): a planning session reached "host readiness
gates passed" with a ten-node plan and two green previews; the review
was built after the loop had sealed the runtime stream, the builder
refused, the refusal could not be recorded (the store is absorbing),
the ContractError escaped, and the goal settled on a Python error with
its whole grant unspent and the reason lost. execution_review_refused
had never been written in 991 campaign streams.
"""

from __future__ import annotations

import json
from types import SimpleNamespace

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.loop import ToolLoopRunner
from chemsmart.agent.runtime.alibaba import (
    Qwen38MaxConfigV1,
    Qwen38MaxToolSession,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from tests.agent.provider_fakes import _DispatchSpyHost, _run_contracts

_STOP = {
    "id": "turn-1",
    "model": "qwen3.8-max",
    "choices": [
        {
            "finish_reason": "stop",
            "message": {
                "role": "assistant",
                "content": "Planned and previewed; ready for review.",
                "reasoning_content": "",
            },
        }
    ],
    "usage": {"prompt_tokens": 10, "completion_tokens": 5},
}


class _ReadyHost(_DispatchSpyHost):
    def __init__(self, review=None, refusal=""):
        super().__init__()
        self._review = review
        self._refusal = refusal
        self.execution_review_refusal = {}
        self.prepared_execution_review = None

    def unapproved_workflow_summary(self):
        return None

    def completion_receipts_for_latest_preflight(self):
        return ()

    def execution_review_wanted(self):
        return True

    def prepare_execution_review(self):
        if self._refusal:
            self.execution_review_refusal = {
                "workflow_id": "ni-thiolate-redox",
                "reason": self._refusal,
            }
            raise ContractError(self._refusal)
        self.prepared_execution_review = self._review
        return self._review


def _run(tmp_path, host):
    def transport(_payload):
        return _STOP

    config = Qwen38MaxConfigV1()
    session = Qwen38MaxToolSession(
        transport=transport,
        messages=[{"role": "user", "content": "Plan the workflow."}],
        config=config,
    )
    stream = tmp_path / "events" / "runtime.jsonl"
    store = RuntimeEventStore(stream, session_id="protocol-session")
    envelope, request_context, network = _run_contracts(host, config)
    result = ToolLoopRunner(host=host, event_store=store).run(
        session=session,
        envelope=envelope,
        request_context=request_context,
        provider_budget=network,
        should_stop=lambda: False,
    )
    events = [
        json.loads(line)
        for line in stream.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    return result, [event["kind"] for event in events], events


@pytest.mark.capability("rule:wake.refusal_is_a_deliverable")
def test_a_refused_review_is_an_event_before_the_terminal_one(tmp_path):
    refusal = (
        "scientific workflow exceeds bounded engine-call budget: 10 nodes "
        "for 4 calls"
    )
    result, kinds, events = _run(tmp_path, _ReadyHost(refusal=refusal))
    assert "execution_review_refused" in kinds
    assert kinds.index("execution_review_refused") < kinds.index(
        "runtime_terminated"
    )
    (refused,) = [e for e in events if e["kind"] == "execution_review_refused"]
    assert refused["payload"]["reason"] == refusal
    assert refused["payload"]["workflow_id"] == "ni-thiolate-redox"
    # With a draft receipt in hand the loop ends planned; the spy host has
    # none, so it ends blocked -- never complete, never an exception.
    assert result.terminal_state in {"planned", "blocked"}
    (terminal,) = [e for e in events if e["kind"] == "runtime_terminated"]
    assert terminal["payload"]["reason"].startswith(
        "execution review refused: "
    )


def test_a_prepared_review_ends_the_loop_waiting_for_approval(tmp_path):
    review = SimpleNamespace(
        scientific_plan=SimpleNamespace(workflow_id="ni-thiolate-redox"),
        review_sha256="b" * 64,
    )
    host = _ReadyHost(review=review)
    result, kinds, events = _run(tmp_path, host)
    assert result.terminal_state == "waiting_for_approval"
    assert "execution_review_prepared" in kinds
    assert kinds.index("execution_review_prepared") < kinds.index(
        "runtime_terminated"
    )
    assert host.prepared_execution_review is review
