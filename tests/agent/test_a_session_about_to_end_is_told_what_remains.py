"""A session about to end with declared observables undelivered is
told, once, what is undelivered and what remains -- and nothing more.

Owner ruling (2026-09-06): never force a session to spend the remaining
budget; when it is about to terminate, provide the remaining budget
purely as informational context. NOVEL-3 po3 stopped with 28 of 30
engine calls and 7 revisions left, its declared observable undelivered
and its diagnosis right; the wake had said the budgets once, in the
attention trough.
"""

from __future__ import annotations

import json
from types import SimpleNamespace

import pytest

from chemsmart.agent.loop import ToolLoopRunner
from chemsmart.agent.rules import rules_by_id
from chemsmart.agent.runtime.alibaba import (
    Qwen38MaxConfigV1,
    Qwen38MaxToolSession,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from tests.agent.provider_fakes import _DispatchSpyHost, _run_contracts

_NOTICE = rules_by_id()["wake.termination_notice"].text


def _stop_turn(ordinal, content):
    return {
        "id": f"turn-{ordinal}",
        "model": "qwen3.8-max",
        "choices": [
            {
                "finish_reason": "stop",
                "message": {
                    "role": "assistant",
                    "content": content,
                    "reasoning_content": "",
                },
            }
        ],
        "usage": {"prompt_tokens": 10, "completion_tokens": 5},
    }


class _NoticingHost(_DispatchSpyHost):
    def __init__(self, notice):
        super().__init__()
        self._notice = notice
        self.asked = 0

    def termination_notice(self):
        self.asked += 1
        return self._notice


def _run(tmp_path, host, turns):
    served = iter(turns)
    calls = {"turns": 0}

    def transport(_payload):
        calls["turns"] += 1
        return next(served)

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
    return result, calls["turns"], events, session.public_history()


@pytest.mark.capability("rule:wake.termination_notice")
def test_the_notice_is_delivered_once_and_a_second_silence_ends(tmp_path):
    notice = {
        "undelivered_declared_observable_ids": ("ddg-c4-c5",),
        "budgets": {
            "engine_calls_remaining": 28,
            "excursion_calls_remaining": 2,
            "wall_seconds_remaining": 6000.0,
            "revisions_remaining": 7,
        },
        "text": _NOTICE + " Undelivered declared observables: ddg-c4-c5.",
    }
    host = _NoticingHost(notice)
    result, turns, events, history = _run(
        tmp_path,
        host,
        (_stop_turn(1, "I will stop here."), _stop_turn(2, "Still done.")),
    )
    assert turns == 2
    assert host.asked == 1
    delivered = [
        event
        for event in events
        if event["kind"] == "termination_notice_delivered"
    ]
    assert len(delivered) == 1
    assert delivered[0]["payload"]["undelivered_declared_observable_ids"] == [
        "ddg-c4-c5"
    ]
    assert delivered[0]["payload"]["budgets"]["engine_calls_remaining"] == 28
    host_messages = [
        message
        for message in history
        if message.get("role") == "user" and _NOTICE in message["content"]
    ]
    assert len(host_messages) == 1
    assert result.final_text == "Still done."


def test_no_notice_when_the_host_has_nothing_to_say(tmp_path):
    host = _NoticingHost(None)
    _result, turns, events, _history = _run(
        tmp_path, host, (_stop_turn(1, "Done."),)
    )
    assert turns == 1
    assert host.asked == 1
    assert not [
        event
        for event in events
        if event["kind"] == "termination_notice_delivered"
    ]


@pytest.mark.capability("rule:wake.execution_wave_decision_pending")
def test_a_pending_execution_boundary_reenters_the_agent_with_typed_context(
    tmp_path,
):
    """The producer notice reaches the loop consumer before it can end."""

    from chemsmart.agent.cohort import build_execution_wave_decision

    decision = build_execution_wave_decision(
        state="undecided",
        workflow_id="w1",
        ready_node_ids=("conf-a-opt", "conf-b-opt"),
    )
    text = (
        rules_by_id()["wake.execution_wave_decision_pending"].text
        + " Workflow w1 currently reports ready: conf-a-opt, conf-b-opt."
    )
    host = _NoticingHost(
        {
            "kind": "execution_wave_decision_pending",
            "execution_wave_decision": decision.public_record(),
            "budgets": {
                "engine_calls_remaining": 4,
                "excursion_calls_remaining": 0,
                "wall_seconds_remaining": 600.0,
                "revisions_remaining": 2,
            },
            "text": text,
        }
    )
    result, turns, events, history = _run(
        tmp_path,
        host,
        (_stop_turn(1, "I am done."), _stop_turn(2, "Still done.")),
    )
    assert turns == 2, "the Agent was not given a decision turn"
    delivered = [
        event
        for event in events
        if event["kind"] == "execution_wave_decision_pending"
    ]
    assert len(delivered) == 1
    assert delivered[0]["payload"]["execution_wave_decision"]["state"] == (
        "undecided"
    )
    assert delivered[0]["payload"]["execution_wave_decision"][
        "ready_node_ids"
    ] == ["conf-a-opt", "conf-b-opt"]
    assert any(
        message.get("role") == "user"
        and "host will not infer a wave from silence" in message["content"]
        for message in history
    )
    assert result.final_text == "Still done."


def _host(tmp_path, **kwargs):
    from tests.agent.test_a_guide_opens_when_something_asks import _host

    return _host(tmp_path, **kwargs)


_J = {
    "observable_id": "j_ohcl",
    "unit": "cm^-1",
    "dimension": (0, 0, 0, 0, 1, 0),
    "meaning": "exchange coupling",
    "expectation_basis": "superexchange",
}


def test_the_host_notice_names_undelivered_ids_and_remaining_lines(
    tmp_path,
):
    host = _host(
        tmp_path,
        approved_requested_observable_declarations=[_J],
        engine_calls_remaining=28,
        excursion_calls_remaining=2,
        wall_seconds_remaining=6000.0,
        revisions_remaining=7,
    )
    notice = host.termination_notice()
    assert notice["undelivered_declared_observable_ids"] == ("j_ohcl",)
    assert notice["budgets"]["engine_calls_remaining"] == 28
    assert notice["text"].startswith(_NOTICE)
    assert "Undelivered declared observables: j_ohcl." in notice["text"]
    assert "engine wall 6000 s." in notice["text"]

    host.analysis_claim_records["r1"] = SimpleNamespace(
        task_spec_sha256="a" * 64,
        claims=(
            SimpleNamespace(
                claim_id="j_ohcl",
                dimension=(0, 0, 0, 0, 1, 0),
                display_value=-12.0,
                display_unit="cm^-1",
            ),
        ),
    )
    assert host.termination_notice() is None


def test_no_notice_outside_a_goal_or_without_budget_or_when_delivered_earlier(
    tmp_path,
):
    assert (
        _host(
            tmp_path, approved_requested_observable_declarations=[_J]
        ).termination_notice()
        is None
    )
    exhausted = _host(
        tmp_path,
        approved_requested_observable_declarations=[_J],
        engine_calls_remaining=0,
        excursion_calls_remaining=0,
        wall_seconds_remaining=0.0,
        revisions_remaining=0,
    )
    assert exhausted.termination_notice() is None
    earlier = _host(
        tmp_path,
        approved_requested_observable_declarations=[_J],
        engine_calls_remaining=3,
        goal_delivered_declared_ids=("j_ohcl",),
    )
    assert earlier.termination_notice() is None
