"""A refusal reaches the durable stream as a failure report.

Over NOVEL-3 (2026-09-05) the agent met 31 refusals and the stream
carried ``failure_report`` on none of them: the loop had the wire, one
raise site used it, and the refusals sessions actually met -- a unit
outside the typed vocabulary, a required output with no producer, an
impossible electronic state, a re-used artifact id, a derivation on an
unbound parent -- were plain messages. A refusal is where the model is
taught; a report that never renders teaches nothing.
"""

from __future__ import annotations

import json

import pytest

from chemsmart.agent._contracts import RoutedContractError
from chemsmart.agent.identity import refuse_impossible_electronic_state
from chemsmart.agent.loop import ToolLoopRunner
from chemsmart.agent.runtime.alibaba import (
    Qwen38MaxConfigV1,
    Qwen38MaxToolSession,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.scientific_toolchain import (
    TOOLCHAIN_REFUSAL_CAUSES,
    TOOLCHAIN_REFUSAL_INVARIANTS,
    ScientificToolchainContractError,
)
from tests.agent.provider_fakes import _DispatchSpyHost, _run_contracts

_TOOL_TURN = {
    "id": "first-turn",
    "model": "qwen3.8-max",
    "choices": [
        {
            "finish_reason": "tool_calls",
            "message": {
                "role": "assistant",
                "content": "",
                "reasoning_content": "",
                "tool_calls": [
                    {
                        "id": "call-1",
                        "type": "function",
                        "function": {
                            "name": "inspect_program_capability",
                            "arguments": '{"program": "orca"}',
                        },
                    }
                ],
            },
        }
    ],
    "usage": {"prompt_tokens": 10, "completion_tokens": 5},
}
_STOP_TURN = {
    "id": "second-turn",
    "model": "qwen3.8-max",
    "choices": [
        {
            "finish_reason": "stop",
            "message": {
                "role": "assistant",
                "content": "The unit was refused; I will redeclare.",
                "reasoning_content": "",
            },
        }
    ],
    "usage": {"prompt_tokens": 12, "completion_tokens": 8},
}


class _RefusingHost(_DispatchSpyHost):
    def dispatch(self, *, tool_name, arguments, **_kwargs):
        self.dispatched.append((tool_name, arguments))
        raise RoutedContractError(
            gate="declaration.unit_is_in_the_typed_vocabulary",
            invariant="a declared unit is one the typed plane measures.",
            diagnosis="declared unit 'furlong' is not in the vocabulary.",
            route="declare a unit from the typed vocabulary.",
        )


@pytest.mark.capability("gate:tool.dispatch.rejected")
def test_a_routed_refusal_rides_the_tool_failed_event(tmp_path):
    """Production path: provider turn, host refusal, durable event."""

    turns = iter((_TOOL_TURN, _STOP_TURN))

    def transport(_payload):
        return next(turns)

    config = Qwen38MaxConfigV1()
    session = Qwen38MaxToolSession(
        transport=transport,
        messages=[{"role": "user", "content": "Plan the workflow."}],
        config=config,
    )
    stream = tmp_path / "events" / "runtime.jsonl"
    store = RuntimeEventStore(stream, session_id="protocol-session")
    host = _RefusingHost()
    envelope, request_context, network = _run_contracts(host, config)

    result = ToolLoopRunner(host=host, event_store=store).run(
        session=session,
        envelope=envelope,
        request_context=request_context,
        provider_budget=network,
        should_stop=lambda: False,
    )
    assert result.failed_tool_calls == 1

    events = [
        json.loads(line)
        for line in stream.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    (failed,) = [event for event in events if event["kind"] == "tool_failed"]
    report = failed["payload"]["failure_report"]
    assert report["gate"] == "declaration.unit_is_in_the_typed_vocabulary"
    assert set(report) == {"gate", "invariant", "diagnosis", "route", "cost"}
    assert report["cost"] == "no engine call"


def test_the_refusals_novel_3_met_are_failure_reports(tmp_path):
    from tests.agent.test_a_guide_opens_when_something_asks import _host

    host = _host(tmp_path, approved_requested_observable_declarations=[])
    with pytest.raises(RoutedContractError) as caught:
        host.dispatch(
            turn_id="t1",
            tool_name="declare_requested_observable",
            arguments={
                "observables": [
                    {
                        "observable_id": "j",
                        "unit": "furlong",
                        "meaning": "an exchange coupling",
                    }
                ]
            },
        )
    report = caught.value.failure_report
    assert report["gate"] == "declaration.unit_is_in_the_typed_vocabulary"
    assert "'furlong'" in report["diagnosis"]
    assert "'cm^-1'" in report["route"]

    with pytest.raises(RoutedContractError) as caught:
        refuse_impossible_electronic_state(
            ("C", "H", "H", "H", "H"), 0, 2, context="scientific identity"
        )
    report = caught.value.failure_report
    assert report["gate"] == "identity.electronic_state_is_possible"
    assert "10 electrons cannot leave 1 unpaired" in report["diagnosis"]
    assert "no state is preferred" in report["route"]

    from chemsmart.agent.tool_runtime import _artifact_id_taken

    report = _artifact_id_taken("sulfone", ["sulfone"]).failure_report
    assert report["gate"] == "artifact.id_is_unused"
    assert "fresh id" in report["route"]


def test_a_typed_toolchain_cause_is_a_failure_report():
    assert set(TOOLCHAIN_REFUSAL_INVARIANTS) == set(TOOLCHAIN_REFUSAL_CAUSES)
    routed = ScientificToolchainContractError(
        "required output(s) have no producer: ['j']",
        cause="required_output_has_no_producer",
        next_legal_route="declare it as an output of the node computing it",
    )
    assert routed.failure_report["gate"] == "required_output_has_no_producer"
    assert routed.failure_report["route"].startswith("declare it")
    assert "producer" in routed.failure_report["invariant"]
    bare = ScientificToolchainContractError("no cause typed")
    # A cause-less refusal still reports, under the one general gate.
    assert bare.failure_report["gate"] == "plan.admission"


def test_the_refusals_reach_1_met_are_now_reports(tmp_path):
    """Twenty-one of REACH-1's twenty-six refusals were bare: an unknown
    artifact id (a digest or a node name off the host's own record), an
    unroutable run reference, a cause-less plan refusal, a periodic
    dihedral refused on its range."""

    from tests.agent.test_a_guide_opens_when_something_asks import _host

    host = _host(tmp_path, approved_requested_observable_declarations=[])
    digest = "0ac07d4977dfbe94" + "a" * 48
    with pytest.raises(RoutedContractError) as caught:
        host._artifact(digest)
    report = caught.value.failure_report
    assert report["gate"] == "artifact.id_is_registered"
    assert "content digest" in report["diagnosis"]
    assert f"<program>-result-{digest[:16]}" in report["diagnosis"]
    assert "inspect_run" in report["route"]
    with pytest.raises(RoutedContractError) as caught:
        host._artifact("ph3-cat-opt2")
    assert (
        "not a registered artifact id"
        in caught.value.failure_report["diagnosis"]
    )

    from chemsmart.agent.tool_runtime import _collect_json_violations

    findings: list[str] = []
    _collect_json_violations(
        "artifact_id",
        digest,
        {"type": "string", "pattern": "^[a-z][a-z0-9_.-]*$"},
        findings,
    )
    assert any("content digest" in item for item in findings)

    bare = ScientificToolchainContractError("expression outputs do not name")
    assert bare.failure_report["gate"] == "plan.admission"
    assert bare.failure_report["diagnosis"].startswith("expression outputs")

    from chemsmart.agent.tool_runtime import _periodic_degrees

    assert _periodic_degrees(285.0) == pytest.approx(-75.0)
    assert _periodic_degrees(-190.0) == pytest.approx(170.0)
