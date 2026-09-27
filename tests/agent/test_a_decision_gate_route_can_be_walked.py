"""The decision gate's route, walked as written, ends in a citation it accepts.

``record_scientific_decision`` refuses a postprocessing citation that is not
a receipt the host minted, and names a route: cite a receipt a tool returned
in this session, or "one inspect_run shows on a recorded run of this goal".
What a decision may cite from a recorded run is an extraction,
thermochemistry, expression, validation or claim receipt
(``_recorded_run_receipt``); what inspect_run shows of a run is each node's
event hashes, artifact digests and anomaly receipts -- none of those five --
so a session that followed the route cited an anomaly receipt and was
refused again (R11 truth, census 3 and 6: 12 anomaly and 10 verification
citations of this shape). R11 truth-3, item 4.

Driven through the goal loop with real hosts: cycle 1's run judges the
PySCF O2 singlet's stability criterion (it fails) and records an anomaly;
cycle 2's woken session cites a digest nothing minted, reads the refusal,
and walks its route as a model would -- the tool the route names and the
receipts that tool shows, or else the receipts of the kinds the route
names, as the wake hands them to it.
"""

from __future__ import annotations

import json
from types import SimpleNamespace

import pytest

from chemsmart.agent.runtime.event_store import RuntimeEventStore

from .test_a_decision_gate_names_what_the_host_minted import (
    _GATE,
    _decision,
    _failed_rows,
    _goal,
    _real_turns,
    _session_stream,
)
from .test_a_failed_criterion_is_a_finding_the_goal_can_deliver import (
    _TASK,
    _engine_prefix,
    _run_turns,
    _stream_rows,
)
from .test_a_partial_delivery_ends_its_session import (
    _call,
    _o2r_turns,
    _turn,
)
from .test_the_goal_loop_recovers_or_returns import (
    _planning_session,
    _review_payload,
)

pytestmark = pytest.mark.capability("tool:record_scientific_decision")

_ANOMALY = "a7" * 32


def _run_with_an_anomaly(tmp_path):
    """Cycle 1's run: an anomaly recorded on its engine node, then the
    approved chain over the archived result, whose stability criterion
    fails; a run's stream holds no decision."""

    def step(run_directory):
        from chemsmart.agent.terminal_states import (
            derive_run_outcome,
            read_run_events,
        )

        _engine_prefix(tmp_path, run_directory)
        events = run_directory / "events.jsonl"
        node = derive_run_outcome(read_run_events(events)).nodes[0]
        session_id = json.loads(events.read_text().splitlines()[0])[
            "session_id"
        ]
        RuntimeEventStore(events, session_id=session_id).append(
            turn_id="exec-anomaly",
            kind="anomaly_observed",
            payload={
                "receipt_sha256": _ANOMALY,
                "status": "unreplicated",
                "node_id": node.node_id,
                "signal_id": "stationary_point.unexpected_order",
                "record": {
                    "signal_id": "stationary_point.unexpected_order",
                    "status": "unreplicated",
                    "values": {"observed_imaginary_modes": 1},
                },
            },
        )
        _run_turns(
            events,
            lambda artifact_id: [
                turn
                for index, turn in enumerate(
                    _o2r_turns(artifact_id, cite_verdict=False)
                )
                if index != 3
            ],
            session_id=session_id,
            scratch=tmp_path,
        )
        return SimpleNamespace(status="completed", analysis_status="completed")

    return step


def _replies(payload):
    return [
        json.loads(message.get("content") or "{}")
        for message in payload.get("messages") or ()
        if message.get("role") == "tool"
    ]


def _digests_named_receipt(value):
    """Every value a reply carries under a key named ``receipt_sha256``."""

    found = []
    if isinstance(value, dict):
        for key, item in value.items():
            if key == "receipt_sha256" and isinstance(item, str):
                found.append(item)
            else:
                found.extend(_digests_named_receipt(item))
    elif isinstance(value, list):
        for item in value:
            found.extend(_digests_named_receipt(item))
    return found


def _walker(run_reference, wake):
    """A woken session that is refused once and then follows the route."""

    def turns(_artifact_id):
        def walked(payload):
            replies = _replies(payload)
            route = next(
                reply["failure_report"]["route"]
                for reply in replies
                if (reply.get("failure_report") or {}).get("gate") == _GATE
            )
            if "inspect_run" in route:
                # The receipts the named tool shows of the recorded run.
                candidates = _digests_named_receipt(replies[-1])
            else:
                # The receipts of the kinds the route names, as the wake
                # hands a failed verdict's receipts to a woken session.
                candidates = [
                    digest
                    for verdict in (
                        wake["context"]["deliverables"][
                            "unanswered_failed_verdicts"
                        ]
                    )
                    for digest in verdict["receipt_sha256s"]
                ]
            assert candidates, (route, replies[-1])
            return _turn(
                3,
                "Following the route.",
                (_decision(3, postprocessing=candidates[:1]),),
            )

        return [
            lambda payload: _turn(
                1,
                "Citing what the earlier cycle found.",
                (_decision(1, postprocessing=("e" * 64,)),),
            ),
            lambda payload: _turn(
                2,
                "Reading the recorded run.",
                (_call(2, "inspect_run", {"run": run_reference}),),
            ),
            walked,
            lambda payload: _turn(4, "Recorded."),
        ]

    return turns


def test_the_refusal_route_ends_in_a_citation_the_gate_accepts(tmp_path):
    goal_id = "goal-route-walked"
    woken = "live-20260928T010000000000Z-truth3-walk"
    workspace = tmp_path / "ws"
    wake = {}

    def woken_session(workspace, kwargs):
        wake["context"] = kwargs.get("goal_context") or {}
        _real_turns(
            _session_stream(workspace, woken),
            _walker(f"goals/{goal_id}/runs/cycle-1", wake),
            workspace=workspace,
            scratch=workspace.parent / f"scratch-{woken}",
            prior_anomalies=wake["context"].get("anomalies") or (),
        )
        return SimpleNamespace(
            terminal_state="complete",
            run_id=woken,
            task_spec_sha256=_TASK,
            selected_execution_wave=(),
        )

    (tmp_path / "first").mkdir()
    _goal(
        tmp_path,
        goal_id=goal_id,
        sessions=[
            _planning_session(
                "live-20260928T000000000000Z-truth3-plan",
                review=_review_payload(),
            ),
            woken_session,
        ],
        executes=[_run_with_an_anomaly(tmp_path / "first")],
    )

    rows = _stream_rows(_session_stream(workspace, woken))
    # The first citation is refused, once, at this gate: nothing minted it.
    refused = _failed_rows(rows)
    assert len(refused) == 1, [row["failure_report"] for row in refused]
    assert "e" * 8 in refused[0]["failure_report"]["diagnosis"]
    # Walked as written, the route ends in a decision the host recorded.
    decisions = [
        row["payload"]["record"]
        for row in rows
        if row["kind"] == "scientific_decision_recorded"
    ]
    assert len(decisions) == 1, [row["kind"] for row in rows]
