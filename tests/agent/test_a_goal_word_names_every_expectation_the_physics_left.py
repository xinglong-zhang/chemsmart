"""A goal's word names every pre-registered expectation the physics left.

The one word a human reads first must not hide what the run found (owner
ruling, 2026-09-03), and a pre-registered expectation the physics left is
one of the things it names: the delivery is ``achieved_with_observations``.
The settlement read those expectations from the completion its delivery
stood on -- one stream -- while it certified numbers an earlier cycle had
delivered and that cycle's own completion had scored. R11 `truth`'s census
of 254 archived goals found the round pin signing plain ``achieved`` over
CUHK r10/q7 g2-scan-modred (cis-barrier 8.30 kcal/mol against a
pre-registered 2-8, oo160-torsion 141.5 deg against 100-135, both scored
diverged by cycle 1's run completion). Here the host declares, claims and
scores in every cycle, and the goal loop settles.
"""

from __future__ import annotations

import json
import shutil
from types import SimpleNamespace

import pytest

from chemsmart.agent.execution import build_program_execution_receipt
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

from .test_runtime_v2_launch_fence import _reserve
from .test_the_goal_loop_recovers_or_returns import (
    _READ_OUTCOME_ROWS,
    _loop,
    _planning_session,
    _review_payload,
)

pytestmark = pytest.mark.capability("rule:declare.diagnostic_has_standing")

_TASK = "a" * 64
_POINTS = "irc-forward-points"


def _rows(path):
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def _named(step, name):
    """A live session names its stream by session id, so the goal records
    it on its own spine; the shared stand-in names none."""

    def named(workspace, kwargs):
        result = step(workspace, kwargs)
        result.session_id = name
        return result

    return named


def _declaring_session(tmp_path, target, role):
    """Cycle 1's session: its own host declares the band and the count."""

    build = tmp_path / "session-live-1"
    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            build / "events.jsonl", session_id="live-1"
        ),
        artifacts={},
        task_spec_sha256s=(_TASK,),
        approved_workspace=build / "ws",
    )
    band = {
        "observable_id": target,
        "unit": "kcal/mol",
        "meaning": "E(TS) - E(reactant)",
        "expected_low": 2.0,
        "expected_high": 8.0,
        "expectation_basis": "a torsional barrier of a few kcal/mol",
    }
    if role == "diagnostic":
        band = {
            **band,
            "role": "diagnostic",
            "failure_update_rule": "a barrier above 8 means the search "
            "left the torsional ridge",
        }
    reply = host.dispatch(
        turn_id="t1",
        tool_name="declare_requested_observable",
        arguments={
            "observables": [
                band,
                {
                    "observable_id": _POINTS,
                    "unit": "1",
                    "meaning": "points on the forward IRC branch",
                },
            ]
        },
    )
    assert reply["status"] == "ok", reply
    declared = tuple(host.requested_observable_declarations.values())
    step = _planning_session(
        "live-1",
        review=_review_payload(),
        wake_rows=_rows(build / "events.jsonl"),
    )
    return _named(step, "live-1"), declared


def _claims_and_completion(host, claims, *, status, findings=()):
    """The host renders each claim from a literal and scores the delivery."""

    receipts = []
    for index, (claim_id, value, unit) in enumerate(claims):
        literal = host.dispatch(
            turn_id="t1",
            tool_name="evaluate_quantity_expression",
            arguments={
                "expression_id": f"literal-{index}",
                "inputs": [],
                "nodes": [
                    {
                        "node_id": "n1",
                        "operation": "literal",
                        "literal_value": value,
                        "literal_unit": unit,
                    }
                ],
                "output_node_ids": ["n1"],
            },
        )
        receipt = literal["result"]["receipt_sha256"]
        claimed = host.dispatch(
            turn_id="t1",
            tool_name="record_analysis_claims",
            arguments={
                "task_spec_sha256": _TASK,
                "claims": [
                    {
                        "claim_id": claim_id,
                        "receipt_sha256": receipt,
                        "quantity_id": "n1",
                        "display_unit": unit,
                    }
                ],
            },
        )
        assert claimed["status"] == "ok", claimed
        receipts.append(receipt)
    # The walk's own completion: the host scores every declared
    # expectation against the claims this stream rendered.
    host._record_toolchain_completion(
        "b" * 64,
        task_spec_sha256=_TASK,
        source_receipt_sha256s=tuple(receipts),
        status=status,
        findings=findings,
    )


def _run(tmp_path, declared, claims, *, validated):
    """One run: the engine receipt, then the chain's claims and completion."""

    def step(run_directory):
        build = tmp_path / f"build-{run_directory.name}"
        store = RuntimeEventStore(
            build / "events.jsonl", session_id="water-session"
        )
        _, plan, _m, _a, invocation = _reserve(store, build)
        store.record_program_execution_receipt(
            turn_id="turn-1",
            workflow_id=plan.workflow_id,
            run_id="run.water-approval",
            receipt=build_program_execution_receipt(
                invocation,
                execution_state="validated" if validated else "failed",
                exit_status=0 if validated else 1,
                child_exit_status=0 if validated else 1,
                engine_complete=validated,
                validated=validated,
                findings=() if validated else ("execution.process.timeout",),
                **(
                    {
                        "validator_receipt_sha256s": ("e" * 64,),
                        "result_validation_receipt_sha256": "e" * 64,
                    }
                    if validated
                    else {}
                ),
                started_at="2026-08-04T00:00:00+00:00",
                finished_at="2026-08-04T00:00:05+00:00",
            ),
        )
        host = CommandCompiledToolHostV1(
            event_store=store,
            artifacts={},
            task_spec_sha256s=(_TASK,),
            approved_workspace=build / "ws",
            approved_requested_observable_declarations=declared,
        )
        _claims_and_completion(
            host,
            claims,
            status="passed" if validated else "partial",
            findings=(
                ()
                if validated
                else (
                    "ext-irc-forward: expected exactly one registered result "
                    "for producer 'irc-forward'; found 0",
                )
            ),
        )
        run_directory.mkdir(parents=True, exist_ok=True)
        shutil.copy(build / "events.jsonl", run_directory / "events.jsonl")
        return SimpleNamespace(
            status="completed" if validated else "partial",
            analysis_status="completed" if validated else "partial",
        )

    return step


def _delivering_session(tmp_path, declared, claims):
    """A woken session that claims from results already in hand, passes
    its completion, and stops: the goal settles on its stream."""

    build = tmp_path / "session-live-2"
    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            build / "events.jsonl", session_id="live-2"
        ),
        artifacts={},
        task_spec_sha256s=(_TASK,),
        approved_workspace=build / "ws",
        approved_requested_observable_declarations=declared,
    )
    _claims_and_completion(host, claims, status="passed")
    step = _planning_session(
        "live-2", terminal="complete", wake_rows=_rows(build / "events.jsonl")
    )
    return _named(step, "live-2")


@pytest.mark.parametrize(
    "role, path, reclaimed",
    [
        ("requested", "planning", None),
        ("requested", "run", None),
        ("diagnostic", "planning", None),
        # The latest typed word about an id governs it: re-claimed inside
        # its band, the number delivered is the new one, and it agrees.
        ("requested", "planning", 5.0),
    ],
    ids=["requested-planning", "requested-run", "diagnostic", "reclaimed"],
)
def test_a_goal_word_names_the_expectations_an_earlier_cycle_scored(
    tmp_path, role, path, reclaimed
):
    """Cycle 1 declares a band of 2-8 kcal/mol and its run claims 11.2
    kcal/mol, so the run's completion scores the expectation diverged, and
    leaves a declared count undelivered; cycle 2 delivers the count and
    settles without claiming the barrier again -- from its session, or
    from a second run. The barrier is delivered from cycle 1, and the word
    names the expectation the physics left, with the receipt that scored
    it; re-claimed inside the band, it names none."""

    target = "barrier-forward" if role == "requested" else "barrier-guess"
    session_one, declared = _declaring_session(tmp_path, target, role)
    second_claims = [(_POINTS, 42, "1")]
    if reclaimed is not None:
        second_claims.append((target, reclaimed, "kcal/mol"))
    if path == "planning":
        sessions = [
            session_one,
            _delivering_session(tmp_path, declared, second_claims),
        ]
        executes = [
            _run(
                tmp_path,
                declared,
                [(target, 11.2, "kcal/mol")],
                validated=False,
            )
        ]
    else:
        sessions = [
            session_one,
            _named(
                _planning_session(
                    "live-2",
                    review=_review_payload(),
                    wake_rows=_READ_OUTCOME_ROWS,
                ),
                "live-2",
            ),
        ]
        executes = [
            _run(
                tmp_path,
                declared,
                [(target, 11.2, "kcal/mol")],
                validated=False,
            ),
            _run(tmp_path, declared, second_claims, validated=True),
        ]

    result = _loop(tmp_path, sessions=sessions, executes=executes)

    agent = tmp_path / "ws" / ".chemsmart-agent"
    ledger = _rows(agent / "goals" / "goal-t1" / "ledger.jsonl")
    (settled,) = [row for row in ledger if row["kind"] == "goal_settled"]
    reasons = " ".join(settled["payload"]["reasons"])
    scoring = [
        row["payload"]
        for row in _rows(
            agent / "goals" / "goal-t1" / "runs" / "cycle-1" / "events.jsonl"
        )
        if row["kind"] == "analysis_completion_evaluated"
    ]
    assert (
        f"falsified_expectation:{target}" in scoring[-1]["anomaly_output_ids"]
    )
    if reclaimed is None:
        assert result.settlement == "achieved_with_observations", reasons
        assert f"falsified_expectation:{target}" in reasons
        assert scoring[-1]["receipt_sha256"] in (
            settled["payload"]["evidence"].get("receipt_sha256s") or ()
        )
    else:
        assert result.settlement == "achieved", reasons
        assert "falsified_expectation" not in reasons
