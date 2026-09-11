"""A goal whose first cycle computed nothing could never compute later.

An analysis-only first cycle creates the goal record from durable
evidence with ``identity=""``, empty conditions, and no initial review
digest -- not because the human approved "no molecule and no solvent",
but because no plan with a molecule in it had been displayed yet.
Revision admission then compared the first executable plan's identity
against that empty string and returned the goal to the human for
"binding different molecular identities than the goal approved".

So the shape OPEN-2's ino3-qwen actually ran -- deliver from registered
results at cycle 1, then be woken to compute -- could deliver its first
calculation only by failing. Found in review 2026-09-09, before it cost
a window.

The repair is not to waive the identity check. An unbound scope is
*established* by the first executable revision and enforced against
every revision after it, and because the goal record is digest-bound
and never rewritten, the binding lives on the ledger.
"""

from __future__ import annotations

import pytest

from chemsmart.agent.driver import _goal_bound_scope
from chemsmart.agent.goal import (
    GOAL_SCHEMA_VERSION,
    GoalLedger,
    GoalRecordV1,
    admit_revision,
    goal_scope_is_unbound,
)

pytestmark = pytest.mark.capability("rule:wake.goal_authority")

_RUN = "goals/goal-scope-01/runs/cycle-1"


def _goal(**overrides):
    body = {
        "schema_version": GOAL_SCHEMA_VERSION,
        "goal_id": "goal-scope-01",
        "task_spec_sha256": "a" * 64,
        "scientific_identity_sha256": "",
        "conditions": {"solvents": (), "thermochemistry": ()},
        "envelope": {
            "allowed_program_engines": (("orca", ("cpu",)),),
            "max_engine_calls": 6,
            "episode_wall_time_seconds": 10800.0,
        },
        "max_revisions": 5,
        "granted_by": "claude-owner-delegated-reviewer",
        "initial_review_sha256": "",
        "created_at": "2026-09-09T00:00:00+00:00",
    }
    body.update(overrides)
    return GoalRecordV1(**body)


def _review(*, solvent=None, temperature=298.15):
    settings = "orca:\n  gas:\n    functional: b3lyp\n"
    if solvent:
        settings = (
            "orca:\n  solv:\n    solvent: %s\n    model: smd\n" % solvent
        )
    return {
        "execution_envelope": {
            "allowed_program_engines": (("orca", ("cpu",)),),
        },
        "node_reviews": ({"project_settings_text": settings},),
        "scientific_toolchain_plan": {
            "analysis_nodes": (
                {
                    "analysis_kind": "thermochemistry",
                    "temperature_k": temperature,
                    "pressure_atm": 1.0,
                    "concentration_mol_l": None,
                },
            )
        },
    }


def _wake_stream(tmp_path):
    path = tmp_path / "wake-events.jsonl"
    path.write_text(
        '{"kind": "run_outcome_inspected", "payload": {"run": "%s"}}\n' % _RUN,
        encoding="utf-8",
    )
    return path


def _admit(tmp_path, goal, review, identity, *, bound=(None, None)):
    return admit_revision(
        goal=goal,
        budgets=GoalLedger(tmp_path / "g").budgets(goal),
        revision_review=review,
        revision_scientific_identity_sha256=identity,
        session_events_path=_wake_stream(tmp_path),
        previous_run_reference=_RUN,
        prior_outcome_evidence_hashes=("e" * 64,),
        bound_scientific_identity_sha256=bound[0],
        bound_conditions=bound[1],
    )


def test_a_goal_with_no_displayed_review_has_no_scope_to_preserve():
    assert goal_scope_is_unbound(_goal())
    assert not goal_scope_is_unbound(_goal(initial_review_sha256="c" * 64))


def test_the_first_executable_revision_establishes_the_scope(tmp_path):
    verdict = _admit(tmp_path, _goal(), _review(), "b" * 64)
    assert verdict.admitted, verdict.reasons
    assert verdict.checks["scope_bound_here"]
    assert verdict.bound_scientific_identity_sha256 == "b" * 64
    assert verdict.bound_conditions is not None


def test_a_later_revision_is_held_to_what_the_first_one_bound(tmp_path):
    """Establishing is not waiving: the check still bites, one cycle on."""

    first = _admit(tmp_path, _goal(), _review(), "b" * 64)
    bound = (
        first.bound_scientific_identity_sha256,
        first.bound_conditions,
    )
    same = _admit(tmp_path, _goal(), _review(), "b" * 64, bound=bound)
    assert same.admitted, same.reasons
    assert not same.checks.get("scope_bound_here")
    assert not same.bound_scientific_identity_sha256

    drifted = _admit(tmp_path, _goal(), _review(), "d" * 64, bound=bound)
    assert not drifted.admitted
    assert not drifted.checks["identity_preserved"]

    resolvated = _admit(
        tmp_path,
        _goal(),
        _review(solvent="acetonitrile"),
        "b" * 64,
        bound=bound,
    )
    assert not resolvated.admitted
    assert not resolvated.checks["conditions_preserved"]


def test_a_goal_that_displayed_a_review_still_compares(tmp_path):
    """The ordinary path is untouched: a bound scope is never re-bound."""

    goal = _goal(
        scientific_identity_sha256="b" * 64,
        initial_review_sha256="c" * 64,
        conditions={
            "solvents": (),
            "thermochemistry": ((298.15, 1.0, None),),
        },
    )
    assert _admit(tmp_path, goal, _review(), "b" * 64).admitted
    changed = _admit(tmp_path, goal, _review(), "d" * 64)
    assert not changed.admitted
    assert not changed.checks["identity_preserved"]


def test_the_binding_is_read_back_from_the_ledger(tmp_path):
    ledger = GoalLedger(tmp_path / "goal")
    ledger.create(_goal())
    assert _goal_bound_scope(ledger) == (None, None)
    ledger.append(
        "goal_scope_bound",
        {
            "cycle": 2,
            "scientific_identity_sha256": "b" * 64,
            "conditions": {"solvents": (), "thermochemistry": ()},
        },
    )
    identity, conditions = _goal_bound_scope(ledger)
    assert identity == "b" * 64
    assert conditions is not None


def test_the_driver_carries_its_own_analysis_evidence_to_admission(tmp_path):
    """The whole wire, driven rather than supplied.

    Every other case in this file hands `admit_revision` a reference and
    a `wake_embedded_run` of its own making. That fabrication hid two
    real breaks at once: the driver wrote the reference relative to the
    workspace while every consumer resolves it under `.chemsmart-agent`,
    and `_wake_context` gated `previous_run` on there having been an
    engine run at all -- so an analysis-only cycle embedded nothing and
    admission compared a real reference against an empty one.

    An independent audit found both. This drives the driver: cycle 1
    delivers from registered results, cycle 2 plans a calculation, and
    the assertion is that the goal reaches execution rather than being
    returned for never having read an outcome.
    """

    from .test_the_goal_loop_recovers_or_returns import (
        _delivery_rows,
        _execute,
        _loop,
        _planning_session,
        _review_payload,
    )

    contexts = []

    def capture(inner):
        def step(workspace, kwargs):
            contexts.append(kwargs["goal_context"])
            return inner(workspace, kwargs)

        return step

    result = _loop(
        tmp_path,
        sessions=[
            # Cycle 1: claims and a decision over registered results,
            # no workflow. This is the shape SUFFICIENCY-1 actually ran.
            capture(
                _planning_session(
                    "live-1",
                    terminal="complete",
                    wake_rows=(
                        [
                            {
                                "kind": "requested_observable_declared",
                                "payload": {
                                    "observables": [
                                        {
                                            "observable_id": "gap",
                                            "unit": "kJ/mol",
                                            "dimension": (1, 0, 0, 0, 0, 0),
                                            "meaning": "the gap",
                                        }
                                    ]
                                },
                            }
                        ]
                        + list(_delivery_rows())
                    ),
                )
            ),
            # Cycle 2: the first executable plan of the goal.
            capture(_planning_session("live-2", review=_review_payload())),
        ],
        executes=[_execute(tmp_path, failed=False, status="completed")],
        max_revisions=3,
    )

    # The wake carried the evidence, in the base the consumer resolves.
    assert contexts[1]["previous_run"], "the wake embedded no evidence"
    assert not contexts[1]["previous_run"].startswith(".chemsmart-agent")
    # And the calculation was admitted rather than returned for never
    # having read an outcome.
    assert "never held the typed outcome" not in " ".join(result.reasons)
    kinds = [entry["kind"] for entry in _ledger_entries(tmp_path)]
    assert "analysis_evidence_recorded" in kinds
    assert "run_recorded" in kinds


def _ledger_entries(tmp_path):
    import json

    path = (
        tmp_path / "ws" / ".chemsmart-agent" / "goals" / "goal-t1"
    ) / "ledger.jsonl"
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
    ]


def test_the_reference_the_gate_needs_comes_from_the_ledger(tmp_path):
    """The connection the tests above assumed, and production lacked.

    Every test in this file hands `admit_revision` a prior run
    reference of its own making. Production reads it from the ledger
    through `_previous_run_reference`, which knew only `run_recorded` --
    and an analysis-only cycle records no run. So the first executable
    revision after such a cycle failed `evidence_read` against an
    outcome named "unrecorded", and the scope repair these tests cover
    could never be reached in a real goal. A fabricated fixture hid the
    missing wire, which is what a fixture is for and why it is not
    evidence.

    A cycle that claimed from registered results and recorded a decision
    left durable typed evidence; it simply launched no engine. The
    reference names that evidence.
    """

    from chemsmart.agent.driver import _previous_run_reference

    ledger = GoalLedger(tmp_path / "goal")
    ledger.create(_goal())
    assert _previous_run_reference(ledger) == ""

    ledger.append(
        "analysis_evidence_recorded",
        {"cycle": 1, "evidence": ".chemsmart-agent/runs/live-abc"},
    )
    reference = _previous_run_reference(ledger)
    assert reference == ".chemsmart-agent/runs/live-abc"

    # And with the wake embedding it, the gate that refused the first
    # calculation is satisfied by evidence that actually exists.
    verdict = admit_revision(
        goal=_goal(),
        budgets=ledger.budgets(_goal()),
        revision_review=_review(),
        revision_scientific_identity_sha256="b" * 64,
        session_events_path=_wake_stream(tmp_path),
        previous_run_reference=reference,
        prior_outcome_evidence_hashes=("e" * 64,),
        wake_embedded_run=reference,
    )
    assert verdict.admitted, verdict.reasons
    assert verdict.checks["evidence_read"]

    # A later engine run still wins as the more recent evidence.
    ledger.append(
        "run_recorded", {"cycle": 2, "run": "goals/goal-scope-01/runs/cycle-2"}
    )
    assert (
        _previous_run_reference(ledger) == "goals/goal-scope-01/runs/cycle-2"
    )
