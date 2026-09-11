"""The task said how good the answer had to be, and nothing held it.

OPEN-2's ino3-qwen (2026-09-07) was asked for a potential versus
ferrocene "good to about plus or minus 0.2 V". It delivered eleven of
eleven declared observables, in the right units, all pre-registered,
with zero engine calls of forty. Its own recorded uncertainties said the
method carries 0.2-0.4 V of systematic error -- "at or beyond the edge
of what this method can promise" -- and named the experiment that would
narrow it, priced at two of the forty calls it had not spent. The host
quoted that sentence into the settlement's own reasons and settled
`achieved`.

Nothing in the goal could hold a requested precision: ``GoalRecordV1``
carried identity, conditions, envelope and revisions, and no tolerance.
``method_resolution`` was the only precision-shaped field and it is a
deadband on the delivered magnitude, which for -0.473 V against 0.2
never fires.

The obligation these tests pin is to *resolve* a requirement, never to
meet one: compute the term that narrows it, re-claim it with the
uncertainty you measured, or refuse it with receipts. A tolerance no method in the envelope can reach is a real
answer and often a better one than a number, so refusing settles
``unreachable_from_evidence`` and the agent keeps the judgement.
"""

from __future__ import annotations

import json

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.delivery import judge_sufficiency
from chemsmart.agent.driver import (
    SUFFICIENCY_SHORT_ROUTE,
    _AnalysisDelivery,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

from .test_the_goal_loop_recovers_or_returns import (
    _loop,
    _planning_session,
)

_KJ_PER_MOL_IN_HARTREE = 2625.4996394798254


def _host(tmp_path):
    return CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="s1"
        ),
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "workspace",
    )


pytestmark = pytest.mark.capability("rule:wake.restate_observable")

#: The couple, its tolerance, and the two legs OPEN-2 delivered.
_TOLERANCE = 0.2
_UNCERTAINTY = 0.30
_PME3 = -0.473
_PH3 = -0.276
_FERROCENIUM = 0.0


def _declaration(**overrides):
    record = {
        "observable_id": "e-couple-vs-fc",
        "unit": "V",
        "dimension": (1, 0, 0, 0, 0, 0, 0, 0, -1),
        "meaning": "the +1/0 couple versus ferrocene in acetonitrile",
        "required_tolerance": _TOLERANCE,
        "tolerance_basis": "good to about plus or minus 0.2 V",
    }
    record.update(overrides)
    return record


def test_the_arithmetic_reproduces_the_case_that_earned_it():
    """The window's own numbers, under the rule they earned.

    All three of its requirements closed `met` on uncertainty 0.2
    against a tolerance of 0.2, basis asserted. They are now `attested`
    and open: the number stands, the word does not.

    The decision question the session was really answering -- does the
    couple clear ferrocenium by more than its own uncertainty -- needs
    no special field. The session showed the route itself: compose the
    margin as a quantity and declare it as its own observable with the
    tolerance the decision needs. That travels through the ordinary
    machinery and carries real provenance, which a hand-typed boundary
    never did.
    """

    couple = _declaration()
    asserted = {"uncertainty": _UNCERTAINTY, "uncertainty_basis": "asserted"}

    # The couple itself: 0.30 against a 0.2 tolerance, on the session's
    # own word. Short, and open.
    pme3 = judge_sufficiency(couple, {**asserted, "display_value": _PME3})
    assert pme3["state"] == "short"

    # The margin, declared as its own observable. The PMe3 leg clears
    # ferrocenium by 0.473 V; a session that needs the decision to a
    # tenth of a volt has an answer, and one that needs it to 0.2 V on
    # an asserted uncertainty does not get to say so.
    margin = {
        "observable_id": "oxidant-margin-fcplus-pme3",
        "unit": "V",
        "required_tolerance": 0.5,
        "tolerance_basis": "the margin only has to separate three "
        "oxidants that sit 0.27 V apart",
    }
    generous = judge_sufficiency(
        margin,
        {
            "uncertainty": _UNCERTAINTY,
            "uncertainty_basis": "measured",
            "uncertainty_evidence_backed": True,
            "display_value": abs(_PME3 - _FERROCENIUM),
        },
    )
    assert generous["state"] == "met"

    # And the PH3 leg, whose margin is 0.276 V, cannot be established to
    # that precision by the same uncertainty.
    tight = judge_sufficiency(
        {**margin, "required_tolerance": 0.2},
        {
            "uncertainty": _UNCERTAINTY,
            "uncertainty_basis": "measured",
            "uncertainty_evidence_backed": True,
            "display_value": abs(_PH3 - _FERROCENIUM),
        },
    )
    assert tight["state"] == "short"

    # Within tolerance and resting on evidence the host resolved: met,
    # and the requirement is discharged.
    met = judge_sufficiency(
        couple,
        {
            "uncertainty": 0.05,
            "uncertainty_basis": "measured",
            "uncertainty_evidence_backed": True,
        },
    )
    assert met["state"] == "met"

    # The same number on the session's own word is attested and stays
    # open. Not a penalty on assertion: that session declined a measured
    # 0.013 V spread to assert an honest 0.2 V, and was right to.
    attested = judge_sufficiency(
        couple, {"uncertainty": 0.05, "uncertainty_basis": "asserted"}
    )
    assert attested["state"] == "attested"

    # A term the session could not quantify holds it open even when
    # everything else is evidence-backed.
    open_term = judge_sufficiency(
        couple,
        {
            "uncertainty": 0.05,
            "uncertainty_basis": "measured",
            "uncertainty_evidence_backed": True,
            "unquantified_components": ["basin coverage of the cation"],
        },
    )
    assert open_term["state"] == "attested"

    # Silence is not sufficiency.
    assert (
        judge_sufficiency(couple, {"display_value": _PME3})["state"]
        == "unstated"
    )

    # An observable that asked for no precision is owed no answer here.
    assert judge_sufficiency(_declaration(required_tolerance=None), {}) is None
    assert judge_sufficiency(None, {"uncertainty": 1.0}) is None


def test_a_tolerance_states_where_it_came_from(tmp_path):
    host = _host(tmp_path)
    with pytest.raises(ContractError, match="tolerance_basis"):
        host._declare_requested_observable(
            "t1", {"observables": [_declaration(tolerance_basis="")]}
        )
    with pytest.raises(ContractError, match="cannot be negative"):
        host._declare_requested_observable(
            "t1", {"observables": [_declaration(required_tolerance=-1.0)]}
        )
    # A diagnostic is the session's own prediction and is never owed, so
    # no tolerance is owed on it either.
    with pytest.raises(ContractError, match="never owed"):
        host._declare_requested_observable(
            "t1",
            {
                "observables": [
                    _declaration(
                        role="diagnostic",
                        expected_sign="negative",
                        expectation_basis="strong sigma donors",
                        failure_update_rule="re-seed the couple",
                    )
                ]
            },
        )
    ok = host._declare_requested_observable(
        "t1", {"observables": [_declaration()]}
    )
    assert ok
    stored = host.requested_observable_declarations["e-couple-vs-fc"]
    assert stored["required_tolerance"] == _TOLERANCE
    assert stored["tolerance_basis"]


def test_an_uncertainty_says_which_kind_it_is(tmp_path):
    """measured, inferred and asserted all count in full.

    The word is not a discount. It exists because a reader cannot tell a
    spread this run computed from a number recalled out of the
    literature, and OPEN-2 declined a measurement on the strength of an
    unsourced one.
    """

    host = _host(tmp_path)
    with pytest.raises(ContractError, match="uncertainty_basis"):
        host._record_analysis_claims(
            "t1",
            {
                "task_spec_sha256": "a" * 64,
                "claims": [
                    {
                        "claim_id": "e-couple-vs-fc",
                        "receipt_sha256": "b" * 64,
                        "quantity_id": "q",
                        "display_unit": "V",
                        "uncertainty": 0.3,
                    }
                ],
            },
        )


def _claims_with(state_rows):
    """A whole delivery: the claims, the decision, and a green gate.

    The point of these cases is that everything else about the run is
    right -- which is what OPEN-2 was.
    """

    return [
        {
            "kind": "result_quantities_extracted",
            "payload": {"receipt_sha256": "e" * 64},
        },
        {
            "kind": "analysis_claims_recorded",
            "payload": {
                "receipt_sha256": "3" * 64,
                "record": {
                    "claims": [
                        {
                            "claim_id": row["observable_id"],
                            "quantity_id": row["observable_id"],
                            "display_value": 1.0,
                            "display_unit": "kJ/mol",
                            "dimension": (1, 0, 0, 0, 0, 0),
                            "source_receipt_sha256": "4" * 64,
                        }
                        for row in state_rows
                    ]
                },
                "sufficiency": state_rows,
            },
        },
        {"kind": "scientific_decision_recorded", "payload": {}},
        {
            "kind": "analysis_completion_evaluated",
            "payload": {
                "receipt_sha256": "c1" + "c" * 62,
                "status": "passed",
                "limitation_output_ids": [],
            },
        },
    ]


def _short_row(observable_id="dg_solv"):
    return {
        "observable_id": observable_id,
        "unit": "kJ/mol",
        "required_tolerance": 2.0,
        "uncertainty": 6.0,
        "uncertainty_basis": "asserted",
        "state": "short",
    }


def _met_row(observable_id="dg_solv"):
    return {
        "observable_id": observable_id,
        "unit": "kJ/mol",
        "required_tolerance": 2.0,
        "uncertainty": 1.0,
        "uncertainty_basis": "measured",
        "state": "met",
    }


def test_the_projection_names_only_the_open_ones():
    delivery = _AnalysisDelivery(
        completion_status="passed",
        limitation_output_ids=(),
        claims=1,
        decisions=1,
        receipt_sha256s=(),
        sufficiency=(
            _short_row("a"),
            _met_row("b"),
            # A state this vocabulary does not know. Rows are read back
            # out of a durable workspace record no version stamps, and
            # this branch itself wrote "separated" two commits ago, so
            # the projection meets one. Open is the only safe reading:
            # the fail-open half discharges a contract nothing assessed.
            {"observable_id": "c", "state": "separated"},
            {"observable_id": "d", "state": "unstated"},
        ),
    )
    assert delivery.unresolved_requirement_ids == ("a", "c", "d")


def _declared(observable_id="dg_solv"):
    return [
        {
            "kind": "requested_observable_declared",
            "payload": {
                "observables": [
                    {
                        "observable_id": observable_id,
                        "unit": "kJ/mol",
                        "dimension": (1, 0, 0, 0, 0, 0),
                        "meaning": observable_id,
                        "required_tolerance": 2.0,
                        "tolerance_basis": "the author's design threshold",
                    }
                ]
            },
        }
    ]


def _ledger(tmp_path):
    path = (
        tmp_path / "ws" / ".chemsmart-agent" / "goals" / "goal-t1"
    ) / "ledger.jsonl"
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
    ]


def test_a_number_short_of_its_tolerance_opens_one_more_cycle(tmp_path):
    """The lever is control flow, not another sentence.

    A host notice would traverse the same branch: `_rewake` returned
    before it ever read the budget, so nothing that only *told* the
    session could reach it.
    """

    contexts = []

    def capture(inner):
        def step(workspace, kwargs):
            contexts.append(kwargs["goal_context"])
            return inner(workspace, kwargs)

        return step

    result = _loop(
        tmp_path,
        sessions=[
            capture(
                _planning_session(
                    "live-1",
                    terminal="planned",
                    wake_rows=_declared() + _claims_with([_short_row()]),
                )
            ),
            capture(
                _planning_session(
                    "live-2",
                    terminal="planned",
                    wake_rows=_declared() + _claims_with([_met_row()]),
                )
            ),
            capture(
                _planning_session(
                    "live-3",
                    terminal="planned",
                    wake_rows=_declared() + _claims_with([_met_row()]),
                )
            ),
        ],
        executes=[],
        max_revisions=3,
    )
    report = contexts[1]["failure_report"]
    assert report["gate"] == "goal.requirement_is_resolved"
    assert "dg_solv" in report["diagnosis"]
    assert report["route"] == SUFFICIENCY_SHORT_ROUTE
    # All three routes are named, and two of them cost nothing.
    assert "declare the margin" in report["route"]
    assert "unreachable_observable_ids" in report["route"]
    assert "no engine call" in report["cost"]
    kinds = [entry["kind"] for entry in _ledger(tmp_path)]
    assert kinds.count("rewake_opened") == 1
    # Resolved on the second pass, so the goal is achieved and the
    # numbers were never withheld.
    assert result.settlement == "achieved"


def test_a_requirement_left_open_returns_to_the_human(tmp_path):
    """No fifth settlement word, and no false `achieved` either."""

    result = _loop(
        tmp_path,
        sessions=[
            _planning_session(
                "live-1",
                terminal="planned",
                wake_rows=_declared() + _claims_with([_short_row()]),
            ),
            _planning_session(
                "live-2",
                terminal="planned",
                wake_rows=_declared() + _claims_with([_short_row()]),
            ),
        ],
        executes=[],
        max_revisions=3,
    )
    assert result.settlement == "returned_to_human"
    # The invariant is that the settlement names the requirement it is
    # handing back and what that requirement stands at -- not the
    # sentence it uses. The wording changed when "the precision the task
    # asked for" turned out to be false of a tolerance the session
    # declared for a decision question of its own.
    joined = " | ".join(result.reasons)
    assert "dg_solv" in joined, result.reasons
    assert "stands short" in joined, result.reasons


def test_a_goal_that_asks_no_tolerance_is_untouched(tmp_path):
    result = _loop(
        tmp_path,
        sessions=[
            _planning_session(
                "live-1",
                terminal="planned",
                wake_rows=[
                    {
                        "kind": "requested_observable_declared",
                        "payload": {
                            "observables": [
                                {
                                    "observable_id": "dg_solv",
                                    "unit": "kJ/mol",
                                    "dimension": (1, 0, 0, 0, 0, 0),
                                    "meaning": "dg_solv",
                                }
                            ]
                        },
                    }
                ]
                + _claims_with([_met_row()]),
            ),
            _planning_session(
                "live-2",
                terminal="planned",
                wake_rows=_claims_with([_met_row()]),
            ),
        ],
        executes=[],
        max_revisions=3,
    )
    assert result.settlement == "achieved"


def test_a_tolerance_and_an_uncertainty_are_compared_in_one_unit(tmp_path):
    """kJ/mol and kcal/mol are both molar energies and not the same number.

    The declaration states its tolerance in the observable's own unit;
    a claim states its uncertainty in the claim's display unit; nothing
    converted. A 2 kJ/mol tolerance certified a 1 kcal/mol uncertainty
    as `met` -- 4.184 kJ/mol against a 2 kJ/mol requirement -- and
    printed it back as "1.0 kJ/mol", so the row was wrong twice.
    """

    host = _host(tmp_path)
    host._declare_requested_observable(
        "t1",
        {
            "observables": [
                {
                    "observable_id": "dg-solv",
                    "unit": "kJ/mol",
                    "meaning": "solvation free energy",
                    "required_tolerance": 2.0,
                    "tolerance_basis": "the author's design threshold",
                }
            ]
        },
    )
    declaration = host.requested_observable_declarations["dg-solv"]

    # Restated into the declared unit, 1 kcal/mol exceeds 2 kJ/mol.
    row = host._restated_sufficiency_row(
        declaration,
        display_unit="kcal/mol",
        display_value=-10.0,
        uncertainty=1.0,
        uncertainty_basis="measured",
    )
    assert row["uncertainty"] == pytest.approx(4.184, rel=1e-3)
    assert judge_sufficiency(declaration, row)["state"] == "short"

    # In the declared unit nothing is touched.
    same = host._restated_sufficiency_row(
        declaration,
        display_unit="kJ/mol",
        display_value=-10.0,
        uncertainty=1.0,
        uncertainty_basis="measured",
    )
    assert same["uncertainty"] == 1.0
    same["uncertainty_evidence_backed"] = True
    assert judge_sufficiency(declaration, same)["state"] == "met"

    # A unit of another dimension says nothing about this tolerance, and
    # the host refuses rather than comparing.
    with pytest.raises(ContractError, match="no conversion"):
        host._restated_sufficiency_row(
            declaration,
            display_unit="angstrom",
            display_value=1.0,
            uncertainty=1.0,
            uncertainty_basis="measured",
        )


def _seeded_host(tmp_path, tolerance=2.0):
    """A host holding one declaration and one extraction to claim from."""

    from types import SimpleNamespace

    host = _host(tmp_path)
    host._declare_requested_observable(
        "t1",
        {
            "observables": [
                {
                    "observable_id": "dg",
                    "unit": "kJ/mol",
                    "meaning": "solvation free energy",
                    "required_tolerance": tolerance,
                    "tolerance_basis": "the author's design threshold",
                }
            ]
        },
    )
    host.quantity_extractions["b" * 64] = SimpleNamespace(
        quantities=(
            SimpleNamespace(
                quantity_id="q1",
                value=-10.0,
                unit="kJ/mol",
                dimension=(1, 0, 0, 0, 0, 0),
                data_kind="scalar",
                value_sha256="9" * 64,
            ),
            # The narrowed spread the re-claim cites: "measured" means
            # the host reads the number where the session got it.
            SimpleNamespace(
                quantity_id="q1_spread",
                # Canonical, as a real receipt carries it: the host
                # converts to the claim's display unit exactly as it
                # does for the value being claimed.
                value=1.0 / _KJ_PER_MOL_IN_HARTREE,
                unit="hartree",
                dimension=(1, 0, 0, 0, 0, 0),
                data_kind="scalar",
                value_sha256="7" * 64,
            ),
        )
    )
    return host


def test_the_menus_own_route_can_be_walked(tmp_path):
    """Re-claiming with a corrected uncertainty was refused outright.

    The uncertainty lived beside the claim rather than inside it, so a
    corrected one produced the same record digest, the same idempotency
    key and a different payload -- and the event store refused it with
    "idempotency key conflicts with persisted action". The sufficiency
    wake offers "claim the observable again with the uncertainty you
    measured" as one of its three routes, and that route could not be
    walked at all. Third strike for a named route the host cannot walk.
    """

    from chemsmart.agent.delivery import unresolved_requirement_ids

    host = _seeded_host(tmp_path)

    def claim(uncertainty, **extra):
        return host._record_analysis_claims(
            "t1",
            {
                "task_spec_sha256": "a" * 64,
                "claims": [
                    {
                        "claim_id": "dg",
                        "receipt_sha256": "b" * 64,
                        "quantity_id": "q1",
                        "display_unit": "kJ/mol",
                        "uncertainty": uncertainty,
                        "uncertainty_basis": "asserted",
                        **extra,
                    }
                ],
            },
        )

    first = claim(6.0)
    # The re-claim cites the receipt its narrowed uncertainty came from,
    # which is what turns the route into a discharge rather than another
    # assertion.
    second = claim(
        1.0,
        uncertainty_basis="measured",
        uncertainty_reference="b" * 64 + ":q1_spread",
    )
    # A number and the uncertainty its author gives it are one
    # statement, so a corrected uncertainty is a different claim.
    assert first.receipt_sha256 != second.receipt_sha256
    assert first.claims[0].uncertainty == 6.0
    assert second.claims[0].uncertainty == 1.0
    assert second.claims[0].uncertainty_basis == "measured"

    rows = [
        row
        for event in _events(tmp_path)
        if event["kind"] == "analysis_claims_recorded"
        for row in event["payload"].get("sufficiency") or ()
    ]
    assert [row["state"] for row in rows] == ["short", "met"]
    assert rows[1]["uncertainty_evidence_backed"] is True
    # And with the current assessment winning, the wake it opened closes.
    assert unresolved_requirement_ids(rows) == ()


def _events(tmp_path):
    return [
        json.loads(line)
        for line in (tmp_path / "events.jsonl")
        .read_text(encoding="utf-8")
        .splitlines()
        if line.strip()
    ]


def test_an_executed_claim_is_assessed_after_the_numbers_exist(tmp_path):
    """The executor's own claim shape, and the join it needs.

    A planned claim node is written while the calculation is still a
    plan, so it carries no uncertainty and cannot: an uncertainty is a
    judgement about a result. The executor therefore supplies four
    fields -- claim_id, receipt, quantity_id, display_unit -- and the
    assessment lands `unstated` by construction. That is the design
    (owner ruling 2026-09-09): the cost is one further cycle and no
    engine call, and the alternative is a prediction wearing an
    assessment's clothes.

    It also needs the join the completion gate uses. A plan carries the
    host-minted quantity_id and its own short input label in claim_id;
    reading claim_id alone here, while the gate reads both, meant a
    legitimate alternate label delivered the observable and produced no
    assessment of it at all.
    """

    host = _seeded_host(tmp_path)
    # The plan's own label differs from the declared observable's id,
    # which is ordinary: the declaration is scientific, the node input
    # is structural.
    host._record_analysis_claims(
        "t1",
        {
            "task_spec_sha256": "a" * 64,
            "claims": [
                {
                    "claim_id": "dg",
                    "receipt_sha256": "b" * 64,
                    "quantity_id": "q1",
                    "display_unit": "kJ/mol",
                }
            ],
        },
    )
    rows = [
        row
        for event in _events(tmp_path)
        if event["kind"] == "analysis_claims_recorded"
        for row in event["payload"].get("sufficiency") or ()
    ]
    assert [row["state"] for row in rows] == ["unstated"]
    assert rows[0]["required_tolerance"] == 2.0


def test_the_assessment_joins_on_the_quantity_id_too(tmp_path):
    """The declared id arriving in quantity_id still gets assessed."""

    from types import SimpleNamespace

    host = _host(tmp_path)
    host._declare_requested_observable(
        "t1",
        {
            "observables": [
                {
                    "observable_id": "dg",
                    "unit": "kJ/mol",
                    "meaning": "solvation free energy",
                    "required_tolerance": 2.0,
                    "tolerance_basis": "the author's design threshold",
                }
            ]
        },
    )
    host.quantity_extractions["b" * 64] = SimpleNamespace(
        quantities=(
            SimpleNamespace(
                quantity_id="dg",
                value=-10.0,
                unit="kJ/mol",
                dimension=(1, 0, 0, 0, 0, 0),
                data_kind="scalar",
                value_sha256="9" * 64,
            ),
        )
    )
    host._record_analysis_claims(
        "t1",
        {
            "task_spec_sha256": "a" * 64,
            "claims": [
                {
                    "claim_id": "expr-out-3",
                    "receipt_sha256": "b" * 64,
                    "quantity_id": "dg",
                    "display_unit": "kJ/mol",
                }
            ],
        },
    )
    rows = [
        row
        for event in _events(tmp_path)
        if event["kind"] == "analysis_claims_recorded"
        for row in event["payload"].get("sufficiency") or ()
    ]
    assert [row["observable_id"] for row in rows] == ["dg"]
    assert rows[0]["state"] == "unstated"


def test_a_correction_keeps_the_obligation_it_corrects(tmp_path):
    """Supersession corrects the observable; it never retires the task's
    requirement.

    A replacement needed only to name an existing id, so a
    tolerance-free declaration -- or a diagnostic, which is never owed
    at all -- could make a requested precision requirement vanish by
    relabelling. That is the cheapest possible escape from a contract
    and it was open.

    What must stay possible is the correction itself: a wrong unit is
    repaired by declaring the right one, and this must not obstruct it.
    """

    host = _host(tmp_path)
    first = {
        "observable_id": "gap-ev",
        "unit": "eV",
        "meaning": "the singlet-triplet gap",
        "required_tolerance": 0.05,
        "tolerance_basis": "the author asked for 0.05 eV",
    }
    host._declare_requested_observable("t1", {"observables": [first]})

    # A diagnostic is never owed, so it cannot absorb an obligation.
    with pytest.raises(ContractError, match="never owed"):
        host._declare_requested_observable(
            "t1",
            {
                "observables": [
                    {
                        "observable_id": "gap-guess",
                        "unit": "eV",
                        "meaning": "a guess at the gap",
                        "role": "diagnostic",
                        "expected_sign": "positive",
                        "expectation_basis": "a strong ligand field",
                        "failure_update_rule": "re-seed the triplet",
                        "supersedes_observable_id": "gap-ev",
                    }
                ]
            },
        )

    # Neither can a replacement that simply states no tolerance.
    with pytest.raises(ContractError, match="states none"):
        host._declare_requested_observable(
            "t1",
            {
                "observables": [
                    {
                        "observable_id": "gap-kcal",
                        "unit": "kcal/mol",
                        "meaning": "the singlet-triplet gap, corrected unit",
                        "supersedes_observable_id": "gap-ev",
                    }
                ]
            },
        )

    # The correction itself is ordinary, and the tolerance is restated
    # in the corrected unit rather than converted by the host, because a
    # supersession is often exactly a unit correction.
    host._declare_requested_observable(
        "t1",
        {
            "observables": [
                {
                    "observable_id": "gap-kcal",
                    "unit": "kcal/mol",
                    "meaning": "the singlet-triplet gap, corrected unit",
                    "required_tolerance": 1.15,
                    "tolerance_basis": "0.05 eV restated in kcal/mol",
                    "supersedes_observable_id": "gap-ev",
                }
            ]
        },
    )
    assert (
        host.requested_observable_declarations["gap-kcal"][
            "required_tolerance"
        ]
        == 1.15
    )


def test_a_retired_requirement_stops_holding_the_goal_open(tmp_path):
    """The mirror: an assessment of a replaced observable is history."""

    from chemsmart.agent.driver import _AnalysisDelivery

    delivery = _AnalysisDelivery(
        completion_status="passed",
        limitation_output_ids=(),
        claims=1,
        decisions=1,
        receipt_sha256s=(),
        declared_observables=(
            {"observable_id": "gap-ev", "unit": "eV"},
            {
                "observable_id": "gap-kcal",
                "unit": "kcal/mol",
                "supersedes_observable_id": "gap-ev",
            },
        ),
        sufficiency=(
            {"observable_id": "gap-ev", "state": "short"},
            {"observable_id": "gap-kcal", "state": "met"},
        ),
    )
    assert delivery.unresolved_requirement_ids == ()


def _refusal_rows(observable_id="dg_solv", *, verified=True):
    """A decision that refuses the requirement, not the observable."""

    return [
        {
            "kind": "scientific_decision_recorded",
            "payload": {
                "receipt_sha256": "d" * 64,
                "record": {
                    "evidence_refs": ["receipt:" + "e" * 64],
                    "uncertainties": [],
                },
                "unreachable_observables": [
                    {
                        "observable_id": observable_id,
                        "statement": "no method in this envelope reaches "
                        "2 kJ/mol",
                        "selector": "",
                        "jobtype": "",
                        "blocked_node_id": "",
                        "receipt_sha256s": ["e" * 64],
                        "verified": verified,
                        "basis": "the requirement stands open",
                    }
                ],
            },
        }
    ]


def test_a_refused_requirement_is_a_deliverable(tmp_path):
    """Route three of the menu reaches the word the charter promises.

    The projections subtract a host-verified refusal, which is right --
    a refused requirement is not owed -- but nothing read what they
    removed. The observable was claimed, so it appears in no
    undelivered list, and the goal settled `achieved`: the refusal
    bought the best word in the vocabulary while the tool's own reply
    promised `unreachable_from_evidence` for it.
    """

    rows = _declared() + _claims_with([_short_row()]) + _refusal_rows()
    result = _loop(
        tmp_path,
        sessions=[
            _planning_session("live-1", terminal="planned", wake_rows=rows)
        ],
        executes=[],
        max_revisions=3,
    )
    assert result.settlement == "unreachable_from_evidence"
    assert any("dg_solv" in reason for reason in result.reasons)


def test_an_unverified_refusal_of_a_requirement_returns_to_the_human(
    tmp_path,
):
    rows = (
        _declared()
        + _claims_with([_short_row()])
        + _refusal_rows(verified=False)
    )
    result = _loop(
        tmp_path,
        sessions=[
            _planning_session("live-1", terminal="planned", wake_rows=rows),
            _planning_session("live-2", terminal="planned", wake_rows=rows),
        ],
        executes=[],
        max_revisions=3,
    )
    assert result.settlement == "returned_to_human"


def test_an_unverified_refusal_does_not_close_the_runs_requirement(
    tmp_path,
):
    """An explanation string never confers authority.

    The cross-stream inheritance passed the *bases* mapping, which
    explains every refusal the session wrote -- including the ones the
    host declined to verify -- and the receiving projection promoted
    every key it received. An unverified refusal arrived carrying
    authority it had been denied, removed the run's open requirement,
    and left a settlement whose refusal branch could not fire: the
    previous audit's own pattern, inside its repair.
    """

    from chemsmart.agent.driver import _analysis_delivery

    session_events = tmp_path / "session.jsonl"
    session_events.write_text(
        "\n".join(json.dumps(row) for row in _refusal_rows(verified=False))
        + "\n",
        encoding="utf-8",
    )
    run_events = tmp_path / "run.jsonl"
    run_events.write_text(
        "\n".join(json.dumps(row) for row in _claims_with([_short_row()]))
        + "\n",
        encoding="utf-8",
    )
    session = _analysis_delivery(session_events)
    assert session.verified_unreachable_ids == ()
    assert session.unverified_unreachable_ids == ("dg_solv",)

    inherited = {
        observable_id: session.unreachable_bases.get(observable_id, "")
        for observable_id in session.verified_unreachable_ids
    }
    joined = _analysis_delivery(run_events, inherited_unreachable=inherited)
    assert joined.unresolved_requirement_ids == ("dg_solv",)
    assert joined.refused_requirement_ids == ()
