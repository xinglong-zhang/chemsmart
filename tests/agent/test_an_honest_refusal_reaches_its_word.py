"""A typed refusal the host can verify settles unreachable_from_evidence
from either delivery route; one it cannot verify returns to the human
naming it; a session cannot refuse its way out.

po3 said "the goal settles as a typed refusal" in prose and recorded
nothing typed; ino3 built five blocked nodes and the settlement joined
blocked-ness to required output ids, not observable ids (NOVEL-3,
2026-09-05). Owner ruling 2026-09-06: host-verified refusal.
"""

from __future__ import annotations

import json

import pytest

from chemsmart.agent._contracts import RoutedContractError
from chemsmart.agent.driver import _analysis_delivery, _settle_from_delivery
from chemsmart.agent.goal import GoalLedger, GoalRecordV1

from .test_a_guide_opens_when_something_asks import _host

pytestmark = pytest.mark.capability("tool:record_scientific_decision")

_SPIN = {
    "observable_id": "spin-nickel",
    "unit": "1",
    "dimension": (0, 0, 0, 0, 0, 0),
    "meaning": "spin population on nickel in the cation",
}
_GAP = {
    "observable_id": "quartet-doublet-gap",
    "unit": "eV",
    "dimension": (1, 0, 0, 0, 0, 0),
    "meaning": "quartet minus doublet",
}


def _decision(host, receipt, **entry):
    base = {
        "decision_id": "d-" + entry.get("observable_id", "x"),
        "assumptions": ["a"],
        "method_rationale": "r",
        "alternatives": ["b"],
        "uncertainties": ["u"],
        "diagnostics": ["g"],
        "stage_order": ["s"],
        "evidence_refs": [],
        "unreachable_observable_ids": [
            {
                "statement": "no atomic spin-population selector",
                "receipt_sha256s": [receipt],
                **entry,
            }
        ],
    }
    return host.dispatch(
        turn_id="t-" + entry.get("observable_id", "x"),
        tool_name="record_scientific_decision",
        arguments=base,
    )


def _extraction_receipt(host):
    from pathlib import Path

    from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256

    resolved = Path("tests/data/ORCATests/outputs/CO2.out").resolve()
    host.artifacts["result.co2"] = TrustedArtifactRefV1(
        artifact_id="result.co2",
        kind="orca_output",
        sha256=file_sha256(resolved),
        size_bytes=resolved.stat().st_size,
        path=str(resolved),
        cli_value=str(resolved),
    )
    reply = host.dispatch(
        turn_id="probe",
        tool_name="extract_result_quantities",
        arguments={
            "artifact_id": "result.co2",
            "program": "orca",
            "selectors": [{"quantity_id": "e", "selector": "energy"}],
        },
    )
    return reply["result"]["receipt_sha256"]


def test_a_selector_no_program_declares_verifies_the_refusal(tmp_path):
    host = _host(
        tmp_path, approved_requested_observable_declarations=[_SPIN, _GAP]
    )
    receipt = _extraction_receipt(host)
    reply = _decision(
        host,
        receipt,
        observable_id="spin-nickel",
        selector="hirshfeld_atomic_spin_populations",
        jobtype="opt",
    )
    (entry,) = reply["result"]["unreachable_observables"]
    assert entry["verified"] is True
    assert "declares selector" in entry["basis"]
    assert (
        "unreachable_from_evidence"
        in reply["result"]["settlement_consequence"]
    )
    events = [
        e
        for e in host.event_store.read_events()
        if e.kind == "scientific_decision_recorded"
    ]
    assert events[-1].payload["unreachable_observables"][0]["verified"]


def test_a_declared_selector_leaves_the_refusal_unverified(tmp_path):
    host = _host(
        tmp_path, approved_requested_observable_declarations=[_SPIN, _GAP]
    )
    receipt = _extraction_receipt(host)
    reply = _decision(
        host,
        receipt,
        observable_id="quartet-doublet-gap",
        selector="energy",
        jobtype="opt",
    )
    (entry,) = reply["result"]["unreachable_observables"]
    assert entry["verified"] is False
    assert "is declared by" in entry["basis"] and "orca/opt" in entry["basis"]


def test_an_undeclared_id_or_a_missing_receipt_is_refused(tmp_path):
    host = _host(tmp_path, approved_requested_observable_declarations=[_SPIN])
    receipt = _extraction_receipt(host)
    with pytest.raises(RoutedContractError) as refused:
        _decision(host, receipt, observable_id="not-declared")
    assert refused.value.failure_report["gate"] == (
        "decision.unreachable_id_is_declared"
    )
    with pytest.raises(RoutedContractError) as refused:
        _decision(host, "9" * 64, observable_id="spin-nickel")
    assert refused.value.failure_report["gate"] == (
        "decision.receipt_is_one_the_host_minted"
    )


def _stream(path, *, unreachable, completion, verified=True, claims=()):
    rows = []
    if claims:
        rows.append(
            {
                "kind": "analysis_claims_recorded",
                "payload": {
                    "receipt_sha256": "5" * 64,
                    "record": {
                        "claims": [
                            {
                                "claim_id": c,
                                "quantity_id": c,
                                "display_value": 0.35,
                                "display_unit": "eV",
                                "source_receipt_sha256": "6" * 64,
                            }
                            for c in claims
                        ]
                    },
                },
            }
        )
    if completion:
        rows.append(
            {
                "kind": "analysis_completion_evaluated",
                "payload": {
                    "receipt_sha256": "7" * 64,
                    "status": "passed",
                    "limitation_output_ids": [
                        f"declared_observable:{i}" for i in unreachable
                    ],
                    "declared_observable_misses": [],
                },
            }
        )
    rows.append(
        {
            "kind": "scientific_decision_recorded",
            "payload": {
                "receipt_sha256": "8" * 64,
                "record": {"decision_id": "d", "evidence_refs": []},
                "unreachable_observables": [
                    {
                        "observable_id": i,
                        "statement": "no spin selector",
                        "verified": verified,
                        "basis": "no program declares it",
                        "receipt_sha256s": ["9" * 64],
                    }
                    for i in unreachable
                ],
            },
        }
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(json.dumps(r) for r in rows) + "\n", encoding="utf-8"
    )
    return path


def _ledger(workspace):
    ledger = GoalLedger(workspace / ".chemsmart-agent" / "goals" / "goal-r")
    ledger.create(
        GoalRecordV1(
            schema_version="chemsmart.goal.v1",
            goal_id="goal-r",
            task_spec_sha256="a" * 64,
            scientific_identity_sha256="",
            conditions={"solvents": (), "thermochemistry": ()},
            envelope={
                "allowed_program_engines": (),
                "max_engine_calls": 30,
                "episode_wall_time_seconds": 21600.0,
                "max_excursion_calls": 0,
            },
            max_revisions=8,
            granted_by="claude-owner-delegated-reviewer",
            initial_review_sha256="",
            created_at="2026-09-05T00:00:00+00:00",
        )
    )
    return ledger


def test_a_verified_refusal_settles_the_word_from_a_certified_chain(
    tmp_path,
):
    stream = _stream(
        tmp_path / "a" / "events.jsonl",
        unreachable=("spin-nickel",),
        completion=True,
        claims=("quartet-doublet-gap",),
    )
    delivery = _analysis_delivery(stream)
    assert delivery.verified_unreachable_ids == ("spin-nickel",)
    result = _settle_from_delivery(
        _ledger(tmp_path / "a"),
        goal_id="goal-r",
        cycles=2,
        revisions_admitted=1,
        events_path=stream,
        terminal="complete",
        workspace=tmp_path / "a",
    )
    assert result.settlement == "unreachable_from_evidence"
    assert "spin-nickel -- no spin selector [no program declares it]" in (
        result.reasons[0]
    )


def test_an_unverified_refusal_returns_naming_it(tmp_path):
    stream = _stream(
        tmp_path / "b" / "events.jsonl",
        unreachable=("spin-nickel",),
        completion=True,
        verified=False,
        claims=("quartet-doublet-gap",),
    )
    result = _settle_from_delivery(
        _ledger(tmp_path / "b"),
        goal_id="goal-r",
        cycles=2,
        revisions_admitted=1,
        events_path=stream,
        terminal="complete",
        workspace=tmp_path / "b",
    )
    assert result.settlement == "returned_to_human"
    assert "could not verify" in result.reasons[0]
    assert "spin-nickel" in result.reasons[0]


def test_a_refusal_without_a_completion_still_reaches_the_word(tmp_path):
    """po3's route: claims made in-session, no plan, no completion. The
    settlement keys the refusal on the goal's declared ids, not on a
    completion's limitation list."""

    workspace = tmp_path / "c"
    ledger = _ledger(workspace)
    ledger.append(
        "observables_declared",
        {
            "cycle": 1,
            "observables": [
                {
                    "observable_id": "ddg-activation-353k",
                    "unit": "kcal/mol",
                    "dimension": (1, 0, 0, 0, 0, 0),
                    "meaning": "ddG",
                }
            ],
        },
    )
    stream = _stream(
        workspace / "events.jsonl",
        unreachable=("ddg-activation-353k",),
        completion=False,
    )
    result = _settle_from_delivery(
        ledger,
        goal_id="goal-r",
        cycles=2,
        revisions_admitted=0,
        events_path=stream,
        terminal="blocked",
        workspace=workspace,
    )
    assert result.settlement == "unreachable_from_evidence"
    assert "ddg-activation-353k" in result.reasons[0]
