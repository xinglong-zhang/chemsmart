"""Seventy claims and a decision ended `blocked` for want of a ceremony.

A session may answer a task from results already in the workspace:
extract, derive, evaluate, claim, decide. The charter names that route
and SUFFICIENCY-1 took it. But finalisation knew only two ways to mint
a completion certificate -- a task-owned analysis policy, or a
preflighted workflow -- so a session with neither ended `blocked`
however much it had delivered, and the goal returned to the human
saying the completion gate had not passed. That is the proximate reason
that window produced no settlement anyone could read.

The certificate is minted from what such a session actually has: the
requirements it declared and the receipts its claims stand on. No plan
is invented and no further model act is asked for.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

pytestmark = pytest.mark.capability("tool:record_analysis_claims")

_TASK = "a" * 64


def _host(tmp_path):
    workspace = tmp_path / "workspace"
    workspace.mkdir(exist_ok=True)
    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="s1"
        ),
        task_spec_sha256s=(_TASK,),
        approved_workspace=workspace,
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
    return host


def _claim(host):
    return host._record_analysis_claims(
        "t1",
        {
            "task_spec_sha256": _TASK,
            "claims": [
                {
                    "claim_id": "dg",
                    "receipt_sha256": "b" * 64,
                    "quantity_id": "dg",
                    "display_unit": "kJ/mol",
                    "uncertainty": 1.0,
                    "uncertainty_basis": "asserted",
                }
            ],
        },
    )


def test_nothing_delivered_is_nothing_to_certify(tmp_path):
    host = _host(tmp_path)
    with pytest.raises(ContractError, match="no analysis claim"):
        host.completion_receipts_for_delivered_claims()


def test_claims_without_a_decision_are_not_a_finished_delivery(tmp_path):
    """The decision is where the science is said; a bare number is not
    a delivery."""

    host = _host(tmp_path)
    _claim(host)
    with pytest.raises(ContractError, match="no scientific decision"):
        host.completion_receipts_for_delivered_claims()


def test_claims_and_a_decision_are_certified_from_the_declarations(tmp_path):
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
    claim_record = _claim(host)
    host._record_scientific_decision(
        "t1",
        {
            "decision_id": "d1",
            "task_spec_sha256": _TASK,
            "method_rationale": "read from results already registered",
            "uncertainties": ["the basis set is the dominant term"],
            "assumptions": [],
            "alternatives": [],
            "diagnostics": [],
            "stage_order": ["extract", "claim", "decide"],
            "evidence_refs": [],
        },
    )

    (digest,) = host.completion_receipts_for_delivered_claims()
    completion = host.analysis_completion_receipts[digest]
    assert completion.status == "passed"
    # The receipts the claims actually stand on, and nothing else.
    assert completion.source_receipt_sha256s == ("b" * 64,)
    assert completion.task_spec_sha256 == _TASK
    # The policy identity is the contract this delivery answers, so two
    # deliveries against the same declarations share it.
    again = host.completion_receipts_for_delivered_claims()
    assert again == (digest,)
    assert claim_record.receipt_sha256
