"""A number standing on an optimisation that printed no frequencies says
its stationary point is uncharacterised.

NOVEL-1's iron runs were opt-only; the trans triplet delivered as an
ordered gap was, under NOVEL-2's Freq on the bit-identical geometry, a
second-order saddle. The owner ruled (2026-09-05) that such a claim
carries a provenance word and is never refused.
"""

from __future__ import annotations

import json

import pytest

from chemsmart.agent.driver import _achieved_word, _analysis_delivery

pytestmark = pytest.mark.capability("rule:stationary_point_order")


def _stream(tmp_path, *, mode_count, jobtype="opt"):
    rows = [
        {
            "kind": "program_result_verified",
            "payload": {
                "record": {
                    "node_id": "trans-t",
                    "state": "valid",
                    "jobtype": jobtype,
                    "observations": {
                        "jobtype": jobtype,
                        "orca": {
                            "vibrational_mode_count": mode_count,
                            "optimization_converged": True,
                        },
                    },
                    "output_artifacts": [{"sha256": "f" * 64}],
                }
            },
        },
        {
            "kind": "result_quantities_extracted",
            "payload": {
                "receipt_sha256": "1" * 64,
                "artifact_sha256": "f" * 64,
                "record": {},
            },
        },
        {
            "kind": "analysis_claims_recorded",
            "payload": {
                "receipt_sha256": "2" * 64,
                "record": {
                    "claims": [
                        {
                            "claim_id": "trans_gap_kjmol",
                            "quantity_id": "trans_gap_kjmol",
                            "source_receipt_sha256": "1" * 64,
                        }
                    ]
                },
            },
        },
    ]
    tmp_path.mkdir(parents=True, exist_ok=True)
    stream = tmp_path / "events.jsonl"
    stream.write_text(
        "\n".join(json.dumps(row) for row in rows) + "\n", encoding="utf-8"
    )
    return stream


def test_a_claim_on_an_opt_without_frequencies_is_uncharacterised(tmp_path):
    delivery = _analysis_delivery(_stream(tmp_path, mode_count=0))
    assert delivery.uncharacterised_source_quantity_ids == ("trans_gap_kjmol",)
    assert delivery.failed_source_quantity_ids == ()
    word, reasons = _achieved_word(delivery, ())
    assert word == "achieved"
    assert any(
        "uncharacterised (no frequencies printed): trans_gap_kjmol" in r
        for r in reasons
    )


def test_a_claim_on_an_opt_with_frequencies_is_not(tmp_path):
    delivery = _analysis_delivery(_stream(tmp_path / "a", mode_count=57))
    assert delivery.uncharacterised_source_quantity_ids == ()
    single_point = _analysis_delivery(
        _stream(tmp_path / "b", mode_count=0, jobtype="sp")
    )
    assert single_point.uncharacterised_source_quantity_ids == ()


def test_the_word_crosses_cycles_through_the_workspace_record(tmp_path):
    """A result verified in cycle 1 and claimed in cycle 3 carries no
    verified record in cycle 3's stream; the workspace record does
    (NOVEL-1 ino1: the trans gap claimed in cycle 3 stood on an
    opt-only result from cycle 1)."""

    from chemsmart.agent.workspace_record import (
        record_run,
        uncharacterised_artifacts,
    )

    workspace = tmp_path / "ws"
    first = _stream(tmp_path / "c1", mode_count=0)
    record_run(
        workspace,
        goal_id="g",
        cycle=1,
        run_events_path=first,
        run="goals/g/runs/cycle-1",
    )
    assert uncharacterised_artifacts(workspace) == ("f" * 64,)
    later = tmp_path / "c3" / "events.jsonl"
    later.parent.mkdir(parents=True)
    rows = [
        json.loads(line)
        for line in first.read_text(encoding="utf-8").splitlines()
        if json.loads(line)["kind"] != "program_result_verified"
    ]
    later.write_text(
        "\n".join(json.dumps(row) for row in rows) + "\n", encoding="utf-8"
    )
    assert _analysis_delivery(later).uncharacterised_source_quantity_ids == ()
    delivery = _analysis_delivery(
        later,
        uncharacterised_artifact_sha256s=uncharacterised_artifacts(workspace),
    )
    assert delivery.uncharacterised_source_quantity_ids == ("trans_gap_kjmol",)
