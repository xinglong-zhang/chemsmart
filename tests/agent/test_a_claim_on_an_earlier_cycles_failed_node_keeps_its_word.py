"""A number claimed on a failed node keeps its word across cycles.

A settlement says, beside a delivered number, that it stands on a node
that did not meet its promise, or on one the session had the host
characterise. Those words were computed from the failed records in the
stream being settled. A goal that types a saddle in cycle one and
claims its characterised numbers from cycle two's session settled with
neither word (PySCF round g5, 2026-09-12: six claims on a characterised
inversion saddle, the settlement silent about all of them). The
workspace record holds the failed results, and both settle-time
deliveries now read them.
"""

from __future__ import annotations

import json

import pytest

from chemsmart.agent.driver import _settle_from_delivery
from chemsmart.agent.goal import GoalLedger
from chemsmart.agent.workspace_record import record_run

pytestmark = pytest.mark.capability("tool:characterise_stationary_point")

SADDLE = "f" * 64
EXTRACTION = "c" * 64


def _write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(json.dumps(row) for row in rows) + "\n", encoding="utf-8"
    )
    return path


def test_the_settlement_names_a_characterised_saddle_from_an_earlier_cycle(
    tmp_path,
):
    workspace = tmp_path / "ws"
    record_run(
        workspace,
        goal_id="g",
        cycle=1,
        run_events_path=_write(
            tmp_path / "cycle-1.jsonl",
            [
                {
                    "kind": "program_result_verified",
                    "payload": {
                        "node_id": "planar_hess",
                        "status": "invalid",
                        "receipt_sha256": "1" * 64,
                        "record": {
                            "node_id": "planar_hess",
                            "program": "pyscf",
                            "jobtype": "hess",
                            "state": "invalid",
                            "findings": ["result.stationary_point_order"],
                            "observations": {
                                "jobtype": "hess",
                                "program": "pyscf",
                                "pyscf": {"vibrational_mode_count": 6},
                            },
                            "output_artifacts": [
                                {"kind": "pyscf_hdf5", "sha256": SADDLE}
                            ],
                        },
                    },
                }
            ],
        ),
        run="goals/g/runs/cycle-1",
    )
    session = _write(
        tmp_path / "cycle-2.jsonl",
        [
            {
                "kind": "stationary_point_characterised",
                "payload": {
                    "node_id": "planar_hess",
                    "receipt_sha256": "2" * 64,
                    "record": {
                        "node_id": "planar_hess",
                        "order_claimed": 1,
                        "observed_imaginary_modes": 1,
                        "result_artifact_sha256": SADDLE,
                    },
                },
            },
            {
                "kind": "result_quantities_extracted",
                "payload": {
                    "receipt_sha256": EXTRACTION,
                    "artifact_sha256": SADDLE,
                    "record": {
                        "artifact_id": "result.planar_hess",
                        "program": "pyscf",
                    },
                },
            },
            {
                "kind": "analysis_claims_recorded",
                "payload": {
                    "receipt_sha256": "5" * 64,
                    "record": {
                        "claims": [
                            {
                                "claim_id": "asgiven_freqs",
                                "quantity_id": "asgiven_freqs",
                                "display_value": -1050.8,
                                "display_unit": "cm^-1",
                                "dimension": [0, 0, 0, 0, 1, 0],
                                "source_receipt_sha256": EXTRACTION,
                            }
                        ]
                    },
                },
            },
            {
                "kind": "analysis_completion_evaluated",
                "payload": {
                    "receipt_sha256": "7" * 64,
                    "status": "passed",
                    "limitation_output_ids": [],
                    "declared_observable_misses": [],
                },
            },
        ],
    )
    goal_dir = workspace / ".chemsmart-agent" / "goals" / "g"
    goal_dir.mkdir(parents=True, exist_ok=True)
    ledger = GoalLedger(goal_dir)
    result = _settle_from_delivery(
        ledger,
        goal_id="g",
        cycles=2,
        revisions_admitted=1,
        events_path=session,
        terminal="complete",
        workspace=workspace,
    )
    assert result.settlement in {"achieved", "achieved_with_observations"}
    text = " ".join(result.reasons)
    assert (
        "delivered from a characterised result: asgiven_freqs" in text
    ), result.reasons


def test_the_record_names_the_results_the_validator_did_not_pass(tmp_path):
    from chemsmart.agent.workspace_record import failed_artifacts

    workspace = tmp_path / "ws"
    record_run(
        workspace,
        goal_id="g",
        cycle=1,
        run_events_path=_write(
            tmp_path / "cycle-1.jsonl",
            [
                {
                    "kind": "program_result_verified",
                    "payload": {
                        "node_id": "n",
                        "status": "invalid",
                        "receipt_sha256": "1" * 64,
                        "record": {
                            "node_id": "n",
                            "program": "pyscf",
                            "jobtype": "hess",
                            "state": "invalid",
                            "observations": {"program": "pyscf"},
                            "output_artifacts": [
                                {"kind": "pyscf_hdf5", "sha256": SADDLE}
                            ],
                        },
                    },
                }
            ],
        ),
        run="goals/g/runs/cycle-1",
    )
    assert failed_artifacts(workspace) == (SADDLE,)
