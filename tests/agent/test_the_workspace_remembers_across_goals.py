"""A host-written record spans goals in one workspace, and a divergence
between two delivered values of one claim is stated in the wake.

Two sealed windows delivered the sulfone's gauche preference in
acetonitrile at -2.23 and -0.91 kJ/mol, at two basis sets, on either
side of the author's 2 kJ/mol cutoff; no session could see both
(NOVEL-1/2 po2, 2026-09-04). The owner ruled the record built and a
host-seeded prior run admissible when written from receipts.
"""

from __future__ import annotations

import json

import pytest

from chemsmart.agent.workspace_record import (
    record_run,
    render_workspace_record,
    workspace_record_path,
)

from .test_the_goal_loop_recovers_or_returns import (
    _execute,
    _loop,
    _planning_session,
    _review_payload,
)

pytestmark = pytest.mark.capability("rule:wake.workspace_record")


def _review(path, basis, **level):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "workflow_execution_review": {
                    "node_reviews": [
                        {
                            "node_id": "sulfone_gauche_mecn",
                            "project_settings_text": json.dumps(
                                {
                                    "functional": "wb97x-d3bj",
                                    "basis": basis,
                                    "solvent_model": "cpcm",
                                    "solvent": "acetonitrile",
                                    **level,
                                }
                            ),
                            "project_settings_text_sha256": basis * 8,
                            "molecular_identity": {
                                "formula": "C8H9FO2S",
                                "charge": 0,
                                "multiplicity": 1,
                                "coordinate_identity": {
                                    "geometry_artifact_sha256": "6" * 64
                                },
                            },
                        }
                    ]
                }
            }
        ),
        encoding="utf-8",
    )
    return path


def _stream(path, *, energy, dg):
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        {
            "kind": "program_result_verified",
            "payload": {
                "record": {
                    "node_id": "sulfone_gauche_mecn",
                    "program": "orca",
                    "state": "valid",
                    "input_artifact_sha256": "2" * 64,
                    "receipt_sha256": "9" * 64,
                    "observations": {
                        "jobtype": "opt",
                        "orca": {
                            "charge": 0,
                            "multiplicity": 1,
                            "energy_hartree": energy,
                            "vibrational_mode_count": 57,
                        },
                    },
                }
            },
        },
        {
            "kind": "analysis_claims_recorded",
            "payload": {
                "receipt_sha256": "5" * 64,
                "record": {
                    "claims": [
                        {
                            "claim_id": "dg_gauche_anti_sulfone_mecn",
                            "quantity_id": "dg_gauche_anti_sulfone_mecn",
                            "display_value": dg,
                            "display_unit": "kJ/mol",
                            "source_receipt_sha256": "1" * 64,
                        }
                    ]
                },
            },
        },
    ]
    path.write_text(
        "\n".join(json.dumps(row) for row in rows) + "\n", encoding="utf-8"
    )
    return path


def test_a_level_names_its_frozen_core_and_the_root_it_followed(tmp_path):
    """Two results at one functional and basis are different calculations
    when one froze its core or optimised on an excited root, so the level
    line carries the convention and the response beside the method: a
    divergence between them is then named beside its cause."""

    workspace = tmp_path / "ws"
    record_run(
        workspace,
        goal_id="goal-fc",
        cycle=1,
        run_events_path=_stream(tmp_path / "fc.jsonl", energy=-958.3, dg=-2.0),
        run="goals/goal-fc/runs/cycle-1",
        review_file=_review(
            tmp_path / "rfc.json",
            "def2svp0",
            frozen_core=1,
            response_method="tda",
            state_manifold="singlet",
            nstates=3,
            excited_state_root=1,
        ),
    )
    (level,) = render_workspace_record(workspace)["levels"].values()
    assert level["basis"] == "def2svp0"
    assert level["frozen_core"] == 1
    assert level["response_method"] == "tda"
    assert level["state_manifold"] == "singlet"
    assert level["nstates"] == 3
    assert level["excited_state_root"] == 1


def test_two_goals_at_two_levels_leave_a_divergence(tmp_path):
    workspace = tmp_path / "ws"
    first = record_run(
        workspace,
        goal_id="goal-w1",
        cycle=1,
        run_events_path=_stream(
            tmp_path / "w1.jsonl", energy=-958.3559, dg=-2.23
        ),
        run="goals/goal-w1/runs/cycle-1",
        review_file=_review(tmp_path / "r1.json", "def2svp0"),
    )
    second = record_run(
        workspace,
        goal_id="goal-w2",
        cycle=1,
        run_events_path=_stream(
            tmp_path / "w2.jsonl", energy=-959.1827, dg=-0.91
        ),
        run="goals/goal-w2/runs/cycle-1",
        review_file=_review(tmp_path / "r2.json", "deftzvp0"),
    )
    assert (first, second) == (2, 2)
    assert workspace_record_path(workspace).exists()
    rendered = render_workspace_record(workspace)
    assert len(rendered["results"]) == 2
    assert {level["basis"] for level in rendered["levels"].values()} == {
        "def2svp0",
        "deftzvp0",
    }
    (divergence,) = rendered["divergences"]
    assert divergence["claim_id"] == "dggaucheantisulfonemecn"
    assert {v["claim_id"] for v in divergence["values"]} == {
        "dg_gauche_anti_sulfone_mecn"
    }
    assert divergence["unit"] == "kJ/mol"
    assert {item["value"] for item in divergence["values"]} == {-2.23, -0.91}
    assert "never explained" in rendered["meaning"]


def test_the_same_run_twice_is_no_divergence(tmp_path):
    workspace = tmp_path / "ws"
    for _ in range(2):
        record_run(
            workspace,
            goal_id="goal-w1",
            cycle=1,
            run_events_path=_stream(
                tmp_path / "w1.jsonl", energy=-958.3559, dg=-2.23
            ),
            run="goals/goal-w1/runs/cycle-1",
            review_file=_review(tmp_path / "r1.json", "def2svp0"),
        )
    assert render_workspace_record(workspace)["divergences"] == ()
    assert render_workspace_record(tmp_path / "fresh") == {}


def test_every_cycle_reads_the_record_and_the_rule(tmp_path):
    contexts = []

    def capture(inner):
        def step(workspace, kwargs):
            contexts.append(kwargs["goal_context"])
            return inner(workspace, kwargs)

        return step

    _loop(
        tmp_path,
        sessions=[
            capture(_planning_session("live-1", review=_review_payload())),
            capture(_planning_session("live-2", review=_review_payload())),
        ],
        executes=[
            _execute(
                tmp_path, failed=True, status="partial", analysis="partial"
            ),
            _execute(
                tmp_path, failed=True, status="partial", analysis="partial"
            ),
        ],
        max_revisions=1,
    )
    first, second = contexts
    assert first["workspace_record"] == {}
    assert "workspace_record" in second
    assert "workspace_record lists" in first["authority"]
    assert "workspace_record lists" in second["authority"]


def test_two_sessions_names_for_one_claim_join(tmp_path):
    """Two windows named the same delivered number dg-sulfone-mecn and
    dg_sulfone_mecn; a join on the raw id saw two claims and no
    disagreement. Case and separators are not identity."""

    from chemsmart.agent.workspace_record import divergences

    first = {
        "kind": "claim",
        "goal_id": "g1",
        "run": "r1",
        "claim_id": "dg-sulfone-mecn",
        "quantity_id": "dg-sulfone-mecn",
        "value": -2.23,
        "unit": "kJ/mol",
        "level_sha256s": ("a",),
    }
    second = {
        **first,
        "goal_id": "g2",
        "run": "r2",
        "claim_id": "dg_sulfone_mecn",
        "quantity_id": "dg_gauche_anti_sulfone_mecn",
        "value": -0.91,
        "level_sha256s": ("b",),
    }
    (found,) = divergences((first, second))
    assert found["claim_id"] == "dgsulfonemecn"
    assert {v["claim_id"] for v in found["values"]} == {
        "dg-sulfone-mecn",
        "dg_sulfone_mecn",
    }
