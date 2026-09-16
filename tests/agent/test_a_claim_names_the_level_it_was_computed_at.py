"""A claim carries the level of the node it stands on, not the run's.

``record_run`` collected every node's ``project_settings_text_sha256``
into one ``levels_in_run`` set and stamped the whole set onto every claim
row. One run was usually one level, so it was harmless and nobody looked.

A wave cohort is the first design that deliberately puts several
independent calculations in one run, and a wave at several levels is
ordinary science -- three conformers at one functional, then the winner
at a better one. Under that, every claim in the cycle would carry every
level in the cycle, and ``divergences()`` uses exactly that field to
decide whether two values of one claim are comparable: it would compare a
number computed at B3LYP against one computed at CCSD(T) and call them a
disagreement about the same quantity.

The rows are durable and append-only, so this cannot be repaired after a
multi-level cohort has written them. It lands first.
"""

from __future__ import annotations

import json

from chemsmart.agent.workspace_record import read_workspace_record, record_run

# No capability marker: no registered gate covers "a claim names the level
# it was computed at" yet. Inventing one would be a declaration with no
# reader, which is the defect class this round exists to remove.

_CHEAP = "a" * 64
_COSTLY = "b" * 64


def _review(tmp_path):
    """Two nodes, two levels -- the ordinary two-stage protocol."""

    path = tmp_path / "review.json"
    path.write_text(
        json.dumps(
            {
                "workflow_execution_review": {
                    "node_reviews": [
                        {
                            "node_id": "conformer-a-opt",
                            "project_settings_text": "B3LYP/def2-SVP",
                            "project_settings_text_sha256": _CHEAP,
                        },
                        {
                            "node_id": "conformer-a-sp",
                            "project_settings_text": "DLPNO-CCSD(T)/def2-TZVP",
                            "project_settings_text_sha256": _COSTLY,
                        },
                    ]
                }
            }
        ),
        encoding="utf-8",
    )
    return path


def _node(node_id, artifact, jobtype):
    return {
        "kind": "program_result_verified",
        "payload": {
            "record": {
                "node_id": node_id,
                "state": "valid",
                "jobtype": jobtype,
                "observations": {"jobtype": jobtype, "orca": {}},
                "output_artifacts": [{"sha256": artifact}],
            }
        },
    }


def _extraction(receipt, artifact):
    return {
        "kind": "result_quantities_extracted",
        "payload": {
            "receipt_sha256": receipt,
            "artifact_sha256": artifact,
            "record": {},
        },
    }


def _claim(claim_id, receipt, value=1.0):
    return {
        "kind": "analysis_claims_recorded",
        "payload": {
            "receipt_sha256": "9" * 64,
            "record": {
                "claims": [
                    {
                        "claim_id": claim_id,
                        "quantity_id": claim_id,
                        "source_receipt_sha256": receipt,
                        "display_value": value,
                        "display_unit": "kcal/mol",
                        "dimension": [0, 0, 0, 0, 0, 0, 0],
                    }
                ]
            },
        },
    }


def _stream(tmp_path):
    rows = [
        _node("conformer-a-opt", "c" * 64, "opt"),
        _node("conformer-a-sp", "d" * 64, "sp"),
        _extraction("1" * 64, "c" * 64),
        _extraction("2" * 64, "d" * 64),
        _claim("geometry_rmsd", "1" * 64),
        _claim("electronic_energy", "2" * 64),
    ]
    stream = tmp_path / "events.jsonl"
    stream.write_text(
        "\n".join(json.dumps(row) for row in rows) + "\n", encoding="utf-8"
    )
    return stream


def _claims(workspace):
    return {
        str(entry.get("claim_id")): entry
        for entry in read_workspace_record(workspace)
        if entry.get("kind") == "claim"
    }


def test_each_claim_carries_only_its_own_node_s_level(tmp_path):
    workspace = tmp_path / "ws"
    record_run(
        workspace,
        goal_id="g",
        cycle=1,
        run_events_path=_stream(tmp_path),
        run="goals/g/runs/cycle-1",
        review_file=_review(tmp_path),
    )
    claims = _claims(workspace)
    assert set(claims) == {"geometry_rmsd", "electronic_energy"}

    assert tuple(claims["geometry_rmsd"]["level_sha256s"]) == (_CHEAP,), (
        "a claim standing on the optimisation carries the whole run's "
        "levels, so the record cannot say what it was computed at: "
        f"{claims['geometry_rmsd']['level_sha256s']}"
    )
    assert tuple(claims["electronic_energy"]["level_sha256s"]) == (_COSTLY,)


def test_two_levels_of_one_quantity_are_seen_rather_than_suppressed(
    tmp_path,
):
    """The consumer, not the field -- and the bug hides rather than invents.

    `divergences()` skips a pair when it is `same_run and same_level`.
    With the run's whole level set stamped on every claim, two values of
    one quantity computed at two levels in one cycle are *identical* in
    that field, so the pair is silently skipped: the host sees a cheap
    number and an expensive number for the same quantity and says
    nothing. Per-node levels make them different, which is what the
    function was written to report -- with each value carrying the level
    it was actually computed at, so a reader can tell which is which.
    """

    from chemsmart.agent.workspace_record import divergences

    workspace = tmp_path / "ws"
    rows = [
        _node("conformer-a-opt", "c" * 64, "opt"),
        _node("conformer-a-sp", "d" * 64, "sp"),
        _extraction("1" * 64, "c" * 64),
        _extraction("2" * 64, "d" * 64),
        # One quantity at two levels -- an ordinary two-stage protocol.
        _claim("binding_energy", "1" * 64, value=10.0),
        _claim("binding_energy", "2" * 64, value=14.0),
    ]
    stream = tmp_path / "two" / "events.jsonl"
    stream.parent.mkdir(parents=True)
    stream.write_text(
        "\n".join(json.dumps(row) for row in rows) + "\n", encoding="utf-8"
    )
    record_run(
        workspace,
        goal_id="g",
        cycle=1,
        run_events_path=stream,
        run="goals/g/runs/cycle-1",
        review_file=_review(tmp_path),
    )
    found = divergences(read_workspace_record(workspace))
    assert len(found) == 1, (
        "two levels of one quantity in one cycle were suppressed, because "
        "both rows carried the same run-wide level set and the pair read "
        f"as same-run-same-level: {found}"
    )
    levels = {tuple(row["level_sha256s"]) for row in found[0]["values"]}
    assert levels == {
        (_CHEAP,),
        (_COSTLY,),
    }, f"the reader cannot tell which value is which level: {levels}"
