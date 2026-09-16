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


# --- What the fallback claimed, and could not know -------------------
#
# Where the receipt join did not resolve, the run's own level was used
# "when the run has exactly one, because then it is unambiguous". Three
# cases make that false. Two of them have an exact answer the record
# already carries and nobody read; the third has none, and answering it
# anyway is how a number acquires a level it was never computed at.


def _thermochemistry(receipt, artifact):
    return {
        "kind": "thermochemistry_derived",
        "payload": {
            "receipt_sha256": receipt,
            "artifact_sha256": artifact,
            "record": {},
        },
    }


def _expression(receipt, source_receipts):
    """An expression receipt, in the shape the host actually emits.

    This invented an ``inputs`` key. ``inputs`` is on the expression
    *request*, which is never emitted; the receipt carries
    ``output_dependencies``, whose rows name the receipts each output
    stands on. The reader read the invented key, so the whole ancestry
    walk was dead code and this test was the reason nobody noticed --
    it is the defect class the round exists to remove, committed by the
    witness rather than by the product.

    Built here from the real dataclass so the shape cannot drift again.
    """

    from chemsmart.agent._contracts import canonical_data
    from chemsmart.analysis.quantity_expressions import (
        QuantityExpressionOutputDependencyV1,
    )

    dependency = QuantityExpressionOutputDependencyV1(
        output_id="composed",
        source_receipt_sha256s=tuple(sorted(set(source_receipts))),
        model_authored_constants=(),
        convention_operations=(),
        arithmetic_node_count=1,
    )
    return {
        "kind": "quantity_expression_evaluated",
        "payload": {
            "receipt_sha256": receipt,
            "record": {
                "output_dependencies": [canonical_data(dependency)],
            },
        },
    }


def _record(tmp_path, rows, *, review=None):
    workspace = tmp_path / "ws"
    stream = tmp_path / "events.jsonl"
    stream.write_text(
        "\n".join(json.dumps(row) for row in rows) + "\n", encoding="utf-8"
    )
    record_run(
        workspace,
        goal_id="g",
        cycle=1,
        run_events_path=stream,
        run="goals/g/runs/cycle-1",
        review_file=review or _review(tmp_path),
    )
    return _claims(workspace)


def test_a_claim_on_a_result_this_run_never_computed_has_no_level(tmp_path):
    """An imported result is not this run's level, and not any level here.

    A goal's later cycle claims from a result an earlier cycle computed,
    or from a registered result the workspace already held. The receipt
    resolves to an artifact; the artifact is simply not one this run
    produced, so this run cannot say what level it was computed at. The
    fallback answered anyway, with whatever single level this run
    happened to have -- a wrong level, which reads as fact where an
    absent one reads as unknown.
    """

    claims = _record(
        tmp_path,
        [
            _node("conformer-a-opt", "c" * 64, "opt"),
            # The extraction reads an artifact no node here produced.
            _extraction("1" * 64, "e" * 64),
            _claim("imported_energy", "1" * 64),
        ],
    )
    assert tuple(claims["imported_energy"]["level_sha256s"]) == (), (
        "a claim on a result this run never computed was stamped with "
        "this run's level: "
        f"{claims['imported_energy']['level_sha256s']}"
    )


def test_a_thermochemistry_claim_names_the_result_it_stands_on(tmp_path):
    """The exact answer was in the receipt all along.

    A thermochemistry receipt carries the artifact it derived from. The
    join skipped it and fell to the run's single level, which is right
    only by luck: in a two-level cohort it is a coin toss.
    """

    claims = _record(
        tmp_path,
        [
            _node("conformer-a-opt", "c" * 64, "opt"),
            _node("conformer-a-sp", "d" * 64, "sp"),
            _thermochemistry("3" * 64, "c" * 64),
            _claim("gibbs_free_energy", "3" * 64),
        ],
    )
    assert tuple(claims["gibbs_free_energy"]["level_sha256s"]) == (_CHEAP,), (
        "the thermochemistry receipt names its own artifact and the "
        "record read the run instead: "
        f"{claims['gibbs_free_energy']['level_sha256s']}"
    )


def test_an_expression_claim_names_every_level_it_composed(tmp_path):
    """A composed number stands on every level underneath it.

    An expression's own receipt names, per input, the receipt it read.
    Following that to the extraction and on to the artifact gives the
    levels the value was actually built from -- which for the ordinary
    high-level-single-point-on-a-cheap-geometry protocol is two, and for
    the run's single-level fallback was one or, in a mixed cohort,
    whichever one the run happened to carry.
    """

    claims = _record(
        tmp_path,
        [
            _node("conformer-a-opt", "c" * 64, "opt"),
            _node("conformer-a-sp", "d" * 64, "sp"),
            _extraction("1" * 64, "c" * 64),
            _extraction("2" * 64, "d" * 64),
            _expression("4" * 64, ["1" * 64, "2" * 64]),
            _claim("reaction_energy", "4" * 64),
        ],
    )
    assert tuple(claims["reaction_energy"]["level_sha256s"]) == tuple(
        sorted((_CHEAP, _COSTLY))
    ), (
        "an expression composing two levels reported "
        f"{claims['reaction_energy']['level_sha256s']}"
    )


def test_a_composed_number_with_one_imported_term_claims_no_level(tmp_path):
    """Partial ancestry is not a level; it is a level this run cannot say.

    A value built from one result this run computed and one it imported
    has two levels and the run knows one. Reporting that one would tell
    `divergences` that two such claims are at the same level when they
    are not -- the suppression per-node attribution exists to prevent,
    reintroduced one hop further down.
    """

    claims = _record(
        tmp_path,
        [
            _node("conformer-a-opt", "c" * 64, "opt"),
            _extraction("1" * 64, "c" * 64),
            # The other term reads a result this run never produced.
            _extraction("2" * 64, "e" * 64),
            _expression("4" * 64, ["1" * 64, "2" * 64]),
            _claim("reaction_energy", "4" * 64),
        ],
    )
    assert tuple(claims["reaction_energy"]["level_sha256s"]) == (), (
        "a half-known ancestry was reported as a level: "
        f"{claims['reaction_energy']['level_sha256s']}"
    )
