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

    tmp_path.mkdir(parents=True, exist_ok=True)
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
    tmp_path.mkdir(parents=True, exist_ok=True)
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


def test_a_fallback_level_is_never_a_pair(tmp_path):
    """The property a composed claim's qualification rests on.

    A claim whose ancestry walk *worked* across two levels carries a
    2-tuple. A claim that fell back carries `single_level`, which is a
    1-tuple when the run has exactly one level and `()` otherwise --
    never two. That asymmetry is what makes a composed claim the only
    test in this machinery that cannot pass by coincidence: no run shape
    produces a pair by accident.

    It is asserted here because the live qualification of the expression
    walk depends on it, and a later change that made the fallback return
    the run's whole level set would silently turn that qualification
    into a tautology without failing anything.
    """

    two = _record(
        tmp_path / "two",
        [
            _node("a", "c" * 64, "opt"),
            _node("b", "d" * 64, "sp"),
            _extraction("1" * 64, "c" * 64),
            _extraction("2" * 64, "d" * 64),
            # No extraction resolves this claim's receipt, so the walk
            # finds nothing and the fallback decides -- on a run that
            # carries two levels.
            _claim("unresolvable", "9" * 64),
        ],
    )
    assert tuple(two["unresolvable"]["level_sha256s"]) == (), (
        "the fallback produced a level set on a two-level run, so a "
        "composed claim's 2-tuple would no longer prove the walk ran"
    )

    one = _record(
        tmp_path / "one",
        [
            _node("a", "c" * 64, "opt"),
            _extraction("1" * 64, "c" * 64),
            _claim("unresolvable", "9" * 64),
        ],
        review=_single_level_review(tmp_path / "one"),
    )
    assert len(one["unresolvable"]["level_sha256s"]) == 1


def _single_level_review(tmp_path):
    """One node, one level -- so the fallback has exactly one answer."""

    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / "review.json"
    path.write_text(
        json.dumps(
            {
                "workflow_execution_review": {
                    "node_reviews": [
                        {
                            "node_id": "a",
                            "project_settings_text": "B3LYP/def2-SVP",
                            "project_settings_text_sha256": _CHEAP,
                        }
                    ]
                }
            }
        ),
        encoding="utf-8",
    )
    return path


def _expression_per_output(receipt, outputs):
    """An expression receipt with several outputs, each on its own pair.

    The shape the host actually emits, measured live on CUHK
    (`water-levels-1`, one water geometry at three basis sets):

        out d-svp-321g   sources [18cb2025, 6ae3f89c]
        out d-tzvp-svp   sources [18cb2025, 3bec2300]
        out d-tzvp-321g  sources [3bec2300, 6ae3f89c]

    Three outputs, two sources each, and the pairs differ.
    """

    from chemsmart.agent._contracts import canonical_data
    from chemsmart.analysis.quantity_expressions import (
        QuantityExpressionOutputDependencyV1,
    )

    return {
        "kind": "quantity_expression_evaluated",
        "payload": {
            "receipt_sha256": receipt,
            "record": {
                "output_dependencies": [
                    canonical_data(
                        QuantityExpressionOutputDependencyV1(
                            output_id=output_id,
                            source_receipt_sha256s=tuple(sorted(sources)),
                            model_authored_constants=(),
                            convention_operations=(),
                            arithmetic_node_count=1,
                        )
                    )
                    for output_id, sources in outputs
                ]
            },
        },
    }


def _claim_of(claim_id, quantity_id, receipt):
    return {
        "kind": "analysis_claims_recorded",
        "payload": {
            "receipt_sha256": "9" * 64,
            "record": {
                "claims": [
                    {
                        "claim_id": claim_id,
                        "quantity_id": quantity_id,
                        "source_receipt_sha256": receipt,
                        "display_value": 1.0,
                        "display_unit": "kcal/mol",
                        "dimension": [0, 0, 0, 0, 0, 0, 0],
                    }
                ]
            },
        },
    }


def _three_level_review(tmp_path):
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / "review.json"
    path.write_text(
        json.dumps(
            {
                "workflow_execution_review": {
                    "node_reviews": [
                        {
                            "node_id": f"calc-{n}",
                            "project_settings_text": f"B3LYP/{n}",
                            "project_settings_text_sha256": sha,
                        }
                        for n, sha in (
                            ("a", "a" * 64),
                            ("b", "b" * 64),
                            ("c", "c" * 64),
                        )
                    ]
                }
            }
        ),
        encoding="utf-8",
    )
    return path


def test_a_composed_claim_names_its_own_terms_not_the_expressions(tmp_path):
    """One expression, three outputs, three different pairs of ancestors.

    The wrong-key repair read `output_dependencies` and then **flattened
    every row into one list keyed by the expression receipt**, so every
    claim standing on that receipt inherited the union of all its
    outputs' sources. Measured live: three basis-set differences, each
    genuinely standing on two of three calculations, every one of them
    recorded with all three levels.

    That is A1's own defect one layer up -- a claim naming levels it does
    not stand on -- and it is worse than the run-level version it
    replaced, because three claims carrying an identical 3-tuple compare
    as `same_level` in `divergences()` and the pair is skipped. The
    suppression the per-node level exists to prevent, reintroduced by
    its own repair.
    """

    claims = _record(
        tmp_path,
        [
            _node("calc-a", "1" * 64, "sp"),
            _node("calc-b", "2" * 64, "sp"),
            _node("calc-c", "3" * 64, "sp"),
            _extraction("a1" + "0" * 62, "1" * 64),
            _extraction("b1" + "0" * 62, "2" * 64),
            _extraction("c1" + "0" * 62, "3" * 64),
            _expression_per_output(
                "e1" + "0" * 62,
                (
                    ("d-ab", ("a1" + "0" * 62, "b1" + "0" * 62)),
                    ("d-bc", ("b1" + "0" * 62, "c1" + "0" * 62)),
                    ("d-ac", ("a1" + "0" * 62, "c1" + "0" * 62)),
                ),
            ),
            _claim_of("delta-ab", "d-ab", "e1" + "0" * 62),
            _claim_of("delta-bc", "d-bc", "e1" + "0" * 62),
        ],
        review=_three_level_review(tmp_path / "review"),
    )

    assert tuple(claims["delta-ab"]["level_sha256s"]) == (
        "a" * 64,
        "b" * 64,
    ), claims["delta-ab"]["level_sha256s"]
    assert tuple(claims["delta-bc"]["level_sha256s"]) == (
        "b" * 64,
        "c" * 64,
    ), claims["delta-bc"]["level_sha256s"]
    assert claims["delta-ab"]["level_sha256s"] != (
        claims["delta-bc"]["level_sha256s"]
    ), (
        "two differences over different pairs of levels compare as the "
        "same level, so divergences() skips them -- which is the "
        "suppression the per-node level exists to prevent"
    )


def test_a_later_analysis_only_run_names_prior_results_own_levels(tmp_path):
    """A later run can resolve the results it reads from the record.

    The production provenance recheck has no new
    ``program_result_verified`` event: it extracts three registered PySCF
    results, composes three two-level differences, and records claims.  The
    initial per-output repair was confined to its own event stream, so all
    three later claims silently lost their levels.  This drives the durable
    result rows into a second, analysis-only projection and checks the
    exact consumer path rather than seeding a private map.
    """

    workspace = tmp_path / "ws"
    first = tmp_path / "first-events.jsonl"
    first.write_text(
        "\n".join(
            json.dumps(row)
            for row in (
                _node("calc-a", "1" * 64, "sp"),
                _node("calc-b", "2" * 64, "sp"),
                _node("calc-c", "3" * 64, "sp"),
            )
        )
        + "\n",
        encoding="utf-8",
    )
    record_run(
        workspace,
        goal_id="water-levels-1",
        cycle=1,
        run_events_path=first,
        run="goals/water-levels-1/runs/cycle-1",
        review_file=_three_level_review(tmp_path / "review"),
    )

    current = tmp_path / "analysis-only-events.jsonl"
    current.write_text(
        "\n".join(
            json.dumps(row)
            for row in (
                _extraction("a1" + "0" * 62, "1" * 64),
                _extraction("b1" + "0" * 62, "2" * 64),
                _extraction("c1" + "0" * 62, "3" * 64),
                _expression_per_output(
                    "e1" + "0" * 62,
                    (
                        ("d-ab", ("a1" + "0" * 62, "b1" + "0" * 62)),
                        ("d-bc", ("b1" + "0" * 62, "c1" + "0" * 62)),
                        ("d-ac", ("a1" + "0" * 62, "c1" + "0" * 62)),
                    ),
                ),
                _claim_of("delta-ab", "d-ab", "e1" + "0" * 62),
                _claim_of("delta-bc", "d-bc", "e1" + "0" * 62),
                _claim_of("delta-ac", "d-ac", "e1" + "0" * 62),
            )
        )
        + "\n",
        encoding="utf-8",
    )
    record_run(
        workspace,
        goal_id="water-levels-provenance-a1-r2",
        cycle=1,
        run_events_path=current,
        run="goals/water-levels-provenance-a1-r2/runs/live",
    )

    claims = {
        str(entry["claim_id"]): entry
        for entry in read_workspace_record(workspace)
        if entry.get("goal_id") == "water-levels-provenance-a1-r2"
    }
    assert {
        claim_id: tuple(claim["level_sha256s"])
        for claim_id, claim in claims.items()
    } == {
        "delta-ab": ("a" * 64, "b" * 64),
        "delta-bc": ("b" * 64, "c" * 64),
        "delta-ac": ("a" * 64, "c" * 64),
    }
