"""The surface is a tree: a stem every session reads, and leaves the host
opens from four signals or the model opens itself. Opening a leaf changes
what the model can express and how much it reads, never what the host
approves; every activation is an event carrying the new schema digest.
"""

from __future__ import annotations

import json

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.guides import (
    GUIDES,
    LEAF_OPERATIONS,
    LEAF_TOOLS,
    guides_from_plan,
    guides_from_states,
    guides_from_text,
    guides_from_workspace,
)
from chemsmart.agent.rules import POLICY_RULES
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.runtime.events import EventKind
from chemsmart.agent.runtime.reducer import replay_events
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.agent.tool_specs import build_command_compiled_tool_surface

pytestmark = pytest.mark.capability("guide:*")


def _names(surface):
    return {item["function"]["name"] for item in surface.tool_definitions}


def _operations(surface):
    evaluator = next(
        item["function"]
        for item in surface.tool_definitions
        if item["function"]["name"] == "evaluate_quantity_expression"
    )
    return set(
        evaluator["parameters"]["properties"]["nodes"]["items"]["properties"][
            "operation"
        ]["enum"]
    )


def test_the_stem_hides_every_leaf_tool_and_operation():
    stem = build_command_compiled_tool_surface()
    assert not (set(LEAF_TOOLS) & _names(stem))
    assert not (set(LEAF_OPERATIONS) & _operations(stem))
    assert "open_guide" in _names(stem)
    # Moved once, 90_000 -> 91_000 (2026-09-06): two owner-ruled
    # point-of-use affordances joined the stem -- the diagnostic
    # declaration (R3) and the input-check probe's sentence (R2) --
    # after their texts were cut to the bone. The next growth trims
    # something else, and did (2026-09-09): fetch_pubchem_geometry
    # joined the stem and the cap held, because the identifier spelling
    # rule stopped being repeated in prose on 45 fields whose JSON
    # pattern already stated it.
    # Not raised again. Stem headroom had fallen to 7 bytes, and the
    # cause was leaf-family guidance sitting in the stem: the CBS
    # exponent's ~500-character description lived on the shared
    # expression-node schema every session reads, while the four CBS
    # operations it serves are leaf operations the stem does not expose,
    # and the `cbs` guide already said the exponent comes from the
    # protocol. Moving it there recovered 834 bytes. What belongs in the
    # stem is what a field IS and what the host records about it; when
    # to reach for it belongs to the guide that opens for that family.
    # Raised 91_000 -> 92_000 (2026-09-16, Round A) for select_execution_wave,
    # measured at 868 bytes on this surface: 275 of them are the placed
    # rule that says an undispatchable wave is a per-member verdict
    # rather than an error, and the rest is the schema. Nothing was
    # trimmed to pay for it, and the record says so rather than
    # implying a trade that did not happen: the stem held 179 bytes of
    # headroom before this change, and the wave model is the contract
    # the round exists for -- the only cheaper home for the tool is a
    # guide, and a session that never opens it would plan waves it
    # cannot name. The three sentences teaching the model *why* it has
    # waves sit in the system prompt, not here, so this pays for the
    # affordance and not for the lesson.
    # Raised 92_000 -> 93_000 (2026-09-17, Round A) for the explicit
    # continue_execution_reasoning affordance. It is a narrow second
    # execution-boundary choice, not a scheduler surface.
    assert len(json.dumps(stem.tool_definitions)) < 93_000


def test_every_guide_adds_exactly_its_tools_and_operations():
    stem = build_command_compiled_tool_surface()
    for guide in GUIDES:
        opened = build_command_compiled_tool_surface(guides=(guide.guide_id,))
        assert _names(opened) - _names(stem) == set(
            guide.tools
        ), guide.guide_id
        assert _operations(opened) - _operations(stem) == set(
            guide.operations
        ), guide.guide_id
        assert opened.tool_schema_sha256 != stem.tool_schema_sha256 or (
            not guide.tools and not guide.operations
        )


def test_the_four_signals():
    assert "scan" in guides_from_text("Run a relaxed torsional scan of butane")
    assert "saddle" in guides_from_text("locate the transition state and IRC")
    assert guides_from_workspace(("chemsmart_db",)) == ("database",)
    assert guides_from_plan(jobtypes=("irc",)) == ("saddle",)
    assert guides_from_plan(operations=("gibbs_to_pka",)) == ("constants",)
    assert guides_from_plan(tools=("edit_molecular_geometry",)) == (
        "structure",
    )
    assert set(guides_from_states(("failed_wrong_stationary_point",))) == {
        "recovery",
        "saddle",
        "structure",
    }
    assert guides_from_states(("timeout_terminated",)) == ("recovery",)
    assert guides_from_states(("validated",)) == ()
    # Two programs in one DAG is the crossprogram signal; one program
    # spelled two ways is not.
    assert guides_from_plan(programs=("xtb", "orca")) == ("crossprogram",)
    assert guides_from_plan(programs=("orca", "ORCA")) == ()
    assert "crossprogram" in guides_from_text(
        "GFN2 optimisation, then B3LYP single points on those geometries"
    )


def _host(tmp_path, **kwargs):
    return CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="guide-session"
        ),
        artifacts={},
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "workspace",
        **kwargs,
    )


def test_opening_a_guide_extends_the_surface_and_records_the_digest(tmp_path):
    host = _host(tmp_path)
    before = host.surface.tool_schema_sha256
    assert "compose_molecular_arrangement" not in _names(host.surface)
    reply = host.dispatch(
        turn_id="t1",
        tool_name="open_guide",
        arguments={"guide_id": "structure"},
    )
    record = reply["result"]
    assert record["guide_id"] == "structure"
    assert "compose_molecular_arrangement" in record["tools_now_available"]
    assert "compose_molecular_arrangement" in _names(host.surface)
    assert host.surface.tool_schema_sha256 != before
    events = [
        json.loads(line)
        for line in (tmp_path / "events.jsonl").read_text().splitlines()
    ]
    activated = [
        e for e in events if e["kind"] == EventKind.GUIDE_ACTIVATED.value
    ]
    assert len(activated) == 1
    assert activated[0]["payload"]["signal"] == "model"
    assert activated[0]["payload"]["tool_schema_sha256"] == (
        host.surface.tool_schema_sha256
    )
    assert replay_events(host.event_store.read_events()).active_guides == [
        "structure"
    ]
    # Opening it again is idempotent: no second event, same digest.
    host.dispatch(
        turn_id="t2",
        tool_name="open_guide",
        arguments={"guide_id": "structure"},
    )
    events = [
        json.loads(line)
        for line in (tmp_path / "events.jsonl").read_text().splitlines()
    ]
    assert (
        len(
            [e for e in events if e["kind"] == EventKind.GUIDE_ACTIVATED.value]
        )
        == 1
    )


def test_a_leaf_tool_called_by_name_opens_its_guide_first(tmp_path):
    host = _host(tmp_path)
    with pytest.raises(ContractError) as failure:
        host.dispatch(
            turn_id="t1",
            tool_name="inspect_database_records",
            arguments={"artifact_id": "db.missing"},
        )
    # The guide opened on the model's own call; the refusal is the
    # handler's (no such artifact), not "tool is not exposed".
    assert "not exposed" not in str(failure.value)
    assert "database" in host.active_guides


def test_a_session_started_with_guides_reads_them(tmp_path):
    host = _host(tmp_path, active_guides=("scan",))
    assert "bind_scan_point_geometry" in _names(host.surface)
    assert "coordinate_at_minimum" in _operations(host.surface)


def test_an_unknown_guide_names_what_exists(tmp_path):
    host = _host(tmp_path)
    with pytest.raises(ContractError, match="guides:.*skills:"):
        host._open_guide("t1", {"guide_id": "kinetics"})


def test_a_skill_is_still_reachable_through_open_guide(tmp_path):
    host = _host(tmp_path)
    reply = host._open_guide("t1", {"guide_id": "method-adequacy"})
    assert reply["skill_id"] == "method-adequacy"


def test_no_rule_is_placed_on_a_leaf_that_does_not_exist():
    guide_ids = {guide.guide_id for guide in GUIDES}
    for rule in POLICY_RULES:
        if rule.placement.startswith("leaf:"):
            assert rule.placement.split(":", 1)[1] in guide_ids, rule.rule_id


def test_a_point_expectation_is_a_band_with_equal_ends(tmp_path):
    """Observed live (W1c, R2c): every session declared is_minimum 1..1
    or imaginary count 0..0 and was refused; one retried three times."""

    host = _host(tmp_path)
    reply = host.dispatch(
        turn_id="t1",
        tool_name="declare_requested_observable",
        arguments={
            "observables": [
                {
                    "observable_id": "is-minimum",
                    "unit": "1",
                    "meaning": "one when every frequency is real",
                    "expected_low": 1,
                    "expected_high": 1,
                    "expectation_basis": "the task asks for a minimum",
                }
            ]
        },
    )
    assert reply["status"] == "ok"
    with pytest.raises(ContractError, match="must not exceed"):
        host.dispatch(
            turn_id="t2",
            tool_name="declare_requested_observable",
            arguments={
                "observables": [
                    {
                        "observable_id": "zpe",
                        "unit": "kcal/mol",
                        "meaning": "harmonic zero-point energy",
                        "expected_low": 35,
                        "expected_high": 30,
                        "expectation_basis": "a typo",
                    }
                ]
            },
        )


def test_a_sign_the_band_excludes_is_refused(tmp_path):
    """Five correct zero imaginary-mode counts printed "diverged"
    because their expectation carried expected_sign positive with a
    0..0 band; a zero has no sign."""

    host = _host(tmp_path)
    with pytest.raises(ContractError, match="a zero has no sign"):
        host.dispatch(
            turn_id="t1",
            tool_name="declare_requested_observable",
            arguments={
                "observables": [
                    {
                        "observable_id": "n-imag",
                        "unit": "1",
                        "meaning": "imaginary modes below -20 cm^-1",
                        "expected_sign": "positive",
                        "expected_low": 0,
                        "expected_high": 0,
                        "expectation_basis": "a minimum has none",
                    }
                ]
            },
        )


def test_the_goals_first_declaration_stands(tmp_path):
    """A woken session re-declared its expectations with a flipped
    sign and wider bands and the completion row printed agreed over a
    falsified first prior. The host is seeded with the goal's first
    declarations; a re-declaration keeps them and the reply says so."""

    first = {
        "observable_id": "cis-barrier",
        "unit": "kcal/mol",
        "dimension": (1, 0, 0, 0, 0, 0),
        "meaning": "syn barrier above anti",
        "expectation_basis": "torsional barriers of chloroethanes",
        "expected_sign": "positive",
        "expected_low": 3.0,
        "expected_high": 8.0,
    }
    host = _host(tmp_path, approved_requested_observable_declarations=[first])
    reply = host.dispatch(
        turn_id="t1",
        tool_name="declare_requested_observable",
        arguments={
            "observables": [
                {
                    "observable_id": "cis-barrier",
                    "unit": "kcal/mol",
                    "meaning": "syn barrier above anti",
                    "expected_sign": "positive",
                    "expected_low": 5.0,
                    "expected_high": 9.0,
                    "expectation_basis": "widened after the fact",
                }
            ]
        },
    )["result"]
    assert list(reply["kept_prior"]) == ["cis-barrier"]
    assert list(reply["declared"]) == []
    kept = host.requested_observable_declarations["cis-barrier"]
    assert (kept["expected_low"], kept["expected_high"]) == (3.0, 8.0)


def test_an_activation_term_matches_whole_words_only():
    """ "base" inside "database" opened the constants guide on a task
    with no constant in it."""

    from chemsmart.agent.guides import guides_from_text

    opened = guides_from_text(
        "five radicals from a workspace database; opt+freq each record"
    )
    assert "database" in opened
    assert "constants" not in opened
    assert "constants" in guides_from_text("the reference acid's pKa")


def test_a_declared_observable_is_answered_only_by_a_claim_carrying_its_id(
    tmp_path,
):
    """The completion gate certified a delivery whose declared endo:exo
    ratio had no claim, because two imaginary-mode counts share its
    dimension, and the expectation row printed agreed on a number the
    session had relabelled. A dimension is not an identity."""

    from types import SimpleNamespace

    ratio = {
        "observable_id": "endo-exo-ratio",
        "unit": "1",
        "dimension": (0, 0, 0, 0, 0, 0),
        "meaning": "endo over exo at 298 K",
        "expectation_basis": "the endo rule",
        "expected_low": 1.2,
        "expected_high": 150.0,
    }
    host = _host(tmp_path, approved_requested_observable_declarations=[ratio])
    claim = SimpleNamespace(
        claim_id="n-imag-ts-a",
        dimension=(0, 0, 0, 0, 0, 0),
        display_value=1.0,
        display_unit="1",
    )
    host.analysis_claim_records["r1"] = SimpleNamespace(
        task_spec_sha256="a" * 64, claims=(claim,)
    )
    misses, limitations = host._declared_observable_completion(
        task_spec_sha256="a" * 64
    )
    assert limitations == ("declared_observable:endo-exo-ratio",)
    assert "no delivered claim named 'endo-exo-ratio'" in misses[0]
    (row,) = host._declared_observable_predictions(task_spec_sha256="a" * 64)
    assert row["agreement"] == "not_comparable"
    assert row["delivered_claim_id"] == ""

    named = SimpleNamespace(
        claim_id="endo-exo-ratio",
        dimension=(0, 0, 0, 0, 0, 0),
        display_value=37.0,
        display_unit="1",
    )
    host.analysis_claim_records["r2"] = SimpleNamespace(
        task_spec_sha256="a" * 64, claims=(named,)
    )
    assert host._declared_observable_completion(task_spec_sha256="a" * 64) == (
        (),
        (),
    )
    (row,) = host._declared_observable_predictions(task_spec_sha256="a" * 64)
    assert row["agreement"] == "agreed"


def test_a_falsified_expectation_is_an_observation_never_a_limitation(
    tmp_path,
):
    """A pre-registered band the physics left is a result: the completion
    carries it under its own prefix in the observation list, stays passed,
    and names no limitation, so the settlement word carries it."""

    from types import SimpleNamespace

    barrier = {
        "observable_id": "cis-barrier",
        "unit": "kcal/mol",
        "dimension": (1, 0, 0, 0, 0, 0),
        "meaning": "syn barrier above anti",
        "expectation_basis": "torsional barriers of chloroethanes",
        "expected_sign": "positive",
        "expected_low": 3.0,
        "expected_high": 8.0,
    }
    host = _host(
        tmp_path, approved_requested_observable_declarations=[barrier]
    )
    host.analysis_claim_records["r1"] = SimpleNamespace(
        task_spec_sha256="a" * 64,
        claims=(
            SimpleNamespace(
                claim_id="cis-barrier",
                dimension=(1, 0, 0, 0, 0, 0),
                display_value=12.0,
                display_unit="kcal/mol",
            ),
        ),
    )
    (row,) = host._declared_observable_predictions(task_spec_sha256="a" * 64)
    assert row["agreement"] == "diverged"
    # The policy identity is a digest, not a plan object: a delivery made
    # directly from registered results has no plan and is certified from
    # its own declarations instead.
    (digest,) = host._record_toolchain_completion(
        "b" * 64,
        task_spec_sha256="a" * 64,
        source_receipt_sha256s=("c" * 64,),
    )
    completion = host.analysis_completion_receipts[digest]
    assert completion.status == "passed"
    assert completion.limitation_output_ids == ()
    assert completion.anomaly_output_ids == (
        "falsified_expectation:cis-barrier",
    )


def test_a_host_opened_guide_delivers_its_body_once(tmp_path):
    """36 host activations in one day's sessions exposed their tools and
    delivered no body; every body that reached the model was a manual
    re-read. The session-start helper returns each opened guide's record,
    body included, exactly once."""

    from chemsmart.agent.live_session import (
        _open_session_guides,
        _public_context,
    )

    host = _host(tmp_path)
    records = _open_session_guides(
        host, {"task": ("saddle", "structure"), "states": ("saddle",)}
    )
    assert [item["guide_id"] for item in records] == ["saddle", "structure"]
    assert all(item["body"] for item in records)
    assert _open_session_guides(host, {"task": ("saddle",)}) == ()
    context = _public_context(
        task="t",
        task_spec_sha256="a" * 64,
        observations=(),
        conformance_records=(),
        registry_sha256="b" * 64,
        live_schema_sha256="c" * 64,
        execution_requested=False,
        execution_available=False,
        open_guides=records,
    )
    assert [item["guide_id"] for item in context["open_guides"]] == [
        "saddle",
        "structure",
    ]
    bare = _public_context(
        task="t",
        task_spec_sha256="a" * 64,
        observations=(),
        conformance_records=(),
        registry_sha256="b" * 64,
        live_schema_sha256="c" * 64,
        execution_requested=False,
        execution_available=False,
    )
    assert "open_guides" not in bare


def test_a_leaf_tool_called_by_name_returns_its_guide_with_the_body(
    tmp_path,
):
    host = _host(tmp_path)
    host._inspect_database_records = lambda turn_id, values: {"records": []}
    assert "database" not in host.active_guides
    reply = host.dispatch(
        turn_id="t1",
        tool_name="inspect_database_records",
        arguments={"database_artifact_id": "db-1"},
    )
    (opened,) = reply["guides_opened"]
    assert opened["guide_id"] == "database" and opened["body"]
    assert "database" in host.active_guides
