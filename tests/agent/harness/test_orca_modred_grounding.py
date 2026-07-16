from __future__ import annotations

from chemsmart.agent.harness.command_semantics import _input_excerpt
from chemsmart.agent.harness.generated_invariants import (
    check_generated_input_invariants,
    electron_multiplicity_evidence,
)
from chemsmart.agent.harness.intent import IntentSpec, evaluate_intent


def test_freeze_bond_intent_maps_to_orca_modred_coordinates() -> None:
    request = (
        "Use ORCA to optimize probe.xyz while keeping bond 2-5 frozen, "
        "charge 0 multiplicity 1."
    )
    expected = IntentSpec.from_request(request)

    assert expected.kind == "orca.modred"
    assert expected.chemistry == {"coordinates": [[2, 5]]}
    assert (
        evaluate_intent(
            "chemsmart run orca -f probe.xyz -c 0 -m 1 "
            "modred --coordinates '[[2,5]]'",
            expected,
        ).verdict
        == "ok"
    )


def test_large_orca_modred_excerpt_preserves_header_constraint() -> None:
    content = (
        "! Opt B3LYP def2-SVP\n"
        "%geom\n  Constraints\n  {B 1 4 C}\n  end\nend\n"
        "* xyz 0 1\n"
        + "C 0.0 0.0 0.0\n" * 300
        + "*\n"
    )
    excerpt = _input_excerpt(content)

    assert "%geom" in excerpt
    assert "{B 1 4 C}" in excerpt
    assert excerpt.endswith("*\n")
    assert (
        check_generated_input_invariants(
            "chemsmart run orca -p demo -f probe.xyz -c 0 -m 1 "
            "modred --coordinates '[[2,5]]'",
            [
                {
                    "path": "probe_modred.inp",
                    "route": "! Opt B3LYP def2-SVP",
                    "content_tail": excerpt,
                }
            ],
        )
        == ()
    )


def test_intent_rejects_regrouped_modred_coordinates() -> None:
    expected = {
        "action": "run",
        "program": "orca",
        "kind": "orca.modred",
        "input_path": "probe.xyz",
        "charge": 0,
        "multiplicity": 1,
        "chemistry": {"coordinates": [[1, 2], [3, 4]]},
    }

    result = evaluate_intent(
        "chemsmart run orca -f probe.xyz -c 0 -m 1 "
        "modred --coordinates '[[1,2,3,4]]'",
        expected,
    )

    assert result.verdict == "reject"
    assert "intent.chemistry.coordinates" in result.failed_rule_ids


def test_modred_generated_input_rejects_wrong_route_and_extra_constraint() -> None:
    issues = check_generated_input_invariants(
        "chemsmart run orca -p demo -f probe.xyz -c 0 -m 1 "
        "modred --coordinates '[[2,5]]'",
        [
            {
                "path": "probe_modred.inp",
                "route": "! OptTS B3LYP def2-SVP",
                "content_tail": (
                    "! OptTS B3LYP def2-SVP\n"
                    "%geom\n  Constraints\n  {B 1 4 C}\n  {A 0 1 2 C}\n"
                    "  end\nend\n* xyz 0 1\nC 0 0 0\n*\n"
                ),
            }
        ],
    )
    rule_ids = {issue.rule_id for issue in issues}

    assert "input.orca.modred.route_opt" in rule_ids
    assert "input.orca.modred.forbidden_route" in rule_ids
    assert "input.orca.modred.constraint_set" in rule_ids


def test_modred_generated_input_rejects_non_internal_constraint_rows() -> None:
    for extra in ("{ C 7 C }", "{ B 1 4 1.25 C }"):
        issues = check_generated_input_invariants(
            "chemsmart run orca -p demo -f probe.xyz -c 0 -m 1 "
            "modred --coordinates '[[2,5]]'",
            [
                {
                    "path": "probe_modred.inp",
                    "route": "! Opt B3LYP def2-SVP",
                    "content_tail": (
                        "! Opt B3LYP def2-SVP\n"
                        "%geom\n  Constraints\n  {B 1 4 C}\n"
                        f"  {extra}\n  end\nend\n"
                        "* xyz 0 1\nC 0 0 0\n*\n"
                    ),
                }
            ],
        )

        assert "input.orca.modred.constraint_set" in {
            issue.rule_id for issue in issues
        }


def test_generated_input_rejects_electron_multiplicity_parity_mismatch() -> None:
    generated = {
        "path": "heme_sp.inp",
        "route": "! B3LYP def2-SVP",
        "content_tail": "! B3LYP def2-SVP\n* xyz 0 3\n*\n",
        "charge": 0,
        "multiplicity": 3,
        "element_counts": {"C": 23, "H": 16, "Fe": 1, "N": 7, "O": 1},
    }

    issues = check_generated_input_invariants(
        "chemsmart run orca -f heme.xyz -c 0 -m 3 sp",
        [generated],
    )

    assert "input.state.electron_multiplicity_parity" in {
        issue.rule_id for issue in issues
    }


def test_generated_input_accepts_electron_multiplicity_parity_match() -> None:
    generated = {
        "path": "heme_sp.inp",
        "route": "! B3LYP def2-SVP",
        "content_tail": "! B3LYP def2-SVP\n* xyz 0 2\n*\n",
        "charge": 0,
        "multiplicity": 2,
        "element_counts": {"C": 23, "H": 16, "Fe": 1, "N": 7, "O": 1},
    }

    assert (
        check_generated_input_invariants(
            "chemsmart run orca -f heme.xyz -c 0 -m 2 sp",
            [generated],
        )
        == ()
    )
    assert electron_multiplicity_evidence(generated) == {
        "charge": 0,
        "multiplicity": 2,
        "element_counts": {"C": 23, "Fe": 1, "H": 16, "N": 7, "O": 1},
        "electron_count": 237,
        "valid": True,
    }
