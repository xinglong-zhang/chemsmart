"""A program family is a leaf, and the program itself is its signal.

PySCF's affordances and limits -- one structure per result, a hess that
moves nothing and an opt that prints no frequencies, an all-real
spectrum that is not stationarity, a functional name that is not a
functional -- lived nowhere the model reads. A leaf cannot key on a
jobtype (``hess`` is xTB's too), so ``GuideV1`` gains ``programs`` and
``guides_from_plan`` opens the leaf when the plan names the program.
"""

from __future__ import annotations

import pytest

from chemsmart.agent.guides import (
    GUIDES_BY_ID,
    guides_from_plan,
    guides_from_text,
    guides_from_workspace,
)
from chemsmart.agent.rules import rules_for

pytestmark = pytest.mark.capability("guide:pyscf")


def test_the_leaf_opens_on_the_program_the_plan_names():
    assert "pyscf" in guides_from_plan(programs=("pyscf",))
    assert "pyscf" in guides_from_plan(programs=("orca", "pyscf"))
    assert "pyscf" not in guides_from_plan(programs=("orca",))
    # An xTB Hessian must not open the PySCF leaf.
    assert "pyscf" not in guides_from_plan(
        jobtypes=("hess",), programs=("xtb",)
    )


def test_the_leaf_opens_on_the_text_and_the_workspace():
    assert "pyscf" in guides_from_text("run a PySCF single point")
    assert "pyscf" in guides_from_workspace(("pyscf_hdf5",))
    assert "pyscf" not in guides_from_workspace(("orca_output",))


@pytest.mark.capability("rule:leaf.pyscf.two_nodes_make_a_minimum")
@pytest.mark.capability("rule:leaf.pyscf.one_structure_per_result")
@pytest.mark.capability(
    "rule:leaf.pyscf.a_matching_name_is_not_a_matching_functional"
)
@pytest.mark.capability(
    "rule:leaf.pyscf.no_imaginary_mode_is_not_a_stationary_point"
)
@pytest.mark.capability("rule:leaf.pyscf.a_root_is_an_index_not_an_identity")
@pytest.mark.capability(
    "rule:leaf.pyscf.an_excited_minimum_has_no_hessian_here"
)
@pytest.mark.capability(
    "rule:leaf.pyscf.correlated_methods_are_ab_initio_values"
)
def test_the_leaf_says_what_pyscf_can_do_and_what_it_cannot():
    guide = GUIDES_BY_ID["pyscf"]
    body = guide.body
    for phrase in (
        "one node each: sp, opt, hess",
        "supplied_positions",
        "reached_positions (opt only)",
        "failed_wrong_stationary_point",
        "No imaginary mode means none was found",
        "b3lyp and b3lypg are one libxc functional",
        "Not available here",
        "only the geometric optimiser is installed",
        "finite differences",
        "isotope-averaged masses",
    ):
        assert phrase in body, phrase
    assert (
        guide.tools == () and guide.operations == ()
    ), "the leaf adds words, never a PySCF-shaped tool"
    placed = {rule.rule_id for rule in rules_for("leaf:pyscf")}
    assert placed == {
        "leaf.pyscf.two_nodes_make_a_minimum",
        "leaf.pyscf.one_structure_per_result",
        "leaf.pyscf.a_matching_name_is_not_a_matching_functional",
        "leaf.pyscf.no_imaginary_mode_is_not_a_stationary_point",
        "leaf.pyscf.a_root_is_an_index_not_an_identity",
        "leaf.pyscf.an_excited_minimum_has_no_hessian_here",
        "leaf.pyscf.correlated_methods_are_ab_initio_values",
    }


@pytest.mark.capability("rule:leaf.crossprogram.frozen_core_is_a_convention")
def test_the_crossprogram_leaf_names_the_frozen_core_convention():
    """A convention that makes two matching strings two calculations is
    placed where a plan naming two programs reads it, never in the stem."""

    placed = {rule.rule_id for rule in rules_for("leaf:crossprogram")}
    assert "leaf.crossprogram.frozen_core_is_a_convention" in placed
    (rule,) = [
        rule
        for rule in rules_for("leaf:crossprogram")
        if rule.rule_id == "leaf.crossprogram.frozen_core_is_a_convention"
    ]
    assert "frozen_core" in rule.text and "auto" in rule.text
    assert "5e-8 Eh" in rule.provenance
