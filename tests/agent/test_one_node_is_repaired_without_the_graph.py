"""Repairing how one node is written must not cost the other eight.

A DAG arrives in a single payload -- median about thirteen kilobytes over nine
nodes -- and the validator used to reject it whole. Roughly three quarters of
recorded rejections were single-node repairs: an identifier's case, a missing
unit, a declared kind the operation does not derive, a selector the result does
not resolve. Each of those was priced at re-emitting the entire graph.

These tests execute the amendment against a real planned toolchain: a repair
lands on the node it names, every other node survives byte for byte, the plan is
appended under a new digest rather than mutated, and an amendment that would
change the science rather than its expression is refused.
"""

from __future__ import annotations

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.scientific_toolchain import (
    AnalysisNodeIntentV1,
    AnalysisOutputIntentV1,
    AnalysisSelectorIntentV1,
    RegisteredResultInputIntentV1,
    build_scientific_toolchain_plan,
)
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1


def _extraction_node(node_id, *, output_id, unit="hartree"):
    return AnalysisNodeIntentV1(
        node_id=node_id,
        analysis_kind="result_extraction",
        dependencies=(),
        inputs=(
            RegisteredResultInputIntentV1(
                input_id="registered-result",
                artifact_id=f"artifact-{node_id}",
            ),
        ),
        selectors=(
            AnalysisSelectorIntentV1(
                quantity_id=f"{node_id}-e", selector="energy"
            ),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id=output_id,
                quantity_kind="electronic_energy",
                unit=unit,
            ),
        ),
        expression_nodes=(),
        expression_output_node_ids=(),
        temperature_k=None,
        pressure_atm=None,
        support_state="planned",
        blocked_reason="",
    )


@pytest.fixture()
def host_with_plan():
    """A host holding one planned analysis-only toolchain."""

    host = CommandCompiledToolHostV1.__new__(CommandCompiledToolHostV1)
    nodes = (
        _extraction_node("extract-a", output_id="e-a"),
        _extraction_node("extract-b", output_id="e-b"),
        _extraction_node("extract-c", output_id="e-c"),
    )
    plan = build_scientific_toolchain_plan(
        plan_id="branching",
        workflow_id="branching",
        command_workflow_draft_sha256="d" * 64,
        calculation_nodes=(),
        calculation_observables={},
        analysis_nodes=nodes,
        required_output_ids=("e-a", "e-b", "e-c"),
    )
    host.scientific_toolchain_plans = {plan.plan_sha256: plan}
    host._scientific_toolchain_command_results = {
        plan.plan_sha256: {
            "workflow_draft": _EmptyDraft(),
            "actionable_node_ids": (),
            "unresolved_node_ids": (),
        }
    }
    host.frozen_workflow_approval = None
    host._bind_program_toolchain = lambda *args, **kwargs: None
    return host, plan


class _EmptyDraft:
    workflow_id = "branching"
    task_spec_id = ""
    nodes = ()
    draft_sha256 = "d" * 64


def _amend(host, **values):
    return host._amend_scientific_workflow(
        "turn-1", {"workflow_id": "branching", **values}
    )


def test_a_unit_repair_lands_on_its_node_and_leaves_the_rest_alone(
    host_with_plan,
):
    host, original = host_with_plan

    result = _amend(
        host,
        analysis_repairs=[
            {"node_id": "extract-b", "outputs": [{"output_id": "e-b", "unit": "1"}]}
        ],
    )

    revised = result["scientific_toolchain_plan"]
    by_id = {node.node_id: node for node in revised.analysis_nodes}
    assert by_id["extract-b"].outputs[0].unit == "1"
    # The nodes not named survive exactly.
    for untouched in ("extract-a", "extract-c"):
        before = next(
            node
            for node in original.analysis_nodes
            if node.node_id == untouched
        )
        assert by_id[untouched] == before


def test_the_revision_is_appended_not_mutated(host_with_plan):
    """Provenance: the previous revision stays addressable by its own digest."""

    host, original = host_with_plan

    result = _amend(
        host,
        analysis_repairs=[
            {
                "node_id": "extract-a",
                "outputs": [{"output_id": "e-a", "new_output_id": "e-alpha"}],
            }
        ],
    )

    revised = result["scientific_toolchain_plan"]
    assert revised.plan_sha256 != original.plan_sha256
    assert host.scientific_toolchain_plans[original.plan_sha256] is original
    assert host.scientific_toolchain_plans[revised.plan_sha256] is revised


def test_a_selector_repair_replaces_only_the_selector(host_with_plan):
    host, _ = host_with_plan

    result = _amend(
        host,
        analysis_repairs=[
            {
                "node_id": "extract-c",
                "selectors": [
                    {"quantity_id": "extract-c-e", "selector": "scf_energy"}
                ],
            }
        ],
    )

    node = next(
        item
        for item in result["scientific_toolchain_plan"].analysis_nodes
        if item.node_id == "extract-c"
    )
    assert node.selectors[0].selector == "scf_energy"
    assert node.selectors[0].quantity_id == "extract-c-e"


def test_an_empty_amendment_is_refused(host_with_plan):
    host, _ = host_with_plan
    with pytest.raises(ContractError, match="changes something"):
        _amend(host)


def test_addressing_an_absent_node_names_what_exists(host_with_plan):
    host, _ = host_with_plan
    with pytest.raises(ContractError, match="extract-a"):
        _amend(host, analysis_repairs=[{"node_id": "extract-z"}])


def test_addressing_an_absent_output_names_what_the_node_declares(
    host_with_plan,
):
    host, _ = host_with_plan
    with pytest.raises(ContractError, match="e-a"):
        _amend(
            host,
            analysis_repairs=[
                {
                    "node_id": "extract-a",
                    "outputs": [{"output_id": "not-there", "unit": "1"}],
                }
            ],
        )


def test_a_registered_result_input_cannot_be_repointed(host_with_plan):
    """Swapping which result is read is a different question, not a repair."""

    host, _ = host_with_plan
    with pytest.raises(ContractError, match="new workflow"):
        _amend(
            host,
            analysis_repairs=[
                {
                    "node_id": "extract-a",
                    "inputs": [
                        {
                            "input_id": "registered-result",
                            "producer_output_id": "something-else",
                        }
                    ],
                }
            ],
        )


def test_an_approved_workflow_cannot_be_amended_in_place(host_with_plan):
    """An amended plan is a new workflow and needs its own review."""

    host, original = host_with_plan

    class _Frozen:
        plan_sha256 = original.plan_sha256

    host.frozen_workflow_approval = _Frozen()

    with pytest.raises(ContractError, match="its own review"):
        _amend(
            host,
            analysis_repairs=[
                {
                    "node_id": "extract-a",
                    "outputs": [{"output_id": "e-a", "unit": "1"}],
                }
            ],
        )
