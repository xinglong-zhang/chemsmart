"""The single /approve now covers the analysis chain it displays.

Through the whole real-ability protocol the typed analysis chain lived only
in session RAM: the human approved review packets that never displayed the
analysis nodes, the approval bundle never carried them, and the executor
had nothing to walk -- every validation rule a session wrote stayed intent,
and every benchmark number was extracted outside the product. The review
and the bundle now carry the digest-bound toolchain plan, additively: a
calculation-only packet keeps its historical bytes.
"""

from __future__ import annotations

import json

from chemsmart.agent._contracts import canonical_data
from chemsmart.agent.execution import (
    approve_workflow_execution_review,
    build_workflow_execution_review,
    workflow_execution_approval_bundle_json,
    workflow_execution_review_json,
)
from chemsmart.agent.live_session import (
    load_workflow_execution_approval_bundle,
    load_workflow_execution_review,
)
from chemsmart.agent.scientific_toolchain import (
    AnalysisInputIntentV1,
    AnalysisNodeIntentV1,
    AnalysisOutputIntentV1,
    AnalysisSelectorIntentV1,
    build_scientific_toolchain_plan,
)
from chemsmart.agent.workflows import (
    ArtifactInputIntentV1,
    ArtifactOutputIntentV1,
    CommandNodeIntentV1,
)

from tests.agent.test_exact_execution_approval_chain import _review


def _toolchain(selector="energy"):
    calculation = CommandNodeIntentV1(
        node_id="sp",
        program="pyscf",
        jobtype="sp",
        project_role="water-sp",
        dependencies=(),
        inputs=(
            ArtifactInputIntentV1(
                binding_id="geometry",
                artifact_class="geometry_xyz",
                artifact_id="water",
                producer_node_id="",
                producer_output_id="",
            ),
        ),
        expected_outputs=(
            ArtifactOutputIntentV1(
                output_id="sp-out", artifact_class="pyscf_hdf5"
            ),
        ),
        unresolved_fields=(),
    )
    analysis = AnalysisNodeIntentV1(
        node_id="extract-sp",
        analysis_kind="result_extraction",
        dependencies=("sp",),
        inputs=(
            AnalysisInputIntentV1(
                input_id="raw",
                source_kind="program_output",
                producer_node_id="sp",
                producer_output_id="sp-out",
            ),
        ),
        selectors=(
            AnalysisSelectorIntentV1(quantity_id="e", selector=selector),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="e", quantity_kind="energy", unit="hartree"
            ),
        ),
        expression_nodes=(),
        expression_output_node_ids=(),
        temperature_k=None,
        pressure_atm=None,
        support_state="planned",
        blocked_reason="",
    )
    return build_scientific_toolchain_plan(
        plan_id="p",
        workflow_id="water-sp",
        command_workflow_draft_sha256="9" * 64,
        calculation_nodes=(calculation,),
        calculation_observables={"sp": ("sp-out",)},
        analysis_nodes=(analysis,),
        required_output_ids=("e",),
    )


def _review_with_chain(base, toolchain):
    return build_workflow_execution_review(
        request=base.request,
        scientific_plan=base.scientific_plan,
        materialized_workflow=base.materialized_workflow,
        execution_resources=base.execution_resources,
        execution_envelope=base.execution_envelope,
        environment_bindings=base.environment_bindings,
        node_reviews=base.node_reviews,
        stationary_point_policy=base.stationary_point_policy,
        non_executable_node_ids=base.non_executable_node_ids,
        scientific_toolchain_plan=toolchain,
    )


def test_the_chain_is_in_the_displayed_bytes_and_the_bundle(tmp_path):
    toolchain = _toolchain()
    base = _review(tmp_path)
    review = _review_with_chain(base, toolchain)

    # Displayed bytes: the packet a human reads carries the analysis node.
    packet = workflow_execution_review_json(review)
    assert "extract-sp" in packet
    path = tmp_path / "review.json"
    path.write_text(packet, encoding="utf-8")
    assert load_workflow_execution_review(path) == review

    # A different analysis chain is a different reviewed digest.
    other = _review_with_chain(base, _toolchain(selector="scf_energy"))
    assert other.review_sha256 != review.review_sha256

    bundle = approve_workflow_execution_review(
        review,
        approval_id="approval-water-sp",
        approved_review_sha256=review.review_sha256,
        actor="human",
        resolution_id="resolution-water-sp",
    )
    assert bundle.scientific_toolchain_plan is not None
    assert (
        bundle.scientific_toolchain_plan.plan_sha256 == toolchain.plan_sha256
    )
    assert (
        bundle.frozen_workflow_approval.scientific_toolchain_plan_sha256
        == toolchain.plan_sha256
    )
    bundle_path = tmp_path / "bundle.json"
    bundle_path.write_text(
        workflow_execution_approval_bundle_json(bundle), encoding="utf-8"
    )
    assert load_workflow_execution_approval_bundle(bundle_path) == bundle
