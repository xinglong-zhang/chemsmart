"""Two host receipts disagreed about the same stage, and the model noticed.

A live session planned a complex-construction stage as declared
non-executable intent; approval readiness excluded it correctly, while the
workflow frontier -- which never consulted support state -- told the same
session "actionable / compile_and_preview" for it. The session flagged the
contradiction itself and named readiness as authoritative. The frontier now
states the same truth the readiness computes, from the same single source:
the declared set plus its data-edge cascade.
"""

from __future__ import annotations

from types import SimpleNamespace

from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.scientific_toolchain import (
    AnalysisInputIntentV1,
    AnalysisNodeIntentV1,
    AnalysisOutputIntentV1,
    AnalysisSelectorIntentV1,
    build_scientific_toolchain_plan,
    project_scientific_toolchain_frontier,
)
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.agent.workflows import (
    ArtifactInputIntentV1,
    ArtifactOutputIntentV1,
    CommandNodeIntentV1,
    ScientificWorkflowNodeV2,
)


def _calculation(node_id, artifact_id="", producer=""):
    return CommandNodeIntentV1(
        node_id=node_id,
        program="orca",
        jobtype="opt",
        project_role="r",
        dependencies=(producer,) if producer else (),
        inputs=(
            ArtifactInputIntentV1(
                binding_id="geometry",
                artifact_class="geometry_xyz",
                artifact_id=artifact_id,
                producer_node_id=producer,
                producer_output_id=f"{producer}-out" if producer else "",
            ),
        ),
        expected_outputs=(
            ArtifactOutputIntentV1(
                output_id=f"{node_id}-out", artifact_class="orca_output"
            ),
        ),
        unresolved_fields=(),
    )


def _extraction(node_id, producer):
    return AnalysisNodeIntentV1(
        node_id=node_id,
        analysis_kind="result_extraction",
        dependencies=(producer,),
        inputs=(
            AnalysisInputIntentV1(
                input_id="raw",
                source_kind="program_output",
                producer_node_id=producer,
                producer_output_id=f"{producer}-out",
            ),
        ),
        selectors=(
            AnalysisSelectorIntentV1(quantity_id="e", selector="energy"),
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


def _v2_node(node_id, support_state="resolvable", blocked_reason=""):
    return ScientificWorkflowNodeV2(
        node_id=node_id,
        stage="opt",
        requested_program="orca",
        program="orca",
        engine="cpu",
        project_role=f"{node_id}-project",
        unresolved_fields=(),
        produces_observables=("energy",),
        support_state=support_state,
        blocked_reason=blocked_reason,
    )


def test_the_frontier_states_the_same_blockage_as_the_readiness(tmp_path):
    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="s1"
        ),
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "workspace",
    )
    v2 = SimpleNamespace(
        nodes=(
            _v2_node(
                "complex-build",
                support_state="blocked_unsupported",
                blocked_reason="no affordance builds a two-molecule complex",
            ),
            _v2_node("complex-opt"),
        ),
        edges=(
            SimpleNamespace(
                edge_kind="data",
                source_node_id="complex-build",
                target_node_id="complex-opt",
                artifact_class="geometry_xyz",
            ),
        ),
    )

    reasons = host._non_executable_reasons(v2)
    assert set(reasons) == set(host._release_non_executable_node_ids(v2))
    assert (
        reasons["complex-build"]
        == "no affordance builds a two-molecule complex"
    )
    assert "complex-build" in reasons["complex-opt"]

    plan = build_scientific_toolchain_plan(
        plan_id="p",
        workflow_id="w",
        command_workflow_draft_sha256="a" * 64,
        calculation_nodes=(
            _calculation("complex-build", artifact_id="start.xyz"),
            _calculation("complex-opt", producer="complex-build"),
        ),
        calculation_observables={
            "complex-build": ("complex-build-out",),
            "complex-opt": ("complex-opt-out",),
        },
        analysis_nodes=(_extraction("extract", "complex-opt"),),
        required_output_ids=("e",),
    )
    frontier = project_scientific_toolchain_frontier(
        plan,
        actionable_calculation_node_ids=("complex-build", "complex-opt"),
        non_executable_calculation_node_ids=reasons,
    )
    by_id = {item["node_id"]: item for item in frontier["nodes"]}

    assert by_id["complex-build"]["state"] == "non_executable"
    assert (
        by_id["complex-build"]["next_action"]
        == "no affordance builds a two-molecule complex"
    )
    assert by_id["complex-opt"]["state"] == "non_executable"
    assert "complex-build" in by_id["complex-opt"]["next_action"]
    assert by_id["extract"]["state"] == "blocked_upstream"
