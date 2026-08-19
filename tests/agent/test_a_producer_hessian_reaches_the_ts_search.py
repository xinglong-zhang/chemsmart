"""A cheap Hessian could not start a TS search; now it can.

The producer-edge vocabulary held exactly two rules, and the only Hessian
rule required a converged final transition state with exactly one
imaginary mode -- useful for feeding an IRC, useless for STARTING a TS
search from a minimum's curvature. The live CLI already exposes
run/orca/ts --inhess-filename; what was missing was the typed edge. The
third selection rule carries any frequency-bearing ORCA producer's
validated Hessian to a TS consumer's inhess_filename role, recording the
observed imaginary-mode count as a fact instead of enforcing it as a
requirement.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from chemsmart.agent._contracts import (
    ContractError,
    TrustedArtifactRefV1,
    file_sha256,
)
from chemsmart.agent.execution import (
    ORCAProducerHessianHandoffV1,
    _frozen_producer_edge_rule,
    build_producer_edge_rule,
    handoff_final_orca_ts_hessian,
    handoff_validated_orca_producer_hessian,
    is_validated_producer_orca_hessian_edge,
)
from chemsmart.agent.workflows import (
    ScientificWorkflowEdgeV2,
    ScientificWorkflowNodeV2,
    build_scientific_workflow_plan,
)

_OUT = Path("tests/data/ORCATests/outputs/water_opt.out").resolve()
_HESS = Path("tests/data/ORCATests/outputs/water_opt.hess").resolve()


def _node(node_id, stage):
    return ScientificWorkflowNodeV2(
        node_id=node_id,
        stage=stage,
        requested_program="orca",
        program="orca",
        engine="cpu",
        project_role=f"{node_id}-project",
        unresolved_fields=(),
    )


def _plan_with_inhess_edge():
    edge = ScientificWorkflowEdgeV2(
        edge_id="freq-to-ts",
        source_node_id="freq-min",
        target_node_id="ts-search",
        edge_kind="data",
        artifact_class="orca_hessian",
        producer_output_id="hessian",
        consumer_input_id="inhess_filename",
    )
    plan = build_scientific_workflow_plan(
        workflow_id="ts-from-min-hessian",
        task_spec_sha256="a" * 64,
        scientific_identity_sha256="b" * 64,
        nodes=(_node("freq-min", "opt"), _node("ts-search", "ts")),
        edges=(edge,),
    )
    return plan, edge


def test_the_hessian_edge_names_its_own_selection_rule():
    plan, edge = _plan_with_inhess_edge()

    assert is_validated_producer_orca_hessian_edge(plan, edge)
    rule = _frozen_producer_edge_rule(
        plan, edge, environment_receipt_sha256="c" * 64
    )
    assert rule.selection_rule == "validated_producer_orca_hessian"

    stray = ScientificWorkflowEdgeV2(
        edge_id="stray",
        source_node_id="freq-min",
        target_node_id="ts-search",
        edge_kind="data",
        artifact_class="orca_hessian",
        producer_output_id="hessian",
        consumer_input_id="mystery_option",
    )
    with pytest.raises(ContractError) as refusal:
        _frozen_producer_edge_rule(
            plan, stray, environment_receipt_sha256="c" * 64
        )
    message = str(refusal.value)
    assert "inhess_filename" in message
    assert "inhess_filename" in message


def _artifact(path, artifact_id, kind):
    return TrustedArtifactRefV1(
        artifact_id=artifact_id,
        kind=kind,
        sha256=file_sha256(path),
        size_bytes=path.stat().st_size,
        path=str(path),
        cli_value=str(path),
    )


def _receipt(result, hessian):
    return SimpleNamespace(
        validated=True,
        node_id="freq-min",
        receipt_sha256="d" * 64,
        output_artifacts=(
            SimpleNamespace(
                artifact_id=result.artifact_id, sha256=result.sha256
            ),
            SimpleNamespace(
                artifact_id=hessian.artifact_id, sha256=hessian.sha256
            ),
        ),
    )


def test_a_minimum_hessian_is_not_required_to_be_a_ts(tmp_path):
    result = _artifact(_OUT, "result.freq-min.1", "orca_output")
    hessian = _artifact(_HESS, "result.freq-min.2", "orca_hessian")
    receipt = _receipt(result, hessian)
    producer_edge = build_producer_edge_rule(
        producer_node_id="freq-min",
        consumer_node_id="ts-search",
        artifact_kind="orca_hessian",
        selection_rule="validated_producer_orca_hessian",
    )
    workspace = tmp_path / "workspace"
    workspace.mkdir()

    artifact, observed = handoff_validated_orca_producer_hessian(
        producer_receipt=receipt,
        result_artifact=result,
        hessian_candidates=(hessian,),
        producer_edge=producer_edge,
        approved_workspace=workspace,
        hessian_artifact_id="hessian.freq-min-to-ts-search",
        expected_charge=0,
        expected_multiplicity=1,
    )

    assert isinstance(observed, ORCAProducerHessianHandoffV1)
    assert observed.observed_imaginary_mode_count == 0
    assert observed.status == "validated_handoff"
    assert Path(artifact.path).exists()
    assert "any imaginary-mode count" in (
        handoff_validated_orca_producer_hessian.__doc__ or ""
    )

    # The final-TS handoff still refuses a minimum: different question.
    ts_edge = build_producer_edge_rule(
        producer_node_id="freq-min",
        consumer_node_id="irc-node",
        artifact_kind="orca_hessian",
        selection_rule="validated_final_orca_ts_hessian",
    )
    with pytest.raises(ContractError):
        handoff_final_orca_ts_hessian(
            producer_receipt=receipt,
            result_artifact=result,
            hessian_candidates=(hessian,),
            producer_edge=ts_edge,
            approved_workspace=tmp_path / "other",
            hessian_artifact_id="hessian.freq-min-to-irc-node",
            expected_charge=0,
            expected_multiplicity=1,
        )
