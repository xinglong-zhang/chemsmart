from __future__ import annotations

import pytest

from chemsmart.agent._contracts import (
    ContractError,
    TrustedArtifactRefV1,
    file_sha256,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.analysis.result_quantities import (
    ENERGY,
    ThermochemistryReceiptV1,
    canonical_quantity_sha256,
    make_quantity_value,
)


def _host(tmp_path, *, artifacts=None):
    return CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="neutral-session"
        ),
        artifacts=artifacts or {},
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "workspace",
    )


def _calculation_only_plan():
    return {
        "workflow_id": "workflow-neutral-energy",
        "task_spec_id": "a" * 64,
        "nodes": [
            {
                "node_id": "energy-node",
                "program": "xtb",
                "jobtype": "sp",
                "project_role": "project-neutral",
                "dependencies": [],
                "inputs": [],
                "expected_outputs": [
                    {
                        "output_id": "structured-result",
                        "artifact_class": "structured_result",
                    }
                ],
                "unresolved_fields": [],
            }
        ],
    }


def test_calculation_only_plan_supports_frontier_and_prepare(tmp_path):
    host = _host(tmp_path)
    arguments = _calculation_only_plan()
    planned = host.dispatch(
        turn_id="turn-1",
        tool_name="plan_command_workflow",
        arguments=arguments,
    )["result"]

    frontier = host.dispatch(
        turn_id="turn-1",
        tool_name="inspect_workflow_frontier",
        arguments={"workflow_id": arguments["workflow_id"]},
    )["result"]
    assert set(frontier) == {
        "scientific_workflow_plan",
        "workflow_context",
    }
    assert (
        frontier["scientific_workflow_plan"]
        == planned["scientific_workflow_plan"]
    )

    prepared = host.dispatch(
        turn_id="turn-1",
        tool_name="prepare_program_node",
        arguments={
            "workflow_id": arguments["workflow_id"],
            "node_id": "energy-node",
        },
    )["result"]
    assert prepared["workflow_id"] == arguments["workflow_id"]
    assert prepared["node_id"] == "energy-node"
    assert prepared["status"] == "needs_clarification"
    assert prepared["finding"] == "program node declares no molecular input"


def _plan_registered_xyz_analysis(host, artifact_id):
    return host.dispatch(
        turn_id="turn-analysis",
        tool_name="plan_scientific_workflow",
        arguments={
            "plan_id": "registered-geometry-analysis",
            "workflow_id": "workflow-registered-geometry",
            "task_spec_id": "a" * 64,
            "calculation_nodes": [],
            "analysis_nodes": [
                {
                    "node_id": "extract-positions",
                    "analysis_kind": "result_extraction",
                    "dependencies": [],
                    "artifact_id": artifact_id,
                    "inputs": [],
                    "selectors": [
                        {
                            "quantity_id": "positions",
                            "selector": "positions",
                        }
                    ],
                    "outputs": [
                        {
                            "output_id": "positions",
                            "quantity_kind": "positions",
                            "unit": "angstrom",
                        }
                    ],
                    "expression_nodes": [],
                    "expression_output_node_ids": [],
                    "support_state": "planned",
                    "blocked_reason": "",
                    "validation_rules": [],
                }
            ],
            "required_output_ids": ["positions"],
        },
    )["result"]


def _plan_registered_thermochemistry(host, artifact_id):
    return host.dispatch(
        turn_id="turn-thermochemistry",
        tool_name="plan_scientific_workflow",
        arguments={
            "plan_id": "registered-thermochemistry",
            "workflow_id": "workflow-registered-thermochemistry",
            "task_spec_id": "a" * 64,
            "calculation_nodes": [],
            "analysis_nodes": [
                {
                    "node_id": "derive-thermochemistry",
                    "analysis_kind": "thermochemistry",
                    "dependencies": [],
                    "artifact_id": artifact_id,
                    "inputs": [],
                    "selectors": [],
                    "outputs": [
                        {
                            "output_id": "gibbs-free-energy",
                            "quantity_kind": "gibbs_free_energy",
                            "unit": "Eh",
                        }
                    ],
                    "expression_nodes": [],
                    "expression_output_node_ids": [],
                    "temperature_k": 298.15,
                    "pressure_atm": 1.0,
                    "support_state": "planned",
                    "blocked_reason": "",
                    "validation_rules": [],
                }
            ],
            "required_output_ids": ["gibbs-free-energy"],
        },
    )["result"]


def test_registered_result_analysis_preserves_analysis_only_workflow(tmp_path):
    xyz = tmp_path / "water.xyz"
    xyz.write_text(
        "3\nwater\nO 0 0 0\nH 0.76 0.58 0\nH -0.76 0.58 0\n",
        encoding="utf-8",
    )
    artifact = TrustedArtifactRefV1(
        artifact_id="registered-water-geometry",
        kind="geometry_xyz",
        sha256=file_sha256(xyz),
        size_bytes=xyz.stat().st_size,
        path=str(xyz.resolve()),
        cli_value=str(xyz.resolve()),
    )
    host = _host(tmp_path, artifacts={artifact.artifact_id: artifact})

    planned = _plan_registered_xyz_analysis(host, artifact.artifact_id)
    inspected = host.dispatch(
        turn_id="turn-analysis",
        tool_name="inspect_workflow_frontier",
        arguments={"workflow_id": "workflow-registered-geometry"},
    )["result"]

    assert planned["calculation_plan"]["scientific_workflow_plan"] is None
    assert planned["scientific_toolchain_plan"]["calculation_node_ids"] == []
    assert inspected["scientific_toolchain_plan"] == (
        planned["scientific_toolchain_plan"]
    )
    assert inspected["workflow_frontier"]["actionable_node_ids"] == [
        "extract-positions"
    ]


def test_workflow_resolver_rejects_equal_but_unbound_registry_copy(tmp_path):
    """A copied result cannot become an alternative registry authority."""

    xyz = tmp_path / "water.xyz"
    xyz.write_text("1\noxygen\nO 0 0 0\n", encoding="utf-8")
    artifact = TrustedArtifactRefV1(
        artifact_id="registered-oxygen-geometry",
        kind="geometry_xyz",
        sha256=file_sha256(xyz),
        size_bytes=xyz.stat().st_size,
        path=str(xyz.resolve()),
        cli_value=str(xyz.resolve()),
    )
    host = _host(tmp_path, artifacts={artifact.artifact_id: artifact})
    planned = _plan_registered_xyz_analysis(host, artifact.artifact_id)
    plan_sha256 = planned["scientific_toolchain_plan"]["plan_sha256"]
    host._scientific_toolchain_command_results[plan_sha256] = dict(
        host._scientific_toolchain_command_results[plan_sha256]
    )

    with pytest.raises(
        ContractError,
        match="lost its exact scientific toolchain binding",
    ):
        host.dispatch(
            turn_id="turn-analysis",
            tool_name="inspect_workflow_frontier",
            arguments={"workflow_id": "workflow-registered-geometry"},
        )


def test_registered_program_result_is_authoritative_thermochemistry_root(
    tmp_path,
):
    output = tmp_path / "frequency.out"
    output.write_text("registered ORCA result\n", encoding="utf-8")
    artifact = TrustedArtifactRefV1(
        artifact_id="registered-orca-frequency",
        kind="orca_output",
        sha256=file_sha256(output),
        size_bytes=output.stat().st_size,
        path=str(output.resolve()),
        cli_value=str(output.resolve()),
    )
    host = _host(tmp_path, artifacts={artifact.artifact_id: artifact})
    planned = _plan_registered_thermochemistry(host, artifact.artifact_id)

    plan = planned["scientific_toolchain_plan"]
    assert plan["analysis_nodes"][0]["inputs"] == [
        {
            "input_id": "registered-result",
            "artifact_id": artifact.artifact_id,
            "source_kind": "registered_result",
        }
    ]
    assert planned["workflow_frontier"]["actionable_node_ids"] == [
        "derive-thermochemistry"
    ]

    quantity = make_quantity_value(
        quantity_id="gibbs_free_energy",
        source_value=-75.0,
        source_unit="Eh",
        value=-75.0,
        unit="Eh",
        dimension=ENERGY,
        evidence_ref=f"artifact:{artifact.artifact_id}#{artifact.sha256}",
    )
    body = {
        "schema_version": "chemsmart.thermochemistry-receipt.v1",
        "artifact_id": artifact.artifact_id,
        "artifact_sha256": artifact.sha256,
        "program": "orca",
        "engine_id": "orca-6.1.1",
        "temperature_k": 298.15,
        "pressure_atm": 1.0,
        "quantities": (quantity,),
        "assumptions": ("harmonic RRHO",),
        "status": "derived",
        "concentration_mol_l": None,
        "entropy_method": "rrho",
        "entropy_cutoff_cm1": None,
        "enthalpy_cutoff_cm1": None,
        "alpha": 4,
        "use_weighted_mass": False,
        "frequency_scale_factor": 1.0,
    }
    receipt = ThermochemistryReceiptV1(
        **body,
        receipt_sha256=canonical_quantity_sha256(body),
    )
    host.thermochemistry_receipts[receipt.receipt_sha256] = receipt

    inspected = host.dispatch(
        turn_id="turn-thermochemistry",
        tool_name="inspect_workflow_frontier",
        arguments={"workflow_id": "workflow-registered-thermochemistry"},
    )["result"]
    assert inspected["workflow_frontier"]["completed_node_ids"] == [
        "derive-thermochemistry"
    ]


def test_registered_geometry_cannot_be_a_thermochemistry_root(tmp_path):
    xyz = tmp_path / "water.xyz"
    xyz.write_text(
        "3\nwater\nO 0 0 0\nH 0.76 0.58 0\nH -0.76 0.58 0\n",
        encoding="utf-8",
    )
    artifact = TrustedArtifactRefV1(
        artifact_id="registered-water-geometry",
        kind="geometry_xyz",
        sha256=file_sha256(xyz),
        size_bytes=xyz.stat().st_size,
        path=str(xyz.resolve()),
        cli_value=str(xyz.resolve()),
    )
    host = _host(tmp_path, artifacts={artifact.artifact_id: artifact})

    with pytest.raises(
        ContractError,
        match="geometry-only registered artifact",
    ):
        _plan_registered_thermochemistry(host, artifact.artifact_id)
