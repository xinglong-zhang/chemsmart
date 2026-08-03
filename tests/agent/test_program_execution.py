from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from chemsmart.agent._contracts import (
    ContractError,
    TrustedArtifactRefV1,
    file_sha256,
)
from chemsmart.agent.execution import (
    build_execution_resource_spec,
    build_locked_pyscf_sp_opt_hess_approval,
    build_program_execution_invocation,
    build_program_execution_receipt,
    handoff_optimized_pyscf_geometry,
    workflow_execution_approval_json,
)
from chemsmart.agent.live_session import _scan_xyz_artifacts, _system_prompt
from chemsmart.agent.tool_runtime import (
    CommandCompiledToolHostV1,
    _prepare_execution_node_workspace,
)
from chemsmart.jobs.pyscf.writer import write_pyscf_h5


def _artifact(path, *, artifact_id, kind):
    return TrustedArtifactRefV1(
        artifact_id=artifact_id,
        kind=kind,
        sha256=file_sha256(path),
        size_bytes=path.stat().st_size,
        path=str(path.resolve()),
        cli_value=str(path.resolve()),
    )


def _locked_approval(tmp_path):
    geometry_path = tmp_path / "water.xyz"
    geometry_path.write_text(
        "3\nwater\nO 0.0 0.0 0.1174\n"
        "H -0.757 0.0 -0.4696\nH 0.757 0.0 -0.4696\n",
        encoding="utf-8",
    )
    project_path = tmp_path / "water-pyscf.yaml"
    shared = {
        "functional": "b3lyp",
        "basis": "def2-svp",
        "defgrid": "defgrid2",
        "density_fit": False,
        "scf_tol": 1e-9,
        "scf_maxiter": 100,
    }
    import yaml

    project_path.write_text(
        yaml.safe_dump(
            {
                "sp": dict(shared),
                "opt": {
                    **shared,
                    "opt_solver": "geometric",
                    "opt_maxsteps": 100,
                },
                "hess": dict(shared),
            },
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    resources = build_execution_resource_spec(
        execution_target="run",
        cores=4,
        memory_gb=4,
        gpu_count=0,
        scratch_policy="none",
        node_timeout_seconds=5400,
    )
    approval = build_locked_pyscf_sp_opt_hess_approval(
        approval_id="water-approval",
        workflow_id="water-pyscf-workflow",
        task_spec_sha256="a" * 64,
        approved_workspace=tmp_path,
        resources=resources,
        initial_artifact=_artifact(
            geometry_path,
            artifact_id="geometry.water.initial",
            kind="geometry_xyz",
        ),
        project_artifact=_artifact(
            project_path,
            artifact_id="project.water.pyscf",
            kind="project_yaml",
        ),
    )
    return approval


def test_locked_pyscf_approval_binds_order_and_optimized_handoff(tmp_path):
    approval = _locked_approval(tmp_path)

    assert tuple(item.node_id for item in approval.node_bindings) == (
        "sp-initial",
        "opt-initial",
        "hess-optimized",
    )
    assert approval.node_bindings[2].input_mode == "producer"
    assert approval.producer_edges[0].producer_node_id == "opt-initial"
    assert approval.producer_edges[0].consumer_node_id == "hess-optimized"
    envelope = json.loads(workflow_execution_approval_json(approval))
    assert envelope["workflow_approval"]["approval_sha256"] == (
        approval.approval_sha256
    )


def test_execution_host_rejects_a_node_before_approved_predecessors(tmp_path):
    host = object.__new__(CommandCompiledToolHostV1)
    host.surface = SimpleNamespace(
        profile="command_compiled_approved_execution"
    )
    host.execution_receipts = {}
    host.workflow_execution_approval = _locked_approval(tmp_path)

    with pytest.raises(ContractError, match="approved predecessors"):
        host._execute_approved_program_node(
            "turn-1", {"node_id": "opt-initial"}
        )


def test_workspace_scan_excludes_host_generated_geometry(tmp_path):
    source = tmp_path / "source.xyz"
    source.write_text("1\nhydrogen\nH 0 0 0\n", encoding="utf-8")
    for directory in (tmp_path / "nodes", tmp_path / "artifacts"):
        directory.mkdir()
        (directory / "derived.xyz").write_text(
            "1\nhydrogen\nH 0 0 0\n", encoding="utf-8"
        )

    observations = _scan_xyz_artifacts(tmp_path)

    assert len(observations) == 1
    assert observations[0].artifact.sha256 == file_sha256(source)


def test_preexisting_node_outputs_block_a_second_launch(tmp_path):
    node_workspace = tmp_path / "nodes" / "sp"
    node_workspace.mkdir(parents=True)
    (node_workspace / "result.out").write_text("existing", encoding="utf-8")

    with pytest.raises(ContractError, match="already contains outputs"):
        _prepare_execution_node_workspace(node_workspace)


def test_validated_opt_handoff_materializes_exact_final_geometry(tmp_path):
    approval = _locked_approval(tmp_path)
    resources = build_execution_resource_spec(
        execution_target="run",
        cores=4,
        memory_gb=4,
        gpu_count=0,
        scratch_policy="none",
        node_timeout_seconds=5400,
    )
    project = _artifact(
        tmp_path / "water-pyscf.yaml",
        artifact_id="project.water.pyscf",
        kind="project_yaml",
    )
    initial = _artifact(
        tmp_path / "water.xyz",
        artifact_id="geometry.water.initial",
        kind="geometry_xyz",
    )
    opt_node = approval.node("opt-initial")
    invocation = build_program_execution_invocation(
        node_id=opt_node.node_id,
        approval=approval,
        project_artifact=project,
        input_artifact=initial,
        scientific_identity_sha256=opt_node.scientific_identity_sha256,
        environment_receipt_sha256="b" * 64,
        resources=resources,
        argv=("chemsmart", "run", "pyscf", "opt"),
    )
    final_positions = np.asarray(
        [
            [0.0, 0.0, 0.1],
            [-0.76, 0.0, -0.47],
            [0.76, 0.0, -0.47],
        ],
        dtype=float,
    )
    h5_path = tmp_path / "water-opt.h5"
    write_pyscf_h5(
        h5_path,
        spec={
            "program": "pyscf",
            "jobtype": "opt",
            "engine": "cpu",
            "stages": ["scf", "opt"],
            "symbols": ["O", "H", "H"],
            "unit": "Angstrom",
            "charge": 0,
            "multiplicity": 1,
            "requested_settings_sha256": "d" * 64,
            "applied_settings_sha256": "d" * 64,
        },
        provenance={"environment_receipt_sha256": "b" * 64},
        status={
            "normal_termination": True,
            "stages": {
                "scf": {"converged": True},
                "opt": {"converged": True},
            },
        },
        results={
            "energies": np.asarray([-76.35]),
            "positions": final_positions,
        },
    )
    result = _artifact(
        h5_path,
        artifact_id="result.opt.hdf5",
        kind="pyscf_hdf5",
    )
    with pytest.raises(ContractError, match="unresolved findings"):
        build_program_execution_receipt(
            invocation,
            execution_state="engine_complete",
            exit_status=0,
            engine_complete=True,
            validated=True,
            output_artifacts=(result,),
            validator_receipt_sha256s=("e" * 64,),
            findings=("pyscf.result.seeded_failure",),
            started_at="2026-08-03T00:00:00+00:00",
            finished_at="2026-08-03T00:00:01+00:00",
        )
    receipt = build_program_execution_receipt(
        invocation,
        execution_state="engine_complete",
        exit_status=0,
        engine_complete=True,
        validated=True,
        output_artifacts=(result,),
        validator_receipt_sha256s=("e" * 64,),
        started_at="2026-08-03T00:00:00+00:00",
        finished_at="2026-08-03T00:00:01+00:00",
    )

    geometry, handoff = handoff_optimized_pyscf_geometry(
        producer_receipt=receipt,
        result_artifact=result,
        producer_edge=approval.producer_edges[0],
        approved_workspace=tmp_path,
        geometry_artifact_id="geometry.water.optimized",
        expected_charge=0,
        expected_multiplicity=1,
    )

    assert handoff.status == "validated_handoff"
    assert handoff.symbols == ("O", "H", "H")
    assert geometry.sha256 == file_sha256(geometry.path)
    assert "source_sha256=" + result.sha256 in Path(
        geometry.path
    ).read_text(encoding="utf-8")


def test_live_prompt_preserves_order_and_loader_bounded_alternatives():
    prompt = _system_prompt(None)

    assert "Preserve every user-declared stage order" in prompt
    assert "current project schema and loader support it" in prompt
    assert "unsupported architecture alternatives" in prompt
