from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest
from ase import units

from chemsmart.agent._contracts import (
    AuxiliaryArtifactBindingV1,
    ContractError,
    TrustedArtifactRefV1,
    canonical_data,
    canonical_sha256,
    file_sha256,
)
from chemsmart.agent.execution import (
    ApprovedNodeBindingV1,
    _typed_result_validation_findings,
    build_execution_resource_spec,
    build_producer_edge_rule,
    build_program_execution_invocation,
    build_program_execution_receipt,
    build_program_result_validation_receipt,
    handoff_optimized_pyscf_geometry,
    build_workflow_execution_approval,
    workflow_execution_approval_json,
)
from chemsmart.agent.capabilities import (
    EnvironmentCapabilityReceiptV1,
    EnvironmentStatus,
    ProgramEnvironmentQueryV1,
    consume_pyscf_compute_environment_receipt,
)
from chemsmart.agent.live_session import _scan_xyz_artifacts, _system_prompt
from chemsmart.agent.runtime.events import EventKind, RuntimeEvent
from chemsmart.agent.tool_runtime import (
    CommandCompiledToolHostV1,
    _PySCFEngineObservation,
    _inspect_pyscf_engine_observation,
    _pyscf_environment_evidence,
    _prepare_execution_node_workspace,
    _staged_auxiliary_input_findings,
    _runner_defers_hessian_classification,
    _validate_stationary_point_policy_binding,
)
from chemsmart.agent.workflows import StationaryPointValidationPolicyV1
from chemsmart.analysis.result_quantities import (
    QuantitySelectorV1,
    ResultQuantityExtractionRequestV1,
    extract_pyscf_quantities,
    validate_pyscf_analysis_artifact,
)
from chemsmart.jobs.pyscf.writer import (
    APPLIED_SPEC_FIELDS,
    RESULT_CONTRACT_VERSION,
    PySCFScriptWriter,
    write_pyscf_h5,
)


def _artifact(path, *, artifact_id, kind):
    return TrustedArtifactRefV1(
        artifact_id=artifact_id,
        kind=kind,
        sha256=file_sha256(path),
        size_bytes=path.stat().st_size,
        path=str(path.resolve()),
        cli_value=str(path.resolve()),
    )


def _test_resources(*, node_timeout_seconds=600):
    return build_execution_resource_spec(
        execution_target="run",
        cores=4,
        memory_gb=4,
        gpu_count=0,
        scratch_policy="none",
        node_timeout_seconds=node_timeout_seconds,
    )


def _test_approval(tmp_path, *, node_timeout_seconds=600):
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
    resources = _test_resources(
        node_timeout_seconds=node_timeout_seconds
    )
    edge = build_producer_edge_rule(
        producer_node_id="opt-initial",
        consumer_node_id="hess-optimized",
        artifact_kind="geometry_xyz",
        selection_rule="validated_optimized_geometry",
    )
    initial_common = {
        "program": "pyscf",
        "engine": "cpu",
        "project_artifact_sha256": _artifact(
            project_path,
            artifact_id="project.water.pyscf",
            kind="project_yaml",
        ).sha256,
        "charge": 0,
        "multiplicity": 1,
        "input_mode": "initial",
        "initial_artifact_id": "geometry.water.initial",
        "initial_artifact_sha256": file_sha256(geometry_path),
        "scientific_identity_sha256": "d" * 64,
        "producer_edge_sha256": "",
    }
    nodes = (
        ApprovedNodeBindingV1(
            node_id="sp-initial",
            jobtype="sp",
            settings_sha256=canonical_sha256({"stage": "sp"}),
            **initial_common,
        ),
        ApprovedNodeBindingV1(
            node_id="opt-initial",
            jobtype="opt",
            settings_sha256=canonical_sha256({"stage": "opt"}),
            **initial_common,
        ),
        ApprovedNodeBindingV1(
            node_id="hess-optimized",
            program="pyscf",
            engine="cpu",
            jobtype="hess",
            project_artifact_sha256=initial_common["project_artifact_sha256"],
            settings_sha256=canonical_sha256({"stage": "hess"}),
            charge=0,
            multiplicity=1,
            input_mode="producer",
            initial_artifact_id="",
            initial_artifact_sha256="",
            scientific_identity_sha256="",
            producer_edge_sha256=edge.edge_sha256,
        ),
    )
    approval = build_workflow_execution_approval(
        approval_id="test-approval",
        workflow_id="test-pyscf-workflow",
        workflow_sha256=canonical_sha256(
            {"nodes": nodes, "producer_edges": (edge,)}
        ),
        task_spec_sha256="a" * 64,
        approved_workspace=tmp_path,
        resources=resources,
        node_bindings=nodes,
        producer_edges=(edge,),
    )
    return approval


def test_generic_approval_binds_order_and_optimized_handoff(tmp_path):
    approval = _test_approval(tmp_path)

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


def test_generic_approval_accepts_operator_selected_timeout(tmp_path):
    approval = _test_approval(tmp_path, node_timeout_seconds=5400)

    assert approval.resource_sha256 == _test_resources(
        node_timeout_seconds=5400
    ).resource_sha256


def test_execution_host_rejects_legacy_v1_approval_before_invocation_lookup(
    tmp_path,
):
    host = object.__new__(CommandCompiledToolHostV1)
    host.surface = SimpleNamespace(
        profile="command_compiled_approved_execution"
    )
    host.execution_receipts = {}
    host.workflow_execution_approval = _test_approval(tmp_path)

    with pytest.raises(
        ContractError,
        match="legacy V1 approval is preview-only; Runtime V2 frozen approval",
    ):
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


def test_auxiliary_input_is_bound_from_approval_through_execution(tmp_path):
    reactant_path = tmp_path / "reactant.xyz"
    reactant_path.write_text(
        "3\nreactant\nO 0.0 0.0 0.1\n"
        "H -0.8 0.0 -0.5\nH 0.8 0.0 -0.5\n",
        encoding="utf-8",
    )
    product_path = tmp_path / "product.xyz"
    product_path.write_text(
        "3\nproduct\nO 0.0 0.0 0.2\n"
        "H -0.8 0.0 -0.4\nH 0.8 0.0 -0.4\n",
        encoding="utf-8",
    )
    product = _artifact(
        product_path,
        artifact_id="geometry.water.product",
        kind="geometry_xyz",
    )
    reactant = _artifact(
        reactant_path,
        artifact_id="geometry.water.reactant",
        kind="geometry_xyz",
    )
    project_path = tmp_path / "water-neb.yaml"
    project_path.write_text(
        "gas:\n  functional: b3lyp\n  basis: def2-svp\n"
        "neb:\n  joboption: NEB-TS\n  nimages: 6\n",
        encoding="utf-8",
    )
    project = _artifact(
        project_path,
        artifact_id="project.water.neb",
        kind="project_yaml",
    )
    binding = AuxiliaryArtifactBindingV1(
        parameter_name="ending_xyzfile",
        artifact_id=product.artifact_id,
        artifact_sha256=product.sha256,
    )
    node = ApprovedNodeBindingV1(
        node_id="water-neb",
        program="orca",
        engine="cpu",
        jobtype="neb",
        project_artifact_sha256=project.sha256,
        settings_sha256="c" * 64,
        charge=0,
        multiplicity=1,
        input_mode="initial",
        initial_artifact_id=reactant.artifact_id,
        initial_artifact_sha256=reactant.sha256,
        scientific_identity_sha256="d" * 64,
        producer_edge_sha256="",
        auxiliary_input_bindings=(binding,),
    )
    approval = build_workflow_execution_approval(
        approval_id="water-with-product",
        workflow_id="water-neb-workflow",
        workflow_sha256="e" * 64,
        task_spec_sha256="a" * 64,
        approved_workspace=tmp_path,
        resources=_test_resources(),
        node_bindings=(node,),
    )
    node = approval.node("water-neb")
    invocation = build_program_execution_invocation(
        node_id=node.node_id,
        approval=approval,
        project_artifact=project,
        input_artifact=reactant,
        scientific_identity_sha256=node.scientific_identity_sha256,
        environment_receipt_sha256="b" * 64,
        resources=_test_resources(),
        argv=("chemsmart", "run", "orca", "neb"),
        auxiliary_input_artifacts={"ending_xyzfile": product},
    )

    assert invocation.auxiliary_input_bindings == (binding,)
    staged = tmp_path / "nodes" / "water-neb"
    staged.mkdir(parents=True)
    staged_product = staged / product_path.name
    staged_product.write_bytes(product_path.read_bytes())
    assert _staged_auxiliary_input_findings(
        node_workspace=staged,
        job_artifact_options=(("ending_xyzfile", product),),
    ) == ()

    product_path.write_text("changed", encoding="utf-8")
    with pytest.raises(ContractError, match="auxiliary input"):
        build_program_execution_invocation(
            node_id=node.node_id,
            approval=approval,
            project_artifact=project,
            input_artifact=reactant,
            scientific_identity_sha256=node.scientific_identity_sha256,
            environment_receipt_sha256="b" * 64,
            resources=_test_resources(),
            argv=("chemsmart", "run", "orca", "neb"),
            auxiliary_input_artifacts={"ending_xyzfile": product},
        )


def test_validated_opt_handoff_materializes_exact_final_geometry(tmp_path):
    approval = _test_approval(tmp_path)
    resources = _test_resources()
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
            execution_state="validated",
            exit_status=0,
            engine_complete=True,
            validated=True,
            output_artifacts=(result,),
            validator_receipt_sha256s=("e" * 64,),
            result_validation_receipt_sha256="e" * 64,
            findings=("pyscf.result.seeded_failure",),
            started_at="2026-08-03T00:00:00+00:00",
            finished_at="2026-08-03T00:00:01+00:00",
        )
    receipt = build_program_execution_receipt(
        invocation,
        execution_state="validated",
        exit_status=0,
        engine_complete=True,
        validated=True,
        output_artifacts=(result,),
        validator_receipt_sha256s=("e" * 64,),
        result_validation_receipt_sha256="e" * 64,
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


def test_execution_receipt_separates_wrapper_and_child_status(tmp_path):
    approval = _test_approval(tmp_path)
    node = approval.node("sp-initial")
    invocation = build_program_execution_invocation(
        node_id=node.node_id,
        approval=approval,
        project_artifact=_artifact(
            tmp_path / "water-pyscf.yaml",
            artifact_id="project.water.pyscf",
            kind="project_yaml",
        ),
        input_artifact=_artifact(
            tmp_path / "water.xyz",
            artifact_id="geometry.water.initial",
            kind="geometry_xyz",
        ),
        scientific_identity_sha256=node.scientific_identity_sha256,
        environment_receipt_sha256="b" * 64,
        resources=_test_resources(),
        argv=("chemsmart", "run", "pyscf", "sp"),
    )

    receipt = build_program_execution_receipt(
        invocation,
        execution_state="engine_complete",
        exit_status=17,
        child_exit_status=0,
        engine_complete=True,
        validated=False,
        engine_receipt_sha256="c" * 64,
        findings=("pyscf.result.seeded_failure",),
        started_at="2026-08-04T00:00:00+00:00",
        finished_at="2026-08-04T00:00:01+00:00",
    )

    assert receipt.wrapper_exit_status == 17
    assert receipt.exit_status == 17
    assert receipt.child_exit_status == 0
    assert receipt.engine_complete is True
    assert receipt.validated is False

    for contradictory_state in ("engine_complete", "failed"):
        with pytest.raises(
            ContractError,
            match="validated execution requires validated execution state",
        ):
            build_program_execution_receipt(
                invocation,
                execution_state=contradictory_state,
                exit_status=0,
                child_exit_status=0,
                engine_complete=True,
                validated=True,
                started_at="2026-08-04T00:00:00+00:00",
                finished_at="2026-08-04T00:00:01+00:00",
            )

    with pytest.raises(ContractError, match="child exit status zero"):
        build_program_execution_receipt(
            invocation,
            execution_state="engine_complete",
            exit_status=0,
            child_exit_status=2,
            engine_complete=True,
            validated=False,
            started_at="2026-08-04T00:00:00+00:00",
            finished_at="2026-08-04T00:00:01+00:00",
        )


def test_execution_receipt_binds_resolvable_result_validation(tmp_path):
    approval = _test_approval(tmp_path)
    node = approval.node("sp-initial")
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
    invocation = build_program_execution_invocation(
        node_id=node.node_id,
        approval=approval,
        project_artifact=project,
        input_artifact=initial,
        scientific_identity_sha256=node.scientific_identity_sha256,
        environment_receipt_sha256="b" * 64,
        resources=_test_resources(),
        argv=("chemsmart", "run", "pyscf", "sp"),
    )
    output_path = tmp_path / "observable.json"
    output_path.write_text("{}", encoding="utf-8")
    output = _artifact(
        output_path,
        artifact_id="result.sp.observation",
        kind="json",
    )
    validation = build_program_result_validation_receipt(
        invocation=invocation,
        validator_id="pyscf-result-validator",
        validator_schema_version="chemsmart.pyscf-result-validation.v1",
        validator_version="1",
        input_artifact=initial,
        project_artifact=project,
        capability_environment_receipt_sha256="b" * 64,
        output_artifacts=(output,),
        observations={"state": "validated", "energy_finite": True},
        findings=(),
        run_environment_receipt_sha256="c" * 64,
        environment_validation_sha256="d" * 64,
    )
    failed_validation = build_program_result_validation_receipt(
        invocation=invocation,
        validator_id="pyscf-result-validator",
        validator_schema_version="chemsmart.pyscf-result-validation.v1",
        validator_version="1",
        input_artifact=initial,
        project_artifact=project,
        capability_environment_receipt_sha256="b" * 64,
        output_artifacts=(output,),
        observations={
            "result_validation": {
                "state": "failed",
                "findings": (),
            }
        },
        findings=(),
    )
    assert failed_validation.state == "invalid"
    assert failed_validation.findings == (
        "result_validation.state_not_validated",
    )
    execution = build_program_execution_receipt(
        invocation,
        execution_state="validated",
        exit_status=0,
        child_exit_status=0,
        engine_complete=True,
        validated=True,
        result_validation_receipt_sha256=validation.receipt_sha256,
        output_artifacts=(output,),
        validator_receipt_sha256s=(validation.receipt_sha256,),
        started_at="2026-08-04T00:00:00+00:00",
        finished_at="2026-08-04T00:00:01+00:00",
    )

    assert execution.result_validation_receipt_sha256 == (
        validation.receipt_sha256
    )

    event = RuntimeEvent.create(
        sequence=1,
        session_id="result-validation-session",
        turn_id="turn-1",
        kind=EventKind.RESULT_VERIFIED,
        payload={
            "receipt_sha256": validation.receipt_sha256,
            "status": "valid",
            "critical_finding_count": 0,
            "record": canonical_data(validation),
        },
    )
    assert event.payload["record"]["receipt_sha256"] == (
        validation.receipt_sha256
    )

    substituted_record = dict(canonical_data(validation))
    substituted_record["validator_version"] = "substituted"
    with pytest.raises(ContractError, match="record digest mismatch"):
        RuntimeEvent.create(
            sequence=2,
            session_id="result-validation-session",
            turn_id="turn-1",
            kind=EventKind.RESULT_VERIFIED,
            payload={
                "receipt_sha256": validation.receipt_sha256,
                "status": "valid",
                "critical_finding_count": 0,
                "record": substituted_record,
            },
        )


def test_stationary_point_policy_rejects_counterfeit_contracts():
    frozen = SimpleNamespace(stationary_point_policy_sha256="a" * 64)
    approved = SimpleNamespace(policy_sha256="a" * 64)

    with pytest.raises(ContractError, match="real approval"):
        _validate_stationary_point_policy_binding(frozen, approved)
    with pytest.raises(ContractError, match="real approval"):
        _validate_stationary_point_policy_binding(
            frozen, SimpleNamespace(policy_sha256="b" * 64)
        )
    with pytest.raises(ContractError, match="real policy"):
        _validate_stationary_point_policy_binding(None, approved)


def _minimum_stationary_point_policy():
    body = {
        "schema_version": "chemsmart.stationary-point-policy.v1",
        "policy_id": "water-minimum",
        "task_spec_sha256": "a" * 64,
        "hessian_node_id": "hess-optimized",
        "stationary_point_kind": "minimum",
        "expected_imaginary_mode_count": 0,
        "imaginary_mode_cutoff_cm1": 20.0,
        "require_finite_modes": True,
        "require_symmetric_hessian": True,
    }
    return StationaryPointValidationPolicyV1(
        **body, policy_sha256=canonical_sha256(body)
    )


def _unclassified_hessian_run_receipt():
    return {
        "state": "engine_complete",
        "child_returncode": 0,
        "engine_complete": True,
        "scientifically_validated": False,
        "scientific_validation_state": "unclassified",
        "findings": [],
        "result_validation": {
            "state": "unclassified",
            "findings": [],
            "frequency_validation": {
                "stationary_point_classification": "unclassified",
                "findings": [],
            },
        },
    }


def test_agent_admits_complete_unclassified_hessian_for_downstream_analysis():
    policy = _minimum_stationary_point_policy()
    run_receipt = _unclassified_hessian_run_receipt()

    assert _runner_defers_hessian_classification(
        run_receipt=run_receipt,
        jobtype="hess",
        hessian_node_id="hess-optimized",
        engine_complete=True,
        stationary_point_policy=None,
        approved_stationary_point_policy_sha256="",
    ) is True
    assert _runner_defers_hessian_classification(
        run_receipt=run_receipt,
        jobtype="hess",
        hessian_node_id="hess-optimized",
        engine_complete=True,
        stationary_point_policy=policy,
        approved_stationary_point_policy_sha256=policy.policy_sha256,
    ) is True
    assert _runner_defers_hessian_classification(
        run_receipt=run_receipt,
        jobtype="hess",
        hessian_node_id="hess-optimized",
        engine_complete=True,
        stationary_point_policy=policy,
        approved_stationary_point_policy_sha256="b" * 64,
    ) is False
    assert _runner_defers_hessian_classification(
        run_receipt=run_receipt,
        jobtype="hess",
        hessian_node_id="another-hessian",
        engine_complete=True,
        stationary_point_policy=policy,
        approved_stationary_point_policy_sha256=policy.policy_sha256,
    ) is False
    assert _runner_defers_hessian_classification(
        run_receipt=run_receipt,
        jobtype="sp",
        hessian_node_id="hess-optimized",
        engine_complete=True,
        stationary_point_policy=policy,
        approved_stationary_point_policy_sha256=policy.policy_sha256,
    ) is False


def test_agent_rejects_incomplete_unclassified_runner_receipt_even_with_policy():
    policy = _minimum_stationary_point_policy()
    run_receipt = _unclassified_hessian_run_receipt()
    run_receipt["result_validation"] = {
        **run_receipt["result_validation"],
        "state": "validated",
    }

    assert _runner_defers_hessian_classification(
        run_receipt=run_receipt,
        jobtype="hess",
        hessian_node_id="hess-optimized",
        engine_complete=True,
        stationary_point_policy=policy,
        approved_stationary_point_policy_sha256=policy.policy_sha256,
    ) is False


def test_parent_policy_validator_promotes_exact_unclassified_hessian_handoff(
    tmp_path,
):
    policy = _minimum_stationary_point_policy()
    result_path = tmp_path / "water-hess.h5"
    result_path.write_bytes(b"host-bound-hessian-result")
    result = _artifact(
        result_path,
        artifact_id="result.hess.hdf5",
        kind="pyscf_hdf5",
    )
    engine = _PySCFEngineObservation(
        child_exit_status=0,
        engine_complete=True,
        run_receipt_sha256="c" * 64,
        run_receipt=_unclassified_hessian_run_receipt(),
        result_artifact=result,
        findings=(),
    )
    host_validation = {
        "schema_version": "chemsmart.pyscf-result-validation.v1",
        "state": "validated",
        "findings": [],
    }

    with (
        patch(
            "chemsmart.jobs.pyscf.validation.validate_pyscf_result",
            return_value=host_validation,
        ) as validate,
        patch(
            "chemsmart.agent.tool_runtime._pyscf_input_geometry",
            return_value=(("O", "H", "H"), np.zeros((3, 3))),
        ),
        patch(
            "chemsmart.io.pyscf.output.read_pyscf_h5",
            return_value=(
                {"program": "pyscf", "engine": "cpu"},
                {},
                {},
                {},
            ),
        ),
    ):
        evaluation = CommandCompiledToolHostV1._evaluate_execution_outputs(
            program="pyscf",
            jobtype="hess",
            charge=0,
            multiplicity=1,
            expected_settings={"engine": "cpu"},
            output_artifacts=(result,),
            exit_status=0,
            pyscf_engine_observation=engine,
            stationary_point_policy=policy,
            approved_stationary_point_policy_sha256=policy.policy_sha256,
            approved_hessian_node_id="hess-optimized",
        )

    assert evaluation.validated is True
    assert evaluation.findings == ()
    assert evaluation.observations["runner_validation_delegation"] == (
        "approved_stationary_point_policy"
    )
    assert evaluation.observations["stationary_point_policy_sha256"] == (
        policy.policy_sha256
    )
    assert validate.call_args.kwargs["stationary_point_policy"] is policy


def test_unclassified_hessian_is_valid_data_for_downstream_scientific_analysis(
    tmp_path,
):
    result_path = tmp_path / "water-hess.h5"
    result_path.write_bytes(b"host-bound-hessian-result")
    result = _artifact(
        result_path,
        artifact_id="result.hess.hdf5",
        kind="pyscf_hdf5",
    )
    engine = _PySCFEngineObservation(
        child_exit_status=0,
        engine_complete=True,
        run_receipt_sha256="c" * 64,
        run_receipt=_unclassified_hessian_run_receipt(),
        result_artifact=result,
        findings=(),
    )
    host_validation = {
        "schema_version": "chemsmart.pyscf-result-validation.v1",
        "state": "unclassified",
        "findings": [],
    }

    with (
        patch(
            "chemsmart.jobs.pyscf.validation.validate_pyscf_result",
            return_value=host_validation,
        ),
        patch(
            "chemsmart.agent.tool_runtime._pyscf_input_geometry",
            return_value=(("O", "H", "H"), np.zeros((3, 3))),
        ),
        patch(
            "chemsmart.io.pyscf.output.read_pyscf_h5",
            return_value=(
                {"program": "pyscf", "engine": "cpu"},
                {},
                {},
                {},
            ),
        ),
    ):
        evaluation = CommandCompiledToolHostV1._evaluate_execution_outputs(
            program="pyscf",
            jobtype="hess",
            charge=0,
            multiplicity=1,
            expected_settings={"engine": "cpu"},
            output_artifacts=(result,),
            exit_status=0,
            pyscf_engine_observation=engine,
            stationary_point_policy=None,
            approved_stationary_point_policy_sha256="",
            approved_hessian_node_id="hess-optimized",
        )

    assert evaluation.validated is True
    assert evaluation.findings == ()
    assert evaluation.observations["runner_validation_delegation"] == (
        "downstream_scientific_analysis"
    )
    assert _typed_result_validation_findings(evaluation.observations) == ()


def _write_digest_bound_receipt(path, payload):
    body = dict(payload)
    body["receipt_sha256"] = canonical_sha256(body)
    path.write_text(json.dumps(body, sort_keys=True), encoding="utf-8")
    return body


def _pyscf_environment_pair(
    tmp_path, *, pyscf_version="2.14.0", libxc_version="7.0.0"
):
    environment_path = tmp_path / "pyscf-environment.receipt.json"
    raw = _write_digest_bound_receipt(
        environment_path,
        {
            "schema_version": "chemsmart.pyscf-environment.v1",
            "status": "available",
            "probe_returncode": 0,
            "interpreter": "/opt/example/bin/python",
            "interpreter_observed": "/opt/example/bin/python",
            "interpreter_sha256": "1" * 64,
            "pyscf_version": pyscf_version,
            "libxc_version": libxc_version,
            "packages": {
                "h5py": "3.15.1",
                "numpy": "2.2.6",
                "pyscf": pyscf_version,
            },
            "dependencies": {
                "h5py": {"available": True, "version": "3.15.1"},
                "numpy": {"available": True, "version": "2.2.6"},
                "pyscf": {"available": True, "version": pyscf_version},
            },
            "solver_callables": {
                "geometric": {"callable": True},
            },
            "cuda_available": 0,
            "cuda_driver_version": 0,
            "cuda_runtime_version": 0,
            "device_count": 0,
            "device_name": "",
            "device_uuid": "",
            "cutensor_runtime": 0,
            "cutensor_compatible": False,
            "cutensor_runtime_available": False,
        },
    )
    compute = consume_pyscf_compute_environment_receipt(raw, engine="cpu")
    query = ProgramEnvironmentQueryV1(
        capability_receipt_sha256="2" * 64,
        program="pyscf",
        engine="cpu",
    )
    body = {
        "schema_version": "chemsmart.environment-capability-receipt.v1",
        "query": query,
        "capability_receipt_sha256": query.capability_receipt_sha256,
        "program": "pyscf",
        "engine": "cpu",
        "status": EnvironmentStatus.AVAILABLE,
        "target_kind": "compute_receipt",
        "locator": "pyscf-cpu",
        "observed_version": pyscf_version,
        "observed_location_sha256": "",
        "compute_interpreter_sha256": compute.compute_interpreter_sha256,
        "compute_evidence_sha256": compute.evidence_sha256,
        "dependency_versions": compute.dependency_versions,
        "solver_evidence": compute.solver_evidence,
        "gpu_evidence": compute.gpu_evidence,
        "observation_method": "trusted_target_compute_receipt",
        "rule_ids": ("environment.compute_receipt.complete",),
    }
    capability = EnvironmentCapabilityReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )
    return environment_path, raw, capability


def test_run_environment_is_semantically_compared_not_digest_equal(tmp_path):
    path, raw, capability = _pyscf_environment_pair(tmp_path)
    artifact = _artifact(
        path,
        artifact_id="result.sp.environment",
        kind="json",
    )

    observation, findings = _pyscf_environment_evidence(
        output_artifacts=(artifact,),
        run_receipt={"environment_receipt_sha256": raw["receipt_sha256"]},
        capability_environment=capability,
    )

    assert capability.receipt_sha256 != raw["receipt_sha256"]
    assert findings == ()
    assert observation["state"] == "valid"
    assert observation["approved_semantic_fingerprint_sha256"] == (
        observation["observed_semantic_fingerprint_sha256"]
    )


def test_run_environment_allows_additional_observed_facts(tmp_path):
    path, raw, capability = _pyscf_environment_pair(tmp_path)
    body = {
        field: getattr(capability, field)
        for field in capability.__dataclass_fields__
        if field != "receipt_sha256"
    }
    body["dependency_versions"] = tuple(
        item
        for item in capability.dependency_versions
        if item[0] in {"h5py", "numpy", "pyscf"}
    )
    body["solver_evidence"] = ()
    body["gpu_evidence"] = ()
    narrowed = EnvironmentCapabilityReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )
    artifact = _artifact(
        path,
        artifact_id="result.sp.environment.superset",
        kind="json",
    )

    observation, findings = _pyscf_environment_evidence(
        output_artifacts=(artifact,),
        run_receipt={"environment_receipt_sha256": raw["receipt_sha256"]},
        capability_environment=narrowed,
    )

    assert findings == ()
    assert observation["state"] == "valid"


def test_run_environment_substitution_is_rejected_semantically(tmp_path):
    original_path, _raw, capability = _pyscf_environment_pair(tmp_path)
    original_path.unlink()
    substituted_path, substituted, _unused = _pyscf_environment_pair(
        tmp_path, pyscf_version="2.13.0"
    )
    artifact = _artifact(
        substituted_path,
        artifact_id="result.sp.environment",
        kind="json",
    )

    observation, findings = _pyscf_environment_evidence(
        output_artifacts=(artifact,),
        run_receipt={
            "environment_receipt_sha256": substituted["receipt_sha256"]
        },
        capability_environment=capability,
    )

    assert observation["state"] == "invalid"
    assert "pyscf.environment.semantic_mismatch" in findings


def test_run_environment_libxc_substitution_is_rejected_semantically(tmp_path):
    original_path, _raw, capability = _pyscf_environment_pair(tmp_path)
    original_path.unlink()
    substituted_path, substituted, _unused = _pyscf_environment_pair(
        tmp_path, libxc_version="6.2.2"
    )
    artifact = _artifact(
        substituted_path,
        artifact_id="result.sp.environment.libxc",
        kind="json",
    )

    observation, findings = _pyscf_environment_evidence(
        output_artifacts=(artifact,),
        run_receipt={
            "environment_receipt_sha256": substituted["receipt_sha256"]
        },
        capability_environment=capability,
    )

    assert observation["state"] == "invalid"
    assert "pyscf.environment.semantic_mismatch" in findings


@pytest.mark.parametrize("bind_source_artifact", [True, False])
def test_nonzero_wrapper_still_inspects_digest_bound_pyscf_result(
    tmp_path, bind_source_artifact
):
    _test_approval(tmp_path)
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
    h5_path = tmp_path / "water-sp.h5"
    spec = {field: None for field in APPLIED_SPEC_FIELDS}
    spec.update(
        {
            "run_id": "agent-run",
            "run_nonce": "agent-nonce",
            "program": "pyscf",
            "result_contract_version": RESULT_CONTRACT_VERSION,
            "preview_only": False,
            "jobtype": "sp",
            "engine": "cpu",
            "ab_initio": "hf",
            "method": "hf",
            "reference_family": "rhf",
            "stages": ["scf"],
            "symbols": ["O", "H", "H"],
            "positions": [
                [0.0, 0.0, 0.1174],
                [-0.757, 0.0, -0.4696],
                [0.757, 0.0, -0.4696],
            ],
            "unit": "Angstrom",
            "charge": 0,
            "spin": 0,
            "multiplicity": 1,
            "num_electrons": 10,
            "nelec": [5, 5],
            "input_geometry_sha256": canonical_sha256(
                {
                    "symbols": ("O", "H", "H"),
                    "positions": (
                        (0.0, 0.0, 0.1174),
                        (-0.757, 0.0, -0.4696),
                        (0.757, 0.0, -0.4696),
                    ),
                    "unit": "Angstrom",
                    "charge": 0,
                    "multiplicity": 1,
                }
            ),
            "input_artifact_kind": (
                "geometry_xyz" if bind_source_artifact else None
            ),
            "input_artifact_sha256": (
                initial.sha256 if bind_source_artifact else None
            ),
            "requested_settings_sha256": "d" * 64,
            "settings_digest": "4" * 64,
        }
    )
    spec["applied_settings_sha256"] = PySCFScriptWriter.settings_digest(
        {field: spec.get(field) for field in APPLIED_SPEC_FIELDS}
    )
    bound_provenance = {
        "run_id": spec["run_id"],
        "run_nonce": spec["run_nonce"],
        "script_sha256": "5" * 64,
        "input_receipt_sha256": "6" * 64,
        "environment_receipt_sha256": "b" * 64,
        "input_geometry_sha256": spec["input_geometry_sha256"],
        "input_artifact_kind": spec["input_artifact_kind"],
        "input_artifact_sha256": spec["input_artifact_sha256"],
        "requested_settings_sha256": spec[
            "requested_settings_sha256"
        ],
        "applied_settings_sha256": spec["applied_settings_sha256"],
        "project_yaml_digest": project.sha256,
        "settings_digest": spec["settings_digest"],
        "engine": "cpu",
        "runtime": {
            "mean_field_class": "pyscf.scf.hf.RHF",
            "mean_field_mro": ["pyscf.scf.hf.RHF", "object"],
        },
    }
    status_payload = {
        "engine_complete": True,
        "normal_termination": True,
        "stages": {"scf": {"converged": True}},
        "failure": None,
        "properties": {"spin_square": {"status": "ok"}},
    }
    results_payload = {
        "energies": np.asarray([-76.3]),
        "positions": np.asarray(
            [
                [0.0, 0.0, 0.1174],
                [-0.757, 0.0, -0.4696],
                [0.757, 0.0, -0.4696],
            ]
        ),
        "atomic_numbers": np.asarray([8, 1, 1]),
        "mo_energy": np.asarray(
            [-20.0, -1.0, -0.8, -0.6, -0.4, 0.5]
        ),
        "mo_occ": np.asarray([2.0, 2.0, 2.0, 2.0, 2.0, 0.0]),
        "spin_square": np.asarray(0.0),
        "spin_square_effective_multiplicity": np.asarray(1.0),
    }
    write_pyscf_h5(
        h5_path,
        spec=spec,
        provenance=bound_provenance,
        status=status_payload,
        results=results_payload,
    )
    from chemsmart.io.pyscf.output import read_pyscf_h5

    roundtrip_spec, _provenance, _status, _results = read_pyscf_h5(h5_path)
    spec = dict(roundtrip_spec)
    spec["applied_settings_sha256"] = PySCFScriptWriter.settings_digest(
        {field: spec.get(field) for field in APPLIED_SPEC_FIELDS}
    )
    bound_provenance["applied_settings_sha256"] = spec[
        "applied_settings_sha256"
    ]
    write_pyscf_h5(
        h5_path,
        spec=spec,
        provenance=bound_provenance,
        status=status_payload,
        results=results_payload,
    )
    result = _artifact(
        h5_path,
        artifact_id="result.sp.hdf5",
        kind="pyscf_hdf5",
    )
    run_path = tmp_path / "water-sp.receipt.json"
    run_receipt = _write_digest_bound_receipt(
        run_path,
        {
            "schema_version": "chemsmart.pyscf-run.v1",
            "state": "validated",
            "fake": False,
            "child_returncode": 0,
                "engine_complete": True,
                "scientifically_validated": True,
                "scientific_validation_state": "unclassified",
                "run_id": spec["run_id"],
            "run_nonce": spec["run_nonce"],
            "script_sha256": bound_provenance["script_sha256"],
            "input_receipt_sha256": bound_provenance[
                "input_receipt_sha256"
            ],
            "environment_receipt_sha256": "b" * 64,
            "input_geometry_sha256": spec["input_geometry_sha256"],
                "input_artifact_kind": spec["input_artifact_kind"],
                "project_yaml_sha256": project.sha256,
                "input_artifact_sha256": spec["input_artifact_sha256"],
            "requested_settings_sha256": spec[
                "requested_settings_sha256"
            ],
            "applied_settings_sha256": spec[
                "applied_settings_sha256"
            ],
            "result_sha256": result.sha256,
            "findings": [],
        },
    )
    run_artifact = _artifact(
        run_path,
        artifact_id="result.sp.run-receipt",
        kind="json",
    )
    outputs = (result, run_artifact)

    engine = _inspect_pyscf_engine_observation(
        outputs, launch_ambiguous=False
    )
    evaluation = CommandCompiledToolHostV1._evaluate_execution_outputs(
        program="pyscf",
        jobtype="sp",
            charge=0,
            multiplicity=1,
            expected_input_artifact=initial,
            expected_project_artifact=project,
        output_artifacts=outputs,
        exit_status=17,
        expected_environment_receipt_sha256="b" * 64,
        pyscf_engine_observation=engine,
    )
    validated, findings = evaluation.validated, evaluation.findings

    assert engine.engine_complete is True
    assert engine.child_exit_status == 0
    assert engine.run_receipt_sha256 == run_receipt["receipt_sha256"]
    assert validated is True, json.dumps(
        evaluation.observations, sort_keys=True
    )
    assert findings == ()
    analysis_output, analysis_receipt = validate_pyscf_analysis_artifact(
        h5_path,
        expected_sha256=result.sha256,
    )
    assert analysis_output.energies[-1] == pytest.approx(-76.3)
    assert analysis_receipt["result_sha256"] == result.sha256
    extraction = extract_pyscf_quantities(
        request=ResultQuantityExtractionRequestV1(
            schema_version="chemsmart.quantity-extraction-request.v1",
            artifact_id=result.artifact_id,
            artifact_sha256=result.sha256,
            program="pyscf",
            selectors=(
                QuantitySelectorV1(quantity_id="homo", selector="homo"),
            ),
        ),
        artifact_path=h5_path,
    )
    assert extraction.quantities[0].source_value == pytest.approx(
        -0.4 * units.Hartree
    )
    assert extraction.quantities[0].value == pytest.approx(-0.4, abs=1e-15)


def test_substituted_pyscf_run_receipt_cannot_prove_engine_completion(
    tmp_path,
):
    h5_path = tmp_path / "result.h5"
    h5_path.write_bytes(b"not inspected without a valid run receipt")
    result = _artifact(
        h5_path,
        artifact_id="result.sp.hdf5",
        kind="pyscf_hdf5",
    )
    run_path = tmp_path / "result.receipt.json"
    run_path.write_text(
        json.dumps(
            {
                "schema_version": "chemsmart.pyscf-run.v1",
                "fake": False,
                "child_returncode": 0,
                "engine_complete": True,
                "result_sha256": result.sha256,
                "receipt_sha256": "0" * 64,
            }
        ),
        encoding="utf-8",
    )
    run_artifact = _artifact(
        run_path,
        artifact_id="result.sp.run-receipt",
        kind="json",
    )

    engine = _inspect_pyscf_engine_observation(
        (result, run_artifact), launch_ambiguous=False
    )

    assert engine.engine_complete is False
    assert engine.run_receipt is None
    assert engine.findings == ("pyscf.run_receipt.digest_invalid",)


def test_output_is_rehashed_immediately_before_validation(tmp_path):
    output_path = tmp_path / "program.out"
    output_path.write_text("original-bytes", encoding="utf-8")
    artifact = _artifact(
        output_path,
        artifact_id="result.seeded.output",
        kind="program_output",
    )
    output_path.write_text("substitute!!!", encoding="utf-8")

    evaluation = CommandCompiledToolHostV1._evaluate_execution_outputs(
        program="gaussian",
        jobtype="sp",
        charge=0,
        multiplicity=1,
        output_artifacts=(artifact,),
        exit_status=0,
    )

    assert "execution.output.artifact_binding_mismatch" in (
        evaluation.findings
    )


def test_orca_error_termination_is_not_engine_completion():
    output_path = Path(
        "tests/data/ORCATests/error_files/GTOInt_error.out"
    ).resolve()
    output = _artifact(
        output_path,
        artifact_id="result.orca.error",
        kind="orca_output",
    )

    evaluation = CommandCompiledToolHostV1._evaluate_execution_outputs(
        program="orca",
        jobtype="opt",
        charge=0,
        multiplicity=1,
        expected_settings={"freq": True},
        output_artifacts=(output,),
        # The command wrapper can return zero even when ORCA reports an
        # internal error; program termination is the scientific authority.
        exit_status=0,
    )

    assert evaluation.validated is False
    assert "orca.result.normal_termination" in evaluation.findings
    assert evaluation.observations["orca"]["normal_termination"] is False


def test_orca_validator_reports_explicit_dft_d_total_not_scf_component():
    output_path = Path("tests/data/ORCATests/outputs/KOH.out").resolve()
    output = _artifact(
        output_path,
        artifact_id="result.orca.dftd",
        kind="orca_output",
    )

    evaluation = CommandCompiledToolHostV1._evaluate_execution_outputs(
        program="orca",
        jobtype="opt",
        charge=0,
        multiplicity=1,
        expected_settings={"freq": True},
        output_artifacts=(output,),
        exit_status=0,
    )

    assert evaluation.validated is True
    assert evaluation.observations["orca"]["energy_hartree"] == pytest.approx(
        -675.522805891018
    )


@pytest.mark.parametrize(
    (
        "neb_converged",
        "ts_converged",
        "expected_validated",
        "expected_finding",
    ),
    (
        (False, False, False, "orca.result.neb_not_converged"),
        (True, False, False, "orca.result.neb_ts_not_converged"),
        (True, True, True, None),
    ),
)
def test_orca_neb_execution_requires_path_and_ts_convergence(
    tmp_path,
    neb_converged,
    ts_converged,
    expected_validated,
    expected_finding,
):
    output_path = tmp_path / "neb.out"
    output_path.write_text("host-bound ORCA output", encoding="utf-8")
    output_artifact = _artifact(
        output_path,
        artifact_id="result.orca.neb",
        kind="orca_output",
    )
    parsed_output = SimpleNamespace(
        vibrational_frequencies=(),
        normal_termination=True,
        converged=False,
        charge=0,
        multiplicity=1,
        final_energy=-100.0,
        route_object=SimpleNamespace(neb_joboption="NEB-TS"),
        neb_converged=neb_converged,
        ts_converged=ts_converged,
    )

    with patch(
        "chemsmart.io.orca.output.ORCANEBOutput",
        return_value=parsed_output,
    ):
        evaluation = CommandCompiledToolHostV1._evaluate_execution_outputs(
            program="orca",
            jobtype="neb",
            charge=0,
            multiplicity=1,
            expected_settings={"joboption": "NEB-TS", "freq": False},
            output_artifacts=(output_artifact,),
            exit_status=0,
        )

    assert evaluation.validated is expected_validated
    if expected_finding is None:
        assert evaluation.findings == ()
    else:
        assert expected_finding in evaluation.findings
    assert evaluation.observations["orca"]["neb_converged"] is (
        neb_converged
    )
    assert evaluation.observations["orca"]["ts_converged"] is ts_converged


def test_live_prompt_preserves_typed_edges_and_evidence_bounded_alternatives():
    prompt = _system_prompt(None)

    assert "do not convert a presentation sequence into a control edge" in prompt
    assert "SP(initial geometry) and OPT(initial geometry) are siblings" in prompt
    assert "Bind scientific identity only to a geometry_xyz" in prompt
    assert "a workflow_draft alone is not the typed scientific DAG" in prompt
    assert "loader-supported, preview-conformant" in prompt
    assert "do not assert quantitative accuracy, cost" in prompt.lower()
    assert "unmaterialized alternative" in prompt
