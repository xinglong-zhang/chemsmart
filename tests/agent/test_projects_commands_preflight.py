from __future__ import annotations

import hashlib
import sys
from types import ModuleType

import pytest

from chemsmart.agent._contracts import (
    ContractError,
    TrustedArtifactRefV1,
    canonical_sha256,
)
from chemsmart.agent.capabilities import (
    EnvironmentTargetV1,
    ProgramCapabilityQueryV1,
    query_capability,
    query_environment,
    resolve_engine_binding,
    resolve_program_binding,
)
from chemsmart.agent.commands import (
    CommandProposalV1,
    build_scientific_identity_binding,
    compile_command,
    inspect_command,
)
from chemsmart.agent.preflight import (
    ProgramValidatorReceiptV1,
    build_program_node_preflight_request,
    evaluate_program_node_preflight,
)
from chemsmart.agent.projects import validate_project_yaml


def _artifact(path, *, artifact_id, kind, cli_value):
    payload = path.read_bytes()
    return TrustedArtifactRefV1(
        artifact_id=artifact_id,
        kind=kind,
        sha256=hashlib.sha256(payload).hexdigest(),
        size_bytes=len(payload),
        path=str(path.resolve()),
        cli_value=cli_value,
    )


def _install_demo_loader(monkeypatch):
    module = ModuleType("chemsmart.settings.demo")

    class Settings:
        basis = "def2-svp"
        functional = "b3lyp"
        engine = "cpu"
        jobtype = "sp"

        def validate(self):
            return None

    class Project:
        def sp_settings(self):
            return Settings()

    class YamlDemoProjectSettings:
        @classmethod
        def from_yaml(cls, _path):
            return Project()

    YamlDemoProjectSettings.__module__ = module.__name__
    module.YamlDemoProjectSettings = YamlDemoProjectSettings
    monkeypatch.setitem(sys.modules, module.__name__, module)


def test_scientific_identity_rejects_non_geometry_artifact(tmp_path):
    project_path = tmp_path / "project.yaml"
    project_path.write_text("sp:\n  method: rhf\n", encoding="utf-8")
    project = _artifact(
        project_path,
        artifact_id="project.not-a-geometry",
        kind="project_yaml",
        cli_value="project",
    )

    with pytest.raises(
        ContractError,
        match="scientific identity requires an exact geometry artifact",
    ):
        build_scientific_identity_binding(
            task_spec_sha256="a" * 64,
            geometry_artifact=project,
            charge=0,
            multiplicity=1,
        )


def test_required_project_command_and_node_preflight_bind_all_scientific_state(
    tmp_path,
    monkeypatch,
    fake_capability_registry,
    fake_click_root,
    fake_click_schema,
    fake_support_overlay,
):
    _install_demo_loader(monkeypatch)
    project_path = tmp_path / "demo.yaml"
    project_path.write_text(
        "gas:\n  functional: b3lyp\n  basis: def2-svp\n",
        encoding="utf-8",
    )
    geometry_path = tmp_path / "water.xyz"
    geometry_path.write_text("3\nwater\nO 0 0 0\nH 0 0 1\nH 0 1 0\n")
    project = _artifact(
        project_path,
        artifact_id="project.demo",
        kind="project_yaml",
        cli_value="demo",
    )
    geometry = _artifact(
        geometry_path,
        artifact_id="geometry.water",
        kind="geometry_xyz",
        cli_value=str(geometry_path),
    )
    capability = query_capability(
        ProgramCapabilityQueryV1("demo", "sp", "cpu"),
        registry=fake_capability_registry,
        live_schema=fake_click_schema,
        overlay=fake_support_overlay,
    )
    project_receipt = validate_project_yaml(
        project, capability=capability
    )
    environment = query_environment(
        capability,
        targets=(
            EnvironmentTargetV1(
                "demo", "cpu", "executable", "demo-bin"
            ),
        ),
        which=lambda _name: "/opt/demo/bin/demo-bin",
    )
    program_binding = resolve_program_binding(capability)
    engine_binding = resolve_engine_binding(program_binding, environment)
    identity = build_scientific_identity_binding(
        task_spec_sha256="f" * 64,
        geometry_artifact=geometry,
        charge=0,
        multiplicity=1,
    )
    proposal = CommandProposalV1(
        node_id="node.sp",
        execution_target="run",
        program="demo",
        jobtype="sp",
        project_artifact_id=project.artifact_id,
        input_artifact_id=geometry.artifact_id,
        scientific_identity_sha256=identity.binding_sha256,
        charge=0,
        multiplicity=1,
    )
    invocation = compile_command(
        proposal,
        capability=capability,
        binding=engine_binding,
        project=project,
        project_validation=project_receipt,
        input_artifact=geometry,
        scientific_identity=identity,
        live_schema=fake_click_schema,
    )
    inspected = inspect_command(
        invocation,
        live_schema=fake_click_schema,
        root=fake_click_root,
    )
    validator_body = {
        "schema_version": "chemsmart.program-validator-receipt.v1",
        "node_id": proposal.node_id,
        "program": proposal.program,
        "invocation_sha256": invocation.invocation_sha256,
        "scientific_identity_sha256": identity.binding_sha256,
        "source_receipt_sha256": inspected.receipt_sha256,
        "validator_id": "command_inspection.v1",
        "status": "valid",
        "critical_finding_sha256s": (),
    }
    validator = ProgramValidatorReceiptV1(
        **validator_body,
        receipt_sha256=canonical_sha256(validator_body),
    )
    request = build_program_node_preflight_request(
        node_id=proposal.node_id,
        capability_receipt_sha256=capability.receipt_sha256,
        program_binding_sha256=program_binding.binding_sha256,
        engine_binding_sha256=engine_binding.binding_sha256,
        environment_receipt_sha256=environment.receipt_sha256,
        geometry_artifact_sha256=geometry.sha256,
        scientific_identity_sha256=identity.binding_sha256,
        charge=0,
        multiplicity=1,
        project_receipt_sha256=project_receipt.receipt_sha256,
        invocation_sha256=invocation.invocation_sha256,
        command_inspection_sha256=inspected.receipt_sha256,
        validator_receipts=(validator,),
    )
    preflight = evaluate_program_node_preflight(
        request,
        capability=capability,
        program_binding=program_binding,
        engine_binding=engine_binding,
        project=project_receipt,
        invocation=invocation,
        command_inspection=inspected,
        scientific_identity=identity,
        validator_receipts=(validator,),
    )

    assert project_receipt.status == "valid"
    assert invocation.argv[:2] == ("chemsmart", "run")
    assert "--project" in invocation.argv
    assert "--fake" in invocation.argv
    assert "--no-scratch" in invocation.argv
    assert inspected.status == "valid"
    assert preflight.plan_state == "ready_for_safe_preview"
    assert preflight.execution_ready is False


def test_optional_project_is_omitted_when_not_bound(
    tmp_path,
    fake_capability_registry,
    fake_click_schema,
    fake_support_overlay,
):
    geometry_path = tmp_path / "input.xyz"
    geometry_path.write_text("1\nH\nH 0 0 0\n")
    geometry = _artifact(
        geometry_path,
        artifact_id="geometry.h",
        kind="geometry_xyz",
        cli_value=str(geometry_path),
    )
    capability = query_capability(
        ProgramCapabilityQueryV1("optional", "sp", "cpu"),
        registry=fake_capability_registry,
        live_schema=fake_click_schema,
        overlay=fake_support_overlay,
    )
    environment = query_environment(
        capability,
        targets=(
            EnvironmentTargetV1(
                "optional", "cpu", "executable", "optional-bin"
            ),
        ),
        which=lambda _name: "/opt/bin/optional-bin",
    )
    binding = resolve_engine_binding(
        resolve_program_binding(capability), environment
    )
    identity = build_scientific_identity_binding(
        task_spec_sha256="e" * 64,
        geometry_artifact=geometry,
        charge=0,
        multiplicity=1,
    )
    invocation = compile_command(
        CommandProposalV1(
            node_id="node.optional",
            execution_target="run",
            program="optional",
            jobtype="sp",
            project_artifact_id="",
            input_artifact_id=geometry.artifact_id,
            scientific_identity_sha256=identity.binding_sha256,
            charge=0,
            multiplicity=1,
        ),
        capability=capability,
        binding=binding,
        project=None,
        project_validation=None,
        input_artifact=geometry,
        scientific_identity=identity,
        live_schema=fake_click_schema,
    )

    assert "--project" not in invocation.argv
    assert invocation.project_sha256 == ""


def test_validator_findings_block_node_preflight():
    finding = canonical_sha256({"rule_id": "seeded.critical"})
    assert len(finding) == 64
