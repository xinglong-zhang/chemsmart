#!/usr/bin/env python3
"""Exercise the installed image's preview-to-materialization transaction.

This is a non-computing qualification.  The PySCF CLI is invoked in exact
``--fake --no-scratch`` mode against the public water fixture; no chemistry
engine or provider call is made.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from types import SimpleNamespace

from chemsmart.agent._contracts import (
    TrustedArtifactRefV1,
    canonical_sha256,
    file_sha256,
)
from chemsmart.agent.capabilities import (
    CapabilityQueryStatus,
    CapabilityQueryV1,
    ProgramSupportRuleV1,
    SupportLevel,
    build_support_overlay,
    load_program_capabilities,
    query_capability,
    query_environment,
    resolve_engine_binding,
    resolve_program_binding,
)
from chemsmart.agent.cli_schema import build_live_click_schema
from chemsmart.agent.commands import (
    CommandProposalV1,
    build_scientific_identity_binding,
    compile_command,
    inspect_command,
)
from chemsmart.agent.preflight import (
    ProgramNodePreflightReceiptV1,
    validator_receipt_from_safe_preview,
)
from chemsmart.agent.preview import (
    SafePreviewReceiptV1,
    execute_safe_preview,
)
from chemsmart.agent.program_verifiers import build_preview_expectation
from chemsmart.agent.projects import validate_project_yaml
from chemsmart.agent.tool_runtime import (
    CommandCompiledToolHostV1,
    _CommandContext,
)
from chemsmart.agent.workflows import (
    ScientificWorkflowNodeV2,
    build_scientific_workflow_plan,
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def artifact(path: Path, *, artifact_id: str, kind: str) -> TrustedArtifactRefV1:
    resolved = path.resolve(strict=True)
    return TrustedArtifactRefV1(
        artifact_id=artifact_id,
        kind=kind,
        sha256=file_sha256(resolved),
        size_bytes=resolved.stat().st_size,
        path=str(resolved),
        cli_value=str(resolved),
    )


def qualify(*, project_path: Path, geometry_path: Path) -> dict[str, object]:
    live_schema = build_live_click_schema()
    registry = load_program_capabilities()
    evidence = {
        name: canonical_sha256(
            {
                "qualification": "cpu-image-transaction-regression-v1",
                "component": name,
            }
        )
        for name in ("compiler", "preview", "preflight", "verifier")
    }
    overlay = build_support_overlay(
        overlay_id="cpu_image_transaction_qualification",
        registry=registry,
        rules=(
            ProgramSupportRuleV1(
                program="pyscf",
                support_level=SupportLevel.PREVIEW_ONLY,
                allowed_jobtypes=("sp",),
                allowed_engines=("cpu",),
                allowed_engine_job_pairs=(("cpu", "sp"),),
                compiler_evidence_sha256=evidence["compiler"],
                preview_evidence_sha256=evidence["preview"],
                preflight_evidence_sha256=evidence["preflight"],
                verifier_evidence_sha256=evidence["verifier"],
                rule_ids=("qualification.image_transaction",),
            ),
        ),
    )
    capability = query_capability(
        CapabilityQueryV1("pyscf", "sp", "cpu"),
        registry=registry,
        overlay=overlay,
        live_schema=live_schema,
    )
    require(
        capability.status is CapabilityQueryStatus.PREVIEW_ONLY,
        f"Unexpected qualification capability: {capability.status}",
    )

    project = artifact(
        project_path,
        artifact_id="project.pyscf.water",
        kind="project_yaml",
    )
    geometry = artifact(
        geometry_path,
        artifact_id="geometry.water",
        kind="geometry_xyz",
    )
    project_validation = validate_project_yaml(project, capability=capability)
    require(
        project_validation.status == "valid",
        f"Project validation failed: {project_validation.rule_ids}",
    )
    environment = query_environment(capability)
    program_binding = resolve_program_binding(capability)
    engine_binding = resolve_engine_binding(program_binding, environment)

    task_spec_sha256 = canonical_sha256(
        {"task": "CPU image safe-preview transaction qualification"}
    )
    identity = build_scientific_identity_binding(
        task_spec_sha256=task_spec_sha256,
        geometry_artifact=geometry,
        charge=0,
        multiplicity=1,
    )
    proposal = CommandProposalV1(
        node_id="pyscf-water-sp",
        execution_target="run",
        program="pyscf",
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
        project_validation=project_validation,
        input_artifact=geometry,
        scientific_identity=identity,
        live_schema=live_schema,
        server="local",
    )
    inspection = inspect_command(invocation, live_schema=live_schema)
    require(
        inspection.status == "valid",
        f"Command inspection failed: {inspection.rule_ids}",
    )

    expectation = build_preview_expectation(
        program="pyscf",
        jobtype="sp",
        input_artifact=geometry,
        project=project_validation,
        charge=0,
        multiplicity=1,
    )
    safe_preview = execute_safe_preview(
        invocation,
        input_artifact=geometry,
        project_artifact=project,
        expectation=expectation,
    )
    require(isinstance(safe_preview, SafePreviewReceiptV1), "Wrong preview type")
    require(
        safe_preview.status == "previewed"
        and safe_preview.program_validation_status == "valid",
        f"Safe preview failed: {safe_preview.rule_ids}",
    )
    validator = validator_receipt_from_safe_preview(
        node_id=proposal.node_id,
        program="pyscf",
        scientific_identity_sha256=identity.binding_sha256,
        safe_preview=safe_preview,
    )
    plan = build_scientific_workflow_plan(
        workflow_id="cpu-image-transaction",
        task_spec_sha256=task_spec_sha256,
        scientific_identity_sha256=identity.binding_sha256,
        nodes=(
            ScientificWorkflowNodeV2(
                node_id=proposal.node_id,
                stage="sp",
                requested_program="pyscf",
                program="pyscf",
                engine="cpu",
                project_role=project.artifact_id,
                unresolved_fields=(),
            ),
        ),
    )

    # Use the real host transition.  Calling build_materialized_workflow here
    # would bypass _node_is_previewed, the exact old-code failure this image
    # label promises has regressed.
    context = _CommandContext(
        proposal=proposal,
        capability=capability,
        program_binding=program_binding,
        engine_binding=engine_binding,
        project_artifact=project,
        project_validation=project_validation,
        input_artifact=geometry,
        scientific_identity=identity,
    )
    materialized_records = []
    host = CommandCompiledToolHostV1(
        event_store=SimpleNamespace(
            read_events=lambda: (),
            record_materialized_workflow=lambda **values: (
                materialized_records.append(values["workflow"])
            )
        ),
        artifacts={project.artifact_id: project, geometry.artifact_id: geometry},
        scientific_identities={identity.binding_sha256: identity},
        registry=registry,
        live_schema=live_schema,
        task_spec_sha256s=(task_spec_sha256,),
        scientific_workflow_plan=plan,
    )
    host.capabilities[capability.receipt_sha256] = capability
    host.program_bindings[program_binding.binding_sha256] = program_binding
    host.engine_bindings[engine_binding.binding_sha256] = engine_binding
    host.project_validations[project_validation.receipt_sha256] = (
        project_validation
    )
    host.invocations[invocation.invocation_sha256] = invocation
    host._command_contexts[invocation.invocation_sha256] = context
    host.command_inspections[inspection.receipt_sha256] = inspection
    host.safe_previews[safe_preview.receipt_sha256] = safe_preview
    host.validators[validator.receipt_sha256] = validator
    host._emit = lambda *_args, **_kwargs: None

    preflight = host._preflight_program_node(
        "cpu-image-transaction",
        {
            "node_id": proposal.node_id,
            "capability_receipt_sha256": capability.receipt_sha256,
            "program_binding_sha256": program_binding.binding_sha256,
            "engine_binding_sha256": engine_binding.binding_sha256,
            "invocation_sha256": invocation.invocation_sha256,
            "command_inspection_receipt_sha256": inspection.receipt_sha256,
            "scientific_identity_sha256": identity.binding_sha256,
            "geometry_artifact_sha256": geometry.sha256,
            "charge": 0,
            "multiplicity": 1,
        },
    )
    require(isinstance(preflight, ProgramNodePreflightReceiptV1), "Wrong preflight type")
    require(
        preflight.plan_state == "previewed"
        and preflight.safe_preview_receipt_sha256
        == safe_preview.receipt_sha256,
        f"Preflight did not bind the preview: {preflight.rule_ids}",
    )
    require(
        host._node_is_previewed(
            proposal.node_id,
            plan_sha256=plan.plan_sha256,
        ),
        "Host did not recognize the exact previewed invocation",
    )
    require(
        len(materialized_records) == 1,
        "Host preflight did not record exactly one materialized workflow",
    )
    materialized = materialized_records[0]
    node = materialized.nodes[0]
    require(
        materialized.status == "previewed"
        and node.state == "previewed"
        and node.preflight_receipt_sha256 == preflight.receipt_sha256,
        "Materialization did not retain exact preflight evidence",
    )
    return {
        "schema_version": "chemsmart.cpu-image-transaction-qualification.v1",
        "status": "qualified",
        "safe_preview_receipt_sha256": safe_preview.receipt_sha256,
        "program_node_preflight_receipt_sha256": preflight.receipt_sha256,
        "materialized_workflow_sha256": materialized.materialized_sha256,
        "generated_artifacts": tuple(
            item.relative_path for item in safe_preview.artifacts
        ),
        "engine_executed": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path("/opt/chemsmart/qualification/pyscf-water.yaml"),
    )
    parser.add_argument(
        "--geometry",
        type=Path,
        default=Path("/opt/chemsmart/qualification/water.xyz"),
    )
    args = parser.parse_args()
    settings_logger = logging.getLogger("chemsmart.jobs.pyscf.settings")
    prior_level = settings_logger.level
    settings_logger.setLevel(logging.ERROR)
    try:
        result = qualify(
            project_path=args.project,
            geometry_path=args.geometry,
        )
    finally:
        settings_logger.setLevel(prior_level)
    print(
        json.dumps(
            result,
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
