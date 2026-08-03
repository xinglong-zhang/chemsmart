from __future__ import annotations

import hashlib
import json

import pytest

from chemsmart.agent._contracts import ContractError, canonical_sha256
from chemsmart.agent.capabilities import (
    CapabilityQueryStatus,
    EnvironmentStatus,
    ProgramCapabilityQueryV1,
    ProgramSupportRuleV1,
    SupportLevel,
    build_command_compiled_preview_overlay,
    build_program_component_conformance_receipt,
    build_support_overlay,
    build_trusted_compute_environment_receipt,
    consume_pyscf_compute_environment_receipt,
    pyscf_compute_environment_target,
    query_capability,
    query_environment,
    resolve_engine_binding,
    resolve_program_binding,
)


def test_conformance_enables_only_observed_jobtype_and_engine(
    fake_capability_registry, fake_click_schema
):
    receipt = build_program_component_conformance_receipt(
        program="pyscf",
        registry_sha256=fake_capability_registry.registry_sha256,
        live_cli_schema_sha256=fake_click_schema.schema_sha256,
        fixture_bundle_sha256="1" * 64,
        covered_jobtypes=("sp",),
        covered_engines=("cpu",),
        compiler_receipt_sha256="2" * 64,
        preview_receipt_sha256="3" * 64,
        preflight_receipt_sha256="4" * 64,
        verifier_receipt_sha256="5" * 64,
        compiler_status="passed",
        preview_status="passed",
        preflight_status="passed",
        verifier_status="passed",
    )

    overlay = build_command_compiled_preview_overlay(
        fake_capability_registry,
        conformance_receipts=(receipt,),
        live_schema=fake_click_schema,
    )

    rule = overlay.get("pyscf")
    assert rule.support_level is SupportLevel.PREVIEW_ONLY
    assert rule.allowed_jobtypes == ("sp",)
    assert rule.allowed_engines == ("cpu",)


def test_capability_receipt_binds_registry_live_schema_and_join(
    fake_capability_registry, fake_click_schema, fake_support_overlay
):
    receipt = query_capability(
        ProgramCapabilityQueryV1("pyscf", "sp", "gpu"),
        registry=fake_capability_registry,
        live_schema=fake_click_schema,
        overlay=fake_support_overlay,
    )

    assert receipt.status is CapabilityQueryStatus.PREVIEW_ONLY
    assert receipt.registry_sha256 == fake_capability_registry.registry_sha256
    assert receipt.live_cli_schema_sha256 == fake_click_schema.schema_sha256
    assert receipt.joined_capability_sha256 == canonical_sha256(
        {
            "registry_sha256": receipt.registry_sha256,
            "live_cli_schema_sha256": receipt.live_cli_schema_sha256,
        }
    )


def test_live_registry_projects_loader_bounded_parameter_domains():
    from chemsmart.agent.capabilities import load_program_capabilities

    pyscf = load_program_capabilities().get("pyscf")

    assert dict(pyscf.project_parameter_domains)["ab_initio"] == ("hf",)


def test_support_overlay_can_narrow_but_not_broaden(
    fake_capability_registry,
):
    overlay = build_support_overlay(
        overlay_id="preview",
        registry=fake_capability_registry,
        rules=(
            ProgramSupportRuleV1(
                "pyscf",
                support_level=SupportLevel.PREVIEW_ONLY,
                allowed_jobtypes=("sp",),
                allowed_engines=("cpu",),
                compiler_evidence_sha256="1" * 64,
                preview_evidence_sha256="2" * 64,
                preflight_evidence_sha256="3" * 64,
                verifier_evidence_sha256="4" * 64,
            ),
        ),
    )
    assert overlay.get("pyscf").allowed_engines == ("cpu",)

    with pytest.raises(ContractError, match="broadens"):
        build_support_overlay(
            overlay_id="invalid",
            registry=fake_capability_registry,
            rules=(
                ProgramSupportRuleV1(
                    "pyscf",
                    support_level=SupportLevel.REFERENCE_ONLY,
                    allowed_jobtypes=("invented",),
                ),
            ),
        )


def test_gpu_environment_requires_complete_target_interpreter_evidence(
    fake_capability_registry, fake_click_schema, fake_support_overlay
):
    capability = query_capability(
        ProgramCapabilityQueryV1("pyscf", "sp", "gpu"),
        registry=fake_capability_registry,
        live_schema=fake_click_schema,
        overlay=fake_support_overlay,
    )
    target = pyscf_compute_environment_target("gpu")
    incomplete = build_trusted_compute_environment_receipt(
        program="pyscf",
        engine="gpu",
        compute_interpreter_sha256="a" * 64,
        dependency_versions={
            "pyscf": "2.14.0",
            "numpy": "1.26.4",
            "h5py": "3.12.1",
            "gpu4pyscf": "1.8.0",
            "cupy": "13.4.1",
        },
        gpu_evidence={
            "cuda_available": True,
            "device_available": True,
            "gpu4pyscf_distribution": True,
            "cutensor_runtime_available": True,
            "cutensor_compatible": False,
        },
        source_probe_id="pyscf-probe",
    )
    receipt = query_environment(
        capability, targets=(target,), compute_receipts=(incomplete,)
    )

    assert receipt.status is EnvironmentStatus.MISSING
    assert "environment.gpu.cutensor_compatible_missing" in receipt.rule_ids
    binding = resolve_engine_binding(
        resolve_program_binding(capability), receipt
    )
    assert binding.state == "blocked"
    assert binding.execution_ready is False


def test_pyscf_program_receipt_adapter_verifies_source_hash():
    raw = {
        "schema_version": "chemsmart.pyscf-environment.v1",
        "status": "available",
        "probe_returncode": 0,
        "interpreter": "/opt/target/bin/python",
        "interpreter_observed": "/opt/target/bin/python",
        "interpreter_sha256": "b" * 64,
        "packages": {
            "pyscf": "2.14.0",
            "numpy": "1.26.4",
            "h5py": "3.12.1",
        },
        "dependencies": {
            "pyscf": {"available": True, "version": "2.14.0"},
            "numpy": {"available": True, "version": "1.26.4"},
            "h5py": {"available": True, "version": "3.12.1"},
        },
        "solver_callables": {"geometric": {"callable": True}},
        "cuda_available": 0,
        "device_count": 0,
        "cutensor_compatible": False,
    }
    raw["receipt_sha256"] = hashlib.sha256(
        json.dumps(raw, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    adapted = consume_pyscf_compute_environment_receipt(raw, engine="cpu")

    assert adapted.compute_interpreter_sha256 == "b" * 64
    assert ("pyscf", "2.14.0") in adapted.dependency_versions

    raw["status"] = "changed"
    with pytest.raises(ContractError, match="digest mismatch"):
        consume_pyscf_compute_environment_receipt(raw, engine="cpu")
