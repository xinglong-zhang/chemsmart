"""The public v1 agent surface must match the support claim exactly."""

from __future__ import annotations

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.bootstrap import bootstrap_program_conformance
from chemsmart.agent.capabilities import (
    build_approved_execution_overlay,
    build_command_compiled_preview_overlay,
    build_program_component_conformance_receipt,
    load_program_capabilities,
)
from chemsmart.agent.cli_schema import build_live_click_schema
from chemsmart.agent.live_session import _preview_server_profile
from chemsmart.settings.capabilities import PROGRAM_CAPABILITIES


def test_v1_agent_programs_exclude_human_only_programs_and_gpu_engine():
    records = {
        item.program: item for item in load_program_capabilities().programs
    }

    assert set(records) == {"gaussian", "orca", "pyscf", "xtb"}
    assert records["pyscf"].engines == ("cpu",)
    assert "gpu" in PROGRAM_CAPABILITIES["pyscf"].engines
    assert "nciplot" in PROGRAM_CAPABILITIES
    assert records["gaussian"].execution_engine_job_pairs == ()
    assert records["orca"].execution_engine_job_pairs == (
        ("cpu", "opt"),
        ("cpu", "sp"),
    )
    assert records["pyscf"].execution_engine_job_pairs == (
        ("cpu", "hess"),
        ("cpu", "opt"),
        ("cpu", "sp"),
    )
    assert records["xtb"].execution_engine_job_pairs


def test_gaussian_preview_cannot_be_upgraded_to_agent_execution():
    registry = load_program_capabilities()
    schema = build_live_click_schema()
    receipt = build_program_component_conformance_receipt(
        program="gaussian",
        registry_sha256=registry.registry_sha256,
        live_cli_schema_sha256=schema.schema_sha256,
        fixture_bundle_sha256="1" * 64,
        covered_jobtypes=("sp",),
        covered_engines=("cpu",),
        compiler_receipt_sha256="2" * 64,
        preview_receipt_sha256="3" * 64,
        preflight_receipt_sha256="4" * 64,
        verifier_receipt_sha256="",
        compiler_status="passed",
        preview_status="passed",
        preflight_status="passed",
        verifier_status="not_observed",
    )
    preview = build_command_compiled_preview_overlay(
        registry, conformance_receipts=(receipt,), live_schema=schema
    )

    with pytest.raises(ContractError, match="preview-only engine-job pair"):
        build_approved_execution_overlay(
            registry=registry,
            preview_overlay=preview,
            approved_nodes=(("gaussian", "sp", "cpu"),),
            execution_evidence_sha256="5" * 64,
        )


def test_bootstrap_verifier_is_explicitly_unobserved(tmp_path):
    registry = load_program_capabilities()
    schema = build_live_click_schema()
    xyz = tmp_path / "h2.xyz"
    xyz.write_text(
        "2\nH2\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8"
    )
    server = tmp_path / "preview-server.yaml"
    server.write_text(_preview_server_profile(), encoding="utf-8")
    receipt = bootstrap_program_conformance(
        program="xtb",
        engine="cpu",
        jobtypes=("sp",),
        input_path=xyz,
        project_path=None,
        server_path=server,
        charge=0,
        multiplicity=1,
        registry_sha256=registry.registry_sha256,
        live_schema=schema,
    )

    assert receipt.compiler_status == "passed"
    assert receipt.preview_status == "passed"
    assert receipt.preflight_status == "passed"
    assert receipt.verifier_status == "not_observed"
    assert receipt.verifier_receipt_sha256 == ""
    overlay = build_command_compiled_preview_overlay(
        registry, conformance_receipts=(receipt,), live_schema=schema
    )
    assert overlay.get("xtb").verifier_evidence_sha256 == ""


def test_an_observed_failed_verifier_still_blocks_preview_support():
    registry = load_program_capabilities()
    schema = build_live_click_schema()
    receipt = build_program_component_conformance_receipt(
        program="xtb",
        registry_sha256=registry.registry_sha256,
        live_cli_schema_sha256=schema.schema_sha256,
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
        verifier_status="failed",
    )

    overlay = build_command_compiled_preview_overlay(
        registry, conformance_receipts=(receipt,), live_schema=schema
    )
    assert overlay.get("xtb").support_level.value == "reference_only"
