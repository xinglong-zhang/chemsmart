from __future__ import annotations

from chemsmart.agent.knowledge_packs import (
    BUILTIN_PROGRAM_PACKS,
    activate_program_knowledge,
)
from chemsmart.agent.tool_specs import (
    build_command_compiled_tool_surface,
)


def test_program_knowledge_pack_is_advisory_and_triggered():
    receipt = activate_program_knowledge(
        request="Use GPU4PySCF with CUDA and cuTENSOR evidence.",
        program="pyscf",
        engine="gpu",
    )
    gpu_pack = next(
        item
        for item in BUILTIN_PROGRAM_PACKS
        if item.pack_id == "gpu4pyscf-advisory"
    )

    assert gpu_pack.readiness_authority is False
    assert receipt.advisory_only is True
    assert gpu_pack.pack_sha256 in receipt.activated_pack_sha256s


def test_model_tool_surface_uses_registry_programs_and_has_no_execution_tools(
    fake_capability_registry,
):
    surface = build_command_compiled_tool_surface(fake_capability_registry)
    names = {item["function"]["name"] for item in surface.tool_definitions}
    blob = str(surface.tool_definitions)

    assert "synthesize_command" in names
    assert "inspect_program_capability" in names
    assert "run_local" not in names
    assert "submit_hpc" not in names
    assert "native_input" not in blob
    synthesize = next(
        item
        for item in surface.tool_definitions
        if item["function"]["name"] == "synthesize_command"
    )
    preflight = next(
        item
        for item in surface.tool_definitions
        if item["function"]["name"] == "preflight_program_node"
    )
    assert "execution_target" not in (
        synthesize["function"]["parameters"]["properties"]
    )
    assert "validator_receipt_sha256s" not in (
        preflight["function"]["parameters"]["properties"]
    )
    bind_identity = next(
        item
        for item in surface.tool_definitions
        if item["function"]["name"] == "bind_scientific_identity"
    )
    plan_workflow = next(
        item
        for item in surface.tool_definitions
        if item["function"]["name"] == "plan_command_workflow"
    )
    assert (
        "never a project YAML or result"
        in bind_identity["function"]["description"]
    )
    assert "called again" in plan_workflow["function"]["description"]
