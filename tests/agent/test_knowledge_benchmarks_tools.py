from __future__ import annotations

from chemsmart.agent._contracts import canonical_sha256
from chemsmart.agent.benchmarks import (
    BenchmarkCaseV1,
    BenchmarkConfigurationV1,
    BenchmarkOracleV1,
    paired_run_plan,
)
from chemsmart.agent.knowledge_packs import (
    BUILTIN_PROGRAM_PACKS,
    activate_program_knowledge,
)
from chemsmart.agent.tool_specs import (
    build_command_compiled_tool_surface,
    build_single_agent_baseline_tool_surface,
)


def test_program_knowledge_pack_is_advisory_and_triggered():
    receipt = activate_program_knowledge(
        request="Use GPU4PySCF with CUDA and cuTENSOR evidence.",
        program="pyscf",
        engine="gpu",
    )
    gpu_pack = next(
        item for item in BUILTIN_PROGRAM_PACKS if item.pack_id == "gpu4pyscf-advisory"
    )

    assert gpu_pack.readiness_authority is False
    assert receipt.advisory_only is True
    assert gpu_pack.pack_sha256 in receipt.activated_pack_sha256s


def test_paired_benchmark_scheduler_keeps_pairing_key_fixed():
    oracle = BenchmarkOracleV1("oracle", "1", "a" * 64)
    case_body = {
        "case_id": "case",
        "source_artifact_sha256s": ("b" * 64,),
        "task_spec_sha256": "c" * 64,
        "geometry_artifact_sha256": "d" * 64,
        "held_out": True,
        "oracle": oracle,
    }
    case = BenchmarkCaseV1(
        **case_body, case_sha256=canonical_sha256(case_body)
    )
    configs = []
    for name, enabled in (("single", False), ("decomposed", True)):
        body = {
            "configuration_id": name,
            "factor_levels": (("decomposition", enabled),),
            "provider_model": "deepseek-v4-flash",
            "prompt_sha256": "e" * 64,
            "tool_schema_sha256": "f" * 64,
            "joined_capability_sha256": "1" * 64,
            "budget_sha256": "2" * 64,
        }
        configs.append(
            BenchmarkConfigurationV1(
                **body, configuration_sha256=canonical_sha256(body)
            )
        )

    runs = paired_run_plan((case,), configs, repeats=2)

    assert len(runs) == 4
    assert runs[0].pairing_key == runs[1].pairing_key
    assert runs[2].pairing_key == runs[3].pairing_key
    assert runs[0].pairing_key != runs[2].pairing_key


def test_model_tool_surface_uses_registry_programs_and_has_no_execution_tools(
    fake_capability_registry,
):
    surface = build_command_compiled_tool_surface(fake_capability_registry)
    names = {
        item["function"]["name"] for item in surface.tool_definitions
    }
    blob = str(surface.tool_definitions)

    assert "synthesize_command" in names
    assert "inspect_program_capability" in names
    assert "run_local" not in names
    assert "submit_hpc" not in names
    assert "native_input" not in blob


def test_h0_and_h1_have_distinct_frozen_surfaces(fake_capability_registry):
    h0 = build_single_agent_baseline_tool_surface(fake_capability_registry)
    h1 = build_command_compiled_tool_surface(fake_capability_registry)
    h0_names = {item["function"]["name"] for item in h0.tool_definitions}
    candidate = next(
        item
        for item in h1.tool_definitions
        if item["function"]["name"] == "assess_program_candidate"
    )
    properties = candidate["function"]["parameters"]["properties"]

    assert h0.tool_schema_sha256 != h1.tool_schema_sha256
    assert "inspect_program_capability" not in h0_names
    assert "approval_ref" not in properties
    assert "functional_semantics_confirmed" not in properties
