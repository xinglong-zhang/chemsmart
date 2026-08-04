from __future__ import annotations

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.feedback import project_tool_feedback


_DIGEST_A = "a" * 64
_DIGEST_B = "b" * 64
_DIGEST_C = "c" * 64


def _tool_result(tool: str, payload: dict, *, status: str = "ok") -> dict:
    return {
        "schema_version": "chemsmart.tool-result.v1",
        "tool": tool,
        "status": status,
        "result": payload,
    }


def test_capability_projection_keeps_query_and_receipt_not_inventory():
    result = _tool_result(
        "inspect_program_capability",
        {
            "schema_version": "chemsmart.capability-query-receipt.v1",
            "query": {"program": "pyscf", "jobtype": "hess", "engine": "cpu"},
            "status": "supported",
            "receipt_sha256": _DIGEST_A,
            "joined_capability_sha256": _DIGEST_B,
            "live_cli_schema_sha256": _DIGEST_C,
            "effective_jobtypes": ["sp", "opt", "hess", "td"] * 100,
            "effective_engines": ["cpu", "gpu4pyscf"] * 100,
            "capability": {
                "program": "pyscf",
                "jobtypes": ["sp", "opt", "hess", "td"] * 100,
                "engines": ["cpu", "gpu4pyscf"] * 100,
                "project_owned_parameters": [
                    f"setting_{index}" for index in range(300)
                ],
            },
            "rule_ids": ["capability.query.supported"],
        },
    )

    projected = project_tool_feedback(
        tool="inspect_program_capability", result=result, mode="causal-v1"
    )
    payload = projected.content["result"]

    assert payload["query"] == {
        "engine": "cpu",
        "jobtype": "hess",
        "program": "pyscf",
    }
    assert payload["receipt_sha256"] == _DIGEST_A
    assert payload["joined_capability_sha256"] == _DIGEST_B
    assert "effective_jobtypes" not in payload
    assert "jobtypes" not in payload.get("capability", {})
    assert projected.receipt.bytes_reduced > 0
    assert projected.receipt.reduction_basis_points > 5000


def test_project_validation_projection_retains_scientific_settings_and_units():
    settings = [
        ["basis", "def2-svp"],
        ["charge", 0],
        ["functional", "b3lypg"],
        ["multiplicity", 1],
        ["temperature", {"unit": "kelvin", "value": 298.15}],
    ]
    result = _tool_result(
        "validate_project_yaml",
        {
            "schema_version": "chemsmart.project-validation-receipt.v1",
            "project_artifact_id": "project.water",
            "project_sha256": _DIGEST_A,
            "capability_receipt_sha256": _DIGEST_B,
            "program": "pyscf",
            "jobtype": "sp",
            "loader_id": "pyscf-project-v1",
            "settings_sha256": _DIGEST_C,
            "settings": settings,
            "scientific_materializations": [
                {
                    "schema_version": (
                        "chemsmart.pyscf-functional-resolution.v1"
                    ),
                    "requested_literal": "b3lyp",
                    "applied_xc": "b3lypg",
                    "correlation_convention": "vwn3_gaussian",
                    "rule_id": "pyscf.functional.b3lypg_vwn3_gaussian",
                    "receipt_sha256": "e" * 64,
                    "evidence_ref": "functional_resolution:" + "e" * 64,
                }
            ],
            "decision_binding": {
                "schema_version": "chemsmart.scientific-decision-binding.v1",
                "status": "required_if_rendering_implementation_semantics",
                "rule_id": "scientific.functional_resolution.decision_binding",
                "next_tool": "record_scientific_decision",
                "evidence_refs": ["functional_resolution:" + "e" * 64],
                "receipt_sha256": "f" * 64,
            },
            "status": "valid",
            "diagnostic": "loader accepted exact stage settings",
            "rule_ids": ["project.loader.valid"],
            "receipt_sha256": "d" * 64,
            "rendered_yaml": "sp:\n  functional: b3lypg\n" * 1000,
            "native_input": "model must never receive this\n" * 1000,
        },
    )

    projected = project_tool_feedback(
        tool="validate_project_yaml", result=result, mode="causal-v1"
    )
    payload = projected.content["result"]

    assert payload["project_artifact_id"] == "project.water"
    assert payload["settings"] == settings
    assert payload["scientific_materializations"][0]["applied_xc"] == "b3lypg"
    assert payload["scientific_materializations"][0][
        "correlation_convention"
    ] == "vwn3_gaussian"
    assert payload["decision_binding"]["next_tool"] == (
        "record_scientific_decision"
    )
    assert payload["decision_binding"]["evidence_refs"] == [
        "functional_resolution:" + "e" * 64
    ]
    assert payload["status"] == "valid"
    assert "rendered_yaml" not in payload
    assert "native_input" not in payload
    assert projected.receipt.provider_feedback_bytes < projected.receipt.canonical_bytes
    assert projected.receipt.reduction_basis_points > 5000


def test_command_chain_projection_retains_only_next_call_references():
    rendered = project_tool_feedback(
        tool="render_project_yaml",
        result=_tool_result(
            "render_project_yaml",
            {
                "program": "pyscf",
                "status": "candidate_rendered",
                "document_sha256": _DIGEST_A,
                "rendered_sha256": _DIGEST_B,
                "receipt_sha256": _DIGEST_C,
                "rendered_yaml": "sp:\n  basis: def2-svp\n" * 500,
            },
        ),
        mode="causal-v1",
    )
    assert rendered.content["result"]["receipt_sha256"] == _DIGEST_C
    assert "rendered_yaml" not in rendered.content["result"]

    promoted = project_tool_feedback(
        tool="promote_project_yaml",
        result=_tool_result(
            "promote_project_yaml",
            {
                "artifact": {
                    "artifact_id": "project.water",
                    "kind": "project_yaml",
                    "sha256": _DIGEST_A,
                    "size_bytes": 123,
                    "path": "/private/disposable/project.yaml",
                    "cli_value": "/private/disposable/project.yaml",
                },
                "promotion": {
                    "status": "materialized",
                    "validation_status": "pending",
                    "render_receipt_sha256": _DIGEST_C,
                    "project_artifact_sha256": _DIGEST_A,
                    "receipt_sha256": _DIGEST_B,
                },
            },
        ),
        mode="causal-v1",
    )
    artifact = promoted.content["result"]["artifact"]
    assert artifact == {
        "artifact_id": "project.water",
        "kind": "project_yaml",
        "sha256": _DIGEST_A,
        "size_bytes": 123,
    }

    compiled = project_tool_feedback(
        tool="synthesize_command",
        result=_tool_result(
            "synthesize_command",
            {
                "invocation": {
                    "node_id": "sp.initial",
                    "command_path": ["run", "pyscf", "sp"],
                    "input_sha256": _DIGEST_A,
                    "project_sha256": _DIGEST_B,
                    "scientific_identity_sha256": _DIGEST_C,
                    "invocation_sha256": "d" * 64,
                    "status": "compiled",
                    "argv": ["chemsmart", "run", "--fake", "pyscf", "sp"],
                    "display_command": "chemsmart run --fake pyscf sp",
                },
                "inspection": {
                    "status": "valid",
                    "invocation_sha256": "d" * 64,
                    "receipt_sha256": "e" * 64,
                    "rule_ids": ["command.schema.valid"],
                },
            },
        ),
        mode="causal-v1",
    )
    invocation = compiled.content["result"]["invocation"]
    assert invocation["invocation_sha256"] == "d" * 64
    assert invocation["command_path"] == ["run", "pyscf", "sp"]
    assert "argv" not in invocation
    assert "display_command" not in invocation


@pytest.mark.parametrize("nested_status", ["invalid", "blocked", "failed"])
def test_preview_projection_never_promotes_nested_red_status(nested_status: str):
    result = _tool_result(
        "preview_command",
        {
            "safe_preview": {
                "status": nested_status,
                "invocation_sha256": _DIGEST_A,
                "input_sha256": _DIGEST_B,
                "project_sha256": _DIGEST_C,
                "receipt_sha256": "d" * 64,
                "program_validation_status": "invalid",
                "rule_ids": ["preview.program_validator.red"],
                "artifacts": [
                    {
                        "relative_path": f"generated-{index}.py",
                        "sha256": "e" * 64,
                        "size_bytes": 50_000,
                        "raw_output": "x" * 20_000,
                    }
                    for index in range(10)
                ],
            },
            "validator": {
                "node_id": "sp.initial",
                "status": "invalid",
                "source_receipt_sha256": "d" * 64,
                "receipt_sha256": "f" * 64,
                "rule_ids": ["preview.program_validator.red"],
            },
        },
    )

    projected = project_tool_feedback(
        tool="preview_command", result=result, mode="causal-v1"
    )

    expected = "failed" if nested_status == "failed" else nested_status
    assert projected.content["status"] == expected
    assert projected.action_signature is not None
    assert projected.action_signature.outcome_status == expected
    assert "artifacts" not in projected.content["result"]["safe_preview"]
    assert projected.receipt.bytes_reduced > 0
    assert projected.receipt.reduction_basis_points > 5000


def test_registered_tools_project_public_error_envelopes_without_downgrade():
    tools = (
        "inspect_program_capability",
        "inspect_program_environment",
        "assess_program_candidate",
        "render_project_yaml",
        "promote_project_yaml",
        "bind_scientific_identity",
        "read_project_yaml",
        "validate_project_yaml",
        "plan_command_workflow",
        "synthesize_command",
        "repair_command",
        "preview_command",
        "preflight_program_node",
        "record_scientific_decision",
        "inspect_calculation_artifact",
        "execute_approved_program_node",
    )
    for tool in tools:
        error = {
            "schema_version": "chemsmart.tool-result.v1",
            "tool": tool,
            "status": "rejected",
            "error_class": "ContractError",
            "message": "typed prerequisite is absent",
            "rule_ids": ["tool.dispatch.rejected"],
        }
        projected = project_tool_feedback(
            tool=tool, result=error, mode="causal-v1"
        )
        payload = projected.content["result"]
        assert projected.content["status"] == "invalid"
        assert payload["error_class"] == "ContractError"
        assert payload["message"] == "typed prerequisite is absent"
        assert payload["rule_ids"] == ["tool.dispatch.rejected"]


def test_causal_mode_rejects_tool_identity_mismatch_and_unknown_tool():
    with pytest.raises(ContractError, match="identity mismatch"):
        project_tool_feedback(
            tool="preview_command",
            result=_tool_result("synthesize_command", {}),
            mode="causal-v1",
        )
    with pytest.raises(ContractError, match="registered tool projection"):
        project_tool_feedback(
            tool="future_unregistered_tool",
            result=_tool_result("future_unregistered_tool", {}),
            mode="causal-v1",
        )


def test_full_mode_remains_exact_even_for_an_unknown_tool():
    result = _tool_result(
        "future_unregistered_tool",
        {"arbitrary": {"native_input": "preserved in full mode"}},
    )
    projected = project_tool_feedback(
        tool="future_unregistered_tool", result=result, mode="full-v1"
    )

    assert projected.content == result
    assert projected.receipt.verdict == "full_result_preserved"
    assert projected.receipt.bytes_reduced == 0
    assert projected.action_signature is None
