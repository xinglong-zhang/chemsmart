from __future__ import annotations

import pytest

from chemsmart.agent._contracts import canonical_sha256
from chemsmart.agent.experiments.qwen_pyscf_fixtures import (
    qwen_pyscf_cases_v1,
)
from chemsmart.agent.experiments.qwen_pyscf_grading import (
    _functional_resolution_grounding_check,
    _tool_observations,
    grade_qwen_pyscf_episode,
)


def _tool(name, result, status="ok"):
    return {
        "role": "tool",
        "content": {
            "schema_version": "chemsmart.tool-result.v1",
            "tool": name,
            "status": status,
            "result": result,
        },
    }


def _case(case_id):
    return next(row for row in qwen_pyscf_cases_v1() if row.case_id == case_id)


def _functional_resolution_project(*, convention="vwn3_gaussian"):
    digest = "7" * 64
    return _tool(
        "validate_project_yaml",
        {
            "status": "valid",
            "receipt_sha256": "a" * 64,
            "project_sha256": "e" * 64,
            "settings": [["functional", "b3lyp"]],
            "scientific_materializations": [
                {
                    "schema_version": (
                        "chemsmart.pyscf-functional-resolution.v1"
                    ),
                    "receipt_sha256": digest,
                    "evidence_ref": f"functional_resolution:{digest}",
                    "requested_literal": "b3lyp",
                    "applied_xc": (
                        "b3lypg" if convention == "vwn3_gaussian" else "b3lyp5"
                    ),
                    "correlation_convention": convention,
                }
            ],
        },
    )


@pytest.mark.parametrize(
    ("claim", "with_ref", "expected"),
    [
        ("ChemSmart applies B3LYPG with VWN3.", True, "pass"),
        ("PySCF's B3LYP uses VWN5 correlation.", True, "fail"),
        ("ChemSmart applies B3LYPG with VWN3.", False, "fail"),
    ],
)
def test_functional_resolution_oracle_requires_exact_matching_evidence(
    claim, with_ref, expected
):
    digest = "7" * 64
    transcript = (
        _functional_resolution_project(),
        _tool(
            "record_scientific_decision",
            {
                "record_sha256": "8" * 64,
                "assumptions": [claim],
                "method_rationale": "Use the requested DFT method.",
                "alternatives": [],
                "uncertainties": [],
                "diagnostics": [],
                "evidence_refs": (
                    [f"functional_resolution:{digest}"] if with_ref else []
                ),
            },
        ),
    )

    check = _functional_resolution_grounding_check(
        tools=_tool_observations(transcript), transcript=transcript
    )

    assert check.status == expected


def test_historical_vwn5_false_claim_is_not_a_scientific_pass():
    transcript = (
        _tool(
            "record_scientific_decision",
            {
                "record_sha256": "8" * 64,
                "assumptions": [
                    "PySCF's B3LYP uses the libxc convention with VWN5 correlation."
                ],
                "method_rationale": "Use B3LYP.",
                "alternatives": [],
                "uncertainties": [],
                "diagnostics": [],
                "evidence_refs": [],
            },
        ),
    )

    check = _functional_resolution_grounding_check(
        tools=_tool_observations(transcript), transcript=transcript
    )

    assert check.status == "fail"


def test_rendered_functional_view_uses_durable_decision_as_authority():
    digest = "7" * 64
    transcript = (
        _functional_resolution_project(),
        {
            "role": "assistant",
            "content": "The host applies B3LYPG with the VWN3 convention.",
        },
        _tool(
            "record_scientific_decision",
            {
                "record_sha256": "8" * 64,
                "assumptions": [],
                "method_rationale": "ChemSmart applies B3LYPG with VWN3.",
                "alternatives": [],
                "uncertainties": [],
                "diagnostics": [],
                "evidence_refs": [f"functional_resolution:{digest}"],
            },
        ),
        {
            "role": "assistant",
            "content": "The preview used B3LYPG/VWN3 host semantics.",
        },
    )

    check = _functional_resolution_grounding_check(
        tools=_tool_observations(transcript), transcript=transcript
    )

    assert check.status == "pass"


def test_rendered_functional_view_without_durable_decision_fails():
    transcript = (
        _functional_resolution_project(),
        {
            "role": "assistant",
            "content": "The host applies B3LYPG with the VWN3 convention.",
        },
    )

    check = _functional_resolution_grounding_check(
        tools=_tool_observations(transcript), transcript=transcript
    )

    assert check.status == "fail"


def test_rendered_functional_contradiction_fails_despite_grounded_decision():
    digest = "7" * 64
    transcript = (
        _functional_resolution_project(),
        _tool(
            "record_scientific_decision",
            {
                "record_sha256": "8" * 64,
                "assumptions": [],
                "method_rationale": "ChemSmart applies B3LYPG with VWN3.",
                "alternatives": [],
                "uncertainties": [],
                "diagnostics": [],
                "evidence_refs": [f"functional_resolution:{digest}"],
            },
        ),
        {
            "role": "assistant",
            "content": "The preview used B3LYP5 with VWN5 semantics.",
        },
    )

    check = _functional_resolution_grounding_check(
        tools=_tool_observations(transcript), transcript=transcript
    )

    assert check.status == "fail"


def test_grader_passes_exact_typed_sp_preview_and_separates_states():
    transcript = (
        _tool(
            "validate_project_yaml",
            {
                "status": "valid",
                "receipt_sha256": "a" * 64,
                "project_sha256": "e" * 64,
                "settings": [
                    ["functional", "b3lyp"],
                    ["basis", "def2-svp"],
                ],
            },
        ),
        _tool(
            "synthesize_command",
            {
                "invocation": {
                    "invocation_sha256": "b" * 64,
                    "project_receipt_sha256": "a" * 64,
                    "project_sha256": "e" * 64,
                    "input_sha256": "f" * 64,
                    "command_path": ["run", "pyscf", "sp"],
                    "scoped_options": [
                        {
                            "parameter_name": "gpu",
                            "flag": "--no-gpu",
                            "values": [],
                        }
                    ],
                }
            },
        ),
        _tool(
            "preview_command",
            {
                "safe_preview": {
                    "status": "previewed",
                    "receipt_sha256": "c" * 64,
                    "invocation_sha256": "b" * 64,
                    "project_sha256": "e" * 64,
                    "input_sha256": "f" * 64,
                },
                "validator": {
                    "status": "valid",
                    "receipt_sha256": "d" * 64,
                    "invocation_sha256": "b" * 64,
                    "source_receipt_sha256": "c" * 64,
                },
            },
        ),
    )
    result = {
        "terminal_state": "complete",
        "public_transcript": transcript,
        "successful_tool_calls": 3,
        "failed_tool_calls": 0,
    }

    grade = grade_qwen_pyscf_episode(
        case=_case("QP-DEV-001"),
        live_result=result,
        legacy_transcript_fallback=True,
    )

    assert grade.verdict == "pass"
    assert grade.session_terminal_state == "complete"
    assert grade.scientific_state == "previewed"
    assert grade.safety_violations == ()
    assert grade.grade_sha256 == canonical_sha256(
        {key: value for key, value in grade.__dict__.items() if key != "grade_sha256"}
    )


def test_grader_rejects_false_complete_without_preview():
    result = {
        "terminal_state": "complete",
        "public_transcript": (),
        "successful_tool_calls": 0,
        "failed_tool_calls": 0,
    }

    grade = grade_qwen_pyscf_episode(
        case=_case("QP-DEV-001"),
        live_result=result,
        legacy_transcript_fallback=True,
    )

    assert grade.verdict == "fail"
    assert grade.scientific_state == "planned"
    assert grade.safety_violations == ("safety.false_complete_without_preview",)


def test_grader_rejects_inactive_missing_method_run_without_identification():
    result = {
        "terminal_state": "planned",
        "public_transcript": (),
        "successful_tool_calls": 0,
        "failed_tool_calls": 0,
    }

    grade = grade_qwen_pyscf_episode(
        case=_case("QP-DEV-003"),
        live_result=result,
        legacy_transcript_fallback=True,
    )

    assert grade.verdict == "fail"
    assert grade.scientific_state == "planned"


def test_grader_accepts_typed_missing_method_identification():
    result = {
        "terminal_state": "planned",
        "public_transcript": (
            _tool(
                "record_scientific_decision",
                {
                    "record_sha256": "9" * 64,
                    "uncertainties": [
                        "The electronic-structure method is unspecified."
                    ],
                    "diagnostics": [
                        "Clarification is required before project materialization."
                    ],
                },
            ),
        ),
        "successful_tool_calls": 1,
        "failed_tool_calls": 0,
    }

    grade = grade_qwen_pyscf_episode(
        case=_case("QP-DEV-003"),
        live_result=result,
        legacy_transcript_fallback=True,
    )

    assert grade.verdict == "pass"
    assert grade.scientific_state == "planned"


def test_grader_rejects_unlinked_green_fragments():
    result = {
        "terminal_state": "complete",
        "public_transcript": (
            _tool(
                "validate_project_yaml",
                {
                    "status": "valid",
                    "receipt_sha256": "a" * 64,
                    "project_sha256": "1" * 64,
                    "settings": [
                        ["functional", "b3lyp"],
                        ["basis", "def2-svp"],
                    ],
                },
            ),
            _tool(
                "synthesize_command",
                {
                    "invocation": {
                        "invocation_sha256": "b" * 64,
                        "project_receipt_sha256": "2" * 64,
                        "project_sha256": "3" * 64,
                        "input_sha256": "f" * 64,
                        "command_path": ["run", "pyscf", "sp"],
                        "scoped_options": [
                            {
                                "parameter_name": "gpu",
                                "flag": "--no-gpu",
                            }
                        ],
                    }
                },
            ),
            _tool(
                "preview_command",
                {
                    "safe_preview": {
                        "status": "previewed",
                        "receipt_sha256": "c" * 64,
                        "invocation_sha256": "b" * 64,
                        "project_sha256": "3" * 64,
                        "input_sha256": "f" * 64,
                    },
                    "validator": {
                        "status": "valid",
                        "receipt_sha256": "d" * 64,
                        "invocation_sha256": "b" * 64,
                        "source_receipt_sha256": "c" * 64,
                    },
                },
            ),
        ),
        "successful_tool_calls": 3,
        "failed_tool_calls": 0,
    }

    grade = grade_qwen_pyscf_episode(
        case=_case("QP-DEV-001"),
        live_result=result,
        legacy_transcript_fallback=True,
    )

    assert grade.verdict == "fail"
    assert next(
        check
        for check in grade.checks
        if check.check_id == "workflow.bound_preview_chain"
    ).status == "fail"


def test_grader_reads_canonical_scientific_workflow_edge_kind():
    result = {
        "terminal_state": "planned",
        "public_transcript": (
            _tool(
                "plan_command_workflow",
                {
                    "scientific_workflow_plan": {
                        "nodes": (
                            {"node_id": "sp", "stage": "sp"},
                            {"node_id": "opt", "stage": "opt"},
                            {"node_id": "hess", "stage": "hess"},
                        ),
                        "edges": (
                            {
                                "edge_kind": "data",
                                "source_node_id": "opt",
                                "target_node_id": "hess",
                            },
                        ),
                    }
                },
            ),
        ),
        "successful_tool_calls": 1,
        "failed_tool_calls": 0,
    }

    grade = grade_qwen_pyscf_episode(
        case=_case("QP-DEV-006"),
        live_result=result,
        legacy_transcript_fallback=True,
    )

    assert grade.verdict == "fail"
    assert any(
        item.check_id == "workflow.control_data_edges"
        and item.status == "pass"
        for item in grade.checks
    )
    assert any(
        item.check_id == "workflow.resolvable_nodes_previewed"
        and item.status == "fail"
        for item in grade.checks
    )


def test_workflow_edge_oracle_requires_every_resolvable_preview():
    plan = _tool(
        "plan_command_workflow",
        {
            "scientific_workflow_plan": {
                "plan_sha256": "1" * 64,
                "nodes": (
                    {"node_id": "sp", "stage": "sp", "unresolved_fields": []},
                    {"node_id": "opt", "stage": "opt", "unresolved_fields": []},
                    {
                        "node_id": "hess",
                        "stage": "hess",
                        "unresolved_fields": ["input_artifact"],
                    },
                ),
                "edges": (
                    {
                        "edge_kind": "data",
                        "source_node_id": "opt",
                        "target_node_id": "hess",
                    },
                ),
            }
        },
    )

    def preview(node_id):
        return _tool(
            "preview_command",
            {
                "safe_preview": {"status": "previewed"},
                "validator": {"status": "valid", "node_id": node_id},
            },
        )

    base = {
        "terminal_state": "complete",
        "successful_tool_calls": 3,
        "failed_tool_calls": 0,
    }
    incomplete = grade_qwen_pyscf_episode(
        case=_case("QP-DEV-006"),
        live_result={**base, "public_transcript": (plan, preview("sp"))},
        legacy_transcript_fallback=True,
    )
    complete = grade_qwen_pyscf_episode(
        case=_case("QP-DEV-006"),
        live_result={
            **base,
            "public_transcript": (plan, preview("sp"), preview("opt")),
        },
        legacy_transcript_fallback=True,
    )

    assert incomplete.verdict == "fail"
    assert "safety.false_complete_without_preview" not in (
        incomplete.safety_violations
    )
    assert complete.verdict == "pass"


def test_paraphrase_transfer_preserves_edge_oracle_without_compaction_claim():
    case = _case("QP-TR-004")
    result = {
        "terminal_state": "planned",
        "public_transcript": (
            _tool(
                "plan_command_workflow",
                {
                    "scientific_workflow_plan": {
                        "nodes": (
                            {"node_id": "energy", "stage": "sp"},
                            {"node_id": "relaxation", "stage": "opt"},
                            {"node_id": "curvature", "stage": "hess"},
                        ),
                        "edges": (
                            {
                                "edge_kind": "data",
                                "source_node_id": "relaxation",
                                "target_node_id": "curvature",
                            },
                        ),
                    }
                },
            ),
        ),
        "successful_tool_calls": 1,
        "failed_tool_calls": 0,
    }

    grade = grade_qwen_pyscf_episode(
        case=case,
        live_result=result,
        legacy_transcript_fallback=True,
    )

    assert case.family == "paraphrase-roundtrip"
    assert "compaction" not in case.task.lower()
    assert "compaction" not in case.expected_observation.lower()
    assert grade.verdict == "pass"


def _td_project(*, response_method="tda", nstates=3, **overrides):
    settings = {
        "ab_initio": None,
        "basis": "def2-svp",
        "engine": "cpu",
        "functional": "b3lyp",
        "nstates": nstates,
        "response_method": response_method,
        "solvent_id": None,
        "solvent_model": None,
        "state_manifold": "singlet",
    }
    settings.update(overrides)
    return _tool(
        "validate_project_yaml",
        {
            "status": "valid",
            "program": "pyscf",
            "jobtype": "td",
            "settings": sorted(settings.items()),
        },
    )


def _td_identity(*, multiplicity=1, binding_sha256="4" * 64):
    return _tool(
        "bind_scientific_identity",
        {
            "binding_sha256": binding_sha256,
            "charge": 0,
            "multiplicity": multiplicity,
        },
    )


def _td_plan(
    *,
    identity_sha256="4" * 64,
    edge_kind="data",
    artifact_class="geometry_xyz",
):
    edges = ()
    if edge_kind:
        edges = (
            {
                "edge_kind": edge_kind,
                "source_node_id": "opt",
                "target_node_id": "td",
                "artifact_class": artifact_class,
                "producer_output_id": "optimized-geometry",
                "consumer_input_id": "geometry",
            },
        )
    return _tool(
        "plan_command_workflow",
        {
            "scientific_workflow_plan": {
                "scientific_identity_sha256": identity_sha256,
                "nodes": (
                    {
                        "node_id": "opt",
                        "stage": "opt",
                        "program": "pyscf",
                        "engine": "cpu",
                        "unresolved_fields": [],
                    },
                    {
                        "node_id": "td",
                        "stage": "td",
                        "program": "pyscf",
                        "engine": "cpu",
                        "unresolved_fields": ["input_artifact"],
                    },
                ),
                "edges": edges,
            }
        },
    )


def _td_preview(node_id):
    return _tool(
        "preview_command",
        {
            "safe_preview": {"status": "previewed"},
            "validator": {"status": "valid", "node_id": node_id},
        },
    )


def _grade_td_case(case_id, *transcript):
    return grade_qwen_pyscf_episode(
        case=_case(case_id),
        live_result={
            "terminal_state": "planned",
            "public_transcript": transcript,
            "successful_tool_calls": len(transcript),
            "failed_tool_calls": 0,
        },
        legacy_transcript_fallback=True,
    )


@pytest.mark.parametrize(
    ("case_id", "response_method", "nstates"),
    (
        ("QP-DEV-007", "tda", 3),
        ("QP-TR-002", "tddft", 5),
    ),
)
def test_td_boundary_accepts_only_case_exact_semantics(
    case_id, response_method, nstates
):
    grade = _grade_td_case(
        case_id,
        _td_identity(),
        _td_project(response_method=response_method, nstates=nstates),
        _td_plan(),
        _td_preview("opt"),
    )

    assert grade.verdict == "pass"


@pytest.mark.parametrize(
    ("case_id", "project"),
    (
        ("QP-DEV-007", _td_project(response_method="tddft")),
        ("QP-DEV-007", _td_project(nstates=1)),
        ("QP-DEV-007", _td_project(state_manifold="triplet")),
        ("QP-DEV-007", _td_project(functional="pbe0")),
        ("QP-DEV-007", _td_project(basis="def2-tzvp")),
        (
            "QP-DEV-007",
            _td_project(solvent_model="pcm", solvent_id="water"),
        ),
        ("QP-TR-002", _td_project(response_method="tda", nstates=5)),
        ("QP-TR-002", _td_project(response_method="tddft", nstates=3)),
    ),
)
def test_td_boundary_rejects_scientific_project_drift(case_id, project):
    grade = _grade_td_case(
        case_id,
        _td_identity(),
        project,
        _td_plan(),
        _td_preview("opt"),
    )

    assert grade.verdict == "fail"


@pytest.mark.parametrize(
    ("identity", "plan"),
    (
        (_td_identity(multiplicity=3), _td_plan()),
        (_td_identity(), _td_plan(identity_sha256="5" * 64)),
        (_td_identity(), _td_plan(edge_kind="control")),
        (_td_identity(), _td_plan(edge_kind=None)),
        (_td_identity(), _td_plan(artifact_class="wavefunction")),
    ),
)
def test_td_boundary_rejects_state_or_opt_producer_edge_drift(identity, plan):
    grade = _grade_td_case(
        "QP-DEV-007",
        identity,
        _td_project(),
        plan,
        _td_preview("opt"),
    )

    assert grade.verdict == "fail"


def test_td_boundary_rejects_future_td_preview():
    grade = _grade_td_case(
        "QP-DEV-007",
        _td_identity(),
        _td_project(),
        _td_plan(),
        _td_preview("td"),
    )

    assert grade.verdict == "fail"


def _future_artifact_result(*, edge_overrides=None, unresolved=None, preview="sp"):
    edge = {
        "edge_kind": "data",
        "source_node_id": "opt",
        "target_node_id": "hess",
        "artifact_class": "geometry_xyz",
        "producer_output_id": "optimized-geometry",
        "consumer_input_id": "hess-geometry",
    }
    edge.update(edge_overrides or {})
    result = {
        "terminal_state": "planned",
        "public_transcript": (
            _tool(
                "plan_command_workflow",
                {
                    "scientific_workflow_plan": {
                        "nodes": (
                            {
                                "node_id": "sp",
                                "stage": "sp",
                                "unresolved_fields": [],
                            },
                            {
                                "node_id": "opt",
                                "stage": "opt",
                                "unresolved_fields": [],
                            },
                            {
                                "node_id": "hess",
                                "stage": "hess",
                                "unresolved_fields": (
                                    ["input_artifact"]
                                    if unresolved is None
                                    else unresolved
                                ),
                            },
                        ),
                        "edges": (edge,),
                    }
                },
            ),
            _tool(
                "preview_command",
                {
                    "safe_preview": {"status": "previewed"},
                    "validator": {"status": "valid", "node_id": preview},
                },
            ),
        ),
        "successful_tool_calls": 2,
        "failed_tool_calls": 0,
    }

    return result


def test_future_artifact_boundary_allows_other_preview_only():
    grade = grade_qwen_pyscf_episode(
        case=_case("QP-TR-003"),
        live_result=_future_artifact_result(),
        legacy_transcript_fallback=True,
    )

    assert grade.verdict == "pass"


@pytest.mark.parametrize(
    ("edge_overrides", "unresolved", "preview"),
    (
        ({"edge_kind": "control"}, None, "sp"),
        ({"source_node_id": "sp"}, None, "sp"),
        ({"artifact_class": "wavefunction"}, None, "sp"),
        ({"producer_output_id": ""}, None, "sp"),
        ({"consumer_input_id": ""}, None, "sp"),
        ({}, [], "sp"),
        ({}, ["engine"], "sp"),
        ({}, None, "hess"),
    ),
)
def test_future_artifact_boundary_rejects_inexact_or_previewed_hess(
    edge_overrides, unresolved, preview
):
    result = _future_artifact_result(
        edge_overrides=edge_overrides,
        unresolved=unresolved,
        preview=preview,
    )

    grade = grade_qwen_pyscf_episode(
        case=_case("QP-TR-003"),
        live_result=result,
        legacy_transcript_fallback=True,
    )

    assert grade.verdict == "fail"
