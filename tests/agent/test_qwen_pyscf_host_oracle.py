from __future__ import annotations

from copy import deepcopy

from chemsmart.agent._contracts import canonical_json, canonical_sha256
from chemsmart.agent.experiments.host_oracle import (
    HostOracleInputBundleV1,
    build_host_oracle_input_bundle,
)
from chemsmart.agent.experiments.qwen_pyscf_fixtures import qwen_pyscf_cases_v1
from chemsmart.agent.experiments.qwen_pyscf_grading import (
    grade_qwen_pyscf_episode,
)
from chemsmart.agent.runtime.events import EventKind, RuntimeEvent


def _case(case_id: str):
    return next(item for item in qwen_pyscf_cases_v1() if item.case_id == case_id)


def _workflow_result():
    return {
        "scientific_workflow_plan": {
            "plan_sha256": "1" * 64,
            "nodes": (
                {"node_id": "sp", "stage": "sp", "unresolved_fields": ()},
                {"node_id": "opt", "stage": "opt", "unresolved_fields": ()},
                {
                    "node_id": "hess",
                    "stage": "hess",
                    "unresolved_fields": ("input_artifact",),
                },
            ),
            "edges": (
                {
                    "edge_kind": "data",
                    "source_node_id": "opt",
                    "target_node_id": "hess",
                },
            ),
        },
        # These host-only/rendered values must not enter an oracle projection.
        "run_directory": "/tmp/private-run",
        "argv": ("chemsmart", "run", "pyscf", "opt"),
        "rendered_yaml": "opt:\n  functional: b3lyp\n",
        "stdout": "private process output",
        "reasoning_content": "provider-private reasoning",
    }


def _preview_result(node_id: str):
    return {
        "safe_preview": {"status": "previewed", "node_id": node_id},
        "validator": {"status": "valid", "node_id": node_id},
    }


def _events(*, mode: str, session_id: str):
    actions = (
        ("plan_command_workflow", _workflow_result()),
        ("preview_command", _preview_result("sp")),
        ("preview_command", _preview_result("opt")),
    )
    rows = []
    previous = ""
    for sequence, (tool, result) in enumerate(actions, start=1):
        canonical_result = {
            "schema_version": "chemsmart.tool-result.v1",
            "tool": tool,
            "status": "ok",
            "result": result,
        }
        event = RuntimeEvent.create(
            sequence=sequence,
            session_id=session_id,
            turn_id=session_id + ".turn-1",
            kind=EventKind.TOOL_SUCCEEDED.value,
            payload={
                "request_id": f"call-{sequence}",
                "tool": tool,
                "result_sha256": canonical_sha256(canonical_result),
                "verdict": "host_observed",
                "typed_receipt_status": "present",
                "canonical_result": canonical_result,
                "feedback_projection": mode,
                "feedback_equivalence_receipt": {
                    "schema_version": "test.feedback-receipt.v1",
                    "mode": mode,
                },
            },
            previous_hash=previous,
            timestamp=f"2026-08-04T00:00:0{sequence}+00:00",
            event_id=f"{mode}-{sequence}",
            idempotency_key=f"tool-{sequence}",
        )
        rows.append(event)
        previous = event.event_hash
    return tuple(rows)


def _bundle(*, mode: str, session_id: str):
    events = _events(mode=mode, session_id=session_id)
    return build_host_oracle_input_bundle(
        events=events,
        session_id=session_id,
        event_stream_head_sha256=events[-1].event_hash,
        successful_tool_calls=3,
        failed_tool_calls=0,
    )


def _live_result(*, mode: str, transcript=()):
    session_id = f"session-{mode}"
    bundle = _bundle(mode=mode, session_id=session_id)
    observations = {
        "schema_version": "chemsmart.live-harness-experiment-observations.v1",
        "experiment_config_sha256": "2" * 64,
        "feedback_projection": mode,
        "usage": {
            "coordinator": {
                "successful_tool_calls": 3,
                "failed_tool_calls": 0,
            }
        },
        "host_oracle_input_bundle": bundle.public_record(),
    }
    observations["record_sha256"] = canonical_sha256(observations)
    return {
        "schema_version": "chemsmart.live-agent-session-result.v1",
        "session_id": session_id,
        "terminal_state": "complete",
        "public_transcript": transcript,
        "artifact_records": (),
        "experiment_observations": observations,
        "successful_tool_calls": 3,
        "failed_tool_calls": 0,
        "event_stream_head_sha256": bundle.event_stream_head_sha256,
    }


def _causal_transcript_without_workflow_nodes():
    return (
        {
            "role": "tool",
            "content": {
                "schema_version": "chemsmart.causal-action-projection.v1",
                "tool": "plan_command_workflow",
                "status": "valid",
                "result": {"status": "planned"},
            },
        },
    )


def test_full_and_causal_feedback_grade_identical_host_actions_identically():
    full = _live_result(mode="full-v1")
    causal = _live_result(
        mode="causal-v1", transcript=_causal_transcript_without_workflow_nodes()
    )

    full_grade = grade_qwen_pyscf_episode(
        case=_case("QP-DEV-006"), live_result=full
    )
    causal_grade = grade_qwen_pyscf_episode(
        case=_case("QP-DEV-006"), live_result=causal
    )

    full_bundle = HostOracleInputBundleV1.from_record(
        full["experiment_observations"]["host_oracle_input_bundle"]
    )
    causal_bundle = HostOracleInputBundleV1.from_record(
        causal["experiment_observations"]["host_oracle_input_bundle"]
    )
    assert full_bundle.tool_actions_sha256 == causal_bundle.tool_actions_sha256
    assert full_grade.verdict == causal_grade.verdict == "pass"
    assert full_grade.scientific_state == causal_grade.scientific_state == "previewed"
    assert tuple((item.check_id, item.status) for item in full_grade.checks) == tuple(
        (item.check_id, item.status) for item in causal_grade.checks
    )


def test_causal_workflow_oracle_retains_nodes_from_host_events():
    result = _live_result(
        mode="causal-v1", transcript=_causal_transcript_without_workflow_nodes()
    )
    bundle = HostOracleInputBundleV1.from_record(
        result["experiment_observations"]["host_oracle_input_bundle"]
    )
    workflow = next(
        item
        for item in bundle.observations
        if item.tool_name == "plan_command_workflow"
    ).oracle_result["scientific_workflow_plan"]

    assert tuple(item["stage"] for item in workflow["nodes"]) == (
        "sp",
        "opt",
        "hess",
    )
    assert grade_qwen_pyscf_episode(
        case=_case("QP-DEV-006"), live_result=result
    ).verdict == "pass"


def test_oracle_projection_is_path_free_but_keeps_scientific_semantics():
    bundle = _bundle(mode="causal-v1", session_id="projection-session")
    workflow_observation = bundle.observations[0]
    serialized = canonical_json(workflow_observation.oracle_result)

    assert "scientific_workflow_plan" in workflow_observation.oracle_result
    assert "nodes" in workflow_observation.oracle_result["scientific_workflow_plan"]
    assert "/tmp/private-run" not in serialized
    assert "rendered_yaml" not in serialized
    assert '"argv"' not in serialized
    assert '"stdout"' not in serialized
    assert "reasoning_content" not in serialized
    assert workflow_observation.omitted_fields


def test_missing_or_tampered_host_bundle_is_inconclusive_without_transcript_fallback():
    missing = _live_result(mode="causal-v1")
    del missing["experiment_observations"]["host_oracle_input_bundle"]
    missing["experiment_observations"].pop("record_sha256")
    missing["experiment_observations"]["record_sha256"] = canonical_sha256(
        missing["experiment_observations"]
    )

    tampered = deepcopy(_live_result(mode="causal-v1"))
    observations = tampered["experiment_observations"][
        "host_oracle_input_bundle"
    ]["observations"]
    observations[0]["oracle_result"]["scientific_workflow_plan"]["nodes"][0][
        "stage"
    ] = "td"
    tampered["experiment_observations"].pop("record_sha256")
    tampered["experiment_observations"]["record_sha256"] = canonical_sha256(
        tampered["experiment_observations"]
    )

    for result in (missing, tampered):
        grade = grade_qwen_pyscf_episode(
            case=_case("QP-DEV-006"), live_result=result
        )
        assert grade.verdict == "inconclusive"
        assert grade.checks[0].check_id == "experiment.host_oracle_input_bundle"
        assert grade.checks[0].status == "inconclusive"


def test_live_and_coordinator_counts_are_bound_to_host_observations():
    result = _live_result(mode="full-v1")
    result["successful_tool_calls"] = 2

    grade = grade_qwen_pyscf_episode(
        case=_case("QP-DEV-006"), live_result=result
    )

    assert grade.verdict == "inconclusive"
    assert grade.successful_tool_calls == 2
    assert grade.failed_tool_calls == 0
