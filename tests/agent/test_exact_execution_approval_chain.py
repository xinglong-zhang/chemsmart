"""Product approval is exact, human-resolved, one-shot, and provider-free."""

from __future__ import annotations

from dataclasses import replace
import json

import pytest

from chemsmart.agent._contracts import ContractError, canonical_data
from chemsmart.agent.api_access import DEFAULT_KEY_LABELS
from chemsmart.agent.execution import (
    ApprovedNodeBindingV1,
    WorkflowEnvironmentBindingV1,
    approve_workflow_execution_review,
    build_execution_resource_spec,
    build_workflow_approval_request,
    build_workflow_execution_review,
    build_workflow_execution_node_review,
    build_workflow_review_resolution,
    execution_server_profile_sha256,
    execution_path_placeholder,
    project_real_execution_argv,
    workflow_execution_approval_bundle_json,
    workflow_execution_review_json,
)
from chemsmart.agent.executor import (
    ApprovedWorkflowExecutor,
    _provider_secret_environment_labels,
)
from chemsmart.agent.execution_envelope import BoundedExecutionEnvelopeV1
from chemsmart.agent.live_session import (
    claim_workflow_execution_approval_bundle,
    load_workflow_execution_approval_bundle,
    load_workflow_execution_review,
    resolve_workflow_execution_review,
)
from chemsmart.agent.workflows import (
    MaterializedNodeV1,
    ScientificWorkflowNodeV2,
    build_materialized_workflow,
    build_scientific_workflow_plan,
)


def _review(tmp_path):
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    resources = build_execution_resource_spec(
        execution_target="run",
        cores=2,
        memory_gb=4.0,
        gpu_count=0,
        scratch_policy="server",
        node_timeout_seconds=300,
    )
    envelope = BoundedExecutionEnvelopeV1(
        schema_version="chemsmart.bounded-execution-envelope.v1",
        mode="bounded-local",
        allowed_program_engines=(("pyscf", ("cpu",)),),
        resources=resources,
        episode_wall_time_seconds=900,
        postprocess_reserve_seconds=60,
        max_engine_calls=1,
        scratch_root=str(tmp_path / "scratch"),
    )
    plan = build_scientific_workflow_plan(
        workflow_id="water-sp",
        task_spec_sha256="a" * 64,
        scientific_identity_sha256="b" * 64,
        nodes=(
            ScientificWorkflowNodeV2(
                node_id="sp",
                stage="sp",
                requested_program="pyscf",
                program="pyscf",
                engine="cpu",
                project_role="water-sp",
                unresolved_fields=(),
            ),
        ),
        edges=(),
    )
    materialized = build_materialized_workflow(
        plan=plan,
        live_cli_schema_sha256="c" * 64,
        resource_sha256=resources.resource_sha256,
        nodes=(
            MaterializedNodeV1(
                node_id="sp",
                input_artifact_sha256="d" * 64,
                project_artifact_sha256="e" * 64,
                project_validation_receipt_sha256="f" * 64,
                environment_receipt_sha256="1" * 64,
                invocation_sha256="2" * 64,
                invocation_identity_sha256="3" * 64,
                preflight_receipt_sha256="4" * 64,
                state="previewed",
            ),
        ),
        unresolved_node_ids=(),
        status="previewed",
    )
    binding = ApprovedNodeBindingV1(
        node_id="sp",
        program="pyscf",
        engine="cpu",
        jobtype="sp",
        project_artifact_sha256="e" * 64,
        settings_sha256="5" * 64,
        charge=0,
        multiplicity=1,
        input_mode="initial",
        initial_artifact_id="water",
        initial_artifact_sha256="d" * 64,
        scientific_identity_sha256="b" * 64,
        producer_edge_sha256="",
    )
    request = build_workflow_approval_request(
        request_id="review-water-sp",
        workflow_id=plan.workflow_id,
        workflow_sha256=plan.plan_sha256,
        task_spec_sha256=plan.task_spec_sha256,
        workspace=workspace,
        resources=resources,
        node_bindings=(binding,),
    )
    server_profile_sha256 = execution_server_profile_sha256(
        resources=resources,
        scratch_root=envelope.scratch_root,
    )
    node_review = build_workflow_execution_node_review(
        node_id="sp",
        stage="sp",
        program="pyscf",
        engine="cpu",
        molecular_identity={
            "identity_evidence_status": "task-bound-geometry-only",
            "coordinate_identity": {
                "kind": "exact-input-artifact",
                "geometry_artifact_sha256": "d" * 64,
            },
            "input_binding_sha256": "d" * 64,
            "charge": 0,
            "multiplicity": 1,
            "electronic_state": "charge-and-multiplicity-specified",
            "scientific_identity_sha256": "b" * 64,
        },
        project_artifact_sha256="e" * 64,
        project_settings_sha256="5" * 64,
        project_settings={"basis": "def2-svp", "method": "b3lyp"},
        real_execution_argv=(
            "chemsmart",
            "run",
            "--no-fake",
            "--delete-scratch",
            "--scratch",
            "--num-cores",
            "2",
            "--num-gpus",
            "0",
            "--mem-gb",
            "4",
            "--server",
            execution_path_placeholder(
                "server-profile", server_profile_sha256
            ),
            "pyscf",
            "--project",
            execution_path_placeholder("project-yaml", "e" * 64),
            "--filename",
            execution_path_placeholder("molecular-input", "d" * 64),
            "sp",
        ),
        execution_resources=resources,
        environment_summary={
            "program": "pyscf",
            "engine": "cpu",
            "status": "available",
            "target_kind": "python_module",
            "locator": "pyscf",
            "observed_version": "2.14.0",
        },
        server_profile_sha256=server_profile_sha256,
        environment_receipt_sha256="1" * 64,
        environment_identity_sha256="6" * 64,
    )
    return build_workflow_execution_review(
        request=request,
        scientific_plan=plan,
        materialized_workflow=materialized,
        execution_resources=resources,
        execution_envelope=canonical_data(envelope),
        environment_bindings=(
            WorkflowEnvironmentBindingV1(
                node_id="sp",
                program="pyscf",
                engine="cpu",
                environment_receipt_sha256="1" * 64,
                environment_identity_sha256="6" * 64,
            ),
        ),
        node_reviews=(node_review,),
    )


def _mixed_review(tmp_path):
    base = _review(tmp_path)
    plan = build_scientific_workflow_plan(
        workflow_id=base.scientific_plan.workflow_id,
        task_spec_sha256=base.scientific_plan.task_spec_sha256,
        scientific_identity_sha256=(
            base.scientific_plan.scientific_identity_sha256
        ),
        nodes=(
            *base.scientific_plan.nodes,
            ScientificWorkflowNodeV2(
                node_id="irc-non-executable",
                stage="irc",
                requested_program="orca",
                program="orca",
                engine="cpu",
                project_role="reaction-path",
                unresolved_fields=(),
                support_state="blocked_unsupported",
                blocked_reason="IRC execution is not release-qualified",
            ),
        ),
        edges=(),
    )
    materialized = build_materialized_workflow(
        plan=plan,
        live_cli_schema_sha256=(
            base.materialized_workflow.live_cli_schema_sha256
        ),
        resource_sha256=base.materialized_workflow.resource_sha256,
        nodes=base.materialized_workflow.nodes,
        unresolved_node_ids=("irc-non-executable",),
        status="partial",
    )
    request = build_workflow_approval_request(
        request_id="review-water-sp-mixed",
        workflow_id=plan.workflow_id,
        workflow_sha256=plan.plan_sha256,
        task_spec_sha256=plan.task_spec_sha256,
        workspace=base.request.workspace,
        resources=base.execution_resources,
        node_bindings=base.request.node_bindings,
    )
    return build_workflow_execution_review(
        request=request,
        scientific_plan=plan,
        materialized_workflow=materialized,
        execution_resources=base.execution_resources,
        execution_envelope=base.execution_envelope,
        environment_bindings=base.environment_bindings,
        node_reviews=base.node_reviews,
        non_executable_node_ids=("irc-non-executable",),
    )


def test_review_is_inert_and_exactly_loadable(tmp_path):
    review = _review(tmp_path)
    path = tmp_path / "review.json"
    path.write_text(workflow_execution_review_json(review), encoding="utf-8")

    loaded = load_workflow_execution_review(path)

    assert loaded == review
    assert loaded.status == "unapproved"
    node = loaded.node_reviews[0]
    assert node.molecular_identity["charge"] == 0
    assert '"basis":"def2-svp"' in node.project_settings_text
    assert "--no-fake" in node.real_execution_argv
    assert str(tmp_path) not in " ".join(node.real_execution_argv)


def test_v1_packets_without_non_executable_field_remain_loadable(tmp_path):
    review = _review(tmp_path)
    legacy_review = canonical_data(review)
    legacy_review.pop("non_executable_node_ids")
    review_path = tmp_path / "legacy-review.json"
    review_path.write_text(
        json.dumps({"workflow_execution_review": legacy_review}),
        encoding="utf-8",
    )

    observed_review = load_workflow_execution_review(review_path)

    assert observed_review == review
    bundle = approve_workflow_execution_review(
        review,
        approval_id="approval-water-sp",
        approved_review_sha256=review.review_sha256,
        actor="human",
        resolution_id="resolution-water-sp",
    )
    legacy_bundle = canonical_data(bundle)
    legacy_bundle.pop("non_executable_node_ids")
    bundle_path = tmp_path / "legacy-approval.json"
    bundle_path.write_text(
        json.dumps({"workflow_execution_approval_bundle": legacy_bundle}),
        encoding="utf-8",
    )

    assert load_workflow_execution_approval_bundle(bundle_path) == bundle


def test_non_executable_partition_round_trips_through_review_and_bundle(tmp_path):
    review = _mixed_review(tmp_path)
    review_path = tmp_path / "mixed-review.json"
    review_path.write_text(
        workflow_execution_review_json(review), encoding="utf-8"
    )

    loaded_review = load_workflow_execution_review(review_path)
    bundle = approve_workflow_execution_review(
        loaded_review,
        approval_id="approval-water-sp-mixed",
        approved_review_sha256=loaded_review.review_sha256,
        actor="human",
        resolution_id="resolution-water-sp-mixed",
    )
    bundle_path = tmp_path / "mixed-approval.json"
    bundle_path.write_text(
        workflow_execution_approval_bundle_json(bundle), encoding="utf-8"
    )

    loaded_bundle = load_workflow_execution_approval_bundle(bundle_path)

    assert loaded_review.non_executable_node_ids == (
        "irc-non-executable",
    )
    assert loaded_bundle.non_executable_node_ids == (
        "irc-non-executable",
    )
    assert loaded_bundle.frozen_workflow_approval.approved_node_ids == (
        "sp",
    )


def test_changed_review_digest_cannot_be_approved(tmp_path):
    review = _review(tmp_path)

    with pytest.raises(ContractError, match="reviewed digest"):
        approve_workflow_execution_review(
            review,
            approval_id="approval-water-sp",
            approved_review_sha256="9" * 64,
            actor="human",
            resolution_id="resolution-water-sp",
        )


@pytest.mark.parametrize("decision", ("deny", "revise", "quit"))
def test_nonapproval_resolutions_never_grant_authority(tmp_path, decision):
    review = _review(tmp_path)

    resolution = build_workflow_review_resolution(
        resolution_id=f"resolution-{decision}",
        review_sha256=review.review_sha256,
        decision=decision,
        actor="human",
    )

    assert resolution.approval_id == ""
    assert resolution.one_shot is False


def test_approved_bundle_round_trips_and_is_consumed_once(tmp_path):
    review = _review(tmp_path)
    review_path = tmp_path / "review.json"
    review_path.write_text(workflow_execution_review_json(review), encoding="utf-8")
    bundle_path = tmp_path / "approval.json"

    resolution, bundle = resolve_workflow_execution_review(
        review_file=review_path,
        reviewed_sha256=review.review_sha256,
        decision="approve",
        actor="human",
        output_file=bundle_path,
    )

    assert resolution.decision == "approve"
    assert bundle is not None and bundle.one_shot is True
    assert load_workflow_execution_approval_bundle(bundle_path) == bundle
    claim_workflow_execution_approval_bundle(
        bundle, workspace=review.request.workspace
    )
    with pytest.raises(ContractError, match="already been consumed"):
        claim_workflow_execution_approval_bundle(
            bundle, workspace=review.request.workspace
        )


def test_tampered_bundle_fails_closed_before_execution(tmp_path):
    review = _review(tmp_path)
    bundle = approve_workflow_execution_review(
        review,
        approval_id="approval-water-sp",
        approved_review_sha256=review.review_sha256,
        actor="human",
        resolution_id="resolution-water-sp",
    )
    payload = json.loads(
        json.dumps(
            {"workflow_execution_approval_bundle": canonical_data(bundle)}
        )
    )
    payload["workflow_execution_approval_bundle"]["execution_resources"][
        "cores"
    ] = 8
    path = tmp_path / "tampered.json"
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ContractError):
        load_workflow_execution_approval_bundle(path)


def test_path_roles_are_digest_bound_and_cannot_be_swapped(tmp_path):
    project = str(tmp_path / "project.yaml")
    geometry = str(tmp_path / "water.xyz")
    argv = ("chemsmart", "run", "pyscf", "--project", project,
            "--filename", geometry, "sp")

    projected = project_real_execution_argv(
        argv,
        path_bindings={
            project: ("project-yaml", "e" * 64),
            geometry: ("molecular-input", "d" * 64),
        },
    )
    swapped = project_real_execution_argv(
        argv,
        path_bindings={
            project: ("molecular-input", "d" * 64),
            geometry: ("project-yaml", "e" * 64),
        },
    )

    assert projected != swapped
    assert str(tmp_path) not in " ".join(projected)
    with pytest.raises(ContractError, match="unreviewed absolute path"):
        project_real_execution_argv(
            (*argv[:-2], str(tmp_path / "changed.xyz"), argv[-1]),
            path_bindings={
                project: ("project-yaml", "e" * 64),
                geometry: ("molecular-input", "d" * 64),
            },
        )


def test_post_review_command_or_identity_mutation_is_rejected(tmp_path):
    node = _review(tmp_path).node_reviews[0]

    with pytest.raises(ContractError, match="argv review digest"):
        replace(node, real_execution_argv=(*node.real_execution_argv, "--changed"))
    with pytest.raises(ContractError, match="molecular identity review digest"):
        replace(
            node,
            molecular_identity={**node.molecular_identity, "charge": 1},
        )


def test_launch_is_compared_before_claim_and_again_before_reuse(
    tmp_path, monkeypatch
):
    calls = []

    class _Bundle:
        def node_review(self, node_id):
            calls.append(("review", node_id))
            return object()

    class _Host:
        reject = True

        def verify_reviewed_real_execution_argv(self, **values):
            calls.append(("compare", values["node_id"]))
            if self.reject:
                raise ContractError("command mismatch")

    def _claim(_bundle, *, workspace):
        calls.append(("claim", Path(workspace).name))

    from pathlib import Path
    from chemsmart.agent import live_session

    monkeypatch.setattr(
        live_session,
        "claim_workflow_execution_approval_bundle",
        _claim,
    )
    executor = object.__new__(ApprovedWorkflowExecutor)
    executor.execution_bundle = _Bundle()
    executor.host = _Host()
    executor.approval_workspace = tmp_path
    executor.claim_workspace_bundle = True
    executor._bundle_claimed = False

    with pytest.raises(ContractError, match="command mismatch"):
        executor._verify_launch_and_claim_once(
            node_id="sp", invocation_sha256="a" * 64
        )
    assert not any(kind == "claim" for kind, _value in calls)

    executor.host.reject = False
    executor._verify_launch_and_claim_once(
        node_id="sp", invocation_sha256="a" * 64
    )
    executor._verify_launch_and_claim_once(
        node_id="sp", invocation_sha256="a" * 64
    )
    assert [kind for kind, _value in calls].count("compare") == 3
    assert [kind for kind, _value in calls].count("claim") == 1


def test_provider_credentials_are_removed_from_engine_environment():
    expected = {
        label for labels in DEFAULT_KEY_LABELS.values() for label in labels
    }

    assert set(_provider_secret_environment_labels()) == expected
