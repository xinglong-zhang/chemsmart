"""A deferred node's deferred producer resolves; the frontier agrees.

REACH-1 ino3 (2026-09-06) planned neutral opt -> cation opt -> TZVP
single point, twice. The frontier deferred every consumer by its edge
shape and said approvable; the review then refused the two second-hop
nodes with "node has no compiled command invocation", because a
deferred node's context was taken from its producer's compiled
invocation and a deferred producer has none. The refusal could not be
recorded and the goal returned with its whole grant unspent. Replayed
from the transcript on the sealed commit, not inferred. Two repairs:
the review resolves a deferred producer as it would the node itself,
and the frontier asks the review's resolver before calling a node
deferred, so the two organs cannot disagree.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.tool_runtime import (
    CommandCompiledToolHostV1,
    ScientificWorkflowEdgeV2,
    ScientificWorkflowNodeV2,
    _CommandContext,
    build_scientific_workflow_plan,
)


def _capability(digest, jobtype):
    return SimpleNamespace(
        receipt_sha256=digest,
        status=SimpleNamespace(value="supported"),
        query=SimpleNamespace(program="orca", jobtype=jobtype, engine="cpu"),
    )


def _engine(capability, binding, environment):
    return SimpleNamespace(
        program="orca",
        selected_engine="cpu",
        state="resolved",
        execution_ready=True,
        capability_receipt_sha256=capability.receipt_sha256,
        program_binding_sha256=binding,
        environment_receipt_sha256=environment,
    )


def _node(node_id, stage, project_role, future=False):
    return ScientificWorkflowNodeV2(
        node_id=node_id,
        stage=stage,
        requested_program="orca",
        program="orca",
        engine="cpu",
        project_role=project_role,
        unresolved_fields=(),
        **({"support_state": "unresolved_future"} if future else {}),
    )


def _edge(source, target, order):
    return ScientificWorkflowEdgeV2(
        edge_id=f"e{order}-{source}-to-{target}",
        source_node_id=source,
        target_node_id=target,
        edge_kind="data",
        artifact_class="geometry_xyz",
        producer_output_id="optimized-geometry",
        consumer_input_id="filename",
    )


def _draft_node(node_id, project_role, producer):
    return SimpleNamespace(
        node_id=node_id,
        project_role=project_role,
        inputs=(
            SimpleNamespace(
                binding_id="filename",
                artifact_class="geometry_xyz",
                producer_node_id=producer,
                producer_output_id="optimized-geometry",
            ),
        ),
    )


def _two_hop_host():
    identity = SimpleNamespace(
        binding_sha256="a" * 64, charge=0, multiplicity=1
    )
    geometry = SimpleNamespace(artifact_id="geometry.root", sha256="b" * 64)
    root_context = _CommandContext(
        proposal=SimpleNamespace(program="orca", jobtype="opt"),
        capability=SimpleNamespace(),
        program_binding=SimpleNamespace(),
        engine_binding=SimpleNamespace(
            engine="cpu", environment_receipt_sha256="c" * 64
        ),
        project_artifact=SimpleNamespace(
            artifact_id="project.opt", sha256="d" * 64
        ),
        project_validation=SimpleNamespace(
            status="valid", settings_sha256="e" * 64
        ),
        input_artifact=geometry,
        scientific_identity=identity,
    )
    plan = build_scientific_workflow_plan(
        workflow_id="neutral-cation-sp",
        task_spec_sha256="f" * 64,
        scientific_identity_sha256=identity.binding_sha256,
        nodes=(
            _node("n0-opt", "opt", "project.opt"),
            _node("cat-opt", "opt", "project.opt", future=True),
            _node("cat-sp", "sp", "project.sp", future=True),
        ),
        edges=(_edge("n0-opt", "cat-opt", 1), _edge("cat-opt", "cat-sp", 2)),
    )
    opt_capability = _capability("2" * 64, "opt")
    sp_capability = _capability("3" * 64, "sp")
    opt_engine = _engine(opt_capability, "5" * 64, "6" * 64)
    sp_engine = _engine(sp_capability, "7" * 64, "8" * 64)
    projects = {
        "project.opt": SimpleNamespace(
            artifact_id="project.opt", sha256="d" * 64
        ),
        "project.sp": SimpleNamespace(
            artifact_id="project.sp", sha256="1" * 64
        ),
    }
    host = object.__new__(CommandCompiledToolHostV1)
    host.workflow_drafts = {
        "draft": SimpleNamespace(
            workflow_id=plan.workflow_id,
            nodes=(
                _draft_node("cat-opt", "project.opt", "n0-opt"),
                _draft_node("cat-sp", "project.sp", "cat-opt"),
            ),
        )
    }
    host.artifacts = dict(projects)
    host.project_validations = {
        "opt": SimpleNamespace(
            project_artifact_id="project.opt",
            project_sha256="d" * 64,
            capability_receipt_sha256=opt_capability.receipt_sha256,
            program="orca",
            jobtype="opt",
            status="valid",
        ),
        "sp": SimpleNamespace(
            project_artifact_id="project.sp",
            project_sha256="1" * 64,
            capability_receipt_sha256=sp_capability.receipt_sha256,
            program="orca",
            jobtype="sp",
            status="valid",
        ),
    }
    host.engine_bindings = {"opt": opt_engine, "sp": sp_engine}
    host.capabilities = {
        c.receipt_sha256: c for c in (opt_capability, sp_capability)
    }
    host.program_bindings = {
        b.program_binding_sha256: SimpleNamespace()
        for b in (opt_engine, sp_engine)
    }
    host._invocation_workflow_plan_sha256s = {}

    def latest(node_id, **_kwargs):
        if node_id == "n0-opt":
            return SimpleNamespace(invocation_sha256="9" * 64), root_context
        raise ContractError("node has no compiled command invocation")

    host._latest_invocation_for_node = latest
    return host, plan, sp_capability, sp_engine, identity


@pytest.mark.capability("gate:plan.deferred_producer_resolves")
def test_a_second_hop_node_resolves_through_its_deferred_producer():
    host, plan, sp_capability, sp_engine, identity = _two_hop_host()
    context = host._bounded_node_context(
        plan=plan,
        planned_node=plan.nodes[2],
        data_target_ids={"cat-opt", "cat-sp"},
    )
    assert context.capability is sp_capability
    assert context.engine_binding is sp_engine
    assert context.scientific_identity is identity
    assert context.proposal.input_artifact_id == "geometry.root"


@pytest.mark.capability("gate:plan.deferred_producer_resolves")
def test_the_frontier_blocks_what_the_review_would_refuse():
    host, plan, *_ = _two_hop_host()
    host.bounded_execution_envelope = SimpleNamespace()
    host._release_non_executable_node_ids = lambda _plan: frozenset()
    host._bounded_deferred_target_ids = lambda _plan: {"cat-opt", "cat-sp"}
    host._node_is_previewed = lambda node_id, **_k: node_id == "n0-opt"
    ready = host._approval_readiness(plan)
    assert ready["approvable"] is True
    assert ready["deferred_node_ids"] == ("cat-opt", "cat-sp")
    assert host._invocation_workflow_plan_sha256s == {}

    del host.project_validations["sp"]
    refused = host._approval_readiness(plan)
    assert refused["approvable"] is False
    assert refused["blocking_node_ids"] == ("cat-sp",)
    (row,) = [n for n in refused["nodes"] if n["node_id"] == "cat-sp"]
    assert row["approval_state"] == "preview_required"
    assert (
        "lacks unique project/environment evidence" in row["deferral_refusal"]
    )


@pytest.mark.capability("gate:plan.deferred_producer_resolves")
def test_a_node_with_a_hessian_edge_beside_its_geometry_is_a_candidate():
    """po3's IRC nodes carried a geometry edge and a Hessian edge from one
    transition-state producer; the frontier counted two incoming edges
    and called them blocking while the review resolved and ran them."""

    plan = build_scientific_workflow_plan(
        workflow_id="ts-irc",
        task_spec_sha256="f" * 64,
        scientific_identity_sha256="a" * 64,
        nodes=(
            _node("ts", "ts", "project.ts"),
            _node("irc", "irc", "project.irc", future=True),
        ),
        edges=(
            _edge("ts", "irc", 1),
            ScientificWorkflowEdgeV2(
                edge_id="e2-ts-hessian-to-irc",
                source_node_id="ts",
                target_node_id="irc",
                edge_kind="data",
                artifact_class="orca_hessian",
                producer_output_id="hessian",
                consumer_input_id="hess_filename",
            ),
        ),
    )
    host = object.__new__(CommandCompiledToolHostV1)
    host.bounded_execution_envelope = SimpleNamespace(
        allows=lambda _program, _engine: True
    )
    assert host._bounded_deferred_target_ids(plan) == {"irc"}
