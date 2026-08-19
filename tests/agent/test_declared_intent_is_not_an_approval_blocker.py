"""A stage the scientist declares unrunnable must not block the runnable rest.

A real session planned the paper's third functional, watched ORCA's validator
refuse it -- the functional is not implemented there -- and correctly said the
stages should remain as declared non-executable scientific intent. Two gaps
then trapped it. The release predicate honoured a `blocked_unsupported`
declaration only for job families outside the executable matrix, so five nodes
riding on plain `orca/opt` and `orca/scan` kept blocking the approval of four
green-previewed stages. And `amend_scientific_workflow` had no channel for the
declaration, so evidence learned mid-session -- the red preview -- could not
become declared intent at all.

Both halves are narrowings. A declared stage stays displayed with the workflow,
is excluded from the approval, and is never launched; the only thing the
declaration can do is shrink what runs. The reverse direction would widen the
executable partition and stays a new workflow for human review.
"""

from __future__ import annotations

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.agent.workflows import ScientificWorkflowNodeV2


def _host(tmp_path):
    return CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="s1"
        ),
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "workspace",
    )


def _node(node_id, stage, support_state="resolvable", blocked_reason=""):
    return ScientificWorkflowNodeV2(
        node_id=node_id,
        stage=stage,
        requested_program="orca",
        program="orca",
        engine="cpu",
        project_role=f"{node_id}-project",
        unresolved_fields=(),
        produces_observables=("energy",),
        support_state=support_state,
        blocked_reason=blocked_reason,
    )


class _Plan:
    """Just the two attributes the predicate reads."""

    plan_sha256 = "b" * 64

    def __init__(self, *nodes):
        self.nodes = tuple(nodes)


def test_a_declared_stage_on_an_executable_family_is_non_executable(tmp_path):
    """The exact trap: `opt` is executable, the functional is not."""

    host = _host(tmp_path)
    plan = _Plan(
        _node("opt-b3lyp", "opt"),
        _node(
            "opt-mn15",
            "opt",
            support_state="blocked_unsupported",
            blocked_reason="the program validator refused this functional",
        ),
    )

    non_executable = host._release_non_executable_node_ids(plan)

    assert non_executable == {"opt-mn15"}


def test_a_declared_stage_no_longer_blocks_approval(tmp_path):
    host = _host(tmp_path)
    # Readiness consults the release partition only under a bounded
    # execution envelope, which every live execution session carries.
    host.bounded_execution_envelope = object()
    plan = _Plan(
        _node("opt-b3lyp", "opt"),
        _node(
            "scan-mn15",
            "scan",
            support_state="blocked_unsupported",
            blocked_reason="functional unavailable in this program",
        ),
    )

    readiness = host._approval_readiness(plan)

    states = {
        item["node_id"]: item["approval_state"] for item in readiness["nodes"]
    }
    assert states["scan-mn15"] == "non_executable"
    assert "scan-mn15" not in readiness["blocking_node_ids"]


def test_an_undeclared_unrunnable_family_still_blocks(tmp_path):
    """The declaration is required; nothing is inferred on the planner's
    behalf. Gaussian execution is unclaimed, so a gaussian node without the
    declaration must keep blocking."""

    host = _host(tmp_path)
    gaussian = ScientificWorkflowNodeV2(
        node_id="g16-opt",
        stage="opt",
        requested_program="gaussian",
        program="gaussian",
        engine="cpu",
        project_role="g16-project",
        unresolved_fields=(),
        produces_observables=("energy",),
    )

    non_executable = host._release_non_executable_node_ids(_Plan(gaussian))

    assert non_executable == frozenset()


def test_non_executability_cascades_to_the_consumer(tmp_path):
    """A consumer of a blocked producer can never receive its input.

    Observed live: a cluster-continuum plan declared its complex-optimisation
    stage blocked (the workspace held only monomers) and fed the TD stage
    from it. Readiness called the workflow approvable while the review
    builder refused it, and the session ended honest-looking with no packet.
    The cascade closes that disagreement as a pure narrowing.
    """

    from types import SimpleNamespace

    host = _host(tmp_path)
    host.bounded_execution_envelope = SimpleNamespace(
        allows=lambda program, engine: True
    )
    producer = _node(
        "opt-complex",
        "opt",
        support_state="blocked_unsupported",
        blocked_reason="no complex arrangement exists in this workspace",
    )
    consumer = _node("td-complex", "td")
    plan = _Plan(producer, consumer)
    plan.edges = (
        SimpleNamespace(
            edge_kind="data",
            artifact_class="geometry_xyz",
            source_node_id="opt-complex",
            target_node_id="td-complex",
        ),
    )

    readiness = host._approval_readiness(plan)

    states = {
        item["node_id"]: item["approval_state"] for item in readiness["nodes"]
    }
    assert states["opt-complex"] == "non_executable"
    assert states["td-complex"] == "non_executable"
    assert readiness["blocking_node_ids"] == ()


def test_a_support_repair_on_an_unknown_node_is_refused(tmp_path):
    host = _host(tmp_path)

    with pytest.raises(ContractError, match="unknown scientific workflow"):
        host._rebind_scientific_workflow_projects(
            "t1",
            {
                "workflow_id": "no-such-workflow",
                "replacements": (),
                "analysis_repairs": (),
                "support_repairs": (
                    {"node_id": "x", "blocked_reason": "because"},
                ),
            },
        )


def test_the_amend_tool_advertises_the_channel():
    from chemsmart.agent.tool_specs import build_command_compiled_tool_surface

    amend = next(
        item
        for item in build_command_compiled_tool_surface().tool_definitions
        if item["function"]["name"] == "amend_scientific_workflow"
    )
    schema = amend["function"]["parameters"]["properties"]["support_repairs"]

    assert schema["items"]["required"] == ["node_id", "blocked_reason"]
    assert "narrow" in schema["description"]


def test_an_empty_amendment_is_still_refused(tmp_path):
    host = _host(tmp_path)

    with pytest.raises(ContractError, match="changes something"):
        host._amend_scientific_workflow(
            "t1", {"workflow_id": "w", "support_repairs": ()}
        )
