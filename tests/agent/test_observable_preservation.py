"""A repair may not clear its findings by deleting the requested stage."""

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.workflows import (
    ScientificWorkflowNodeV2,
    build_scientific_workflow_plan,
)

_TASK = "a" * 64
_IDENTITY = "b" * 64


def _node(node_id, **overrides):
    values = dict(
        node_id=node_id,
        stage="sp",
        requested_program="orca",
        program="orca",
        engine="cpu",
        project_role="sp",
        unresolved_fields=(),
    )
    values.update(overrides)
    return ScientificWorkflowNodeV2(**values)


def _plan(nodes, required_observables=()):
    return build_scientific_workflow_plan(
        workflow_id="wf.probe",
        task_spec_sha256=_TASK,
        scientific_identity_sha256=_IDENTITY,
        nodes=nodes,
        required_observables=required_observables,
    )


def test_plan_records_the_observable_its_nodes_produce():
    plan = _plan(
        (_node("energy", produces_observables=("relative-energy",)),),
        required_observables=("relative-energy",),
    )
    assert plan.required_observables == ("relative-energy",)
    assert plan.nodes[0].produces_observables == ("relative-energy",)


def test_plan_without_a_producer_for_a_required_observable_is_refused():
    with pytest.raises(ContractError, match="no node produces"):
        _plan((_node("energy"),), required_observables=("relative-energy",))


def test_a_blocked_node_still_carries_the_observable():
    """The honest alternative to deletion is an explicit blocked stage."""

    plan = _plan(
        (
            _node(
                "extrapolation",
                produces_observables=("relative-energy",),
                support_state="blocked_unsupported",
                blocked_reason="no complete-basis-set executor is available",
            ),
        ),
        required_observables=("relative-energy",),
    )
    node = plan.nodes[0]
    assert node.support_state == "blocked_unsupported"
    assert node.blocked_reason
    assert node.produces_observables == ("relative-energy",)


def test_declaring_an_observable_changes_the_plan_identity():
    plain = _plan(
        (_node("energy", produces_observables=("relative-energy",)),)
    )
    declared = _plan(
        (_node("energy", produces_observables=("relative-energy",)),),
        required_observables=("relative-energy",),
    )
    assert plain.plan_sha256 != declared.plan_sha256


def _bare_host():
    from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

    host = CommandCompiledToolHostV1.__new__(CommandCompiledToolHostV1)
    host.scientific_plans = {}
    return host


def test_replanning_that_drops_a_stage_is_refused():
    runtime = _bare_host()

    first = _plan(
        (
            _node("energy", produces_observables=("relative-energy",)),
            _node("extrapolation", produces_observables=("cbs-energy",)),
        )
    )
    runtime.scientific_plans[first.workflow_id] = first

    shrunk = _plan(
        (_node("energy", produces_observables=("relative-energy",)),)
    )
    with pytest.raises(ContractError, match="removed workflow stage"):
        runtime._refuse_observable_regression(shrunk)


def test_replanning_that_keeps_every_stage_is_accepted():
    runtime = _bare_host()
    nodes = (
        _node("energy", produces_observables=("relative-energy",)),
        _node("extrapolation", produces_observables=("cbs-energy",)),
    )
    first = _plan(nodes)
    runtime.scientific_plans[first.workflow_id] = first

    repaired = _plan(
        (
            nodes[0],
            _node(
                "extrapolation",
                produces_observables=("cbs-energy",),
                support_state="blocked_unsupported",
                blocked_reason="no complete-basis-set executor is available",
            ),
        )
    )
    runtime._refuse_observable_regression(repaired)
