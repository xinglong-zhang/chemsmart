"""A session required to reproduce an approved plan must be shown it.

an execution-enabled session refuses to execute unless its own plan hashes to the
approved `plan_sha256`:

    planned workflow differs from frozen execution approval

That digest covers nine top-level plan fields, ten fields per node and seven
per edge. The public approval projection disclosed one, four and two of them.
A session therefore had to re-derive `project_role` (`"opt-hf-631gd"`),
`required_observables`, `complexity_factors`
(`["multiple_stages", "producer_artifact_edge"]`) and `edge_id`
(`"data.opt.hess.optimized_geometry"`) -- free-form strings chosen by an
earlier, different session -- and was told only that its plan "differs".

Measured on a real run: the approval loaded, the composition threaded the
frozen body, `execute_approved_program_node` appeared on the surface for the
first time, and the session was then refused on a digest it had no way to hit.
Reproducing an approved plan was not merely hard, it was unspecified.

Disclosing the approved plan grants nothing. It is the plan the user signed,
execution still requires the digest to match, and the body is disclosed only
when it hashes to the digest the frozen approval already pins -- so a reviewer
cannot show the session one plan while approving another.
"""

import inspect
from types import SimpleNamespace

from chemsmart.agent import live_session

#: Every field `plan_sha256` is computed over, by level.
_PLAN_FIELDS = {
    "schema_version",
    "workflow_id",
    "task_spec_sha256",
    "scientific_identity_sha256",
    "nodes",
    "edges",
    "complexity_factors",
    "status",
    "required_observables",
}
_NODE_FIELDS = {
    "node_id",
    "stage",
    "requested_program",
    "program",
    "engine",
    "project_role",
    "unresolved_fields",
    "produces_observables",
    "support_state",
    "blocked_reason",
    "charge",
    "multiplicity",
}
_EDGE_FIELDS = {
    "edge_id",
    "source_node_id",
    "target_node_id",
    "edge_kind",
    "artifact_class",
    "producer_output_id",
    "consumer_input_id",
}


def test_the_digest_covers_exactly_the_fields_this_gate_names():
    """If the plan contract grows, this gate must be revisited deliberately."""
    import dataclasses as dc

    from chemsmart.agent.workflows import (
        ScientificWorkflowEdgeV2,
        ScientificWorkflowNodeV2,
        ScientificWorkflowPlanV2,
    )

    plan_fields = {f.name for f in dc.fields(ScientificWorkflowPlanV2)}
    assert plan_fields == _PLAN_FIELDS | {"plan_sha256"}
    assert {
        f.name for f in dc.fields(ScientificWorkflowNodeV2)
    } == _NODE_FIELDS
    assert {
        f.name for f in dc.fields(ScientificWorkflowEdgeV2)
    } == _EDGE_FIELDS


def test_the_projection_discloses_the_plan_it_will_be_judged_against():
    source = inspect.getsource(live_session._public_workflow_approval)
    assert "approved_scientific_plan" in source
    assert "approved_plan_sha256" in source
    assert "plan_reproduction_rule" in source, (
        "the session must be told that reproduction is required, not only "
        "given the material to reproduce"
    )


def test_the_disclosed_plan_is_self_verifying():
    """A reviewer cannot show one plan while approving another."""
    source = inspect.getsource(live_session._execution_composition_inputs)
    assert "plan.plan_sha256 != frozen.plan_sha256" in source
    assert "approved scientific plan is not the plan that was frozen" in source


def test_disclosure_stays_optional():
    """An approval with no plan body keeps working, minus the disclosure."""
    source = inspect.getsource(live_session._execution_composition_inputs)
    assert "if raw_plan is not None:" in source

    projection = inspect.signature(live_session._public_workflow_approval)
    assert projection.parameters["approved_plan"].default is None


def test_the_execution_gate_it_serves_still_exists():
    """If this refusal is ever removed, the disclosure loses its reason."""
    from chemsmart.agent import tool_runtime

    source = inspect.getsource(tool_runtime)
    assert "planned workflow differs from frozen execution approval" in source


def test_the_approved_plan_is_registered_for_execution_resolution():
    """The consequence gate, observed twelve times in one run.

    ``_execution_scientific_plan`` resolves the DAG to execute by looking up
    the frozen approval's ``plan_sha256`` in the host's plan registry, which
    holds only plans this session produced. An approved plan comes from an
    earlier session, so the lookup raised "frozen workflow approval has no
    registered scientific plan" and execution could not start even after the
    session was shown what to reproduce.
    """
    from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

    assert (
        "scientific_workflow_plan"
        in inspect.signature(CommandCompiledToolHostV1).parameters
    )
    composition = inspect.getsource(live_session.run_live_agent_session)
    assert 'host_kwargs["scientific_workflow_plan"] = approved_plan' in (
        composition
    ), (
        "the approved plan must be registered on the host, or execution "
        "cannot resolve the DAG the approval authorises"
    )

    lookup = inspect.getsource(
        CommandCompiledToolHostV1._execution_scientific_plan
    )
    assert "frozen_approval.plan_sha256" in lookup, (
        "registration is keyed by plan_sha256; if the lookup key changes, "
        "this wiring must change with it"
    )


def test_public_plan_omits_absent_node_state_but_keeps_explicit_state():
    """The disclosed body must itself be accepted by the model tool schema."""

    from chemsmart.agent.workflows import (
        ScientificWorkflowNodeV2,
        build_scientific_workflow_plan,
    )

    approval = SimpleNamespace(
        workflow_id="state-disclosure",
        node_bindings=(
            SimpleNamespace(
                node_id="vertical-energy",
                program="pyscf",
                engine="cpu",
                jobtype="sp",
                input_mode="initial",
            ),
        ),
        producer_edges=(),
        status="approved",
    )

    def disclosed_node(*, charge=None, multiplicity=None):
        plan = build_scientific_workflow_plan(
            workflow_id="state-disclosure",
            task_spec_sha256="a" * 64,
            scientific_identity_sha256="b" * 64,
            nodes=(
                ScientificWorkflowNodeV2(
                    node_id="vertical-energy",
                    stage="sp",
                    requested_program="pyscf",
                    program="pyscf",
                    engine="cpu",
                    project_role="project.vertical",
                    unresolved_fields=(),
                    charge=charge,
                    multiplicity=multiplicity,
                ),
            ),
        )
        return live_session._public_workflow_approval(
            approval, approved_plan=plan
        )["approved_scientific_plan"]["nodes"][0]

    inherited = disclosed_node()
    assert "charge" not in inherited
    assert "multiplicity" not in inherited
    explicit = disclosed_node(charge=-1, multiplicity=2)
    assert (explicit["charge"], explicit["multiplicity"]) == (-1, 2)
