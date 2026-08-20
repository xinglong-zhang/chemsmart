"""The executor supplies the caller, never the science.

Measured across repeated live sessions: seventeen ``agent plan`` sessions
rejected 42 of 868 tool calls (4.8%); the one execution-enabled session
rejected 20 of 69 (29%) -- twelve execution attempts and eight re-plans
*inside* an execution session -- and completed no nodes. That mode was
``run_live_agent_session(execution_enabled=True)``: the same planning loop
with one extra tool, asked to re-derive a plan it was already handed.

These gates pin the division of labour that replaces it. They are deliberately
about what the executor may *not* do.
"""

import inspect
from types import SimpleNamespace

import pytest

from chemsmart.agent import executor as executor_module
from chemsmart.agent._contracts import canonical_sha256
from chemsmart.agent.execution import (
    build_workflow_run_state,
    derive_ready_node_ids,
    argv_shape,
    invocation_identity_sha256,
)
from chemsmart.agent.workflows import (
    ScientificWorkflowEdgeV2,
    ScientificWorkflowNodeV2,
    build_scientific_workflow_plan,
)


def _code_of(module) -> str:
    """Module source minus its docstring, so prose is not mistaken for code."""

    source = inspect.getsource(module)
    docstring = inspect.getdoc(module) or ""
    return source.replace(docstring, "") if docstring else source


def test_the_executor_never_calls_a_provider():
    """The whole claim is that execution needs no model."""

    code = _code_of(executor_module)
    for forbidden in (
        "load_secret_lease(",
        "load_agent_provider_selection(",
        "UnifiedSessionRunner(",
        "run_live_agent_session(",
    ):
        assert forbidden not in code, (
            f"{forbidden} would give the executor a way to reach a model; "
            "the point is that an approved plan needs no further judgement"
        )


def test_readiness_comes_from_the_plan_not_from_the_executor():
    """A hand-rolled order would be the executor deciding the science."""

    source = inspect.getsource(executor_module.ApprovedWorkflowExecutor.run)
    assert "derive_ready_node_ids" in source
    assert "sort" not in source, (
        "ordering nodes by anything other than the plan's own edges would "
        "make the executor a planner"
    )


def test_a_node_is_only_settled_into_the_state_its_receipt_reports():
    source = inspect.getsource(
        executor_module.ApprovedWorkflowExecutor._settle
    )
    assert "outcome.validated" in source
    assert '"validated"' in source and '"engine_complete"' in source
    assert "result_validation_receipt=" in source, (
        "a consumer may only become ready behind a genuinely validated "
        "producer, and that transition requires the typed receipt"
    )


def test_an_approved_input_is_located_by_content_not_by_path():
    source = inspect.getsource(executor_module._locate_by_digest)
    assert "file_sha256" in source


def test_each_sibling_root_uses_its_own_approved_initial_artifact():
    executor = object.__new__(executor_module.ApprovedWorkflowExecutor)
    executor.initial_artifacts = {
        "neutral": SimpleNamespace(artifact_id="neutral", sha256="a" * 64),
        "anion": SimpleNamespace(artifact_id="anion", sha256="b" * 64),
    }

    neutral = SimpleNamespace(
        input_mode="initial",
        initial_artifact_id="neutral",
        initial_artifact_sha256="a" * 64,
    )
    anion = SimpleNamespace(
        input_mode="initial",
        initial_artifact_id="anion",
        initial_artifact_sha256="b" * 64,
    )

    assert executor._input_artifact_id(neutral) == "neutral"
    assert executor._input_artifact_id(anion) == "anion"


def test_preparation_block_blocks_descendants_but_not_independent_siblings():
    plan = build_scientific_workflow_plan(
        workflow_id="preparation-isolation",
        task_spec_sha256="a" * 64,
        scientific_identity_sha256="b" * 64,
        nodes=(
            ScientificWorkflowNodeV2(
                node_id="producer",
                stage="opt",
                requested_program="pyscf",
                program="pyscf",
                engine="cpu",
                project_role="producer-opt",
                unresolved_fields=(),
            ),
            ScientificWorkflowNodeV2(
                node_id="descendant",
                stage="hess",
                requested_program="pyscf",
                program="pyscf",
                engine="cpu",
                project_role="producer-hess",
                unresolved_fields=(),
            ),
            ScientificWorkflowNodeV2(
                node_id="sibling",
                stage="sp",
                requested_program="pyscf",
                program="pyscf",
                engine="cpu",
                project_role="independent-sp",
                unresolved_fields=(),
            ),
        ),
        edges=(
            ScientificWorkflowEdgeV2(
                edge_id="producer-to-descendant",
                source_node_id="producer",
                target_node_id="descendant",
                edge_kind="data",
                artifact_class="geometry_xyz",
                producer_output_id="optimized-geometry",
                consumer_input_id="geometry",
            ),
        ),
    )
    approval = SimpleNamespace(
        workflow_id=plan.workflow_id,
        plan_sha256=plan.plan_sha256,
        approved_node_ids=tuple(sorted(node.node_id for node in plan.nodes)),
        approval_id="preparation-isolation-approval",
        approval_sha256="c" * 64,
    )
    run = build_workflow_run_state(
        run_id="preparation-isolation-run",
        plan=plan,
        approval=approval,
        approval_consumed=True,
    )
    executor = object.__new__(executor_module.ApprovedWorkflowExecutor)
    executor.plan = plan
    outcome = executor_module.ExecutedNodeV1(
        node_id="producer",
        program="pyscf",
        jobtype="opt",
        state="blocked",
        invocation_identity_sha256=canonical_sha256("producer"),
        execution_receipt_sha256="",
        rule_ids=(),
        failure="preparation rejected",
        invocation_sha256="",
    )

    settled = executor._settle(run, outcome)
    states = {node.node_id: node.state for node in settled.nodes}

    assert states == {
        "descendant": "blocked",
        "producer": "blocked",
        "sibling": "pending",
    }
    assert settled.state == "running"
    assert derive_ready_node_ids(plan, settled) == ("sibling",)


class _Node:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def test_invocation_identity_survives_a_change_of_directory():
    """The defect: an approval could not be run by any later process.

    ``invocation_sha256`` covers argv, and argv holds absolute paths into the
    planning session's own timestamped run directory. Driving one approved
    plan from two run directories produced two digests and one identity.
    """

    common = dict(
        program="pyscf",
        engine="cpu",
        jobtype="opt",
        project_sha256="a" * 64,
        input_sha256="b" * 64,
        scientific_identity_sha256="c" * 64,
    )
    planned = invocation_identity_sha256(
        argv=(
            "chemsmart",
            "run",
            "--server",
            "/runs/live-20260806T0739-abc/s.yaml",
            "--fake",
            "pyscf",
            "--project",
            "/runs/live-20260806T0739-abc/p.yaml",
            "--filename",
            "/ws-a/water.xyz",
            "--charge",
            "0",
            "opt",
        ),
        **common,
    )
    executed = invocation_identity_sha256(
        argv=(
            "chemsmart",
            "run",
            "--server",
            "/elsewhere/execution-server.yaml",
            "--fake",
            "pyscf",
            "--project",
            "/other/p.yaml",
            "--filename",
            "/ws-b/water.xyz",
            "--charge",
            "0",
            "opt",
        ),
        **common,
    )
    assert planned == executed


@pytest.mark.parametrize(
    "argv,reason",
    [
        (
            (
                "chemsmart",
                "run",
                "--filename",
                "/ws/x.xyz",
                "--charge",
                "1",
                "opt",
            ),
            "a different charge is different chemistry",
        ),
        (
            (
                "chemsmart",
                "run",
                "--filename",
                "/ws/x.xyz",
                "--charge",
                "0",
                "sp",
            ),
            "a different job type is a different calculation",
        ),
        (
            (
                "chemsmart",
                "run",
                "--filename",
                "/ws/x.xyz",
                "--charge",
                "0",
                "--no-gpu",
                "opt",
            ),
            "an added flag changes what runs",
        ),
    ],
)
def test_invocation_identity_still_fails_on_a_changed_command(argv, reason):
    common = dict(
        program="pyscf",
        engine="cpu",
        jobtype="opt",
        project_sha256="a" * 64,
        input_sha256="b" * 64,
        scientific_identity_sha256="c" * 64,
    )
    baseline = invocation_identity_sha256(
        argv=(
            "chemsmart",
            "run",
            "--filename",
            "/ws/x.xyz",
            "--charge",
            "0",
            "opt",
        ),
        **common,
    )
    assert invocation_identity_sha256(argv=argv, **common) != baseline, reason


def test_argv_shape_erases_only_absolute_paths():
    assert argv_shape(
        ("chemsmart", "run", "/abs/path.yaml", "--charge", "0")
    ) == ("chemsmart", "run", "<path>", "--charge", "0")


def test_the_frozen_preview_still_pins_content():
    """Dropping the unmatched digests must not drop the real constraints."""

    from chemsmart.agent import execution
    from chemsmart.agent.runtime import records

    source = inspect.getsource(records.build_workflow_node_launch_reservation)
    for clause in (
        "invocation.input_sha256",
        "invocation.project_sha256",
        "invocation_identity_sha256",
    ):
        assert clause in source, clause

    # Molecular identity is node-specific.  It is checked while the host
    # builds the execution invocation from the approved node, not against the
    # workflow's aggregate multi-geometry identity at launch reservation.
    identity_source = inspect.getsource(
        execution.build_program_execution_invocation
    )
    assert (
        "identity_sha256 != node.scientific_identity_sha256" in identity_source
    )
    assert "!= plan.scientific_identity_sha256" not in source
