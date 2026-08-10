"""The executor supplies the caller, never the science.

Measured across the campaign's live sessions: seventeen ``agent plan`` sessions
rejected 42 of 868 tool calls (4.8%); the one ``agent run`` session rejected 20
of 69 (29%) -- twelve execution attempts and eight re-plans *inside* an
execution session -- and completed no nodes. ``agent run`` is
``run_live_agent_session(execution_enabled=True)``: the same planning loop with
one extra tool, asked to re-derive a plan it was already handed.

These gates pin the division of labour that replaces it. They are deliberately
about what the executor may *not* do.
"""

import inspect

import pytest

from chemsmart.agent import executor as executor_module
from chemsmart.agent.execution import (
    argv_shape,
    invocation_identity_sha256,
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
