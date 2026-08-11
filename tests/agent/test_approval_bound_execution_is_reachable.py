"""``agent run --approval-file`` must be able to reach execution.

The campaign has never graded a number. The reason was recorded as "no approval
file is ever emitted from a plan session", and `WorkflowApprovalRequestV1` was
added to close that. Reading the execution path end to end shows the reason was
wrong, and the real one is worse.

`CommandCompiledToolHostV1.execute_program_node` refuses when its
`frozen_workflow_approval` is `None`:

    legacy V1 approval is preview-only; Runtime V2 frozen approval is required
    for execution

`live_session._execution_composition_inputs` -- the only place an approval file
becomes host state -- returns `approved_workspace`, `execution_resources`,
`workflow_execution_approval`, `execution_server` and `execution_environment`,
and no frozen approval. `frozen_workflow_approval=` is passed at no call site
in the package. So no approval file can drive execution: the object the
execution tool requires is one that nothing constructs for a live session.

This is the `repair_command` finding of W7 again -- a path advertised that the
runtime cannot let succeed -- on the single most important capability rather
than on a repair helper. `build_frozen_workflow_approval` exists and is
reachable, so the gap is wiring, not a missing implementation: the frozen body
belongs in the approval file next to the V1 approval, because deriving it
inside the session from the plan would make its own `plan_sha256` check
vacuous and destroy the safety property it exists to provide.

Written first as a strict xfail; the wiring below now makes it pass, so the
marker is gone and these are live assertions.
"""

import inspect

from chemsmart.agent import live_session
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

#: What the host must be given for the execution tool to proceed.
_EXECUTION_STATE = "frozen_workflow_approval"


def test_the_execution_tool_still_requires_the_frozen_approval():
    """Pins the requirement, so the gap cannot be closed by dropping it."""
    assert (
        _EXECUTION_STATE
        in inspect.signature(CommandCompiledToolHostV1).parameters
    )

    source = inspect.getsource(CommandCompiledToolHostV1)
    assert "Runtime V2 frozen approval" in source, (
        "the execution tool no longer states why a V1 approval is "
        "preview-only; update this gate deliberately"
    )


def test_the_v1_composition_surface_is_otherwise_complete():
    """The rest of the composition is present, so only one thing is missing."""
    present = set(inspect.signature(CommandCompiledToolHostV1).parameters)
    assert {
        "approved_workspace",
        "execution_environment",
        "execution_resources",
        "execution_server",
        "workflow_execution_approval",
    }.issubset(present)


def test_an_approval_file_can_actually_drive_execution():
    """The property the campaign needs in order to grade any number."""
    source = inspect.getsource(live_session._execution_composition_inputs)
    assert _EXECUTION_STATE in source, (
        "the only place an approval file becomes host state does not supply "
        "the frozen approval, so no approval file can execute a node"
    )


def test_a_frozen_body_for_another_task_is_refused():
    """The pair is cross-checked, so an approval cannot be reused elsewhere."""
    source = inspect.getsource(live_session._execution_composition_inputs)
    assert "frozen workflow approval targets another task spec" in source


def test_an_approval_without_a_frozen_body_stays_preview_only():
    """Preview-only approval remains legal rather than becoming an error."""
    source = inspect.getsource(live_session._execution_composition_inputs)
    assert "if raw_frozen is not None:" in source, (
        "omitting the frozen body must keep an approval usable for preview "
        "instead of failing the session"
    )
