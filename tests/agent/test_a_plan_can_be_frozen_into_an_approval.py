"""A plan must be freezable, or the documented approval flow is circular.

`agent plan`, `agent review`, then `agent run --approval-file` is the documented flow, and it
could not be walked. Wiring the frozen approval into the session composition
was only the outer half; the inner half is that nobody could *build* the frozen
body from a plan at all.

`build_materialized_workflow` records `resource_sha256` from the host's
execution resources, and a plan session has none -- it makes no engine call --
so it records the preview sentinel instead. `build_frozen_workflow_approval`
then demanded

    materialized_workflow.resource_sha256 == resources.resource_sha256

against the locked *run* profile. A plan's materialization can therefore never
be frozen. The only materialization that satisfied the equality came from a
session that already held execution resources, and holding those requires an
approval, which requires a freeze. Nothing could produce the first one. That is
why `frozen_workflow_approval=` appeared at no call site in the package: not an
oversight, an unreachable precondition.

Binding resources is the reviewer's act, so the freeze is where it belongs. The
sentinel is now accepted and the resources come from the freeze. The safety
chain is unchanged and pinned below: a workflow materialized under a *real*
profile must still match it, so an approval cannot be carried between resource
profiles.
"""

import pytest

from chemsmart.agent.execution import (
    build_execution_resource_spec,
    build_frozen_workflow_approval,
)
from chemsmart.agent.workflows import PREVIEW_RESOURCE_SHA256


def _run_resources(cores=4):
    return build_execution_resource_spec(
        execution_target="run",
        cores=cores,
        memory_gb=4.0,
        gpu_count=0,
        scratch_policy="none",
        node_timeout_seconds=600,
    )


def test_the_preview_sentinel_is_what_a_plan_session_records():
    """Pin the sentinel's identity; it is a contract between two modules."""
    from chemsmart.agent._contracts import canonical_sha256

    assert PREVIEW_RESOURCE_SHA256 == canonical_sha256(
        {
            "schema_version": "chemsmart.preview-resource.v1",
            "chemistry_engine_calls": 0,
        }
    )
    assert PREVIEW_RESOURCE_SHA256 != _run_resources().resource_sha256


def test_the_freeze_accepts_a_plan_materialization():
    """The property that makes the documented flow walkable at all."""
    import inspect

    source = inspect.getsource(build_frozen_workflow_approval)
    assert "PREVIEW_RESOURCE_SHA256" in source, (
        "a plan session materializes under the preview sentinel; refusing it "
        "makes plan -> approve -> run circular and unreachable"
    )


def test_a_real_resource_profile_must_still_match():
    """The safety property: an approval cannot move between profiles."""
    import inspect

    source = inspect.getsource(build_frozen_workflow_approval)
    assert "resources.resource_sha256" in source, (
        "a workflow materialized under a real profile must still be checked "
        "against the profile being frozen"
    )
    assert "materialized workflow resource binding differs" in source


def test_two_run_profiles_are_distinguishable():
    """If they were not, the retained check would be vacuous."""
    assert (
        _run_resources(cores=4).resource_sha256
        != _run_resources(cores=8).resource_sha256
    )


def test_a_frozen_approval_still_refuses_a_foreign_plan():
    """The other bindings the freeze checks are untouched by this change."""
    import inspect

    source = inspect.getsource(build_frozen_workflow_approval)
    for message in (
        "materialized workflow belongs to another plan",
        "materialized workflow plan digest differs",
        "materialized workflow task identity differs",
        "materialized workflow scientific identity differs",
    ):
        assert message in source, message


@pytest.mark.parametrize(
    "field",
    [
        "workflow_id",
        "plan_sha256",
        "task_spec_sha256",
        "scientific_identity_sha256",
    ],
)
def test_every_other_binding_is_still_compared(field):
    import inspect

    source = inspect.getsource(build_frozen_workflow_approval)
    assert f"materialized_workflow.{field}" in source
