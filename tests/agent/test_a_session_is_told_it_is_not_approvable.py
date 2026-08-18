"""A session that cannot get its workflow approved should be told so.

A live session ends on its first assistant message with no tool calls. At that
moment the host classified the stop as "workflow draft recorded; action-grade
gates remain" -- which names nothing -- and the model's own prose was preserved
verbatim. In the cycle-038 paper task that prose said the workflow was "ready
... for approval" while a node still blocked it, and nothing contradicted it.

The host already knew. `_approval_readiness` computes the blocking node ids on
every call, but it was attached to exactly one receipt, the plan's own, which a
session reads *before* it prepares anything. So the only snapshot a session ever
saw predated its own work.

Two repairs, checked here: the readiness travels with the frontier -- the tool
the system prompt sends a session to for host-derived next actions -- and the
stop classification names what still blocks approval.

This is host-owned status, not a grade: the ids are the host's own bookkeeping,
and nothing here judges the chemistry.
"""

from __future__ import annotations

import pytest

from chemsmart.agent.loop import _unapprovable_reason


def test_the_reason_names_the_nodes_that_block():
    reason = _unapprovable_reason(
        {
            "workflow_id": "bal-conformers",
            "blocking_node_ids": ("cal-conf-psimin", "cal-opt631"),
            "deferred_node_ids": ("cal-scan-psi",),
            "workflow_blocked_reason": "",
        }
    )

    assert "not approvable" in reason
    assert "cal-conf-psimin" in reason
    assert "cal-opt631" in reason


def test_a_stated_workflow_reason_wins_over_the_node_list():
    """When the host has an explicit reason, that is the more useful sentence."""

    reason = _unapprovable_reason(
        {
            "blocking_node_ids": ("n1",),
            "workflow_blocked_reason": "no executable stage remains",
        }
    )

    assert "no executable stage remains" in reason


def test_a_workflow_with_nothing_to_say_still_says_it_is_not_approvable():
    assert "not approvable" in _unapprovable_reason({})


def test_an_approvable_workflow_reports_nothing(tmp_path):
    """`None` is the signal that the stop needs no correction at all.

    Checked against the real host rather than a stand-in, because the point is
    that an analysis-only session -- which holds no scientific plan -- must be
    left exactly as it was.
    """

    from chemsmart.agent.runtime.event_store import RuntimeEventStore
    from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="s1"
        ),
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "workspace",
    )

    assert host.unapproved_workflow_summary() is None


@pytest.mark.parametrize(
    "method", ("unapproved_workflow_summary",)
)
def test_the_loop_host_protocol_is_declared_by_the_real_host(method):
    """The loop calls this at every stop, so the host must actually have it."""

    from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

    assert callable(getattr(CommandCompiledToolHostV1, method))
