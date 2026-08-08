"""A plan must be told which of its nodes block its own approval.

Freezing an approval requires a green preview for *every* materialized node,
so one unpreviewable node makes the whole plan unexecutable. Until now that
rule only spoke at freeze time -- after the session had ended.

A live session hit it exactly: it planned Gaussian, ORCA and PySCF for one
task, previewed PySCF green immediately, could not repair the other two, and
left their nodes in place. The session terminated cleanly with 66 of 70 tool
calls succeeding, and the plan was unfreezable:

    FREEZE REFUSED -> every currently materialized approved node requires
                      exact preview

Knowing this is what makes program fallback possible rather than fatal: keep
the program the task named while it can be repaired, and if it truly cannot
preview green, plan again without that node instead of carrying it.
"""

import inspect

from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1


def test_the_readiness_predicate_is_the_materializers_own():
    """Two spellings of "previewed" would let the disclosure lie."""

    materialize = inspect.getsource(
        CommandCompiledToolHostV1._materialize_scientific_workflow
    )
    readiness = inspect.getsource(
        CommandCompiledToolHostV1._approval_readiness
    )
    assert "_node_is_previewed" in materialize
    assert "_node_is_previewed" in readiness


def test_readiness_names_the_blocking_nodes_not_only_a_count():
    source = inspect.getsource(CommandCompiledToolHostV1._approval_readiness)
    assert "blocking_node_ids" in source
    assert "approvable" in source


def test_readiness_states_the_rule_it_is_enforcing():
    """A verdict without its rule cannot be acted on."""

    source = inspect.getsource(CommandCompiledToolHostV1._approval_readiness)
    assert '"rule"' in source
    for phrase in (
        "green preview",
        "plan the workflow again",
        "preview_command",
    ):
        assert phrase in source, phrase


def test_planning_returns_readiness_to_the_model():
    source = inspect.getsource(
        CommandCompiledToolHostV1._plan_command_workflow
    )
    assert '"approval_readiness"' in source


def test_the_prompt_states_program_preference_and_its_fallback():
    """Fallback without a rule to drop the abandoned node repeats the bug."""

    from chemsmart.agent import live_session

    source = inspect.getsource(live_session)
    for phrase in (
        "When the task names a program, plan that program",
        "cannot preview green",
        "plan again without",
        "approval_readiness",
    ):
        assert phrase in source, phrase
