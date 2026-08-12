"""A plan must distinguish preview blockers from causal future stages.

Freezing an exact approval requires a green preview for every materialized
node. Bounded continuous execution is intentionally different: a consumer of
an optimized geometry cannot preview before the producer runs, so its exact
producer-data edge is admitted as deferred rather than presented as a stage
the model should delete.
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
    assert "deferred_node_ids" in source
    assert "deferred_admissible" in source
    assert "approvable" in source


def test_readiness_states_the_rule_it_is_enforcing():
    """A verdict without its rule cannot be acted on."""

    source = inspect.getsource(CommandCompiledToolHostV1._approval_readiness)
    assert '"rule"' in source
    for phrase in (
        "green preview",
        "producer-data target",
        "preview_command",
        "do not delete",
    ):
        assert phrase in source, phrase


def test_planning_returns_readiness_to_the_model():
    source = inspect.getsource(
        CommandCompiledToolHostV1._plan_command_workflow
    )
    assert '"approval_readiness"' in source


def test_the_prompt_preserves_bounded_causal_future_nodes():
    """The system prompt must not contradict bounded readiness feedback."""

    from chemsmart.agent import live_session

    prompt = live_session._system_prompt(
        {"authorization_mode": "bounded_continuous"}
    )
    for phrase in (
        "When the task names a program, plan that program",
        "cannot preview green",
        "deferred_admissible",
        "keep that causal stage",
        "approval_readiness",
    ):
        assert phrase in prompt, phrase
    assert "plan again without" not in prompt
