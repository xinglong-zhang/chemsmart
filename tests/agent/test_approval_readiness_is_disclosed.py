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
    assert "non_executable_node_ids" in source
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

    prompt = live_session._system_prompt({}, bounded_review_requested=True)
    for phrase in (
        "When the task names a program, plan that program",
        "cannot preview green",
        "deferred_admissible",
        "keep that causal stage",
        "non_executable",
        "will not be approved or launched",
        "approval_readiness",
    ):
        assert phrase in prompt, phrase
    assert "inert exact workflow review" in prompt
    assert "Execution is not exposed" in prompt
    assert "execute_approved_program_node" not in prompt
    assert "plan again without" not in prompt

    context = live_session._public_context(
        task="bounded review",
        task_spec_sha256="a" * 64,
        observations=(),
        conformance_records=(),
        registry_sha256="b" * 64,
        live_schema_sha256="c" * 64,
        execution_requested=False,
        execution_available=False,
        execution_review_requested=True,
        bounded_execution_record={
            "schema_version": (
                "chemsmart.public-bounded-execution-envelope.v1"
            ),
            "mode": "bounded-local",
            "resources": {"cores": 1, "memory_gb": 4},
        },
    )

    assert context["execution_requested"] is False
    assert context["execution_tool_available"] is False
    assert context["execution_review_requested"] is True
    assert context["approved_execution_contract"] == {}
    assert context["bounded_execution_envelope"]["mode"] == "bounded-local"
    assert "scratch_root" not in context["bounded_execution_envelope"]
