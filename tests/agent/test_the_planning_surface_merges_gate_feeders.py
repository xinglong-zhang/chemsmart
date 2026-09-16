"""The model reads one tool where it used to read four.

Eleven of the withdrawn tools existed so the model could carry a receipt
digest from one call to the next; the merged tool does the join and
returns every receipt, so a red preview costs one turn. The executor
still drives nodes through the legacy names on its own surface, one step
at a time, and nothing about receipts or events changed.
"""

from __future__ import annotations

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.executor import PROGRAM_NODE_SEQUENCE
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.agent.tool_specs import (
    MERGED_PLANNING_TOOLS,
    build_approved_execution_tool_surface,
    build_command_compiled_tool_surface,
)

pytestmark = pytest.mark.capability("tool:*")


def _names(surface):
    return {item["function"]["name"] for item in surface.tool_definitions}


def test_the_planning_surface_exposes_the_merged_tools_only():
    planning = _names(build_command_compiled_tool_surface())
    withdrawn = {
        legacy for group in MERGED_PLANNING_TOOLS.values() for legacy in group
    }
    assert set(MERGED_PLANNING_TOOLS) <= planning
    assert not (withdrawn & planning), sorted(withdrawn & planning)
    # 16 -> 17 (2026-09-09): fetch_pubchem_geometry joined the stem. It
    # belongs beside bind_scientific_identity rather than in the
    # structure guide, because it is how a molecule *enters* a session
    # and a guide the session never opens is an affordance it never
    # sees. Its bytes were paid for by de-duplicating the identifier
    # spelling rule, which the JSON pattern already states.
    # 17 -> 18 (2026-09-16, Round A): select_execution_wave joined the
    # stem. Which ready calculations belong in one wave is the Agent's
    # scientific strategy, and a tool it can only reach by opening a
    # guide is an affordance a session that never opens that guide does
    # not have -- so a wave would have been chosen by the host, which is
    # exactly the decision this design moves to the model.
    assert len(planning) == 18, "the stem"
    from chemsmart.agent.guides import GUIDES

    everything = _names(
        build_command_compiled_tool_surface(
            guides=tuple(guide.guide_id for guide in GUIDES)
        )
    )
    assert len(everything) == 29, "the stem with every leaf open"


def test_the_executor_keeps_its_step_by_step_surface():
    execution = _names(build_approved_execution_tool_surface())
    for name in PROGRAM_NODE_SEQUENCE:
        assert name in execution, name
    assert "execute_approved_program_node" in execution


def _host(tmp_path):
    return CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="merge-session"
        ),
        artifacts={},
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "workspace",
    )


def test_project_yaml_names_what_each_action_needs(tmp_path):
    host = _host(tmp_path)
    with pytest.raises(ContractError, match="needs.*missing.*sections"):
        host.dispatch(
            turn_id="t1",
            tool_name="project_yaml",
            arguments={"action": "render", "program": "orca"},
        )
    with pytest.raises(ContractError, match="action must be one of"):
        host._project_yaml("t1", {"action": "materialise"})


def test_inspect_run_needs_the_reader_for_a_result(tmp_path):
    host = _host(tmp_path)
    with pytest.raises(ContractError, match="program beside artifact_id"):
        host._inspect_run("t1", {"artifact_id": "result.x.1"})


def test_the_legacy_names_are_not_callable_from_the_planning_surface(
    tmp_path,
):
    host = _host(tmp_path)
    with pytest.raises(ContractError, match="not exposed"):
        host.dispatch(
            turn_id="t1",
            tool_name="preview_command",
            arguments={"invocation_sha256": "b" * 64},
        )
