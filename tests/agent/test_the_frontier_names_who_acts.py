"""An actionable analysis node names the tool that executes it.

A woken session read "actionable: execute the registered analysis
operation" on fourteen thermochemistry nodes, took it for a stage
awaiting approval, and stopped with every number in hand (NOVEL-2 po2,
2026-09-04). The host knows the tool for each analysis kind; the
frontier says so on every actionable analysis node, and a node the host
already executed under a wake arrives completed.
"""

from __future__ import annotations

import pytest

from chemsmart.agent.scientific_toolchain import (
    ANALYSIS_NEXT_TOOL,
    project_scientific_toolchain_frontier,
)

from .test_an_analysis_only_plan_executes_at_wake import _registered_chain

pytestmark = pytest.mark.capability("tool:inspect_workflow_frontier")


def test_every_analysis_kind_names_its_tool():
    from chemsmart.agent.scientific_toolchain import (
        ANALYSIS_INTENT_KINDS as ANALYSIS_KINDS,
    )

    executable = [
        kind for kind in ANALYSIS_KINDS if kind != "unsupported_external"
    ]
    assert set(executable) == set(ANALYSIS_NEXT_TOOL)


def test_an_actionable_analysis_node_carries_next_tool():
    frontier = project_scientific_toolchain_frontier(
        _registered_chain("orca-result-8ae1cdc683f8eb7d")
    )
    by_id = {node["node_id"]: node for node in frontier["nodes"]}
    root = by_id["extract-registered"]
    assert root["state"] == "actionable"
    assert root["next_tool"] == "extract_result_quantities"
    assert (
        "this session calls extract_result_quantities" in root["next_action"]
    )
    assert by_id["claims"]["state"] == "waiting_for_artifact"


def test_a_host_executed_node_arrives_completed():
    frontier = project_scientific_toolchain_frontier(
        _registered_chain("orca-result-8ae1cdc683f8eb7d"),
        completed_analysis_node_ids=("extract-registered", "claims"),
    )
    assert {node["state"] for node in frontier["nodes"]} == {"completed"}
    assert all("next_tool" not in node for node in frontier["nodes"])
