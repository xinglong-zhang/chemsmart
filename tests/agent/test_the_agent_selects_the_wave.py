"""The Agent selects a wave; the host says what it sees about it.

Readiness is a host/DAG fact with one authority, and this tool does not
compute a second one: it asks ``workflow_context`` -- the same projection
`inspect_workflow_frontier` already serves the model -- and judges the
Agent's proposal against that frontier and the plan's own producer edges.

Which ready calculations belong in one wave is scientific strategy and
stays the Agent's: three conformers together is a decision about what
evidence to see at once, and the host has no basis for it. So this tool
never chooses, never reorders, and never refuses. An undispatchable wave
comes back as a typed per-member verdict the Agent reads and selects from
again, because an exception is what teaches a session to carry
workarounds for something the host should simply have reported.
"""

from __future__ import annotations

from types import SimpleNamespace

from chemsmart.agent.cohort import validate_wave


def _context(ready, waiting=()):
    return SimpleNamespace(
        ready_node_ids=tuple(ready),
        waiting_node_ids=tuple(waiting),
    )


def _edges(*pairs):
    return tuple(pairs)


def test_a_wave_of_independent_ready_nodes_is_ready():
    verdict = validate_wave(
        proposed=("conf-a-opt", "conf-b-opt", "conf-c-opt"),
        ready=("conf-a-opt", "conf-b-opt", "conf-c-opt"),
        edges=(),
    )
    assert all(row.status == "ready" for row in verdict.rows)
    assert verdict.summary == "3 calculations ready to run together"


def test_the_tool_is_declared_and_reaches_a_handler():
    """A tool the model cannot call is a paragraph."""

    from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
    from chemsmart.agent.tool_specs import (
        build_command_compiled_tool_surface,
    )

    names = {
        item["function"]["name"]
        for item in build_command_compiled_tool_surface().tool_definitions
    }
    assert "select_execution_wave" in names, (
        "the Agent has no way to name a wave, so the cohort could only "
        "ever be chosen by the host -- which is the scientific decision "
        "this whole design moves to the model"
    )
    assert "select_execution_wave" in CommandCompiledToolHostV1.TOOL_HANDLERS


def test_the_declared_input_is_a_list_of_nodes_and_nothing_about_machines():
    """Scientific width is the Agent's; the machine is not."""

    from chemsmart.agent.tool_specs import (
        build_command_compiled_tool_surface,
    )

    spec = next(
        item
        for item in build_command_compiled_tool_surface().tool_definitions
        if item["function"]["name"] == "select_execution_wave"
    )
    properties = set((spec["function"]["parameters"].get("properties") or {}))
    assert properties == {"workflow_id", "node_ids"}, properties
    forbidden = {
        "cores",
        "num_cores",
        "threads",
        "num_threads",
        "memory_gb",
        "mem_gb",
        "max_concurrent",
        "concurrency",
        "queue",
        "qos",
        "account",
        "ntasks",
        "array",
    }
    assert not (properties & forbidden)


def test_the_host_reports_an_undispatchable_wave_rather_than_refusing():
    """Typed evidence, not an exception."""

    verdict = validate_wave(
        proposed=("conf-a-opt", "conf-a-hess"),
        ready=("conf-a-opt",),
        edges=_edges(("conf-a-opt", "conf-a-hess")),
    )
    statuses = {row.node_id: row.status for row in verdict.rows}
    assert statuses == {"conf-a-opt": "ready", "conf-a-hess": "depends_on"}
    detail = next(
        row.detail for row in verdict.rows if row.node_id == "conf-a-hess"
    )
    assert "conf-a-opt" in detail


def test_the_handler_asks_the_one_readiness_authority(monkeypatch, tmp_path):
    """No second frontier: the host's own projection is what is judged."""

    from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

    host = CommandCompiledToolHostV1.__new__(CommandCompiledToolHostV1)
    asked: list[str] = []

    resolved = SimpleNamespace(
        draft=SimpleNamespace(
            workflow_id="w1",
            nodes=(
                SimpleNamespace(node_id="a1", inputs=()),
                SimpleNamespace(node_id="a2", inputs=()),
                SimpleNamespace(
                    node_id="b1",
                    inputs=(
                        SimpleNamespace(
                            binding_id="filename",
                            producer_node_id="a1",
                            producer_output_id="geometry",
                        ),
                    ),
                ),
            ),
        ),
        scientific_plan=SimpleNamespace(plan_sha256="d" * 64),
    )

    def _resolve(workflow_id):
        asked.append(str(workflow_id))
        return resolved

    def _context_for(draft, *, scientific_plan_sha256=""):
        asked.append("workflow_context")
        return _context(("a1", "a2"), waiting=("b1",))

    host._resolve_program_workflow = _resolve
    host._workflow_context = _context_for

    reply = host._select_execution_wave(
        "t1", {"workflow_id": "w1", "node_ids": ["a1", "a2"]}
    )
    assert "workflow_context" in asked, (
        "the handler computed readiness itself instead of asking the "
        "authority that already owns it"
    )
    assert reply["status"] == "ready"
    assert tuple(reply["node_ids"]) == ("a1", "a2")

    # b1 consumes a1, which is in the same proposed wave: that is the
    # more specific verdict, and it names the sibling rather than the
    # bare fact that b1 is not on the frontier.
    mixed = host._select_execution_wave(
        "t2", {"workflow_id": "w1", "node_ids": ["a1", "b1"]}
    )
    assert mixed["status"] == "not_dispatchable"
    rows = {row["node_id"]: row["status"] for row in mixed["members"]}
    assert rows == {"a1": "ready", "b1": "depends_on"}

    # And a member that is simply not ready, with no sibling involved.
    alone = host._select_execution_wave(
        "t3", {"workflow_id": "w1", "node_ids": ["b1"]}
    )
    assert {row["node_id"]: row["status"] for row in alone["members"]} == {
        "b1": "not_ready"
    }


def test_a_selected_wave_is_recorded_where_the_dispatcher_reads_it(
    monkeypatch, tmp_path
):
    """The selection must survive to dispatch, or it is only a reply.

    The model selects in the planning session; the dispatcher runs in the
    driver afterwards. A wave that lives only in a tool reply is a wave
    the array never hears about.
    """

    from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

    host = CommandCompiledToolHostV1.__new__(CommandCompiledToolHostV1)
    host._resolve_program_workflow = lambda workflow_id: SimpleNamespace(
        draft=SimpleNamespace(
            workflow_id="w1",
            nodes=(
                SimpleNamespace(node_id="a1", inputs=()),
                SimpleNamespace(node_id="a2", inputs=()),
            ),
        ),
        scientific_plan=SimpleNamespace(plan_sha256="d" * 64),
    )
    host._workflow_context = lambda draft, **_kw: _context(("a1", "a2"))

    host._select_execution_wave(
        "t1", {"workflow_id": "w1", "node_ids": ["a2", "a1"]}
    )
    assert host.selected_execution_wave == ("a2", "a1"), (
        "the order the Agent chose is the order the array elements take, "
        "and the host must not reorder a scientific decision"
    )


def test_an_invalid_wave_is_not_recorded_as_the_selection():
    from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

    host = CommandCompiledToolHostV1.__new__(CommandCompiledToolHostV1)
    host._resolve_program_workflow = lambda workflow_id: SimpleNamespace(
        draft=SimpleNamespace(
            workflow_id="w1",
            nodes=(
                SimpleNamespace(node_id="a1", inputs=()),
                SimpleNamespace(
                    node_id="b1",
                    inputs=(
                        SimpleNamespace(
                            binding_id="filename",
                            producer_node_id="a1",
                            producer_output_id="geometry",
                        ),
                    ),
                ),
            ),
        ),
        scientific_plan=SimpleNamespace(plan_sha256="d" * 64),
    )
    host._workflow_context = lambda draft, **_kw: _context(("a1",))

    reply = host._select_execution_wave(
        "t1", {"workflow_id": "w1", "node_ids": ["a1", "b1"]}
    )
    assert reply["status"] == "not_dispatchable"
    assert getattr(host, "selected_execution_wave", ()) == ()
