from __future__ import annotations

from rich.console import Group

from chemsmart.agent.tui.tool_meta import (
    _TOOL_META,
    render_tool_result_detail,
    render_tool_result_summary,
)
from chemsmart.agent.tui.widgets.cells.tool_call_cell import ToolCallCell


def test_tool_call_cell_renders_read_result_summary_and_detail():
    result = {
        "path": "/tmp/sample.txt",
        "start_line": 10,
        "end_line": 16,
        "total_lines": 200,
        "truncated": True,
        "content": "\n".join(
            [
                "    10\talpha",
                "    11\tbeta",
                "    12\tgamma",
                "    13\tdelta",
                "    14\tepsilon",
                "    15\tzeta",
                "    16\teta",
            ]
        ),
    }

    cell = ToolCallCell(
        tool="read",
        status="ok",
        description="Read local file lines",
        result=result,
    )

    renderables = _renderables(cell)
    assert any(
        item.plain == "Result · L10-16 of 200, truncated"
        for item in renderables
    )
    detail = render_tool_result_detail("read", result)
    assert detail is not None
    assert len(detail) <= 6
    assert detail[-1].plain == "… truncated"


def test_tool_call_cell_renders_ssh_probe_result_summary_and_detail():
    result = {
        "server": "cluster-a",
        "probe": "survey.slurm.scontrol_partition",
        "returncode": 0,
        "stdout_truncated": "PartitionName=debug\nState=UP\nTotalNodes=4",
        "stderr_truncated": "",
        "parsed": {"queues": [{"name": "debug"}]},
        "duration_s": 0.126,
    }

    cell = ToolCallCell(
        tool="ssh_probe",
        status="ok",
        description="Run a remote inspection probe",
        result=result,
    )

    renderables = _renderables(cell)
    assert any(item.plain == "Result · rc=0 in 0.13s" for item in renderables)
    detail = render_tool_result_detail("ssh_probe", result)
    assert detail is not None
    assert len(detail) <= 6
    assert detail[0].plain == "server: cluster-a"
    assert detail[1].plain == "probe: survey.slurm.scontrol_partition"


def test_tool_call_cell_renders_scheduler_query_result_summary_and_detail():
    queue_result = {
        "scheduler": "slurm",
        "partition_or_queue": "debug",
        "nodes": 4,
        "total_cpus": 128,
        "raw": {"queues": [{"name": "debug"}]},
    }
    job_result = {
        "scheduler": "slurm",
        "job_id": "12345",
        "state": "RUNNING",
        "queue": "debug",
        "user": "alice",
        "node": "node001",
    }

    cell = ToolCallCell(
        tool="scheduler_query",
        status="ok",
        description="Inspect remote scheduler state",
        result=queue_result,
    )

    renderables = _renderables(cell)
    assert any(
        item.plain == "Result · slurm: debug 4n/128cpu" for item in renderables
    )
    assert render_tool_result_summary("scheduler_query", job_result) == (
        "slurm: RUNNING"
    )
    detail = render_tool_result_detail("scheduler_query", queue_result)
    assert detail is not None
    assert len(detail) <= 6
    assert [line.plain for line in detail] == [
        "scheduler: slurm",
        "queue: debug",
        "nodes: 4",
        "cpus: 128",
    ]


def test_tool_call_cell_renders_log_tail_result_summary_and_detail():
    result = {
        "server": "cluster-a",
        "path": "/scratch/job.log",
        "lines_requested": 200,
        "lines_returned": 12,
        "content_truncated": "...",
        "errors": [
            {
                "kind": "oom_killed",
                "line": "slurmstepd: error: Detected 1 oom-kill event(s)",
                "line_no": 4,
            },
            {
                "kind": "walltime_exceeded",
                "line": "Job exceeded walltime limit",
                "line_no": 8,
            },
            {
                "kind": "missing_module",
                "line": "module: command not found",
                "line_no": 11,
            },
            {
                "kind": "segfault",
                "line": "Segmentation fault",
                "line_no": 12,
            },
        ],
        "duration_s": 0.2,
    }

    cell = ToolCallCell(
        tool="log_tail",
        status="ok",
        description="Tail remote logs for diagnostics",
        result=result,
    )

    renderables = _renderables(cell)
    assert any(
        item.plain
        == "Result · 12L, 4 errors: oom_killed, walltime_exceeded, missing_module"
        for item in renderables
    )
    detail = render_tool_result_detail("log_tail", result)
    assert detail is not None
    assert len(detail) <= 6
    assert detail[-1].plain == "… truncated"


def test_tool_call_cell_result_none_stays_backward_compatible():
    cell = ToolCallCell(
        tool="build_job",
        status="pending",
        description="Build",
        arguments={"a": 1},
        result=None,
    )

    renderables = _renderables(cell)
    assert [item.plain for item in renderables] == [
        "● PENDING build_job [local-state]",
        "Build",
        "effect: creates a job object handle",
        "",
        '{\n  "a": 1\n}',
    ]


def test_tool_meta_existing_entries_are_unchanged_and_mva_entries_exist():
    expected_existing = {
        "build_molecule": {
            "risk": "read-only",
            "read_only": True,
            "style": "warning",
            "summary": "reads a structure file",
        },
        "recommend_method": {
            "risk": "read-only",
            "read_only": True,
            "style": "warning",
            "summary": "computes an advisory method recommendation",
        },
        "build_gaussian_settings": {
            "risk": "read-only",
            "read_only": True,
            "style": "warning",
            "summary": "builds validated Gaussian settings",
        },
        "build_orca_settings": {
            "risk": "read-only",
            "read_only": True,
            "style": "warning",
            "summary": "builds validated ORCA settings",
        },
        "build_job": {
            "risk": "local-state",
            "read_only": False,
            "style": "warning",
            "summary": "creates a job object handle",
        },
        "dry_run_input": {
            "risk": "mutates-state",
            "read_only": False,
            "style": "warning",
            "summary": "writes a local input file",
        },
        "validate_runtime": {
            "risk": "inspection",
            "read_only": True,
            "style": "warning",
            "summary": "checks local/runtime state",
        },
        "run_local": {
            "risk": "risky",
            "read_only": False,
            "style": "error",
            "summary": "starts a local job",
        },
        "extract_optimized_geometry": {
            "risk": "read-only",
            "read_only": True,
            "style": "warning",
            "summary": "reads output logs to extract geometry",
        },
        "submit_hpc": {
            "risk": "risky",
            "read_only": False,
            "style": "error",
            "summary": "may submit to a remote queue",
        },
        "wizard_probe": {
            "risk": "inspection",
            "read_only": True,
            "style": "warning",
            "summary": "probes local or remote server state",
        },
        "wizard_refresh": {
            "risk": "inspection",
            "read_only": True,
            "style": "warning",
            "summary": "refreshes or reuses the wizard node cache",
        },
        "wizard_verify": {
            "risk": "inspection",
            "read_only": True,
            "style": "warning",
            "summary": "verifies wizard/server transport wiring",
        },
        "wizard_write": {
            "risk": "risky",
            "read_only": False,
            "style": "error",
            "summary": "writes a server YAML config",
        },
    }

    assert len(expected_existing) == 14
    assert {
        key: _TOOL_META[key] for key in expected_existing
    } == expected_existing
    assert _TOOL_META["read"] == {
        "risk": "inspection",
        "read_only": True,
        "style": "warning",
        "summary": "reads local file lines",
    }
    assert _TOOL_META["ssh_probe"] == {
        "risk": "inspection",
        "read_only": True,
        "style": "warning",
        "summary": "runs a remote inspection probe",
    }
    assert _TOOL_META["scheduler_query"] == {
        "risk": "inspection",
        "read_only": True,
        "style": "warning",
        "summary": "inspects remote scheduler state",
    }
    assert _TOOL_META["log_tail"] == {
        "risk": "inspection",
        "read_only": True,
        "style": "warning",
        "summary": "tails remote logs for diagnostics",
    }


def test_render_tool_result_error_uses_compact_error_summary_only():
    result = {
        "error": "invalid_path",
        "path": "/tmp/job.log",
        "message": "Path must be absolute.",
        "stdout": '{"should_not": "show"}',
    }

    cell = ToolCallCell(
        tool="log_tail",
        status="error",
        description="Tail remote logs for diagnostics",
        result=result,
    )

    renderables = _renderables(cell)
    summary_line = next(
        item for item in renderables if item.plain.startswith("Result ·")
    )
    assert summary_line.plain == "Result · error: invalid_path"
    assert summary_line.spans[-1].style == "error"
    assert not any('"should_not"' in item.plain for item in renderables)
    assert render_tool_result_summary("log_tail", result) == (
        "error: invalid_path"
    )


def _renderables(cell: ToolCallCell):
    group = cell.renderable
    assert isinstance(group, Group)
    return list(group.renderables)
