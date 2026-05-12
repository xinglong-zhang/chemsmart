from __future__ import annotations

import pytest

from chemsmart.agent.permissions import RuntimePermissionMode
from chemsmart.agent.registry import ToolRegistry
from chemsmart.agent.tool_protocol import RuntimeToolMetadata
from chemsmart.agent.tools_hpc import log_tail
from chemsmart.agent.transport import ExecResult, MockExecTransport


def test_log_tail_happy_path_returns_content_without_errors(monkeypatch):
    transport = MockExecTransport(
        scripted_results=[
            ExecResult(
                returncode=0,
                stdout="line 1\nline 2\n",
                stderr="",
                command="",
                server="cluster-a",
                duration_s=0.2,
            )
        ]
    )
    monkeypatch.setattr(
        "chemsmart.agent.tools_hpc._TRANSPORT_FACTORY",
        lambda: transport,
    )

    result = log_tail(
        server="cluster-a",
        path="/scratch/job.log",
        lines=5,
    )

    assert result == {
        "server": "cluster-a",
        "path": "/scratch/job.log",
        "lines_requested": 5,
        "lines_returned": 2,
        "content_truncated": "line 1\nline 2",
        "errors": [],
        "duration_s": 0.2,
    }
    assert transport.calls == [
        {
            "command": "tail -n 5 /scratch/job.log",
            "server": "cluster-a",
            "timeout_s": 15,
        }
    ]


def test_log_tail_detects_oom_errors(monkeypatch):
    transport = MockExecTransport(
        scripted_results=[
            ExecResult(
                returncode=0,
                stdout=(
                    "normal line\n"
                    "slurmstepd: error: Detected 1 oom-kill event(s)\n"
                ),
                stderr="",
                command="",
                server="cluster-a",
                duration_s=0.1,
            )
        ]
    )
    monkeypatch.setattr(
        "chemsmart.agent.tools_hpc._TRANSPORT_FACTORY",
        lambda: transport,
    )

    result = log_tail(server="cluster-a", path="/scratch/job.log")

    assert result["errors"][0]["kind"] == "oom_killed"
    assert result["errors"][0]["line_no"] == 2


@pytest.mark.parametrize("lines", [0, 10001])
def test_log_tail_rejects_invalid_line_counts(lines):
    result = log_tail(
        server="cluster-a",
        path="/scratch/job.log",
        lines=lines,
    )

    assert result == {"error": "invalid_lines"}


@pytest.mark.parametrize(
    "grep",
    [
        "x" * 129,
        "oom\nkill",
        "(",
    ],
)
def test_log_tail_rejects_invalid_grep(grep):
    result = log_tail(
        server="cluster-a",
        path="/scratch/job.log",
        grep=grep,
    )

    assert result == {"error": "invalid_grep"}


def test_log_tail_rejects_invalid_path_before_transport_exec(monkeypatch):
    transport = MockExecTransport()
    monkeypatch.setattr(
        "chemsmart.agent.tools_hpc._TRANSPORT_FACTORY",
        lambda: transport,
    )

    result = log_tail(
        server="cluster-a",
        path="relative/job.log",
    )

    assert result["error"].startswith("invalid_")
    assert "path" in result["error"]
    assert transport.calls == []


def test_log_tail_registry_metadata_and_runtime_modes():
    registry = ToolRegistry.default()
    tool = registry.get_tool("log_tail")

    assert tool is not None
    assert tool.metadata == RuntimeToolMetadata(
        read_only=True,
        ui_summary_template="Tail {path} on {server} ({lines}L)",
    )

    expected = {
        RuntimePermissionMode.READ_ONLY: True,
        RuntimePermissionMode.ACCEPT_EDITS: True,
        RuntimePermissionMode.BYPASS: True,
        RuntimePermissionMode.PLAN: False,
    }
    for mode, present in expected.items():
        names = {item.name for item in registry.assemble_tool_pool(mode)}
        assert ("log_tail" in names) is present
