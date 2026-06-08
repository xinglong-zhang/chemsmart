from __future__ import annotations

from pathlib import Path

from chemsmart.agent.permissions import RuntimePermissionMode
from chemsmart.agent.registry import ToolRegistry
from chemsmart.agent.tools_fs import read


def test_read_happy_path_defaults(tmp_path, monkeypatch):
    target = tmp_path / "notes.txt"
    target.write_text("alpha\nbeta\ngamma\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)

    result = read("notes.txt")

    assert result["path"] == str(target.resolve())
    assert result["start_line"] == 1
    assert result["end_line"] == 3
    assert result["total_lines"] == 3
    assert result["truncated"] is False
    assert result["content"].startswith("     1\talpha")
    assert result["content"].splitlines() == [
        "     1\talpha",
        "     2\tbeta",
        "     3\tgamma",
    ]


def test_read_supports_paging(tmp_path, monkeypatch):
    target = tmp_path / "paged.txt"
    target.write_text(
        "\n".join(f"line {index}" for index in range(1, 21)) + "\n",
        encoding="utf-8",
    )
    monkeypatch.chdir(tmp_path)

    result = read("paged.txt", start_line=10, limit=5)

    assert result["start_line"] == 10
    assert result["end_line"] == 14
    assert result["total_lines"] == 20
    assert result["truncated"] is True
    assert result["content"].splitlines() == [
        "    10\tline 10",
        "    11\tline 11",
        "    12\tline 12",
        "    13\tline 13",
        "    14\tline 14",
    ]


def test_read_marks_truncation_when_limit_cuts_file(tmp_path, monkeypatch):
    target = tmp_path / "large.txt"
    target.write_text(
        "\n".join(f"row {index}" for index in range(1, 8)) + "\n",
        encoding="utf-8",
    )
    monkeypatch.chdir(tmp_path)

    result = read("large.txt", limit=3)

    assert result["end_line"] == 3
    assert result["total_lines"] == 7
    assert result["truncated"] is True
    assert result["content"].splitlines()[-1] == "     3\trow 3"


def test_read_rejects_missing_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    result = read("missing.txt")

    assert result["error"] == "file_not_found"
    assert result["path"] == str((tmp_path / "missing.txt").resolve())


def test_read_rejects_binary_files(tmp_path, monkeypatch):
    target = tmp_path / "data.bin"
    target.write_bytes(b"hello\x00world")
    monkeypatch.chdir(tmp_path)

    result = read("data.bin")

    assert result["error"] == "binary_file_refused"
    assert result["path"] == str(target.resolve())


def test_read_rejects_paths_outside_cwd(tmp_path, monkeypatch):
    outside = tmp_path.parent / "outside.txt"
    outside.write_text("secret\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)

    result = read("../outside.txt")

    assert result["error"] == "path_outside_cwd"
    assert Path(result["path"]) == outside.resolve()
    assert result["cwd"] == str(tmp_path.resolve())


def test_read_tool_is_registered_as_read_only():
    registry = ToolRegistry.default()

    tool = registry.get_tool("read")

    assert tool is not None
    assert tool.metadata.read_only is True


def test_read_tool_pool_visibility_matches_runtime_modes():
    registry = ToolRegistry.default()

    read_only_names = {
        tool.name
        for tool in registry.assemble_tool_pool(
            RuntimePermissionMode.READ_ONLY
        )
    }
    edit_names = {
        tool.name
        for tool in registry.assemble_tool_pool(
            RuntimePermissionMode.ACCEPT_EDITS
        )
    }
    bypass_names = {
        tool.name
        for tool in registry.assemble_tool_pool(RuntimePermissionMode.BYPASS)
    }
    plan_names = {
        tool.name
        for tool in registry.assemble_tool_pool(RuntimePermissionMode.PLAN)
    }

    assert "read" in read_only_names
    assert "read" in edit_names
    assert "read" in bypass_names
    assert "read" not in plan_names
