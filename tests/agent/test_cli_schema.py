"""Tests for the ChemSmart CLI schema introspector."""

from __future__ import annotations

import json
import os
from collections.abc import Generator
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner

from chemsmart.agent.cli_schema import (
    build_chemsmart_cli_schema,
    dump_schema_to_json,
)

SNAPSHOT_DIR = Path(__file__).parent / "snapshots"
TS_SNAPSHOT = SNAPSHOT_DIR / "cli_schema_sub_gaussian_ts.json"


@pytest.fixture(autouse=True)
def reset_deferred_entry_point_groups() -> Generator[None, None, None]:
    """Keep schema tests from changing process-global lazy CLI caches."""

    yield
    from chemsmart.cli.main import entry_point

    for command in entry_point.commands.values():
        if hasattr(command, "_loaded"):
            command._loaded = None


def _option_with(options: list[dict[str, Any]], opt: str) -> dict[str, Any]:
    for option in options:
        if opt in option["opts"]:
            return option
    raise AssertionError(f"option {opt!r} not found")


def _path_options(*nodes: dict[str, Any]) -> list[dict[str, Any]]:
    options: list[dict[str, Any]] = []
    for node in nodes:
        options.extend(node["options"])
    return options


def test_top_level_subcommands() -> None:
    schema = build_chemsmart_cli_schema()

    assert set(schema["subcommands"]) == {
        "run",
        "sub",
        "agent",
        "config",
        "update",
    }


def test_sub_gaussian_subcommands_include_expected_commands() -> None:
    schema = build_chemsmart_cli_schema()
    gaussian = schema["subcommands"]["sub"]["subcommands"]["gaussian"]

    assert {
        "opt",
        "ts",
        "irc",
        "com",
        "modred",
        "nci",
        "qrc",
        "resp",
        "wbi",
        "traj",
        "link",
        "crest",
    } <= set(gaussian["subcommands"])
    # Current Click command names for custom Gaussian jobs and TDDFT are the
    # syntactically valid CLI spellings.
    assert {"userjob", "td"} <= set(gaussian["subcommands"])
    assert "qmmm" in gaussian["subcommands"]["ts"]["subcommands"]


def test_sub_gaussian_ts_path_options_include_required_parent_options() -> (
    None
):
    schema = build_chemsmart_cli_schema()
    sub = schema["subcommands"]["sub"]
    gaussian = sub["subcommands"]["gaussian"]
    ts = gaussian["subcommands"]["ts"]
    options = _path_options(sub, gaussian, ts)

    assert _option_with(options, "--server")["opts"] == ["-s", "--server"]
    assert _option_with(options, "--server")["type"] == "str"
    assert _option_with(options, "--project")["opts"] == ["-p", "--project"]
    assert _option_with(options, "--project")["type"] == "str"
    assert _option_with(options, "--basis")["opts"] == ["-b", "--basis"]
    assert _option_with(options, "--basis")["type"] == "str"


def test_choice_and_path_types_round_trip() -> None:
    schema = build_chemsmart_cli_schema()
    dias = schema["subcommands"]["sub"]["subcommands"]["gaussian"][
        "subcommands"
    ]["dias"]
    mode = _option_with(dias["options"], "--mode")
    assert mode["type"] == {"type": "choice", "choices": ["irc", "ts"]}
    assert mode["choices"] == ["irc", "ts"]

    iterate = schema["subcommands"]["sub"]["subcommands"]["iterate"]
    directory = _option_with(iterate["options"], "--directory")
    assert directory["type"] == {
        "type": "path",
        "exists": False,
        "file_okay": False,
        "dir_okay": True,
    }


def test_schema_is_json_serializable() -> None:
    schema = build_chemsmart_cli_schema()

    json.dumps(schema)


def test_dump_schema_to_json_adds_metadata(tmp_path: Path) -> None:
    output = tmp_path / "schema.json"

    document = dump_schema_to_json(output)
    loaded = json.loads(output.read_text())

    assert loaded == document
    assert len(loaded["_meta"]["schema_hash"]) == 64
    assert loaded["_meta"]["chemsmart_version"]


def test_hidden_dump_cli_schema_command_outputs_json(tmp_path: Path) -> None:
    from chemsmart.cli.main import entry_point

    runner = CliRunner()

    result = runner.invoke(entry_point, ["agent", "_dump-cli-schema"])

    assert result.exit_code == 0
    document = json.loads(result.output)
    assert document["_meta"]["schema_hash"]

    output = tmp_path / "cli-schema.json"
    result = runner.invoke(
        entry_point, ["agent", "_dump-cli-schema", "--out", str(output)]
    )

    assert result.exit_code == 0
    assert result.output == ""
    assert json.loads(output.read_text())["_meta"]["schema_hash"]


def test_sub_gaussian_ts_schema_snapshot() -> None:
    schema = build_chemsmart_cli_schema()
    ts_schema = schema["subcommands"]["sub"]["subcommands"]["gaussian"][
        "subcommands"
    ]["ts"]
    rendered = json.dumps(ts_schema, indent=2, sort_keys=True) + "\n"

    if os.environ.get("CHEMSMART_UPDATE_SNAPSHOTS") == "1":
        TS_SNAPSHOT.write_text(rendered)

    assert rendered == TS_SNAPSHOT.read_text()
