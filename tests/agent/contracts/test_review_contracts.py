"""Sanity checks for the frozen pre-refactor agent contracts."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from chemsmart.agent.v8_kind_index import KIND_SETTINGS


REPO_ROOT = Path(__file__).resolve().parents[3]
CONTRACT_DIR = REPO_ROOT / "tests/agent/contracts"


def _load(name: str) -> dict[str, object]:
    return json.loads((CONTRACT_DIR / name).read_text(encoding="utf-8"))


def test_full26_fixture_has_every_canonical_kind_once() -> None:
    fixture = _load("full26_compact_specs.json")
    cases = fixture["cases"]
    assert isinstance(cases, list)
    kinds = [case["kind"] for case in cases]
    assert len(kinds) == 26
    assert len(set(kinds)) == 26
    assert sorted(kinds) == sorted(KIND_SETTINGS)


def test_frozen_contract_has_complete_cli_and_agent_inventory() -> None:
    contract = _load("agent_contract_baseline.json")
    cli_help = contract["cli_help"]
    tools = contract["tools"]
    full26 = contract["full26"]
    assert isinstance(cli_help, dict)
    assert len(cli_help) == 30
    assert {row["exit_code"] for row in cli_help.values()} == {0}
    assert isinstance(tools, list)
    tool_names = [row["name"] for row in tools]
    assert len(tool_names) == len(set(tool_names)) == 30
    assert isinstance(full26, list)
    assert len(full26) == 26
    assert all(row["adapter_valid"] is True for row in full26)


def test_frozen_contract_records_semantic_and_intent_baselines() -> None:
    full26 = _load("agent_contract_baseline.json")["full26"]
    assert isinstance(full26, list)
    semantic = Counter(row["semantic_verdict"] for row in full26)
    intent = Counter(row["intent"]["verdict"] for row in full26)
    assert semantic == {"ok": 18, "reject": 8}
    assert intent == {"ok": 21, "reject": 5}
    assert all(row["intent"]["assertions"] for row in full26)


def test_frozen_contract_contains_no_machine_specific_or_secret_text() -> None:
    text = (CONTRACT_DIR / "agent_contract_baseline.json").read_text(
        encoding="utf-8"
    )
    forbidden = (
        "/Users/",
        "/private/tmp/",
        "/var/folders/",
        "hf_",
        "sk-",
    )
    assert not any(value in text for value in forbidden)


def test_architecture_baseline_records_all_completion_deficits() -> None:
    path = REPO_ROOT / "docs/review/agent-architecture-baseline.json"
    report = json.loads(path.read_text(encoding="utf-8"))
    summary = report["summary"]
    assert summary == {
        "classes": 202,
        "classes_over_500_lines": 4,
        "files_over_800_lines": 14,
        "functions": 1694,
        "functions_over_100_lines": 39,
        "import_cycles": 1,
        "physical_lines": 44691,
        "pylint_duplicate_groups": 8,
        "python_files": 151,
        "ruff_c901_violations": 65,
        "ruff_default_violations": 0,
    }
    assert report["coverage"]["percent_covered"] == 74.63279642375647
