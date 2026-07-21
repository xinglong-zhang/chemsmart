#!/usr/bin/env python3
"""Freeze normalized ChemSmart agent contracts before structural refactors."""

from __future__ import annotations

import argparse
import ast
import difflib
import hashlib
import importlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from dataclasses import asdict, is_dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Iterator

ANSI_RE = re.compile(r"\x1b\[[0-9;?]*[ -/]*[@-~]")
LOG_TIMESTAMP_RE = re.compile(r"\b\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3}\b")
RULE_ID_RE = re.compile(
    r"^(?:cmd|gaussian|harness|input|intent|orca|project|runtime|spec|"
    r"terminal|workflow|yaml)\.[a-z0-9_.-]+$"
)

CLI_HELP_PATHS = (
    ("--help",),
    ("agent", "--help"),
    ("run", "--help"),
    ("sub", "--help"),
    ("run", "gaussian", "--help"),
    ("run", "orca", "--help"),
    (
        "run",
        "gaussian",
        "-p",
        "test",
        "-f",
        "examples/h2o.xyz",
        "opt",
        "--help",
    ),
    (
        "run",
        "gaussian",
        "-p",
        "test",
        "-f",
        "examples/h2o.xyz",
        "sp",
        "--help",
    ),
    (
        "run",
        "gaussian",
        "-p",
        "test",
        "-f",
        "examples/h2o.xyz",
        "ts",
        "--help",
    ),
    (
        "run",
        "gaussian",
        "-p",
        "test",
        "-f",
        "examples/h2o.xyz",
        "irc",
        "--help",
    ),
    (
        "run",
        "gaussian",
        "-p",
        "test",
        "-f",
        "examples/h2o.xyz",
        "scan",
        "--help",
    ),
    (
        "run",
        "gaussian",
        "-p",
        "test",
        "-f",
        "examples/h2o.xyz",
        "modred",
        "--help",
    ),
    (
        "run",
        "gaussian",
        "-p",
        "test",
        "-f",
        "examples/h2o.xyz",
        "nci",
        "--help",
    ),
    (
        "run",
        "gaussian",
        "-p",
        "test",
        "-f",
        "examples/h2o.xyz",
        "resp",
        "--help",
    ),
    (
        "run",
        "gaussian",
        "-p",
        "test",
        "-f",
        "examples/h2o.xyz",
        "td",
        "--help",
    ),
    (
        "run",
        "gaussian",
        "-p",
        "test",
        "-f",
        "examples/h2o.xyz",
        "dias",
        "--help",
    ),
    (
        "run",
        "gaussian",
        "-p",
        "test",
        "-f",
        "examples/h2o.xyz",
        "crest",
        "--help",
    ),
    (
        "run",
        "gaussian",
        "-p",
        "test",
        "-f",
        "examples/h2o.xyz",
        "traj",
        "--help",
    ),
    (
        "run",
        "gaussian",
        "-p",
        "test",
        "-f",
        "examples/h2o.xyz",
        "wbi",
        "--help",
    ),
    (
        "run",
        "gaussian",
        "-p",
        "test",
        "-f",
        "examples/h2o.xyz",
        "qrc",
        "--help",
    ),
    (
        "run",
        "gaussian",
        "-p",
        "test",
        "-f",
        "examples/h2o.xyz",
        "opt",
        "qmmm",
        "--help",
    ),
    ("run", "orca", "-p", "test", "-f", "examples/h2o.xyz", "sp", "--help"),
    ("run", "orca", "-p", "test", "-f", "examples/h2o.xyz", "opt", "--help"),
    ("run", "orca", "-p", "test", "-f", "examples/h2o.xyz", "ts", "--help"),
    ("run", "orca", "-p", "test", "-f", "examples/h2o.xyz", "irc", "--help"),
    ("run", "orca", "-p", "test", "-f", "examples/h2o.xyz", "scan", "--help"),
    (
        "run",
        "orca",
        "-p",
        "test",
        "-f",
        "examples/h2o.xyz",
        "modred",
        "--help",
    ),
    ("run", "orca", "-p", "test", "-f", "examples/h2o.xyz", "neb", "--help"),
    ("run", "orca", "-p", "test", "-f", "examples/h2o.xyz", "qrc", "--help"),
    (
        "run",
        "orca",
        "-p",
        "test",
        "-f",
        "examples/h2o.xyz",
        "opt",
        "qmmm",
        "--help",
    ),
)

SESSION_ARTIFACTS = (
    "decision_log.jsonl",
    "handles.jsonl",
    "harness_result.json",
    "runtime_state.json",
    "session.json",
    "session_metadata.json",
    "state.json",
    "turn_{turn:02d}_step_{step:02d}.json",
)


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def normalize_text(text: str, replacements: dict[str, str]) -> str:
    normalized = ANSI_RE.sub("", str(text)).replace("\r\n", "\n")
    normalized = LOG_TIMESTAMP_RE.sub("<TIMESTAMP>", normalized)
    for source, target in sorted(
        replacements.items(), key=lambda item: len(item[0]), reverse=True
    ):
        if source:
            normalized = normalized.replace(source, target)
    normalized = re.sub(
        r"<FIXTURE>/tmp/[^/\\\s\"']+",
        "<FIXTURE>/tmp/<RUN>",
        normalized,
    )
    return "\n".join(line.rstrip() for line in normalized.splitlines()).strip()


def normalize_value(value: Any, replacements: dict[str, str]) -> Any:
    if is_dataclass(value) and not isinstance(value, type):
        return normalize_value(asdict(value), replacements)
    if isinstance(value, Enum):
        return normalize_value(value.value, replacements)
    if isinstance(value, dict):
        return {
            str(key): normalize_value(item, replacements)
            for key, item in sorted(
                value.items(), key=lambda pair: str(pair[0])
            )
        }
    if isinstance(value, (list, tuple)):
        return [normalize_value(item, replacements) for item in value]
    if isinstance(value, (set, frozenset)):
        return [
            normalize_value(item, replacements)
            for item in sorted(value, key=repr)
        ]
    if isinstance(value, Path):
        return normalize_text(str(value), replacements)
    if isinstance(value, str):
        return normalize_text(value, replacements)
    return value


@contextmanager
def isolated_environment(
    repo_root: Path,
) -> Iterator[tuple[Path, dict[str, str]]]:
    with tempfile.TemporaryDirectory(prefix="chemsmart-contract-") as raw:
        root = Path(raw)
        home = root / "home"
        workspace = root / "workspace"
        runtime_tmp = root / "tmp"
        home.mkdir()
        workspace.mkdir()
        runtime_tmp.mkdir()
        old_home = os.environ.get("HOME")
        old_tmpdir = os.environ.get("TMPDIR")
        old_pythonpath = os.environ.get("PYTHONPATH")
        old_tempfile_dir = tempfile.tempdir
        os.environ["HOME"] = str(home)
        os.environ["TMPDIR"] = str(runtime_tmp)
        os.environ["PYTHONPATH"] = str(repo_root)
        tempfile.tempdir = str(runtime_tmp)
        templates = repo_root / "chemsmart/settings/templates/.chemsmart"
        for source, destination in (
            (
                templates / "server/local.yaml",
                home / ".chemsmart/server/local.yaml",
            ),
            (
                templates / "gaussian/test.yaml",
                home / ".chemsmart/gaussian/test.yaml",
            ),
            (templates / "orca/test.yaml", home / ".chemsmart/orca/test.yaml"),
        ):
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
        replacements = {
            str(repo_root): "<REPO>",
            str(root): "<FIXTURE>",
            str(Path.home()): "<HOME>",
            sys.executable: "<PYTHON>",
        }
        try:
            yield workspace, replacements
        finally:
            tempfile.tempdir = old_tempfile_dir
            if old_home is None:
                os.environ.pop("HOME", None)
            else:
                os.environ["HOME"] = old_home
            if old_tmpdir is None:
                os.environ.pop("TMPDIR", None)
            else:
                os.environ["TMPDIR"] = old_tmpdir
            if old_pythonpath is None:
                os.environ.pop("PYTHONPATH", None)
            else:
                os.environ["PYTHONPATH"] = old_pythonpath


def cli_help_contract(
    repo_root: Path, replacements: dict[str, str]
) -> dict[str, Any]:
    environment = dict(os.environ)
    environment.update(
        {
            "COLUMNS": "120",
            "NO_COLOR": "1",
            "TERM": "dumb",
            "PYTHONPATH": str(repo_root),
        }
    )
    output: dict[str, Any] = {}
    for path in CLI_HELP_PATHS:
        result = subprocess.run(
            [sys.executable, "-m", "chemsmart.cli.main", *path],
            cwd=repo_root,
            env=environment,
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
        text = normalize_text(result.stdout + result.stderr, replacements)
        output[" ".join(path)] = {
            "exit_code": result.returncode,
            "sha256": sha256_text(text),
            "text": text,
        }
    return output


def tui_contract(repo_root: Path) -> dict[str, Any]:
    from chemsmart.agent.tui.bindings import BINDINGS
    from chemsmart.agent.tui.config import DEFAULT_KEYBINDINGS
    from chemsmart.agent.tui.slash_catalog import SLASH_PALETTE_COMMANDS

    safety = [
        {
            "key": binding.key,
            "action": binding.action,
            "description": str(binding.description),
            "priority": bool(binding.priority),
        }
        for binding in BINDINGS
    ]
    return {
        "slash_commands": [
            {"command": command, "description": description}
            for command, description in SLASH_PALETTE_COMMANDS
        ],
        "configurable_keybindings": dict(sorted(DEFAULT_KEYBINDINGS.items())),
        "safety_bindings": safety,
    }


def tool_registry_contract() -> list[dict[str, Any]]:
    from chemsmart.agent.registry import ToolRegistry

    rows: list[dict[str, Any]] = []
    for tool in ToolRegistry.default().list_tools():
        metadata = asdict(tool.metadata) if is_dataclass(tool.metadata) else {}
        rows.append(
            {
                "name": tool.name,
                "description": tool.openai_tool_def()["function"][
                    "description"
                ],
                "parameters": tool.openai_tool_def()["function"]["parameters"],
                "metadata": metadata,
            }
        )
    return rows


def public_import_contract() -> list[dict[str, str]]:
    package = importlib.import_module("chemsmart.agent")
    rows: list[dict[str, str]] = []
    for name in package.__all__:
        value = getattr(package, name)
        rows.append(
            {
                "name": name,
                "module": getattr(value, "__module__", type(value).__module__),
                "qualified_name": getattr(
                    value, "__qualname__", type(value).__qualname__
                ),
                "kind": type(value).__name__,
            }
        )
    return rows


def prompt_contract(repo_root: Path) -> dict[str, Any]:
    from chemsmart.agent.cli_schema import build_chemsmart_cli_schema
    from chemsmart.agent.prompts import load_prompt
    from chemsmart.agent.prompts.identity import build_system_prompt
    from chemsmart.agent.prompts.synthesis import build_synthesis_system_prompt
    from chemsmart.agent.registry import ToolRegistry
    from chemsmart.agent.schema_prune import prune_schema_for_request

    prompt_dir = repo_root / "chemsmart/agent/prompts"
    files: dict[str, Any] = {}
    for path in sorted(prompt_dir.iterdir()):
        if path.suffix not in {".md", ".py"}:
            continue
        text = path.read_text(encoding="utf-8")
        files[path.name] = {"chars": len(text), "sha256": sha256_text(text)}

    registry = ToolRegistry.default()
    stage = load_prompt("unified_agent.md")
    requests = {
        "command": "Prepare an ORCA transition-state search for ts.xyz.",
        "project_yaml": "Build a mixed-basis Gaussian project YAML.",
        "diagnostics": "Inspect job 123 on server cluster-a.",
    }
    rendered: dict[str, Any] = {}
    cli_schema = build_chemsmart_cli_schema()
    synthesis: dict[str, Any] = {}
    for label, request in requests.items():
        system = build_system_prompt(
            registry=registry,
            stage_instructions=stage,
            session_meta={"cwd": "<WORKSPACE>", "project": "demo"},
            conversation_context=None,
            request=request,
            max_chars=4096,
        )
        rendered[label] = {
            "chars": len(system),
            "sha256": sha256_text(system),
            "text": system,
        }
        pruned = prune_schema_for_request(
            cli_schema,
            request,
            workspace_program="gaussian",
        )
        synthesis_prompt = build_synthesis_system_prompt(pruned)
        synthesis[label] = {
            "chars": len(synthesis_prompt),
            "sha256": sha256_text(synthesis_prompt),
            "text": synthesis_prompt,
        }
    cli_body = json.dumps(cli_schema, sort_keys=True, separators=(",", ":"))
    return {
        "files": files,
        "rendered_system_prompts": rendered,
        "rendered_synthesis_prompts": synthesis,
        "full_cli_schema_sha256": sha256_text(cli_body),
    }


def workspace_yaml_contract(
    workspace: Path, replacements: dict[str, str]
) -> dict[str, Any]:
    from chemsmart.settings.workspace_project import resolve_workspace_project

    def project_status(selected_path: Path | None = None) -> dict[str, Any]:
        status = resolve_workspace_project(
            cwd=workspace, selected_path=selected_path
        )
        return normalize_value(asdict(status), replacements)

    statuses = {"empty": project_status()}
    gaussian = workspace / ".chemsmart/gaussian/demo.yaml"
    gaussian.parent.mkdir(parents=True)
    gaussian.write_text(
        "gas:\n  functional: b3lyp\n  basis: def2svp\n", encoding="utf-8"
    )
    statuses["single"] = project_status()
    orca = workspace / ".chemsmart/orca/orca_demo.yaml"
    orca.parent.mkdir(parents=True)
    orca.write_text(
        "gas:\n  functional: pbe0\n  basis: def2-svp\n", encoding="utf-8"
    )
    statuses["multiple_unselected"] = project_status()
    statuses["multiple_selected"] = project_status(orca)
    return statuses


def session_artifact_contract(repo_root: Path) -> dict[str, Any]:
    from chemsmart.agent.core import SessionState
    from chemsmart.agent.services.session_finalizer import (
        SESSION_METADATA_KEYS,
    )

    return {
        "artifact_names": list(SESSION_ARTIFACTS),
        "metadata_keys": list(SESSION_METADATA_KEYS),
        "session_state_schema": SessionState.model_json_schema(),
    }


def rule_id_contract(repo_root: Path) -> list[str]:
    values: set[str] = set()
    for path in sorted((repo_root / "chemsmart/agent").rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                candidate = node.value.strip()
                if RULE_ID_RE.fullmatch(candidate):
                    values.add(candidate)
    return sorted(values)


def _write_xyz(path: Path, *, frames: int = 1) -> None:
    frame = "4\ncontract H4\nH 0 0 0\nH 0 0 0.74\nH 0 1.5 0\nH 0 1.5 0.74\n"
    path.write_text(frame * frames, encoding="utf-8")


def _prepare_full26_workspace(repo_root: Path, workspace: Path) -> None:
    _write_xyz(workspace / "molecule.xyz")
    _write_xyz(workspace / "product.xyz")
    _write_xyz(workspace / "trajectory.xyz", frames=3)
    shutil.copy2(
        repo_root
        / "tests/data/GaussianTests/outputs/link/oxygen_openshell_singlet_ts_link.log",
        workspace / "gaussian_ts.log",
    )
    shutil.copy2(
        repo_root / "tests/data/ORCATests/outputs/sn2_ts.out",
        workspace / "orca_ts.out",
    )
    shutil.copy2(
        repo_root / "tests/data/ORCATests/outputs/sn2_ts.hess",
        workspace / "sn2_ts.hess",
    )
    templates = repo_root / "chemsmart/settings/templates/.chemsmart"
    for program in ("gaussian", "orca"):
        target = workspace / ".chemsmart" / program
        target.mkdir(parents=True, exist_ok=True)
        shutil.copy2(
            templates / program / "test.yaml", target / "contract.yaml"
        )
        shutil.copy2(
            templates / program / "test_qmmm.yaml",
            target / "contract_qmmm.yaml",
        )
    # xTB has no QMMM leaf, so it only needs the base project.
    xtb_target = workspace / ".chemsmart" / "xtb"
    xtb_target.mkdir(parents=True, exist_ok=True)
    shutil.copy2(templates / "xtb" / "test.yaml", xtb_target / "contract.yaml")
    server_dir = Path.home() / ".chemsmart/server"
    server_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(templates / "server/local.yaml", server_dir / "local.yaml")


def _normalized_generated_input(
    item: dict[str, Any], replacements: dict[str, str]
) -> dict[str, Any]:
    output = normalize_value(item, replacements)
    if isinstance(output, dict) and "path" in output:
        output["path"] = Path(str(output["path"])).name
    content = output.get("content") if isinstance(output, dict) else None
    if isinstance(content, str):
        output["content_sha256"] = sha256_text(content)
    return output


def full26_contract(
    repo_root: Path,
    workspace: Path,
    replacements: dict[str, str],
) -> list[dict[str, Any]]:
    from chemsmart.agent.harness.command_semantics import (
        evaluate_command_semantics,
    )
    from chemsmart.agent.harness.intent import IntentSpec, evaluate_intent
    from chemsmart.agent.v8_adapter import adapt
    from chemsmart.agent.v8_kind_index import KIND_SETTINGS

    fixture_path = (
        repo_root / "tests/agent/contracts/full26_compact_specs.json"
    )
    fixture = json.loads(fixture_path.read_text(encoding="utf-8"))
    cases = fixture["cases"]
    kinds = [case["kind"] for case in cases]
    if sorted(kinds) != sorted(KIND_SETTINGS):
        raise ValueError(
            "full26 fixture must contain every canonical kind exactly once"
        )
    _prepare_full26_workspace(repo_root, workspace)
    rows: list[dict[str, Any]] = []
    for case in cases:
        kind = case["kind"]
        job: dict[str, Any] = {
            "id": 1,
            "kind": kind,
            "file": case.get("input", "molecule.xyz"),
            "charge": 0,
            "mult": 1,
        }
        if case.get("settings"):
            job["settings"] = case["settings"]
        if case.get("product_file"):
            job["product_file"] = case["product_file"]
        spec = {"intent": "workflow", "jobs": [job]}
        adapted = adapt(
            spec,
            default_project=case.get("project", "contract"),
        )
        commands = adapted.get("commands") or []
        command = commands[0] if commands else ""
        semantic = (
            evaluate_command_semantics(command, cwd=workspace)
            if command
            else None
        )
        expected_intent = IntentSpec.from_dict(
            {
                "action": "run",
                "program": kind.split(".", 1)[0],
                "kind": kind,
                "project": case.get("project", "contract"),
                "input_path": case.get("input", "molecule.xyz"),
                "charge": 0,
                "multiplicity": 1,
                "execution_mode": "local",
                "chemistry": case.get("intent_chemistry", {}),
            }
        )
        intent = (
            evaluate_intent(command, expected_intent, cwd=str(workspace))
            if command
            else None
        )
        rows.append(
            {
                "kind": kind,
                "spec": spec,
                "adapter_valid": adapted.get("valid"),
                "adapter_errors": normalize_value(
                    adapted.get("errors") or [], replacements
                ),
                "command": normalize_text(command, replacements),
                "semantic_verdict": (
                    semantic.verdict if semantic else "not_run"
                ),
                "failed_rule_ids": (
                    semantic.failed_rule_ids if semantic else []
                ),
                "issues": normalize_value(
                    (
                        [issue.to_dict() for issue in semantic.issues]
                        if semantic
                        else []
                    ),
                    replacements,
                ),
                "generated_inputs": [
                    _normalized_generated_input(item, replacements)
                    for item in (semantic.generated_inputs if semantic else ())
                ],
                "intent": normalize_value(
                    intent.to_dict() if intent is not None else {},
                    replacements,
                ),
            }
        )
    return rows


def build_snapshot(repo_root: Path) -> dict[str, Any]:
    with isolated_environment(repo_root) as (workspace, replacements):
        high_risk = (
            repo_root / "tests/agent/harness/fixtures/high_risk_matrix.json"
        ).read_text(encoding="utf-8")
        snapshot = {
            "schema_version": 1,
            "cli_help": cli_help_contract(repo_root, replacements),
            "tui": tui_contract(repo_root),
            "tools": tool_registry_contract(),
            "public_imports": public_import_contract(),
            "prompts": prompt_contract(repo_root),
            "workspace_yaml": workspace_yaml_contract(workspace, replacements),
            "session_artifacts": session_artifact_contract(repo_root),
            "semantic_rule_ids": rule_id_contract(repo_root),
            "full26": full26_contract(repo_root, workspace, replacements),
            "high_risk48_fixture_sha256": sha256_text(high_risk),
        }
        return normalize_value(snapshot, replacements)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    destination = parser.add_mutually_exclusive_group(required=True)
    destination.add_argument("--output", type=Path)
    destination.add_argument("--check", type=Path)
    args = parser.parse_args()
    repo_root = args.repo.resolve()
    sys.path.insert(0, str(repo_root))
    snapshot = build_snapshot(repo_root)
    payload = json.dumps(snapshot, indent=2, sort_keys=True) + "\n"
    if args.check is not None:
        expected = args.check.read_text(encoding="utf-8")
        if expected != payload:
            diff = difflib.unified_diff(
                expected.splitlines(),
                payload.splitlines(),
                fromfile=str(args.check),
                tofile="current agent contract",
                lineterm="",
            )
            print("\n".join(diff))
            return 1
        print(f"matched {args.check} ({sha256_text(payload)})")
        return 0
    assert args.output is not None
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload, encoding="utf-8")
    print(f"wrote {args.output} ({sha256_text(payload)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
