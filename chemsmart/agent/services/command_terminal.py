"""Deterministic terminal receipts for local and mock-submission execution."""

from __future__ import annotations

import shlex
from pathlib import Path
from typing import Any

from chemsmart.agent.harness.sub_intent import build_sub_intent_assertions
from chemsmart.agent.harness.terminal_state import (
    assertion,
    build_terminal_state,
)
from chemsmart.agent.model_command_parser import parse_model_command

JsonDict = dict[str, Any]


def execution_terminal_state(
    command: str,
    *,
    returncode: int,
    test: bool,
    cwd: Path,
    workflow_scope: str,
    intent_expected: JsonDict | None,
    selected_project: Any = None,
    submit_scripts_before: dict[str, tuple[int, int]] | None = None,
) -> JsonDict:
    """Create observable post-execution assertions for the training ledger."""

    tokens = _command_tokens(command)
    parsed = parse_model_command(command)
    is_submit = len(tokens) > 1 and tokens[:2] == ["chemsmart", "sub"]
    assertions = [
        assertion("command.returncode", expected=0, observed=returncode),
        assertion("command.process_completed", expected=True, observed=True),
    ]
    if not is_submit:
        return build_terminal_state(
            action="execute_command",
            command=command,
            assertions=assertions,
            returncode=returncode,
            server=None,
            scheduler=None,
            execution_mode="local",
            safe_execution_mode="test_fake" if test else "execute",
            required_assertion_ids=[
                "command.returncode",
                "command.process_completed",
            ],
            artifacts=[],
        )
    return _submission_terminal_state(
        command,
        tokens=tokens,
        parsed_project=parsed.project,
        returncode=returncode,
        test=test,
        cwd=cwd,
        workflow_scope=workflow_scope,
        intent_expected=intent_expected,
        selected_project=selected_project,
        submit_scripts_before=submit_scripts_before,
        assertions=assertions,
    )


def _submission_terminal_state(
    command: str,
    *,
    tokens: list[str],
    parsed_project: str | None,
    returncode: int,
    test: bool,
    cwd: Path,
    workflow_scope: str,
    intent_expected: JsonDict | None,
    selected_project: Any,
    submit_scripts_before: dict[str, tuple[int, int]] | None,
    assertions: list[JsonDict],
) -> JsonDict:
    safe_mode = "test_fake" if test else "execute"
    required_ids = [
        "command.returncode",
        "command.process_completed",
        "sub.intent_provenance_present",
    ]
    assertions.append(
        assertion("sub.execution_mode", expected="submit", observed="submit")
    )
    intent_rows = _submission_intent_assertions(
        command,
        cwd=cwd,
        workflow_scope=workflow_scope,
        intent_expected=intent_expected,
    )
    assertions.extend(intent_rows)
    required_ids.extend(str(row["id"]) for row in intent_rows[1:])
    assertions.append(
        assertion(
            "sub.safe_execution_mode",
            expected=safe_mode,
            observed=safe_mode,
        )
    )
    scripts = changed_submit_scripts(cwd, submit_scripts_before or {})
    assertions.append(
        assertion(
            "sub.submit_script_exists",
            expected=True,
            observed=bool(scripts),
            evidence={"cwd": str(cwd)},
        )
    )
    server_name = _server_name(tokens)
    server_path = _server_yaml_path(server_name)
    scheduler = _server_scheduler(server_path)
    marker = {"PBS": "#PBS", "SLURM": "#SBATCH", "SGE": "#$"}.get(
        str(scheduler or "").upper(), ""
    )
    marker_present = bool(marker) and any(
        marker in path.read_text(encoding="utf-8", errors="replace")
        for path in scripts
    )
    assertions.extend(
        _server_assertions(server_path, scheduler, marker, marker_present)
    )
    project_hash = (
        selected_project.sha256
        if selected_project is not None
        and selected_project.name == parsed_project
        else ""
    )
    assertions.append(
        assertion(
            "sub.project_hash_present",
            expected=True,
            observed=bool(project_hash),
            evidence={"project": parsed_project, "sha256": project_hash},
        )
    )
    required_ids.extend(
        [
            "sub.execution_mode",
            "sub.safe_execution_mode",
            "sub.submit_script_exists",
            "sub.server_yaml_exists",
            "sub.server_scheduler_present",
            "sub.scheduler_marker",
            "sub.project_hash_present",
        ]
    )
    artifacts = [
        {"kind": "submit_script", "path": str(path)} for path in scripts
    ]
    if server_path.is_file():
        artifacts.append({"kind": "server_yaml", "path": str(server_path)})
    return build_terminal_state(
        action="submit_job",
        command=command,
        assertions=assertions,
        returncode=returncode,
        server=server_name,
        scheduler=scheduler,
        execution_mode="submit",
        safe_execution_mode=safe_mode,
        required_assertion_ids=required_ids,
        artifacts=artifacts,
    )


def _submission_intent_assertions(
    command: str,
    *,
    cwd: Path,
    workflow_scope: str,
    intent_expected: JsonDict | None,
) -> list[JsonDict]:
    present = intent_expected is not None
    rows = [
        assertion(
            "sub.intent_provenance_present",
            expected=True,
            observed=present,
            evidence={
                "source": "synthesis_intent_contract"
                if present
                else "missing",
                "workflow_scope": workflow_scope,
            },
        )
    ]
    if intent_expected is not None:
        rows.extend(
            build_sub_intent_assertions(
                command,
                intent_expected,
                cwd=str(cwd),
            )
        )
    return rows


def _server_assertions(
    server_path: Path,
    scheduler: str | None,
    marker: str,
    marker_present: bool,
) -> list[JsonDict]:
    return [
        assertion(
            "sub.server_yaml_exists",
            expected=True,
            observed=server_path.is_file(),
            evidence={"server_yaml": str(server_path)},
        ),
        assertion(
            "sub.server_scheduler_present",
            expected=True,
            observed=bool(scheduler),
            evidence={"scheduler": scheduler},
        ),
        assertion(
            "sub.scheduler_marker",
            expected=True,
            observed=marker_present,
            evidence={"marker": marker},
        ),
    ]


def submit_script_fingerprints(cwd: Path) -> dict[str, tuple[int, int]]:
    """Snapshot submit scripts so receipts can identify newly written files."""

    return {
        str(path): fingerprint
        for path in cwd.glob("chemsmart_sub_*.sh")
        if (fingerprint := _file_fingerprint(path)) is not None
    }


def changed_submit_scripts(
    cwd: Path,
    before: dict[str, tuple[int, int]],
) -> list[Path]:
    return [
        path
        for path in sorted(cwd.glob("chemsmart_sub_*.sh"))
        if before.get(str(path)) != _file_fingerprint(path)
    ]


def _command_tokens(command: str) -> list[str]:
    try:
        return shlex.split(command)
    except ValueError:
        return []


def _file_fingerprint(path: Path) -> tuple[int, int] | None:
    try:
        stat = path.stat()
    except OSError:
        return None
    return stat.st_mtime_ns, stat.st_size


def _server_name(tokens: list[str]) -> str | None:
    for index, token in enumerate(tokens[:-1]):
        if token in {"-s", "--server"}:
            value = str(tokens[index + 1]).strip()
            return value or None
    return None


def _server_yaml_path(server: str | None) -> Path:
    name = str(server or "").strip()
    if not name:
        return Path.home() / ".chemsmart" / "server" / "__missing__.yaml"
    suffix = Path(name).suffix.lower()
    filename = name if suffix in {".yaml", ".yml"} else f"{name}.yaml"
    return Path.home() / ".chemsmart" / "server" / filename


def _server_scheduler(path: Path) -> str | None:
    if not path.is_file():
        return None
    try:
        import yaml

        payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        server = payload.get("SERVER") if isinstance(payload, dict) else None
        if isinstance(server, dict):
            value = server.get("SCHEDULER") or server.get("scheduler")
            return str(value).strip() or None
    except Exception:
        return None
    return None


__all__ = ["execution_terminal_state", "submit_script_fingerprints"]
