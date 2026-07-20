#!/usr/bin/env python3
"""Run one bounded, reproducible batch from the frozen agentic matrix.

The runner is deliberately opt-in for provider calls. ``--live`` uses the
already configured ChemSmart provider, executes only the semantic gate's safe
fake mode, and writes result JSONL for ``score_agentic_matrix.py``. A batch is
limited to four cases so a provider failure or quota change cannot corrupt a
large comparison run.
"""

from __future__ import annotations

import argparse
import contextlib
import json
import os
import tempfile
import time
from pathlib import Path
from typing import Any, Iterator

from chemsmart.agent.harness.evaluation import (
    OutcomeClass,
    classify_agent_result,
    load_case_matrix,
)
from chemsmart.agent.harness.terminal_state import (
    terminal_state_is_positive,
    validate_terminal_state,
)
from chemsmart.agent.harness.workflow_state import (
    reset_workflow_state,
    workflow_state_scope,
)
from chemsmart.agent.providers import get_provider
from chemsmart.agent.synthesis import SynthesisSession
from chemsmart.agent.tools_command import (
    execute_chemsmart_command,
    register_command_intent,
)


def planned_cases(
    matrix: Path,
    *,
    offset: int,
    limit: int,
) -> list[dict[str, Any]]:
    cases = load_case_matrix(matrix)[offset : offset + limit]
    return [
        {
            "case_id": case.case_id,
            "family": case.family,
            "expected_outcome": case.expected_outcome,
            "turns": list(case.turns),
            "intent": case.intent.to_dict(),
        }
        for case in cases
    ]


def run_batch(
    matrix: Path,
    *,
    offset: int,
    limit: int,
    repeats: int,
) -> list[dict[str, Any]]:
    provider = get_provider()
    rows: list[dict[str, Any]] = []
    for trial in range(1, repeats + 1):
        for case in load_case_matrix(matrix)[offset : offset + limit]:
            with tempfile.TemporaryDirectory(
                prefix="chemsmart-agentic-matrix-"
            ) as temp:
                workspace = Path(temp) / "workspace"
                home = Path(temp) / "home"
                _materialize_fixture(
                    workspace, home, case.fixture, case.intent.to_dict()
                )
                with _workspace_environment(workspace, home):
                    reset_workflow_state()
                    scope = f"matrix:{case.case_id}:{trial}"
                    with workflow_state_scope(scope):
                        session = SynthesisSession(
                            provider=provider,
                            default_project="",
                            enable_intent_router=False,
                        )
                        started = time.monotonic()
                        result: dict[str, Any] = {}
                        turn_results: list[dict[str, Any]] = []
                        for turn in case.turns:
                            result = session.prepare_command(turn)
                            turn_results.append(_public_turn_result(result))

                        command = str(result.get("command") or "").strip()
                        semantic = (
                            result.get("semantic")
                            if isinstance(result.get("semantic"), dict)
                            else {}
                        )
                        intent = (
                            result.get("intent_assertion")
                            if isinstance(result.get("intent_assertion"), dict)
                            else {}
                        )
                        execution: dict[str, Any] | None = None
                        if (
                            command
                            and str(result.get("status") or "") == "ready"
                            and semantic.get("verdict") in {"ok", "warn"}
                        ):
                            case_intent = register_command_intent(
                                command,
                                case.intent,
                            )
                            intent = case_intent.to_dict()
                            if case_intent.verdict == "ok":
                                execution = execute_chemsmart_command(
                                    command,
                                    test=True,
                                    timeout_s=120,
                                )
                        elapsed = round(time.monotonic() - started, 4)
            status = str(result.get("status") or "")
            outcome = classify_agent_result(
                status=status,
                expected_outcome=case.expected_outcome,
                semantic_rule_ids=semantic.get("failed_rule_ids") or (),
                intent_rule_ids=intent.get("failed_rule_ids") or (),
                repaired=_has_repair(turn_results),
            )
            terminal = (
                execution.get("terminal_state")
                if isinstance(execution, dict)
                and isinstance(execution.get("terminal_state"), dict)
                else None
            )
            execution_stage = (
                "completed" if execution is not None else "not_reached"
            )
            terminal_rule_ids = (
                validate_terminal_state(terminal)
                if execution is not None
                else []
            )
            if case.expected_outcome == "valid_ask":
                passed = outcome is OutcomeClass.VALID_ASK
            else:
                passed = outcome in {
                    OutcomeClass.DIRECT_PASS,
                    OutcomeClass.REPAIR_PASS,
                } and terminal_state_is_positive(terminal)
            rows.append(
                {
                    "schema_version": 1,
                    "case_id": case.case_id,
                    "family": case.family,
                    "provider": str(getattr(provider, "name", "configured")),
                    "model": str(getattr(provider, "model", "configured")),
                    "trial": trial,
                    "status": status,
                    "outcome": outcome.value,
                    "passed": passed,
                    "repaired": outcome is OutcomeClass.REPAIR_PASS,
                    "semantic_rule_ids": semantic.get("failed_rule_ids") or [],
                    "intent_rule_ids": intent.get("failed_rule_ids") or [],
                    "command": str(result.get("command") or ""),
                    "turn_results": turn_results,
                    "execution_stage": execution_stage,
                    "execution": _public_execution_result(execution),
                    "terminal_state": terminal,
                    "terminal_rule_ids": terminal_rule_ids,
                    "token_estimate": {
                        "input": _estimate_tokens("\n".join(case.turns)),
                        "output": _estimate_tokens(session._last_raw_response),
                        "source": "character_estimate",
                    },
                    "latency_s": elapsed,
                }
            )
    return rows


def _materialize_fixture(
    workspace: Path,
    home: Path,
    fixture: dict[str, Any],
    intent: dict[str, Any],
) -> None:
    workspace.mkdir(parents=True, exist_ok=True)
    (home / ".chemsmart" / "server").mkdir(parents=True, exist_ok=True)
    (home / ".chemsmart" / "server" / "hpc1.yaml").write_text(
        "SERVER:\n  SCHEDULER: SLURM\n  NUM_CORES: 1\n"
        "ORCA:\n  EXEFOLDER: /tmp\n  LOCAL_RUN: false\n  SCRATCH: false\n"
        "GAUSSIAN:\n  EXEFOLDER: /tmp\n  LOCAL_RUN: false\n  SCRATCH: false\n",
        encoding="utf-8",
    )
    program = str(intent.get("program") or "gaussian")
    project_names = set(
        str(name) for name in fixture.get("workspace_projects") or []
    )
    if intent.get("project"):
        project_names.add(str(intent["project"]))
    for name in project_names:
        project_dir = workspace / ".chemsmart" / program
        project_dir.mkdir(parents=True, exist_ok=True)
        (project_dir / f"{name}.yaml").write_text(
            _project_yaml(program),
            encoding="utf-8",
        )
    paths = [fixture.get("input"), fixture.get("endpoint")]
    for path in paths:
        if isinstance(path, str) and path.endswith(".xyz"):
            _write_xyz(workspace / path, int(fixture.get("atom_count") or 12))


def _project_yaml(program: str) -> str:
    base = "gas:\n  functional: b3lyp\n  basis: def2svp\nsolv:\n  functional: b3lyp\n  basis: def2svp\n"
    if program == "gaussian":
        return base + "td:\n  functional: cam-b3lyp\n  basis: def2svp\n"
    return (
        base
        + "qmmm:\n  high_level_functional: b3lyp\n  high_level_basis: def2svp\n  low_level_method: AMBER=HardFirst\n"
    )


def _write_xyz(path: Path, atom_count: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    coordinates = [str(max(atom_count, 3)), "matrix fixture"]
    coordinates.extend(
        f"H {index * 0.7:.3f} 0.0 0.0" for index in range(max(atom_count, 3))
    )
    path.write_text("\n".join(coordinates) + "\n", encoding="utf-8")


def _public_turn_result(result: dict[str, Any]) -> dict[str, Any]:
    semantic = (
        result.get("semantic")
        if isinstance(result.get("semantic"), dict)
        else {}
    )
    intent = (
        result.get("intent_assertion")
        if isinstance(result.get("intent_assertion"), dict)
        else {}
    )
    return {
        "status": result.get("status"),
        "command": result.get("command"),
        "semantic_rule_ids": semantic.get("failed_rule_ids") or [],
        "intent_rule_ids": intent.get("failed_rule_ids") or [],
    }


def _has_repair(turns: list[dict[str, Any]]) -> bool:
    return any(turn.get("status") != "ready" for turn in turns[:-1])


def _public_execution_result(
    result: dict[str, Any] | None,
) -> dict[str, Any] | None:
    if not isinstance(result, dict):
        return None
    return {
        "ok": result.get("ok"),
        "status": result.get("status"),
        "test": result.get("test"),
        "returncode": result.get("returncode"),
        "executed_argv": result.get("executed_argv"),
    }


def _estimate_tokens(text: str) -> int:
    return max(0, round(len(str(text or "")) / 4))


@contextlib.contextmanager
def _workspace_environment(workspace: Path, home: Path) -> Iterator[None]:
    previous_cwd = Path.cwd()
    previous_home = os.environ.get("HOME")
    os.environ["HOME"] = str(home)
    os.chdir(workspace)
    try:
        yield
    finally:
        os.chdir(previous_cwd)
        if previous_home is None:
            os.environ.pop("HOME", None)
        else:
            os.environ["HOME"] = previous_home


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--matrix",
        type=Path,
        default=Path("tests/agent/harness/fixtures/high_risk_matrix.json"),
    )
    parser.add_argument("--offset", type=int, default=0)
    parser.add_argument("--limit", type=int, default=4)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--live", action="store_true")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    if not 1 <= args.limit <= 4:
        parser.error("--limit must be between 1 and 4")
    if args.offset < 0 or args.repeats < 1:
        parser.error(
            "--offset must be non-negative and --repeats must be positive"
        )

    if not args.live:
        rendered = json.dumps(
            {
                "mode": "dry_run",
                "cases": planned_cases(
                    args.matrix, offset=args.offset, limit=args.limit
                ),
            },
            ensure_ascii=False,
            indent=2,
        )
    else:
        rendered = "\n".join(
            json.dumps(row, ensure_ascii=False, sort_keys=True)
            for row in run_batch(
                args.matrix,
                offset=args.offset,
                limit=args.limit,
                repeats=args.repeats,
            )
        )
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(rendered + "\n", encoding="utf-8")
    else:
        print(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
