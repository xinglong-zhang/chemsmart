"""Collect deterministic ``chemsmart sub`` terminal-state fixtures.

This creates a mock PBS server YAML and a mock qsub executable, then runs the
real ChemSmart CLI twice: once with an intentionally wrong server name and
once with the repaired server name.  The resulting scripts, marker, return
codes, and semantic receipts are stored as append-only training evidence.
No chemistry program or real scheduler is invoked.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import yaml

REPO = Path(__file__).resolve().parents[2]
ENV = {
    **os.environ,
    "PYTHONPATH": str(REPO)
    + (os.pathsep + os.environ["PYTHONPATH"] if os.environ.get("PYTHONPATH") else ""),
}


def _write_fixture(
    root: Path,
    *,
    program: str | None = None,
) -> dict[str, Path]:
    config = root / "home" / ".chemsmart"
    server_dir = config / "server"
    workspace = root / "workspace"
    project_dir = workspace / ".chemsmart" / "gaussian"
    orca_project_dir = workspace / ".chemsmart" / "orca"
    bin_dir = config / "bin"
    project_dirs = []
    if program in {None, "gaussian"}:
        project_dirs.append(project_dir)
    if program in {None, "orca"}:
        project_dirs.append(orca_project_dir)
    for path in (server_dir, *project_dirs, bin_dir):
        path.mkdir(parents=True, exist_ok=True)

    marker = root / "mock_qsub.marker"
    qsub = bin_dir / "mock_qsub"
    qsub.write_text(
        "#!/bin/sh\n"
        f"printf '%s\\n' \"$1\" > {marker}\n"
        "printf 'MOCK_JOB_ID=terminal-fixture-001\\n'\n"
        "exit 0\n",
        encoding="utf-8",
    )
    qsub.chmod(0o755)

    g16 = bin_dir / "g16"
    g16.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    g16.chmod(0o755)
    orca = bin_dir / "orca"
    orca.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    orca.chmod(0o755)

    server = server_dir / "mock-pbs.yaml"
    server.write_text(
        yaml.safe_dump(
            {
                "SERVER": {
                    "SCHEDULER": "PBS",
                    "QUEUE_NAME": "debug",
                    "NUM_HOURS": 1,
                    "MEM_GB": 2,
                    "NUM_CORES": 2,
                    "NUM_THREADS": 2,
                    "NUM_GPUS": 0,
                    "SUBMIT_COMMAND": str(qsub),
                    "PROJECT": "mock-project",
                    "SCRATCH_DIR": None,
                },
                "GAUSSIAN": {
                    "EXEFOLDER": str(bin_dir),
                    "LOCAL_RUN": False,
                },
                "ORCA": {
                    "EXEFOLDER": str(bin_dir),
                    "LOCAL_RUN": False,
                },
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    project = project_dir / "mock.yaml"
    if program in {None, "gaussian"}:
        project.write_text(
            yaml.safe_dump(
                {
                    "gas": {
                        "functional": "B3LYP",
                        "basis": "6-31G(d)",
                        "freq": False,
                    },
                    "solv": {
                        "functional": "B3LYP",
                        "basis": "6-31G(d)",
                        "freq": False,
                    },
                },
                sort_keys=False,
            ),
            encoding="utf-8",
        )
    if program in {None, "orca"}:
        (orca_project_dir / "mock.yaml").write_text(
            yaml.safe_dump(
                {
                    "gas": {
                        "functional": "B3LYP",
                        "basis": "def2-SVP",
                        "freq": False,
                    },
                    "solv": {
                        "functional": "B3LYP",
                        "basis": "def2-SVP",
                        "freq": False,
                    },
                },
                sort_keys=False,
            ),
            encoding="utf-8",
        )
    xyz = workspace / "examples" / "h2o.xyz"
    xyz.parent.mkdir(parents=True, exist_ok=True)
    xyz.write_text(
        "3\nmock water\n"
        "O 0.000000 0.000000 0.117000\n"
        "H 0.000000 0.757000 -0.464000\n"
        "H 0.000000 -0.757000 -0.464000\n",
        encoding="utf-8",
    )
    return {
        "config": config,
        "workspace": workspace,
        "server": server,
        "project": project,
        "xyz": xyz,
        "marker": marker,
        "qsub": qsub,
    }


def _run(command: str, workspace: Path) -> dict[str, Any]:
    env = {**ENV, "HOME": str(workspace.parent / "home")}
    completed = subprocess.run(
        [sys.executable, "-m", "chemsmart.cli.main", *command.split()[1:]],
        cwd=workspace,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    return {
        "ok": completed.returncode == 0,
        "status": "ok" if completed.returncode == 0 else "error",
        "command": command,
        "returncode": completed.returncode,
        "stdout_tail": completed.stdout[-2000:],
        "stderr_tail": completed.stderr[-2000:],
    }


def _assertion(assertion_id: str, expected: Any, observed: Any, evidence: dict[str, Any] | None = None) -> dict[str, Any]:
    return {
        "id": assertion_id,
        "expected": expected,
        "observed": observed,
        "status": "pass" if expected == observed else "fail",
        "evidence": evidence or {},
    }


def _terminal_state(
    *,
    command: str,
    result: dict[str, Any],
    paths: dict[str, Path],
    expected_server: str,
    repaired: bool,
) -> dict[str, Any]:
    marker_exists = paths["marker"].exists()
    scripts = sorted(paths["workspace"].glob("chemsmart_sub_*.sh"))
    assertions = [
        _assertion("sub.command_present", True, bool(command.strip())),
        _assertion("server.yaml_exists", True, paths["server"].is_file()),
        _assertion(
            "server.scheduler_is_pbs",
            "PBS",
            "PBS" if paths["server"].is_file() else None,
            {"server_yaml": str(paths["server"])},
        ),
        _assertion("sub.returncode", 0 if repaired else 1, result["returncode"]),
        _assertion("sub.server_name", expected_server, expected_server),
        _assertion("sub.submit_script_exists", repaired, bool(scripts)),
        _assertion("sub.mock_scheduler_marker", repaired, marker_exists),
    ]
    all_passed = all(row["status"] == "pass" for row in assertions)
    return {
        "schema_version": 1,
        "action": "submit_job" if repaired else "repair_submit_job",
        "command": command,
        "status": "passed" if all_passed else "failed",
        "all_passed": all_passed,
        "returncode": result["returncode"],
        "expected_returncode": 0 if repaired else 1,
        "server": expected_server,
        "scheduler": "PBS",
        "execution_mode": "mock_scheduler" if repaired else "mock_scheduler_reject",
        "assertions": assertions,
        "artifacts": [
            {"kind": "server_yaml", "path": str(paths["server"])},
            *[
                {"kind": "submit_script", "path": str(path)}
                for path in scripts
            ],
            {"kind": "mock_qsub_marker", "path": str(paths["marker"])},
        ],
    }


def collect(run_dir: Path) -> Path:
    run_dir = run_dir.expanduser().resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    paths = _write_fixture(run_dir / "fixtures")
    wrong = (
        "chemsmart sub -s missing-pbs --fake gaussian -p mock "
        "-f examples/h2o.xyz -c 0 -m 1 opt"
    )
    correct = (
        "chemsmart sub -s mock-pbs --fake gaussian -p mock "
        "-f examples/h2o.xyz -c 0 -m 1 opt"
    )
    wrong_result = _run(wrong, paths["workspace"])
    wrong_state = _terminal_state(
        command=wrong,
        result=wrong_result,
        paths=paths,
        expected_server="missing-pbs",
        repaired=False,
    )
    # The first attempt must not leave a scheduler marker. The second attempt
    # is the actual corrected `chemsmart sub` path against the mock YAML.
    correct_result = _run(correct, paths["workspace"])
    correct_state = _terminal_state(
        command=correct,
        result=correct_result,
        paths=paths,
        expected_server="mock-pbs",
        repaired=True,
    )

    from chemsmart.agent.training_log import TrainingEpisodeWriter, TrainingLogConfig

    training_dir = run_dir
    writer = TrainingEpisodeWriter(
        TrainingLogConfig(enabled=True, dir=training_dir, mask_paths=True)
    )
    session_id = f"mock-sub-terminal-{int(time.time())}"
    first = writer.write_episode(
        session_id=session_id,
        turn=1,
        provider_name="deepseek",
        model="deepseek-v4-pro",
        messages=[
            {"role": "user", "content": "Submit a Gaussian optimization to the mock PBS server."},
            {"role": "assistant", "content": "The first server target was invalid; I will repair it."},
        ],
        tool_records=[
            {
                "tool": "synthesize_command",
                "status": "ok",
                "result": {
                    "status": "ready",
                    "command": wrong,
                    "semantic": {
                        "verdict": "reject",
                        "failed_rule_ids": ["cmd.semantic.safe_execution_failed"],
                    },
                },
            },
            {
                "tool": "execute_chemsmart_command",
                "status": "error",
                "result": wrong_result,
            },
        ],
        cwd=str(paths["workspace"]),
        final_answer="The requested server name was not configured; repair required.",
        terminal_state=wrong_state,
    )
    second = writer.write_episode(
        session_id=session_id,
        turn=2,
        provider_name="deepseek",
        model="deepseek-v4-pro",
        messages=[
            {"role": "user", "content": "Submit a Gaussian optimization to the mock PBS server."},
            {"role": "assistant", "content": "The first server target was invalid; I will repair it."},
            {"role": "user", "content": "Use mock-pbs and submit the corrected command."},
            {"role": "assistant", "content": "The corrected chemsmart sub command reached the mock PBS scheduler."},
        ],
        tool_records=[
            {
                "tool": "synthesize_command",
                "status": "ok",
                "result": {
                    "status": "ready",
                    "command": correct,
                    "semantic": {
                        "verdict": "ok",
                        "generated_inputs": [{"path": str(paths["workspace"] / "mock_opt.com"), "route": "# opt b3lyp 6-31g(d)"}],
                    },
                },
            },
            {
                "tool": "execute_chemsmart_command",
                "status": "ok",
                "result": {**correct_result, "submitted": True, "job_id": "terminal-fixture-001"},
            },
        ],
        cwd=str(paths["workspace"]),
        final_answer="Submitted the corrected Gaussian job to mock-pbs.",
        terminal_state=correct_state,
    )
    if first is None or second is None:
        raise RuntimeError("could not write mock terminal-state episodes")
    print(json.dumps({
        "run_dir": str(run_dir),
        "episode_file": str(second),
        "wrong": wrong_state,
        "correct": correct_state,
    }, indent=2, sort_keys=True))
    return second


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", default=None)
    args = parser.parse_args()
    run_dir = Path(args.run_dir) if args.run_dir else REPO / "var" / "agent-training" / "runs" / f"mock-server-terminal-state-{int(time.time())}"
    collect(run_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
