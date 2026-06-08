import os
import statistics
import subprocess
import time
from pathlib import Path

import pytest
from click.testing import CliRunner

import chemsmart.cli.main as cli_main

BANNER = " ____ _   _ _____ "


def test_agent_help_has_no_banner():
    result = CliRunner().invoke(cli_main.entry_point, ["agent", "--help"])

    assert result.exit_code == 0
    assert BANNER not in result.output


def test_run_help_keeps_banner():
    result = CliRunner().invoke(cli_main.entry_point, ["run", "--help"])

    assert result.exit_code == 0
    assert BANNER in result.output


def test_missing_agent_extras_emits_install_hint(monkeypatch):
    real_load_group = cli_main._load_group

    def fake_load_group(module_path: str, attr_name: str):
        if module_path == "chemsmart.agent.cli":
            raise ImportError("optional deps missing")
        return real_load_group(module_path, attr_name)

    monkeypatch.setattr(cli_main, "_load_group", fake_load_group)
    result = CliRunner(mix_stderr=False).invoke(
        cli_main.entry_point,
        ["agent"],
    )

    assert result.exit_code == 1
    assert 'pip install -e ".[agent-tui]"' in result.stderr


@pytest.mark.slow
def test_agent_help_cold_start_budget():
    if os.environ.get("CI"):
        pytest.xfail("Timing budgets are noisy on shared CI runners.")

    repo_root = Path(__file__).resolve().parents[2]
    command = ["/opt/anaconda3/bin/chemsmart", "agent", "--help"]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(repo_root) + (
        os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else ""
    )
    durations = []
    for _ in range(3):
        started = time.perf_counter()
        completed = subprocess.run(
            command,
            cwd=repo_root,
            env=env,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        durations.append(time.perf_counter() - started)
        assert completed.returncode == 0

    assert statistics.median(durations) <= 0.8, durations
