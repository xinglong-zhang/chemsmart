from __future__ import annotations

import subprocess
import sys

from chemsmart.agent.harness.command_semantics import (
    evaluate_command_semantics,
)


def test_run_command_semantic_gate_uses_safe_fake_execution(
    monkeypatch,
    tmp_path,
) -> None:
    def fake_run(argv, cwd, env, text, stdout, stderr, timeout, check):
        assert text is True
        assert stdout is subprocess.PIPE
        assert stderr is subprocess.PIPE
        assert check is False
        assert timeout == 30.0
        assert "PYTHONPATH" in env
        assert argv[:5] == [
            sys.executable,
            "-m",
            "chemsmart.cli.main",
            "--no-verbose",
            "run",
        ]
        assert "--fake" in argv
        assert "--no-scratch" in argv
        (tmp_path / "water_orca_opt.inp").write_text(
            "! Opt B3LYP def2-SVP\n"
            "* xyz 0 1\n"
            "O 0 0 0\n"
            "H 0 0 1\n"
            "H 0 1 0\n"
            "*\n",
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(argv, 0, "", "")

    monkeypatch.setattr(
        "chemsmart.agent.harness.command_semantics.subprocess.run",
        fake_run,
    )

    result = evaluate_command_semantics(
        "chemsmart run orca -f water.xyz -c 0 -m 1 -l water_orca_opt opt",
        cwd=tmp_path,
    )

    assert result.verdict == "ok"
    assert result.generated_inputs[0]["route"] == "! Opt B3LYP def2-SVP"


def test_safe_execution_failure_reports_missing_runtime_info(
    monkeypatch,
    tmp_path,
) -> None:
    def fake_run(argv, **_kwargs):
        return subprocess.CompletedProcess(
            argv,
            1,
            "",
            "ValueError: No server implemented for local.yaml.\n"
            "Currently available servers: ['colab_orca']\n",
        )

    monkeypatch.setattr(
        "chemsmart.agent.harness.command_semantics.subprocess.run",
        fake_run,
    )

    result = evaluate_command_semantics(
        "chemsmart run orca -f water.xyz -c 0 -m 1 opt",
        cwd=tmp_path,
    )

    assert result.verdict == "reject"
    assert result.failed_rule_ids == ["cmd.semantic.safe_execution_failed"]
    assert "valid chemsmart server configuration" in result.missing_info
    assert "available servers: ['colab_orca']" in result.missing_info


def test_safe_execution_failure_reports_missing_project_info(
    monkeypatch,
    tmp_path,
) -> None:
    def fake_run(argv, **_kwargs):
        return subprocess.CompletedProcess(
            argv,
            1,
            "",
            "FileNotFoundError: No project settings implemented for test.\n"
            "Currently available projects: []\n",
        )

    monkeypatch.setattr(
        "chemsmart.agent.harness.command_semantics.subprocess.run",
        fake_run,
    )

    result = evaluate_command_semantics(
        "chemsmart run gaussian -p test -f water.xyz -c 0 -m 1 opt",
        cwd=tmp_path,
    )

    assert result.verdict == "reject"
    assert "valid chemsmart project configuration" in result.missing_info
    assert "available projects: []" in result.missing_info


def test_wrong_option_order_is_rejected_before_safe_execution(
    monkeypatch,
    tmp_path,
) -> None:
    def should_not_run(*_args, **_kwargs):  # pragma: no cover - defensive
        raise AssertionError("semantic gate should reject before subprocess")

    monkeypatch.setattr(
        "chemsmart.agent.harness.command_semantics.subprocess.run",
        should_not_run,
    )

    result = evaluate_command_semantics(
        "chemsmart run orca opt -f water.xyz -c 0 -m 1",
        cwd=tmp_path,
    )

    assert result.verdict == "reject"
    assert result.failed_rule_ids == ["cmd.semantic.option_order"]


def test_submit_success_without_observed_input_warns(
    monkeypatch,
    tmp_path,
) -> None:
    def fake_run(argv, **_kwargs):
        assert argv[:5] == [
            sys.executable,
            "-m",
            "chemsmart.cli.main",
            "--no-verbose",
            "sub",
        ]
        assert "--test" in argv
        assert "--fake" in argv
        return subprocess.CompletedProcess(argv, 0, "", "")

    monkeypatch.setattr(
        "chemsmart.agent.harness.command_semantics.subprocess.run",
        fake_run,
    )

    result = evaluate_command_semantics(
        "chemsmart sub -s cluster gaussian -p test -f water.xyz -c 0 -m 1 opt",
        cwd=tmp_path,
    )

    assert result.verdict == "warn"
    assert result.failed_rule_ids == [
        "cmd.semantic.submit_generated_input_not_observed"
    ]
