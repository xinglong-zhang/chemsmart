from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from chemsmart.agent.cli import agent

from ._agent_session_helpers import FakeProvider, critic_ok, planner_plan


def test_agent_cli_run_prints_plan_and_writes_decision_log(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    import chemsmart.agent.tools as agent_tools

    provider = FakeProvider(
        [planner_plan(single_molecule_xyz_file, "cli_case"), critic_ok()]
    )

    def fake_get_provider():
        return provider

    def fake_validate_runtime(job, server=None):
        return {
            "ok": "ok",
            "local_ok": True,
            "local_issues": [],
            "remote_unknown": [],
        }

    monkeypatch.setattr(
        "chemsmart.agent.providers.get_provider", fake_get_provider
    )
    monkeypatch.setattr("chemsmart.agent.core.get_provider", fake_get_provider)
    monkeypatch.setattr(agent_tools, "validate_runtime", fake_validate_runtime)
    monkeypatch.setattr(
        "chemsmart.agent.core._default_session_root",
        lambda: str(tmp_path / "sessions"),
    )

    runner = CliRunner()
    result = runner.invoke(
        agent,
        ["run", "--dry-submit", f"optimize {single_molecule_xyz_file}"],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    assert "Plan:" in result.output
    assert "inputfile:" in result.output
    assert "decision log:" in result.output

    session_dir = tmp_path / "sessions"
    decision_logs = list(session_dir.rglob("decision_log.jsonl"))
    assert len(decision_logs) == 1
    assert Path(decision_logs[0]).read_text()


def test_agent_cli_tools_lists_registered_tools():
    runner = CliRunner()
    result = runner.invoke(agent, ["tools"], catch_exceptions=False)

    assert result.exit_code == 0, result.output
    assert "registered tools:" in result.output
    assert "build_molecule" in result.output
