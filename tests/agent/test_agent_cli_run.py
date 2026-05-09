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


def test_agent_cli_run_wraps_runtime_errors_as_click_exceptions(
    monkeypatch,
):
    def fake_run(self, request, **kwargs):
        raise RuntimeError(
            "run_local failed with returncode 1; see /tmp/run.stderr"
        )

    monkeypatch.setattr("chemsmart.agent.cli.AgentSession.run", fake_run)

    runner = CliRunner()
    result = runner.invoke(
        agent,
        ["run", "--dry-submit", "optimize examples/h2o.xyz"],
        catch_exceptions=False,
    )

    assert result.exit_code == 1
    assert "Error: run_local failed with returncode 1" in result.output
    assert "Traceback" not in result.output


def test_agent_cli_run_surfaces_advisory_only_answers(
    monkeypatch,
    tmp_path: Path,
):
    provider = FakeProvider(
        [
            {
                "steps": [],
                "rationale": "Use M06-2X/def2-SVP for the TS search and refine with def2-TZVP after one imaginary frequency check.",
                "estimated_cost": "low",
            }
        ]
    )

    def fake_get_provider():
        return provider

    monkeypatch.setattr(
        "chemsmart.agent.providers.get_provider", fake_get_provider
    )
    monkeypatch.setattr("chemsmart.agent.core.get_provider", fake_get_provider)
    monkeypatch.setattr(
        "chemsmart.agent.core._default_session_root",
        lambda: str(tmp_path / "sessions"),
    )

    runner = CliRunner()
    result = runner.invoke(
        agent,
        [
            "run",
            "--dry-submit",
            "Recommend a TS method and basis set for a Cope rearrangement.",
        ],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    assert "Advice:" in result.output
    assert "M06-2X/def2-SVP" in result.output
    assert "critic verdict:" not in result.output


def test_ask_renders_advisory_plan(
    monkeypatch,
    tmp_path: Path,
):
    provider = FakeProvider(
        [
            {
                "steps": [],
                "rationale": "Advisory text: start with wB97X-D/def2-SVP, confirm one imaginary frequency, then refine with def2-TZVP.",
                "estimated_cost": "low",
            }
        ]
    )

    def fake_get_provider():
        return provider

    monkeypatch.setattr(
        "chemsmart.agent.providers.get_provider", fake_get_provider
    )
    monkeypatch.setattr("chemsmart.agent.core.get_provider", fake_get_provider)
    monkeypatch.setattr(
        "chemsmart.agent.core._default_session_root",
        lambda: str(tmp_path / "sessions"),
    )

    runner = CliRunner()
    result = runner.invoke(
        agent,
        [
            "ask",
            "Recommend method/basis for a Cope rearrangement TS and explain trade-offs.",
        ],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    assert "Advice" in result.output
    assert "Advisory text" in result.output
    assert "Plan" not in result.output
