from __future__ import annotations

from click.testing import CliRunner

from chemsmart.agent.cli import agent, sanitize_inline_cli_output


def test_agent_wizard_verify_raises_for_missing_yaml(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))

    result = CliRunner().invoke(agent, ["wizard-verify", "missing-server"])

    assert result.exit_code == 1
    assert "Server YAML not found" in result.output


def test_agent_wizard_verify_prints_expected_fields(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    target = tmp_path / ".chemsmart" / "server" / "perlmutter.yaml"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        """SERVER:
  SCHEDULER: SLURM
  QUEUE_NAME: debug
  NUM_HOURS: 8
  MEM_GB: 64
  NUM_CORES: 16
  SUBMIT_COMMAND: sbatch
  HOST: localhost
GAUSSIAN:
  EXEFOLDER: /apps/gaussian
  LOCAL_RUN: true
  SCRATCH: true
  CONDA_ENV: ""
  MODULES: ""
  ENVARS: ""
""",
        encoding="utf-8",
    )

    result = CliRunner().invoke(
        agent,
        ["wizard-verify", "perlmutter"],
        catch_exceptions=False,
    )

    output = sanitize_inline_cli_output(result.output)
    assert result.exit_code == 0
    assert "server_name" in output
    assert "perlmutter" in output
    assert "mode" in output
    assert "local" in output
    assert "would_submit_via" in output
    assert "transport_invocation" in output
