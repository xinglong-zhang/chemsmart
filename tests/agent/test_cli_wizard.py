from __future__ import annotations

from click.testing import CliRunner

from chemsmart.agent.cli import agent, sanitize_inline_cli_output
from chemsmart.agent.wizard.orchestrator import WizardOutcome
from chemsmart.agent.wizard.render import ServerYamlPlan
from chemsmart.agent.wizard.validate import ValidationResult

YAML_TEXT = """SERVER:
  SCHEDULER: SLURM
GAUSSIAN: {}
"""


def _wizard_outcome() -> WizardOutcome:
    return WizardOutcome(
        plan=ServerYamlPlan(
            text=YAML_TEXT,
            server_block={"SCHEDULER": "SLURM"},
            program_blocks={"GAUSSIAN": {}},
            notes=[],
        ),
        validation=ValidationResult(
            ok=True,
            errors=[],
            parsed={"SERVER": {"SCHEDULER": "SLURM"}},
        ),
        surveys={},
        written_path=None,
    )


def test_agent_wizard_prints_yaml_stdout(monkeypatch):
    monkeypatch.setattr(
        "chemsmart.agent.cli.run_wizard",
        lambda *args, **kwargs: _wizard_outcome(),
    )

    result = CliRunner().invoke(
        agent,
        ["wizard", "perlmutter"],
        catch_exceptions=False,
    )

    assert result.exit_code == 0
    assert sanitize_inline_cli_output(result.output) == YAML_TEXT.rstrip()


def test_agent_wizard_write_yes_writes_to_home(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "chemsmart.agent.cli.run_wizard",
        lambda *args, **kwargs: _wizard_outcome(),
    )
    monkeypatch.setenv("HOME", str(tmp_path))

    result = CliRunner().invoke(
        agent,
        ["wizard", "perlmutter", "--write", "--yes"],
        catch_exceptions=False,
    )

    expected = tmp_path / ".chemsmart" / "server" / "perlmutter.yaml"
    assert result.exit_code == 0
    assert expected.read_text(encoding="utf-8") == YAML_TEXT
    assert f"Wrote {expected}" in result.output
