import builtins
import importlib

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
    real_import = builtins.__import__

    def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "chemsmart.agent.cli":
            raise ImportError("optional deps missing")
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    try:
        main = importlib.reload(cli_main)
        result = CliRunner(mix_stderr=False).invoke(
            main.entry_point, ["agent"]
        )
    finally:
        monkeypatch.setattr(builtins, "__import__", real_import)
        importlib.reload(cli_main)

    assert result.exit_code == 1
    assert 'pip install -e ".[agent-tui]"' in result.stderr
