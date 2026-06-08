from __future__ import annotations

import builtins

from click.testing import CliRunner

from chemsmart.agent.cli import agent


def test_bare_agent_prints_install_hint_when_textual_is_missing(monkeypatch):
    real_import = builtins.__import__

    def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "chemsmart.agent.tui":
            raise ImportError("textual not installed")
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", fake_import)

    runner = CliRunner()
    result = runner.invoke(agent, [], catch_exceptions=False)

    assert result.exit_code == 1
    assert 'pip install -e ".[agent-tui]"' in result.output
