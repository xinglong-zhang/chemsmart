from __future__ import annotations

import pytest


@pytest.fixture(autouse=True)
def _isolate_tui_provider_config(monkeypatch):
    """Keep TUI tests off the developer's live ``~/.chemsmart/agent/agent.yaml``.

    ``ChatScreen`` resolves its interaction mode from the active provider config
    (a ``local`` provider defaults to ask/synthesis, everything else to
    run/harness). Reading the real on-disk config would make TUI tests depend on
    whatever provider the machine happens to have active. Default to "no active
    config" (→ run mode), and let individual tests monkeypatch
    ``load_active_provider_config`` when they need a specific provider.
    """

    monkeypatch.setattr(
        "chemsmart.agent.tui.screens.chat.load_active_provider_config",
        lambda: None,
        raising=False,
    )
