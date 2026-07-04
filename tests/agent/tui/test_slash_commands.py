from __future__ import annotations

import asyncio

from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.widgets.composer import Composer

from .._agent_session_helpers import FakeProvider
from ._helpers import (
    assert_matches_snapshot,
    normalize_svg,
    write_session_fixture,
)


def _set_composer_text(app: ChemsmartTuiApp, text: str) -> None:
    composer = app.query_one(Composer)
    composer.load_text(text)


def test_phase1_slash_commands_match_snapshots(monkeypatch, tmp_path):
    session_root = tmp_path / "sessions"
    write_session_fixture(session_root)

    def fake_get_provider():
        return FakeProvider([])

    monkeypatch.setattr(
        "chemsmart.agent.providers.get_provider", fake_get_provider
    )
    monkeypatch.setattr("chemsmart.agent.core.get_provider", fake_get_provider)
    monkeypatch.setenv("AI_PROVIDER", "openai")
    monkeypatch.setenv("ai_api_key", "test-key")
    monkeypatch.setenv("CHEMSMART_AGENT_TUI_MODE", "run")

    async def snapshot_for(
        command: str, name: str, *, patch_exit: bool = False
    ):
        app = ChemsmartTuiApp(session_root=session_root)
        exit_calls: list[bool] = []
        if patch_exit:
            app.exit = lambda *args, **kwargs: exit_calls.append(True)  # type: ignore[method-assign]
        async with app.run_test() as pilot:
            await pilot.pause()
            _set_composer_text(app, command)
            await pilot.press("enter")
            await pilot.pause()
            if command == "/sessions":
                await pilot.pause()
            if command.startswith("/resume"):
                await pilot.pause(0.2)
            if patch_exit:
                assert exit_calls == [True]
            assert_matches_snapshot(
                name, normalize_svg(app.export_screenshot())
            )

    async def scenario() -> None:
        await snapshot_for("/help", "slash_help")
        await snapshot_for("/sessions", "slash_sessions")
        await snapshot_for("/resume session-001", "slash_resume")
        await snapshot_for("/clear", "slash_clear")
        await snapshot_for("/quit", "slash_quit", patch_exit=True)

    asyncio.run(scenario())


def test_init_enters_project_yaml_build_mode(monkeypatch, tmp_path):
    session_root = tmp_path / "sessions"
    write_session_fixture(session_root)

    def fake_get_provider():
        return FakeProvider([])

    monkeypatch.setattr(
        "chemsmart.agent.providers.get_provider", fake_get_provider
    )
    monkeypatch.setattr("chemsmart.agent.core.get_provider", fake_get_provider)
    monkeypatch.setenv("AI_PROVIDER", "openai")
    monkeypatch.setenv("ai_api_key", "test-key")
    monkeypatch.setenv("CHEMSMART_AGENT_TUI_MODE", "run")

    class _FakeWorker:
        is_finished = True

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=session_root)
        async with app.run_test() as pilot:
            await pilot.pause()
            screen = app.chat_screen
            captured: list[str] = []
            monkeypatch.setattr(
                screen,
                "run_build_session",
                lambda text: captured.append(text) or _FakeWorker(),
            )

            # Bare /init enters a persistent build mode.
            _set_composer_text(app, "/init")
            await pilot.press("enter")
            await pilot.pause()
            assert screen._build_mode is True

            # A subsequent message routes to the project-YAML build harness.
            _set_composer_text(
                app, "Optimize in water with B3LYP-D3BJ/def2-SVP, call it h2o."
            )
            await pilot.press("enter")
            await pilot.pause()
            assert captured == [
                "Optimize in water with B3LYP-D3BJ/def2-SVP, call it h2o."
            ]
            assert screen._build_mode is True

            # /mode leaves build mode.
            _set_composer_text(app, "/mode ask")
            await pilot.press("enter")
            await pilot.pause()
            assert screen._build_mode is False

    asyncio.run(scenario())


def test_init_slash_command_is_registered_in_palette():
    from chemsmart.agent.tui.screens.chat import _SLASH_PALETTE_COMMANDS

    commands = {name for name, _desc in _SLASH_PALETTE_COMMANDS}
    assert "/init" in commands
