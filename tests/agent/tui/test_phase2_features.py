from __future__ import annotations

import asyncio
from pathlib import Path

from click.testing import CliRunner
from textual import events

from chemsmart.agent.cli import agent
from chemsmart.agent.tui.app import ChemsmartTuiApp, launch_tui
from chemsmart.agent.tui.widgets.composer import Composer
from chemsmart.agent.tui.widgets.popups import ApprovalOverlay, ApprovalResult

from .._agent_session_helpers import FakeProvider
from .._loop_helpers import openai_final_response


def test_approval_overlay_blocks_run_until_yes(
    tmp_path: Path,
):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        results: list[ApprovalResult | None] = []
        async with app.run_test() as pilot:
            await pilot.pause()
            app.push_screen(
                ApprovalOverlay(action="run_local", request="optimize water"),
                results.append,
            )
            await pilot.pause()
            await pilot.press("y")
            await pilot.pause()
            assert results == [ApprovalResult("y")]

    asyncio.run(scenario())


def test_decline_and_revise_restarts_with_correction(
    tmp_path: Path,
):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            submitted: list[str] = []
            app.chat_screen.start_request = submitted.append
            app.chat_screen._handle_inline_approval("revise use def2-tzvp")
            assert submitted == [
                "Corrective instruction: use def2-tzvp\n\nOriginal request:\n"
            ]

    asyncio.run(scenario())


def test_composer_large_paste_placeholder_and_external_editor(
    monkeypatch, tmp_path: Path
):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            composer = app.query_one(Composer)
            await composer._on_paste(events.Paste("x" * 10001))
            assert composer.text == "[Pasted 10001 chars]"
            assert composer.resolve_text() == "x" * 10001

            def fake_run(cmd, check=False):
                path = Path(cmd[-1])
                path.write_text("edited externally", encoding="utf-8")

            monkeypatch.setattr("subprocess.run", fake_run)
            monkeypatch.setenv("EDITOR", "vi")
            composer.action_external_editor()
            assert composer.text == "edited externally"

    asyncio.run(scenario())


def test_composer_submit_is_guarded_until_screen_handles_message(
    tmp_path: Path,
):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        submitted: list[str] = []
        app.chat_screen.start_request = submitted.append

        async with app.run_test() as pilot:
            await pilot.pause()
            composer = app.query_one(Composer)
            composer.load_text("single submit only")
            composer.action_submit()
            composer.action_submit()
            await pilot.pause()
            assert submitted == ["single submit only"]
            composer.clear_text()
            composer.load_text("next request")
            composer.action_submit()
            await pilot.pause()
            assert submitted == ["single submit only", "next request"]

    asyncio.run(scenario())


def test_approval_overlay_stops_s_key_from_reaching_composer(tmp_path: Path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        results: list[ApprovalResult | None] = []
        async with app.run_test() as pilot:
            await pilot.pause()
            composer = app.query_one(Composer)
            composer.focus()
            app.push_screen(
                ApprovalOverlay(action="run_local", request="optimize water"),
                results.append,
            )
            await pilot.pause()
            await pilot.press("s")
            await pilot.pause()
            assert results
            assert results[0] == ApprovalResult("s")
            assert composer.text == ""
            assert app.screen is app.chat_screen

    asyncio.run(scenario())


def test_launch_tui_uses_fullscreen_by_default(monkeypatch, tmp_path: Path):
    captured: list[tuple[bool, dict, str]] = []

    def fake_run(self, **kwargs) -> None:
        captured.append((self.plain, kwargs, self.animation_level))

    monkeypatch.setattr(ChemsmartTuiApp, "run", fake_run)

    launch_tui(session_root=tmp_path / "sessions")

    assert captured == [(False, {}, "full")]


def test_launch_tui_plain_mode_stays_inline(monkeypatch, tmp_path: Path):
    captured: list[tuple[bool, dict, str]] = []

    def fake_run(self, **kwargs) -> None:
        captured.append((self.plain, kwargs, self.animation_level))

    monkeypatch.setattr(ChemsmartTuiApp, "run", fake_run)

    launch_tui(plain=True, session_root=tmp_path / "sessions")

    assert captured == [
        (
            True,
            {"inline": True, "inline_no_clear": True, "mouse": False},
            "none",
        )
    ]


def test_agent_cli_ask_streams_without_tui_import(
    monkeypatch, single_molecule_xyz_file, tmp_path: Path
):
    import builtins

    provider = FakeProvider(
        [{"__raw_response__": openai_final_response("Loop answer.")}]
    )
    real_import = builtins.__import__

    def fake_get_provider():
        return provider

    def fake_validate_runtime(job, server=None):
        return {
            "ok": "ok",
            "local_ok": True,
            "local_issues": [],
            "remote_unknown": [],
        }

    def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "chemsmart.agent.tui":
            raise AssertionError("ask should not import the TUI")
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", fake_import)
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
        ["ask", f"optimize {single_molecule_xyz_file}"],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    assert "Assistant" in result.output
    assert "Loop answer." in result.output
