from __future__ import annotations

import asyncio
from pathlib import Path

from click.testing import CliRunner
from textual import events

from chemsmart.agent.cli import agent
from chemsmart.agent.tui.app import ChemsmartTuiApp, launch_tui
from chemsmart.agent.tui.widgets.composer import Composer
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.popups import ApprovalOverlay, ApprovalResult
from chemsmart.agent.tui.widgets.transcript import Transcript

from .._agent_session_helpers import FakeProvider, critic_ok, planner_plan


def test_approval_overlay_blocks_run_until_yes(
    monkeypatch, single_molecule_xyz_file, tmp_path: Path
):
    import chemsmart.agent.tools as agent_tools

    provider = FakeProvider(
        [
            planner_plan(
                single_molecule_xyz_file, "approval_case", "run_local"
            ),
            critic_ok(),
        ]
    )
    run_local_calls = {"count": 0}

    def fake_get_provider():
        return provider

    def fake_validate_runtime(job, server=None):
        return {
            "ok": "ok",
            "local_ok": True,
            "local_issues": [],
            "remote_unknown": [],
        }

    def fake_run_local(job):
        run_local_calls["count"] += 1
        return {
            "ok": True,
            "returncode": 0,
            "stdout_path": str(Path(job.folder) / "run.stdout"),
            "stderr_path": str(Path(job.folder) / "run.stderr"),
            "output_summary": {},
        }

    monkeypatch.setattr(
        "chemsmart.agent.providers.get_provider", fake_get_provider
    )
    monkeypatch.setattr("chemsmart.agent.core.get_provider", fake_get_provider)
    monkeypatch.setattr(agent_tools, "validate_runtime", fake_validate_runtime)
    monkeypatch.setattr(agent_tools, "run_local", fake_run_local)

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            composer = app.query_one(Composer)
            composer.load_text(f"optimize {single_molecule_xyz_file}")
            await pilot.press("enter")
            for _ in range(20):
                await pilot.pause(0.1)
                if (
                    "approval required"
                    in str(app.query_one(FooterWidget).renderable).lower()
                ):
                    break
            assert run_local_calls["count"] == 0
            app.chat_screen._handle_approval_result(ApprovalResult("y"))
            for _ in range(20):
                await pilot.pause(0.1)
                if run_local_calls["count"] == 1:
                    break
            assert run_local_calls["count"] == 1
            assert (
                "finished"
                in str(app.query_one(FooterWidget).renderable).lower()
            )

    asyncio.run(scenario())


def test_decline_and_revise_restarts_with_correction(
    monkeypatch, single_molecule_xyz_file, tmp_path: Path
):
    import chemsmart.agent.tools as agent_tools

    provider = FakeProvider(
        [
            planner_plan(single_molecule_xyz_file, "revise_case", "run_local"),
            critic_ok(),
            planner_plan(
                single_molecule_xyz_file, "revise_case_2", "run_local"
            ),
            critic_ok(),
        ]
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

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            composer = app.query_one(Composer)
            composer.load_text(f"optimize {single_molecule_xyz_file}")
            await pilot.press("enter")
            for _ in range(20):
                await pilot.pause(0.1)
                if (
                    "approval required"
                    in str(app.query_one(FooterWidget).renderable).lower()
                ):
                    break
            app.chat_screen._handle_approval_result(
                ApprovalResult("r", corrective_text="use def2-tzvp")
            )
            for _ in range(20):
                await pilot.pause(0.1)
                transcript = app.query_one(Transcript).query_one("#cells")
                texts = [
                    str(child.renderable) for child in transcript.children
                ]
                if any("Corrective instruction:" in text for text in texts):
                    break
            transcript = app.query_one(Transcript).query_one("#cells")
            texts = [str(child.renderable) for child in transcript.children]
            assert any(
                "Corrective instruction: use def2-tzvp" in text
                for text in texts
            )

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

    import chemsmart.agent.tools as agent_tools

    provider = FakeProvider(
        [planner_plan(single_molecule_xyz_file, "ask_case"), critic_ok()]
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
    monkeypatch.setattr(agent_tools, "validate_runtime", fake_validate_runtime)
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
    assert "Plan" in result.output
    assert "Critic" in result.output
