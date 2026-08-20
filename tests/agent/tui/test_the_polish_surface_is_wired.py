"""Palette, skill tagging, and export all read the one command registry."""

from __future__ import annotations

import asyncio
from pathlib import Path

import pytest

pytest.importorskip("textual")

from chemsmart.agent.tui.app import ChemSmartAgentApp  # noqa: E402
from chemsmart.agent.tui.commands import COMMANDS, command_for  # noqa: E402
from chemsmart.agent.tui.controller import (  # noqa: E402
    AgentSessionConfigV1,
    AgentTuiController,
)


def _app(tmp_path: Path) -> ChemSmartAgentApp:
    secret = tmp_path / "secret.env"
    secret.write_text("PROVIDER_KEY=not-used\n", encoding="utf-8")
    return ChemSmartAgentApp(
        AgentTuiController(
            AgentSessionConfigV1(workspace=tmp_path, secret_file=secret)
        ),
        plain=True,
    )


def test_every_registered_command_has_a_bound_handler(tmp_path: Path):
    app = _app(tmp_path)

    async def observe() -> None:
        async with app.run_test(size=(120, 35)):
            handlers = app._command_handlers()
            for spec in COMMANDS:
                assert spec.name in handlers, spec.name

    asyncio.run(observe())


def test_a_tagged_skill_prefixes_the_next_request_visibly(tmp_path: Path):
    app = _app(tmp_path)

    async def observe() -> None:
        async with app.run_test(size=(120, 35)) as pilot:
            composer = app.query_one("#composer")
            composer.value = "/skill scientific-conventions"
            await pilot.press("enter")
            await pilot.pause()
            assert app._pending_skill == "scientific-conventions"

            captured = {}

            def fake_begin(task: str) -> str:
                captured["task"] = task
                raise __import__(
                    "chemsmart.agent._contracts", fromlist=["ContractError"]
                ).ContractError("stop before any provider call")

            app.controller.begin_planning = fake_begin  # type: ignore
            composer.value = "optimize water with xtb"
            await pilot.press("enter")
            await pilot.pause()
            assert captured["task"].startswith(
                "The user tagged domain skill 'scientific-conventions'"
            )
            assert "optimize water with xtb" in captured["task"]
            assert app._pending_skill == ""

    asyncio.run(observe())


def test_an_unknown_skill_is_refused_with_the_known_list(tmp_path: Path):
    app = _app(tmp_path)

    async def observe() -> None:
        async with app.run_test(size=(120, 35)) as pilot:
            composer = app.query_one("#composer")
            composer.value = "/skill not-a-skill"
            await pilot.press("enter")
            await pilot.pause()
            assert app._pending_skill == ""

    asyncio.run(observe())


def test_export_writes_the_transcript_into_the_workspace(tmp_path: Path):
    app = _app(tmp_path)

    async def observe() -> None:
        async with app.run_test(size=(120, 35)) as pilot:
            composer = app.query_one("#composer")
            composer.value = "/help"
            await pilot.press("enter")
            await pilot.pause()
            composer.value = "/export"
            await pilot.press("enter")
            await pilot.pause()
            exports = list(tmp_path.glob("chemsmart-transcript-*.txt"))
            assert len(exports) == 1
            content = exports[0].read_text(encoding="utf-8")
            assert "/approve" in content

    asyncio.run(observe())


def test_the_registry_knows_the_full_production_surface():
    for slash in (
        "/approve",
        "/deny",
        "/revise",
        "/status",
        "/capabilities",
        "/skills",
        "/skill",
        "/dag",
        "/report",
        "/runs",
        "/export",
        "/help",
        "/quit",
        "/exit",
    ):
        assert command_for(slash) is not None, slash
