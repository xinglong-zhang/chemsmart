"""Escape must be armed, then confirmed -- and it only reaches planning.

A single stray keypress cannot end a session; the first esc arms with an
in-place footer change, the second sets the cooperative event the loop
checks between provider turns. Outside planning, esc does nothing.
"""

from __future__ import annotations

import asyncio
from pathlib import Path

import pytest

pytest.importorskip("textual")

from chemsmart.agent.tui.app import ChemSmartAgentApp  # noqa: E402
from chemsmart.agent.tui.controller import (  # noqa: E402
    AgentSessionConfigV1,
    AgentTuiController,
    AgentTuiPhase,
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


def test_first_esc_arms_and_second_esc_sets_the_cancel_event(tmp_path: Path):
    app = _app(tmp_path)

    async def observe() -> None:
        async with app.run_test(size=(120, 35)) as pilot:
            app._busy = True
            app.controller.phase = AgentTuiPhase.PLANNING

            await pilot.press("escape")
            assert app._esc_armed is True
            assert not app.controller.cancel_planning.is_set()
            footer = app.query_one("#footer")
            assert "esc again to cancel planning" in str(footer.render())

            await pilot.press("escape")
            assert app.controller.cancel_planning.is_set()
            assert app._esc_armed is False
            app._busy = False

    asyncio.run(observe())


def test_esc_outside_planning_arms_nothing(tmp_path: Path):
    app = _app(tmp_path)

    async def observe() -> None:
        async with app.run_test(size=(120, 35)) as pilot:
            await pilot.press("escape")
            assert app._esc_armed is False
            assert not app.controller.cancel_planning.is_set()

    asyncio.run(observe())


def test_the_disarm_timer_restores_the_footer(tmp_path: Path):
    app = _app(tmp_path)

    async def observe() -> None:
        async with app.run_test(size=(120, 35)) as pilot:
            app._busy = True
            app.controller.phase = AgentTuiPhase.PLANNING
            await pilot.press("escape")
            assert app._esc_armed is True
            app._disarm_interrupt()
            assert app._esc_armed is False
            footer = app.query_one("#footer")
            assert "esc again" not in str(footer.render())
            app._busy = False

    asyncio.run(observe())
