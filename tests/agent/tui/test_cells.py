from __future__ import annotations

import asyncio

from textual.app import App, ComposeResult

from chemsmart.agent.core import CriticVerdict
from chemsmart.agent.tui.widgets.cells import (
    AgentMessageCell,
    CriticVerdictCell,
    DryRunInputCell,
    ErrorCell,
    PlanCell,
    UserMessageCell,
)

from ._helpers import assert_matches_snapshot, normalize_svg


class CellApp(App[None]):
    def __init__(self, cell) -> None:
        super().__init__()
        self.cell = cell

    def compose(self) -> ComposeResult:
        yield self.cell


_CELL_FACTORIES = {
    "user_message": lambda: UserMessageCell("build water + recommend method"),
    "agent_message": lambda: AgentMessageCell("Phase 1 help output."),
    "plan": lambda: PlanCell(
        'Plan:\n1. build_molecule {"filepath": "water.xyz"}\n'
        "   - Load the starting structure."
    ),
    "dry_run_input": lambda: DryRunInputCell(
        "#p b3lyp/6-31g* opt\n\nwater\n\n0 1\nO 0 0 0\n",
        inputfile="/tmp/water.com",
    ),
    "critic_verdict": lambda: CriticVerdictCell(
        CriticVerdict(
            verdict="ok",
            confidence=0.92,
            issues=[],
            rationale="The generated input looks reasonable.",
        )
    ),
    "error": lambda: ErrorCell(
        "RuntimeError",
        "run_local failed with returncode 1",
        {"tool": "run_local"},
    ),
}


def test_phase1_cells_match_snapshots():
    async def scenario() -> None:
        for name, factory in _CELL_FACTORIES.items():
            app = CellApp(factory())
            async with app.run_test() as pilot:
                await pilot.pause()
                assert_matches_snapshot(
                    name, normalize_svg(app.export_screenshot())
                )

    asyncio.run(scenario())
