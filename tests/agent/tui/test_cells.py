from __future__ import annotations

import asyncio

from textual.app import App, ComposeResult

from chemsmart.agent.core import CriticVerdict
from chemsmart.agent.tui.widgets.cells import (
    AgentMessageCell,
    CriticVerdictCell,
    DryRunInputCell,
    ErrorCell,
    GeometryHandoffCell,
    MethodCell,
    PlanCell,
    RuntimeValidationCell,
    SubmissionPreviewCell,
    UserMessageCell,
)
from chemsmart.io.molecules.structure import Molecule

from ._helpers import assert_matches_snapshot, normalize_svg


class CellApp(App[None]):
    def __init__(self, cell) -> None:
        super().__init__()
        self.cell = cell

    def compose(self) -> ComposeResult:
        yield self.cell


def _sample_molecule() -> Molecule:
    return Molecule(
        symbols=["O", "H", "H"],
        positions=[[0, 0, 0], [0, 0, 1], [0, 1, 0]],
        energy=-76.123456,
    )


def _method_cell() -> MethodCell:
    return MethodCell(
        {
            "match": "dft_default",
            "functional": "wb97x-d",
            "basis": "def2-svp",
            "solvent_model": "smd",
            "solvent_id": "water",
            "rationale": "Matched project_hint rule.",
        }
    )


def _geometry_cell(expanded: bool = False) -> GeometryHandoffCell:
    cell = GeometryHandoffCell(_sample_molecule(), session_dir="/tmp")
    if expanded:
        cell.expanded = True
        cell.update(cell._build_renderable())
    return cell


_CELL_FACTORIES = {
    "user_message": lambda: UserMessageCell("build water + recommend method"),
    "agent_message": lambda: AgentMessageCell("Phase 2 help output."),
    "plan": lambda: PlanCell(
        'Plan:\n1. build_molecule {"filepath": "water.xyz"}\n'
        "   - Load the starting structure."
    ),
    "dry_run_input": lambda: DryRunInputCell(
        "#p b3lyp/6-31g* opt\n\nwater\n\n0 1\nO 0 0 0\n",
        inputfile="/tmp/water.com",
    ),
    "dry_run_input_diff": lambda: DryRunInputCell(
        "#p wb97x-d/def2svp opt\n\nwater\n\n0 1\nO 0 0 0\n",
        inputfile="/tmp/water.com",
        previous_content="#p b3lyp/6-31g* opt\n\nwater\n\n0 1\nO 0 0 0\n",
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
    "method": _method_cell,
    "geometry_handoff": _geometry_cell,
    "geometry_handoff_expanded": lambda: _geometry_cell(expanded=True),
    "runtime_validation": lambda: RuntimeValidationCell(
        {
            "ok": "partial",
            "local_ok": True,
            "local_issues": [],
            "remote_unknown": ["ssh login reachable"],
        }
    ),
    "submission_preview": lambda: SubmissionPreviewCell(
        {
            "transport": "LocalDryRunTransport",
            "script_path": "/tmp/submit.sh",
            "command_executed": ["sbatch", "submit.sh"],
            "duplicate_check": {"duplicate": False},
        }
    ),
}


def test_phase2_cells_match_snapshots():
    async def scenario() -> None:
        for name, factory in _CELL_FACTORIES.items():
            app = CellApp(factory())
            async with app.run_test() as pilot:
                await pilot.pause()
                cell = app.cell
                if name == "dry_run_input_diff":
                    cell.action_toggle_diff()
                    await pilot.pause()
                assert_matches_snapshot(
                    name, normalize_svg(app.export_screenshot())
                )

    asyncio.run(scenario())
