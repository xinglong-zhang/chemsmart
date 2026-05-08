from __future__ import annotations

import asyncio

from textual.app import App, ComposeResult

from chemsmart.agent.tui.widgets.cells import (
    JobStatusCell,
    MoleculeCell,
    RunResultCell,
    WorkflowCell,
)
from chemsmart.io.molecules.structure import Molecule

from ._helpers import assert_matches_snapshot, normalize_svg


class CellApp(App[None]):
    def __init__(self, cell) -> None:
        super().__init__()
        self.cell = cell

    def compose(self) -> ComposeResult:
        yield self.cell


def _molecule() -> Molecule:
    return Molecule(
        symbols=["O", "H", "H"],
        positions=[[0.0, 0.0, 0.1], [-0.8, 0.0, -0.5], [0.8, 0.0, -0.5]],
        charge=0,
        multiplicity=1,
    )


_FACTORIES = {
    "job_status": lambda: JobStatusCell(
        "12345.remote",
        {
            "name": "water_opt",
            "scheduler": "slurm",
            "status": "running",
            "started": "2026-05-09 00:10",
            "runtime": "5m 12s",
            "host": "remote-hpc",
        },
    ),
    "run_result": lambda: RunResultCell(
        {
            "energy": -76.432198,
            "delta_g": -76.401234,
            "frequencies": [102.4, 155.2, 324.1],
            "imag_freqs": [],
            "output_path": "/tmp/water_opt.log",
            "input_path": "/tmp/water_opt.com",
            "normal_termination": True,
        }
    ),
    "workflow": lambda: WorkflowCell(
        "opt+freq",
        [
            {"name": "opt", "status": "done", "detail": "geometry converged"},
            {
                "name": "freq",
                "status": "running",
                "detail": "parsing frequencies",
            },
            {"name": "irc", "status": "queued", "detail": "optional"},
        ],
    ),
    "workflow_expanded": lambda: WorkflowCell(
        "opt+freq",
        [
            {"name": "opt", "status": "done", "detail": "geometry converged"},
            {
                "name": "freq",
                "status": "running",
                "detail": "parsing frequencies",
            },
            {"name": "irc", "status": "queued", "detail": "optional"},
        ],
        expanded=True,
    ),
    "molecule": lambda: MoleculeCell(_molecule(), source="examples/h2o.xyz"),
    "molecule_coords": lambda: _coords_cell(),
}


def _coords_cell() -> MoleculeCell:
    cell = MoleculeCell(_molecule(), source="examples/h2o.xyz")
    cell.show_coords = True
    cell.update(cell._build_renderable())
    return cell


def test_phase3_cells_match_snapshots():
    async def scenario() -> None:
        for name, factory in _FACTORIES.items():
            app = CellApp(factory())
            async with app.run_test() as pilot:
                await pilot.pause()
                assert_matches_snapshot(
                    name, normalize_svg(app.export_screenshot())
                )

    asyncio.run(scenario())
