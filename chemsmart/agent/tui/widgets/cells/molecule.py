"""Molecule summary transcript cell."""

from __future__ import annotations

from rich.console import Group
from rich.table import Table
from rich.text import Text
from textual.binding import Binding

from chemsmart.io.molecules.structure import Molecule

from .base import BaseCell


class MoleculeCell(BaseCell):
    BINDINGS = [Binding("c", "toggle_coords", "Coords", show=False)]

    def __init__(
        self, molecule: Molecule, *, source: str | None = None
    ) -> None:
        self.molecule = molecule
        self.source = source
        self.show_coords = False
        super().__init__(
            self._build_renderable(),
            title="Molecule",
            classes="molecule-cell",
        )

    def action_toggle_coords(self) -> None:
        self.show_coords = not self.show_coords
        self.update(self._build_renderable())

    def _build_renderable(self):
        coords_rows = [
            tuple(float(value) for value in coords[:3])
            for coords in self.molecule.positions
            if len(coords) >= 3
        ]
        bbox = (0.0, 0.0, 0.0)
        if coords_rows:
            x_values = [coords[0] for coords in coords_rows]
            y_values = [coords[1] for coords in coords_rows]
            z_values = [coords[2] for coords in coords_rows]
            bbox = (
                max(x_values) - min(x_values),
                max(y_values) - min(y_values),
                max(z_values) - min(z_values),
            )
        header = Text(
            (
                f"{self.molecule.chemical_formula or 'unknown'} · q={self.molecule.charge} "
                f"mult={self.molecule.multiplicity} · {len(self.molecule)} atoms · "
                f"bbox={bbox[0]:.2f}×{bbox[1]:.2f}×{bbox[2]:.2f} Å"
            ),
            style="bold",
        )
        header.append("  [c coords]", style="dim")
        header.append(
            f"\nsource: {self.source or 'no data'}",
            style="dim",
        )
        if not self.show_coords:
            return Group(header)
        table = Table(show_header=True, header_style="bold")
        table.add_column("#", justify="right")
        table.add_column("atom")
        table.add_column("x", justify="right")
        table.add_column("y", justify="right")
        table.add_column("z", justify="right")
        rows_added = 0
        for index, (symbol, coords) in enumerate(
            zip(self.molecule.symbols, self.molecule.positions),
            start=1,
        ):
            if len(coords) < 3:
                continue
            rows_added += 1
            table.add_row(
                str(index),
                str(symbol),
                f"{float(coords[0]):.4f}",
                f"{float(coords[1]):.4f}",
                f"{float(coords[2]):.4f}",
            )
        if rows_added == 0:
            table.add_row("—", "no data", "", "", "")
        return Group(header, table)
