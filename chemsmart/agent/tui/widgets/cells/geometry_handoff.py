"""Geometry handoff transcript cell."""

from __future__ import annotations

import subprocess

from rich.console import Group
from rich.syntax import Syntax
from rich.text import Text
from textual.binding import Binding

from chemsmart.io.molecules.structure import Molecule

from ._common import molecule_xyz, write_temp_xyz
from .base import BaseCell


class GeometryHandoffCell(BaseCell):
    BINDINGS = [
        Binding("enter", "toggle_expand", "Expand", show=False),
        Binding("space", "toggle_expand", "Expand", show=False),
        Binding("o", "open_geometry", "Open", show=False),
    ]

    def __init__(
        self,
        molecule: Molecule,
        *,
        session_dir: str | None = None,
        converged: bool = True,
    ) -> None:
        self.molecule = molecule
        self.session_dir = session_dir
        self.converged = converged
        self.expanded = False
        super().__init__(
            self._build_renderable(),
            title="Geometry handoff",
            classes="geometry-cell",
        )

    def action_toggle_expand(self) -> None:
        self.expanded = not self.expanded
        self.update(self._build_renderable())

    def action_open_geometry(self) -> None:
        xyz_path = write_temp_xyz(self.molecule, directory=self.session_dir)
        subprocess.run(["xdg-open", str(xyz_path)], check=False)

    def _build_renderable(self):
        energy = self.molecule.energy
        energy_text = "unknown"
        if energy is not None:
            energy_text = f"{float(energy):.6f} Eh"
        header = Text(
            f"{self.molecule.chemical_formula} · {len(self.molecule)} atoms · "
            f"{energy_text} · converged={str(self.converged).lower()}"
        )
        header.stylize("bold")
        header.append("  [enter expand] [o open]", style="dim")
        if not self.expanded:
            return Group(header)
        return Group(
            header,
            Text(""),
            Syntax(molecule_xyz(self.molecule), "text", theme="ansi_dark"),
        )
