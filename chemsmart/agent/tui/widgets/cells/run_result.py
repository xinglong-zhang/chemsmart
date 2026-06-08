"""Parsed run-result transcript cell."""

from __future__ import annotations

from rich.console import Group
from rich.table import Table
from rich.text import Text

from .base import BaseCell


class RunResultCell(BaseCell):
    def __init__(self, result: dict | None) -> None:
        self.result = result or {}
        super().__init__(
            self._build_renderable(),
            title="Run result",
            classes="run-result-cell",
        )

    def _build_renderable(self):
        table = Table.grid(expand=True)
        table.add_column(ratio=2)
        table.add_column(ratio=4)
        energy = self.result.get("energy")
        delta_g = self.result.get("delta_g")
        table.add_row(
            "E_SCF",
            f"{float(energy):.6f} Eh" if energy is not None else "no data",
        )
        table.add_row(
            "ΔG",
            f"{float(delta_g):.6f} Eh" if delta_g is not None else "no data",
        )
        imag = list(self.result.get("imag_freqs") or [])
        freqs = list(self.result.get("frequencies") or [])
        if imag:
            freq_text = ", ".join(f"{float(freq):.1f}" for freq in imag)
            table.add_row("imag freqs", Text(freq_text, style="error"))
        else:
            freq_text = ", ".join(f"{float(freq):.1f}" for freq in freqs)
            table.add_row("lowest freqs", freq_text or "no data")
        table.add_row(
            "output",
            str(
                self.result.get("output_path")
                or self.result.get("source")
                or "no data"
            ),
        )
        input_path = self.result.get("input_path")
        if input_path:
            table.add_row("input", str(input_path))
        heading = Text(
            (
                "normal termination"
                if self.result.get("normal_termination")
                else "parsed result"
            ),
            style="bold",
        )
        return Group(heading, table)
