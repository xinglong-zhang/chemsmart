"""Mutable transcript receipt for one calculation lifecycle."""

from __future__ import annotations

from rich.console import Group
from rich.table import Table
from rich.text import Text

from .base import BaseCell


class CalculationReceiptCell(BaseCell):
    def __init__(self, run: dict[str, object]) -> None:
        self.run = dict(run)
        super().__init__(
            self._build_renderable(),
            title="Calculation",
            classes="calculation-receipt-cell",
        )

    @property
    def run_id(self) -> str:
        return str(self.run.get("run_id") or "")

    def update_run(self, run: dict[str, object]) -> None:
        self.run.update(run)
        self.update(self._build_renderable())

    def _build_renderable(self):
        status = str(self.run.get("status") or "unknown")
        table = Table.grid(expand=True)
        table.add_column(style="dim", width=18)
        table.add_column(ratio=1)
        self._add_method_rows(table)
        self._add_progress_rows(table)
        self._add_runtime_rows(table, status)
        return Group(self._heading(status), table)

    def _heading(self, status: str) -> Text:
        style = (
            "success"
            if status == "completed"
            else (
                "accent"
                if status in {"validating", "starting", "running"}
                else "error"
            )
        )
        program = str(self.run.get("program") or "calculation").upper()
        kind = str(self.run.get("kind") or "job").upper()
        label = str(self.run.get("label") or "calculation")
        return Text.assemble(
            (status.upper(), f"bold {style}"),
            (f" · {program} {kind} · ", "dim"),
            (label, "bold"),
        )

    def _add_method_rows(self, table: Table) -> None:
        method = "/".join(
            item
            for item in (
                str(self.run.get("method") or ""),
                str(self.run.get("basis") or ""),
            )
            if item
        )
        if method:
            table.add_row("method", method)
        stage = str(self.run.get("stage") or "")
        if stage:
            table.add_row("stage", stage)
        if self.run.get("reused_output") is True:
            table.add_row("output source", "existing completed calculation")
        energy = self.run.get("energy")
        if isinstance(energy, (int, float)):
            table.add_row("E(SCF)", f"{float(energy):.12f} Eh")
        scf_cycles = self.run.get("scf_cycles")
        if isinstance(scf_cycles, int):
            table.add_row("SCF cycles", str(scf_cycles))

    def _add_progress_rows(self, table: Table) -> None:
        opt_cycles = self.run.get("optimization_cycles")
        if isinstance(opt_cycles, int):
            convergence = self.run.get("optimization_converged")
            suffix = (
                " · converged"
                if convergence is True
                else " · not confirmed" if convergence is False else ""
            )
            table.add_row("optimization", f"{opt_cycles} cycles{suffix}")
        imag = list(self.run.get("imag_freqs") or [])
        imag_count = self.run.get("imaginary_frequency_count")
        if imag or isinstance(imag_count, int):
            table.add_row(
                "imag frequencies",
                f"{imag_count if isinstance(imag_count, int) else len(imag)}"
                + (
                    " · "
                    + ", ".join(f"{float(value):.1f}" for value in imag)
                    + " cm^-1"
                    if imag
                    else ""
                ),
            )
        scan_completed = self.run.get("scan_points_completed")
        if isinstance(scan_completed, int):
            scan_total = self.run.get("scan_points_total")
            table.add_row(
                "scan points",
                (
                    f"{scan_completed}/{scan_total}"
                    if isinstance(scan_total, int)
                    else str(scan_completed)
                ),
            )
        neb_images = self.run.get("neb_images")
        if isinstance(neb_images, int):
            table.add_row("NEB images", str(neb_images))
        qm_atoms = self.run.get("qmmm_qm_atoms")
        mm_atoms = self.run.get("qmmm_mm_atoms")
        if isinstance(qm_atoms, int) or isinstance(mm_atoms, int):
            table.add_row(
                "QMMM regions",
                f"QM {qm_atoms if isinstance(qm_atoms, int) else '?'} · "
                f"MM {mm_atoms if isinstance(mm_atoms, int) else '?'}",
            )
        direction = str(self.run.get("path_direction") or "")
        if direction:
            table.add_row("path direction", direction)

    def _add_runtime_rows(self, table: Table, status: str) -> None:
        output_path = str(self.run.get("output_path") or "")
        if output_path:
            table.add_row("output", output_path)
        elapsed = self.run.get("elapsed_s")
        if isinstance(elapsed, (int, float)):
            table.add_row("command elapsed", f"{float(elapsed):.3f} s")
        chemistry_elapsed = self.run.get("chemistry_elapsed_s")
        if isinstance(chemistry_elapsed, (int, float)):
            table.add_row(
                "program runtime", f"{float(chemistry_elapsed):.3f} s"
            )
        termination = self.run.get("normal_termination")
        if termination is not None:
            table.add_row(
                "chemistry termination",
                "normal" if termination else "not confirmed",
            )
        if status not in {"validating", "starting", "running", "completed"}:
            error = str(self.run.get("error") or "No diagnostic message.")
            table.add_row("error", "\n".join(error.splitlines()[-40:]))


__all__ = ["CalculationReceiptCell"]
