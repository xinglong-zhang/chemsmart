"""Job runner that calculates and reports Sterimol parameters in Python."""

import json
import logging

from chemsmart.jobs.runner import JobRunner
from chemsmart.utils.utils import get_range_from_list

logger = logging.getLogger(__name__)

RADIUS_SOURCE = "RDKit periodic table (van der Waals radii)"
LENGTH_CONVENTION = (
    "Measured along the axis from the first axis atom, without Verloop's "
    "+0.40 Angstrom correction"
)


class SterimolJobRunner(JobRunner):
    """Run a SterimolJob in the current process, with no external program."""

    JOBTYPES = ["sterimol"]
    PROGRAM = "Sterimol"
    FAKE = False
    SCRATCH = False

    @property
    def executable(self):
        """Sterimol jobs do not use an external executable."""
        return None

    def run(self, job, **kwargs):
        """Calculate the parameters, print them, and write the JSON report."""
        report = self._build_report(job)
        # The terminal locally, and the scheduler's output file on a node
        print(self._format_report(report), flush=True)
        with open(job.outputfile, "w") as f:
            json.dump(report, f, indent=2)
        print(f"Report:       {job.outputfile}", flush=True)

    @staticmethod
    def _build_report(job):
        """Calculate the parameters and return them with their metadata."""
        molecule = job.molecule
        keep = {job.atom1, *job.substituent_atoms}
        excluded = [
            idx for idx in range(1, molecule.num_atoms + 1) if idx not in keep
        ]
        result = molecule.calculate_sterimol_parameters(
            atom1=job.atom1, atom2=job.atom2, excluded_atoms=excluded
        )
        symbols = molecule.chemical_symbols
        return {
            "file": job.filename,
            "index": job.index,
            "formula": molecule.chemical_formula,
            "axis": {
                "atoms": [job.atom1, job.atom2],
                "symbols": [
                    str(symbols[job.atom1 - 1]),
                    str(symbols[job.atom2 - 1]),
                ],
            },
            "selection": job.selection,
            "substituent_atoms": job.substituent_atoms,
            "L": float(result["L"]),
            "B1": float(result["B1"]),
            "B5": float(result["B5"]),
            "units": "Angstrom",
            "radius_source": RADIUS_SOURCE,
            "length_convention": LENGTH_CONVENTION,
        }

    @staticmethod
    def _format_report(report):
        """Return the report as text."""
        atom1, atom2 = report["axis"]["atoms"]
        symbol1, symbol2 = report["axis"]["symbols"]
        atoms = ",".join(get_range_from_list(report["substituent_atoms"]))
        return "\n".join(
            [
                f"Sterimol parameters: {report['file']} ({report['formula']})",
                f"Axis:         {symbol1}{atom1} -> {symbol2}{atom2}",
                f"Substituent:  {atoms} ({report['selection']})",
                f"L  = {report['L']:.3f} Å",
                f"B1 = {report['B1']:.3f} Å",
                f"B5 = {report['B5']:.3f} Å",
                f"Radii:        {report['radius_source']}",
                f"L convention: {report['length_convention']}",
            ]
        )
