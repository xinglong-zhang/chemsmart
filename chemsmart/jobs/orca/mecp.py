"""
ORCA Minimum Energy Cross Point (MECP) job implementation.

This module contains the ORCAMECPJob class for running minimum energy
crossing-point optimizations using ORCA's native SurfCrossOpt feature.
"""

import logging
import os
from typing import Type

from chemsmart.jobs.orca.job import ORCAJob
from chemsmart.jobs.orca.settings import ORCAMECPJobSettings
from chemsmart.utils.periodictable import PeriodicTable

logger = logging.getLogger(__name__)


class ORCAMECPJob(ORCAJob):
    """
    ORCA Minimum Energy Cross Point (MECP) job.

    Wraps ORCA's ``SurfCrossOpt`` geometry optimisation on the crossing
    seam between two spin states of the **same charge** and using the
    **same level of theory**.  The two-state multiplicities are specified
    via the ``%mecp`` input block (PES2 multiplicity) and the ``* xyz``
    charge/multiplicity line (PES1 multiplicity).

    Attributes:
        TYPE (str): Job type identifier ('orcamecp').
        molecule: Molecule object used for the MECP optimization.
        settings: ORCAMECPJobSettings configuration for the job.
        label (str): Job identifier used for file naming.
        jobrunner: Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "orcamecp"
    DEFAULT_ENERGY_GAP_TOLERANCE = 1.0e-4

    @classmethod
    def settings_class(cls) -> Type[ORCAMECPJobSettings]:
        return ORCAMECPJobSettings

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        """
        Initialize ORCAMECPJob.

        Args:
            molecule: Molecule object for the MECP optimization.
            settings: ORCAMECPJobSettings instance.
            label: Job label for identification.
            jobrunner: Job runner instance.
            **kwargs: Additional keyword arguments.
        """
        settings = ORCAMECPJobSettings.from_settings(settings)
        settings.validate()
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
        if self.settings.mode == "numfreq" and len(self.molecule) < 3:
            raise ValueError(
                "ORCA SurfCrossNumFreq requires at least 3 atoms; "
                f"the supplied structure contains {len(self.molecule)}."
            )
        self._validate_broken_symmetry()

    def _validate_broken_symmetry(self):
        """Reject broken-symmetry requests incompatible with PES2.

        ORCA's ``brokenSym NA,NB`` describes two antiferromagnetically
        coupled centres carrying ``NA`` and ``NB`` unpaired electrons.  The
        resulting determinant has multiplicity ``abs(NA - NB) + 1``.
        """
        broken_sym = self.settings.broken_sym
        if broken_sym is None:
            return

        unpaired_a, unpaired_b = broken_sym
        broken_sym_multiplicity = abs(unpaired_a - unpaired_b) + 1
        if broken_sym_multiplicity != self.settings.multiplicity2:
            raise ValueError(
                f"brokenSym {unpaired_a},{unpaired_b} generates PES2 "
                "multiplicity "
                f"{broken_sym_multiplicity}, but --multiplicity2 is "
                f"{self.settings.multiplicity2}."
            )

        periodic_table = PeriodicTable()
        electron_count = (
            sum(
                periodic_table.to_atomic_number(symbol)
                for symbol in self.molecule.symbols
            )
            - self.settings.charge
        )
        if (electron_count - (broken_sym_multiplicity - 1)) % 2:
            raise ValueError(
                f"brokenSym {unpaired_a},{unpaired_b} generates "
                "multiplicity "
                f"{broken_sym_multiplicity}, which is incompatible with "
                f"the molecule's {electron_count} electrons."
            )

    @property
    def report_file(self):
        """Path to the concise, user-facing MECP quality report."""
        return os.path.join(self.folder, f"{self.label}_report.log")

    def log_result(self, energy_gap_tolerance=None):
        """Write a concise quality report for a completed ORCA MECP job.

        ORCA remains responsible for the optimization and frequency
        calculation.  This report collects the most important acceptance
        checks in one place so users do not need to inspect the large ORCA
        output and its auxiliary files manually.

        Returns:
            str | None: Report path, or ``None`` if no ORCA output exists.
        """
        output = self._output()
        if output is None:
            return None

        tolerance = (
            self.DEFAULT_ENERGY_GAP_TOLERANCE
            if energy_gap_tolerance is None
            else float(energy_gap_tolerance)
        )
        if tolerance <= 0:
            raise ValueError("MECP energy-gap tolerance must be positive.")

        result = output.mecp_result
        if result is None:
            return None

        issues = []
        if not result.normal_termination:
            issues.append("ORCA did not terminate normally")
        if not result.converged:
            issues.append("MECP geometry optimization did not converge")
        if result.energy_gap is None:
            issues.append("final two-state energy gap was not found")
        elif abs(result.energy_gap) > tolerance:
            issues.append(
                "final two-state energy gap exceeds the acceptance tolerance"
            )
        if result.numfreq_requested and not result.numfreq_completed:
            issues.append(
                "requested SurfCrossNumFreq calculation is incomplete"
            )
        elif result.numfreq_completed and result.is_minimum is False:
            issues.append("imaginary mode detected on the crossing hyperline")

        if not result.normal_termination or not result.converged:
            status = "FAILED"
        elif issues:
            status = "WARNING"
        else:
            status = "PASSED"

        def value_or_na(value, precision=12):
            return "N/A" if value is None else f"{value:.{precision}f}"

        gap_kcal = (
            None
            if result.energy_gap is None
            else abs(result.energy_gap) * 627.509474
        )
        lines = [
            "CHEMSMART ORCA MECP quality report",
            f"job={self.label}",
            f"status={status}",
            "",
            "Completion checks",
            f"normal_termination={result.normal_termination}",
            f"optimization_converged={result.converged}",
            "",
            "Final crossing-point energies (Hartree)",
            f"state_1_energy={value_or_na(result.state_1_energy)}",
            f"state_2_energy={value_or_na(result.state_2_energy)}",
            f"energy_gap={value_or_na(result.energy_gap)}",
            f"absolute_energy_gap={value_or_na(None if result.energy_gap is None else abs(result.energy_gap))}",
            f"absolute_energy_gap_kcal_mol={value_or_na(gap_kcal, 6)}",
            f"energy_gap_tolerance={tolerance:.8f}",
            f"energy_gap_accepted={result.energy_gap is not None and abs(result.energy_gap) <= tolerance}",
            "",
            "Crossing-hyperline frequency check",
            f"numfreq_requested={result.numfreq_requested}",
            f"numfreq_completed={result.numfreq_completed}",
            f"state_1_imaginary_frequencies_cm-1={list(result.state_1_imaginary_frequencies)}",
            f"state_2_imaginary_frequencies_cm-1={list(result.state_2_imaginary_frequencies)}",
            f"is_minimum={result.is_minimum}",
            "",
            "Assessment",
        ]
        if issues:
            lines.extend(f"- {issue}" for issue in issues)
        else:
            lines.append("- All requested MECP quality checks passed.")
        if status == "WARNING":
            lines.append(
                "- Review the ORCA output; refinement from the final geometry "
                "may be appropriate."
            )
        lines.append("")

        with open(self.report_file, "w", encoding="utf-8") as report:
            report.write("\n".join(lines))
        logger.info(f"Wrote ORCA MECP quality report: {self.report_file}")
        return self.report_file

    @property
    def results(self):
        """Return the parsed native SurfCrossOpt result, if available."""
        output = self._output()
        return None if output is None else output.mecp_result
