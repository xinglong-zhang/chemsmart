"""
ORCA Minimum Energy Cross Point (MECP) job implementation.

This module contains the ORCAMECPJob class for running minimum energy
crossing-point optimizations using ORCA's native SurfCrossOpt feature.
"""

import copy
import logging
import os
import re
from typing import Type

import numpy as np
from ase import units

from chemsmart.io.mecp import read_report_fields, state_energy_lines
from chemsmart.jobs.mecp import MECPSeamMixin
from chemsmart.jobs.orca.job import ORCAGeneralJob, ORCAJob
from chemsmart.jobs.orca.settings import ORCAJobSettings, ORCAMECPJobSettings
from chemsmart.utils.periodictable import PeriodicTable

logger = logging.getLogger(__name__)


class ORCAMECPJob(MECPSeamMixin, ORCAJob):
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
    STATE_OUTPUT_EXTENSION = ".out"
    DEFAULT_ENERGY_GAP_TOLERANCE = 1.0e-4

    @classmethod
    def settings_class(cls) -> Type[ORCAMECPJobSettings]:
        return ORCAMECPJobSettings

    def __init__(
        self,
        molecule,
        settings,
        label,
        jobrunner=None,
        is_seam_follow_branch=False,
        **kwargs,
    ):
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
        if settings.uses_seam_workflow:
            if self.molecule.frozen_atoms is not None or settings.modred:
                raise ValueError(
                    "ORCA seam analysis currently requires an unconstrained molecule"
                )
            if settings.input_string:
                raise ValueError(
                    "ORCA seam analysis requires structured settings"
                )
            for route in (
                settings.route_to_be_written,
                settings.additional_route_parameters,
            ):
                if route and any(char in route for char in ("\n", "%", "*")):
                    raise ValueError(
                        "ORCA seam analysis requires simple route keywords, not embedded blocks"
                    )
        self.is_seam_follow_branch = is_seam_follow_branch
        self._reset_execution_state()

    @property
    def seam_report_file(self):
        return self.final_report_file

    def _job_is_complete(self):
        if not self.settings.uses_seam_workflow:
            return super()._job_is_complete()
        report = read_report_fields(self.final_report_file)
        return (
            report.get("mode") == self.settings.mode
            and report.get("hessian_method") == self.hessian_method
            and report.get("status") == "COMPLETED"
            and os.path.isfile(self.frequency_file)
        )

    def _run(self, **kwargs):
        if not self.settings.uses_seam_workflow:
            return super()._run(**kwargs)
        self._reset_execution_state()
        # Remove an earlier success marker before starting a fresh workflow.
        if os.path.isfile(self.seam_report_file):
            os.remove(self.seam_report_file)
        if self.settings.mecp_numfreq_only:
            self._run_frequency_only()
            return
        self._run_native_optimization(**kwargs)
        result = self._run_seam_minimum_check(
            self.molecule.positions / units.Bohr
        )
        self._last_seam_result = result
        self._final_convergence_metrics = {
            "energy_diff": result["energy_diff"]
        }
        self._write_final_report(
            self._final_optimization_steps,
            result["energy_a"],
            result["energy_b"],
            result,
            frequency_only=self.settings.mecp_numfreq_only,
        )

    def _run_native_optimization(self, **kwargs):
        """Let ORCA optimize the crossing, then adopt its final geometry."""
        super()._run(**kwargs)
        native = self.results
        if (
            native is None
            or not native.normal_termination
            or not native.converged
        ):
            raise RuntimeError(
                "ORCA native MECP optimization did not converge"
            )
        if native.final_structure is None:
            raise RuntimeError("ORCA final MECP structure is missing")
        self.molecule.positions = native.final_structure.positions.copy()
        self._final_optimization_steps = len(native.energy_gap_history)

    def _final_report_legacy_fields(self, result, frequency_only):
        return {
            "energy_gap": f"{result['energy_diff']:+.12e}",
            "selected_branch": getattr(
                self, "_selected_seam_follow_branch", "none"
            ),
        }

    def verify_seam_minimum(self, *args, **kwargs):
        # Native convergence alone is not sufficient for projected analysis.
        # Macro steps intentionally keep a nonzero force along the progress mode.
        kwargs.setdefault(
            "require_stationary", kwargs.get("macro_step") is None
        )
        return super().verify_seam_minimum(*args, **kwargs)

    def _remove_checkpoint_set(self, checkpoint_tag):
        # Retain ORCA orbitals on disk for inspection and explicit reuse.
        for state in ("A", "B"):
            self._state_checkpoint_files.pop((checkpoint_tag, state), None)

    def _state_settings(self, state, checkpoint_tag, analytic):
        """Build a fixed-geometry state input and its orbital guess."""
        settings = ORCAJobSettings()
        for key in vars(settings):
            if hasattr(self.settings, key):
                setattr(
                    settings, key, copy.deepcopy(getattr(self.settings, key))
                )
        settings.jobtype = "sp"
        settings.freq = settings.numfreq = settings.forces = False
        settings.multiplicity = getattr(
            self.settings, f"multiplicity_{state.lower()}"
        )
        previous = self._state_checkpoint_files.get((checkpoint_tag, state))
        if analytic and (not previous or not os.path.isfile(previous)):
            raise RuntimeError(
                "ORCA analytic Hessian requires its reference GBW"
            )
        broken_sym = self.settings.broken_sym if state == "B" else None
        if broken_sym is not None:
            settings.scf_hftyp = "UHF"
            settings.scf_broken_sym = None
            if not previous or not os.path.isfile(previous):
                # ORCA's BrokenSym SCF starts from NA+NB unpaired electrons
                # and finishes on the requested BS surface before EnGrad.
                settings.multiplicity = sum(broken_sym) + 1
                settings.scf_broken_sym = list(broken_sym)
        # Remove native crossing/optimization tasks from the single-state job.
        route = settings.route_string
        route = re.sub(
            r"\b(?:SurfCrossOpt|SurfCrossNumFreq|TightOpt|LooseOpt|VeryTightOpt|Opt|Freq|AnFreq|NumFreq|EnGrad|SP)\b",
            "",
            route,
            flags=re.I,
        )
        route = " ".join(route.split()) + (
            " AnFreq" if analytic else " EnGrad"
        )
        if broken_sym is not None:
            # Explicitly keep an unrestricted BS singlet when reading orbitals.
            route = re.sub(r"\b(?:RHF|RKS|ROHF|ROKS)\b", "", route, flags=re.I)
        if previous and os.path.isfile(previous):
            settings.moinp = previous
            route += f' MORead\n%moinp "{os.path.basename(previous)}"'
        settings.route_to_be_written = route
        settings.additional_route_parameters = None
        return settings

    def _state_job(
        self, positions_bohr, step_idx, state, checkpoint_tag, settings
    ):
        """Prepare a fresh state job using shared MECP naming conventions."""
        molecule = self.molecule.copy()
        molecule.positions = np.asarray(positions_bohr) * units.Bohr
        folder = (
            self.frequency_type_folder
            if checkpoint_tag == "check"
            else self.optimization_folder
        )
        os.makedirs(folder, exist_ok=True)
        label = self.mecp_files.state_label(step_idx, state, checkpoint_tag)
        job = ORCAGeneralJob(
            molecule=molecule,
            settings=settings,
            label=label,
            skip_completed=False,
        )
        job.set_folder(folder)
        self._propagate_runner(self.jobrunner, job)
        # Stale results must never substitute for a failed/incomplete new run.
        for suffix in (".out", ".engrad", ".hess", ".gbw"):
            path = os.path.join(folder, label + suffix)
            if os.path.isfile(path):
                os.remove(path)
        return job

    def _read_state_gradient(self, job):
        """Validate ORCA's atomic units and atom ordering before using gradients."""
        from chemsmart.io.orca.output import ORCAEngradFile

        try:
            gradient = ORCAEngradFile(
                os.path.join(job.folder, job.label + ".engrad")
            )
        except (OSError, ValueError, IndexError, TypeError) as error:
            raise RuntimeError(
                f"Cannot read ORCA energy/gradient: {job.label}"
            ) from error
        if (
            gradient.num_atoms != len(job.molecule.symbols)
            or gradient.energy is None
            or gradient.gradient is None
            or not np.isfinite(gradient.energy)
            or not np.isfinite(gradient.gradient).all()
            or gradient.molecule is None
            or list(gradient.molecule.symbols) != list(job.molecule.symbols)
            or not np.allclose(
                gradient.molecule.positions,
                job.molecule.positions,
                atol=1e-6,
                rtol=0,
            )
        ):
            raise RuntimeError(
                f"Invalid ORCA energy/gradient or geometry: {job.label}"
            )
        return gradient.energy, gradient.gradient

    def _run_state(
        self,
        positions_bohr,
        step_idx,
        state,
        checkpoint_tag=None,
        analytic=False,
    ):
        from chemsmart.io.orca.hessian import read_hessian

        settings = self._state_settings(state, checkpoint_tag, analytic)
        job = self._state_job(
            positions_bohr, step_idx, state, checkpoint_tag, settings
        )
        job.run()
        if not job.is_complete():
            raise RuntimeError(
                f"ORCA state {state} did not terminate normally: {job.label}"
            )
        if os.path.isfile(job.gbwfile):
            self._state_checkpoint_files[(checkpoint_tag, state)] = job.gbwfile
        elif state == "B" and self.settings.broken_sym is not None:
            raise RuntimeError(
                f"Missing broken-symmetry state B orbitals: {job.gbwfile}"
            )
        if not analytic:
            return self._read_state_gradient(job)
        try:
            self._analytic_hessians[state] = read_hessian(
                os.path.join(job.folder, job.label + ".hess"),
                len(job.molecule.symbols),
            )
        except (OSError, ValueError) as error:
            raise RuntimeError(
                f"Cannot read ORCA analytic Hessian: {job.label}"
            ) from error
        return None, None

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
            f"label={self.label}",
            "backend=orca",
            f"status={status}",
            *state_energy_lines(
                result.state_1_energy, result.state_2_energy, result.energy_gap
            ),
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
