"""
Gaussian Minimum Energy Cross Point (MECP) job implementation.

This module provides the GaussianMECPJob class for performing
Minimum Energy Cross Point calculations using Gaussian.
"""

import logging
import os
from typing import Type

import numpy as np

from chemsmart.jobs.gaussian.job import GaussianGeneralJob, GaussianJob
from chemsmart.jobs.gaussian.settings import GaussianMECPJobSettings

logger = logging.getLogger(__name__)


class GaussianMECPJob(GaussianJob):
    """
    Gaussian job class for Minimum Energy Cross Point (MECP) calculations.

    Performs iterative molecular geometry optimization to find
    the minimum energy cross point between two spin multiplicities
    for a molecular species.

    Attributes:
        TYPE (str): Job type identifier ('g16mecp').
        molecule (Molecule): Molecular structure to optimize.
        settings (GaussianJobSettings): Calculation configuration options.
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "g16mecp"
    BOHR_TO_ANGSTROM = 0.529177210903
    MIN_DIFF_GRAD_NORM_SQ = 1.0e-20

    def __init__(
        self,
        molecule,
        settings,
        label,
        jobrunner=None,
        **kwargs,
    ):
        settings = GaussianMECPJobSettings.from_settings(settings)
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )

    @classmethod
    def settings_class(cls) -> Type[GaussianMECPJobSettings]:
        return GaussianMECPJobSettings

    @property
    def report_file(self):
        return os.path.join(self.folder, f"{self.label}_report.log")

    @property
    def trajectory_file(self):
        return os.path.join(self.folder, f"{self.label}_traj.xyz")

    def _job_is_complete(self):
        """Check MECP completion by looking for a 'Converged' marker in the report file."""
        if not os.path.isfile(self.report_file):
            return False
        with open(self.report_file, encoding="utf-8") as f:
            return any("Converged" in line for line in f)

    def _state_settings(self, charge, multiplicity, title):
        state_settings = self.settings.copy()
        state_settings.jobtype = "sp"
        state_settings.freq = False
        state_settings.numfreq = False
        state_settings.forces = True
        state_settings.charge = charge
        state_settings.multiplicity = multiplicity
        state_settings.title = title
        return state_settings

    def _run_state(self, positions_bohr, step_idx, state):
        if state == "A":
            charge = self.settings.charge_a
            multiplicity = self.settings.multiplicity_a
            title = self.settings.title_a
        else:
            charge = self.settings.charge_b
            multiplicity = self.settings.multiplicity_b
            title = self.settings.title_b

        mol = self.molecule.copy()
        mol.positions = positions_bohr * self.BOHR_TO_ANGSTROM
        settings = self._state_settings(
            charge=charge,
            multiplicity=multiplicity,
            title=f"{title} step {step_idx}",
        )
        state_label = f"{self.label}_step{step_idx:03d}_{state}"
        job = GaussianGeneralJob(
            molecule=mol,
            settings=settings,
            label=state_label,
            jobrunner=self.jobrunner,
            skip_completed=False,
        )
        job.run()
        output = job._output()
        if output is None or output.energies is None or not output.energies:
            raise RuntimeError(f"Failed to parse energy for {state_label}.")
        if output.forces is None or len(output.forces) == 0:
            raise RuntimeError(
                f"Failed to parse forces for {state_label}. "
                "MECP requires `force` calculations."
            )
        energy = output.energies[-1]
        gradient = -np.array(output.forces[-1], dtype=float)
        return energy, gradient

    def _write_trajectory_frame(self, positions_bohr, step_idx):
        positions = positions_bohr * self.BOHR_TO_ANGSTROM
        mode = "a" if os.path.exists(self.trajectory_file) else "w"
        with open(self.trajectory_file, mode) as f:
            f.write(f"{len(self.molecule.symbols)}\n")
            f.write(f"MECP step {step_idx}\n")
            for symbol, pos in zip(self.molecule.symbols, positions):
                f.write(
                    f"{symbol:>3s} {pos[0]:15.8f} {pos[1]:15.8f} "
                    f"{pos[2]:15.8f}\n"
                )

    @staticmethod
    def _rms(array):
        return float(np.sqrt(np.mean(np.square(array))))

    def _effective_gradient(self, energy_diff, grad_a, grad_b):
        diff_grad = grad_a - grad_b
        diff_norm_sq = float(np.sum(diff_grad * diff_grad))
        if diff_norm_sq < self.MIN_DIFF_GRAD_NORM_SQ:
            raise RuntimeError(
                "Difference gradient is too small; cannot continue MECP step."
            )
        f_component = (energy_diff / diff_norm_sq) * diff_grad
        proj = np.sum(grad_a * diff_grad) / diff_norm_sq
        g_component = grad_a - proj * diff_grad
        return f_component + g_component

    def _apply_trust_radius(self, displacement):
        step = np.array(displacement, dtype=float)
        norms = np.linalg.norm(step, axis=1)
        exceed = norms > self.settings.trust_radius
        if np.any(exceed):
            scale = self.settings.trust_radius / norms[exceed]
            step[exceed] = step[exceed] * scale[:, None]
        return step

    def _is_converged(self, energy_diff, eff_grad, displacement):
        grad_max = float(np.max(np.abs(eff_grad)))
        grad_rms = self._rms(eff_grad)
        disp_max = float(np.max(np.abs(displacement)))
        disp_rms = self._rms(displacement)
        return (
            abs(energy_diff) <= self.settings.energy_diff_tol
            and grad_max <= self.settings.force_max_tol
            and grad_rms <= self.settings.force_rms_tol
            and disp_max <= self.settings.disp_max_tol
            and disp_rms <= self.settings.disp_rms_tol
        )

    def _log_step(self, f, step_idx, ea, eb, eff_grad, displacement):
        grad_max = float(np.max(np.abs(eff_grad)))
        grad_rms = self._rms(eff_grad)
        disp_max = float(np.max(np.abs(displacement)))
        disp_rms = self._rms(displacement)
        f.write(
            f"step={step_idx:03d} E_A={ea:.10f} E_B={eb:.10f} "
            f"dE={ea - eb:+.10e} grad_max={grad_max:.3e} "
            f"grad_rms={grad_rms:.3e} disp_max={disp_max:.3e} "
            f"disp_rms={disp_rms:.3e}\n"
        )

    def _run(self, **kwargs):
        positions_bohr = (
            np.array(self.molecule.positions, dtype=float)
            / self.BOHR_TO_ANGSTROM
        )
        displacement = np.zeros_like(positions_bohr)
        logger.info(
            f"Starting MECP optimization for {self.label} at position: "
            f"{positions_bohr} Bohr and displacement: {displacement}\n"
        )

        with open(self.report_file, "w") as report:
            report.write("CHEMSMART self-contained MECP optimization\n")
            report.write(
                f"max_steps={self.settings.max_steps} step_size={self.settings.step_size} "
                f"trust_radius={self.settings.trust_radius}\n"
            )
            for step_idx in range(self.settings.max_steps):
                self._write_trajectory_frame(positions_bohr, step_idx)
                ea, grad_a = self._run_state(positions_bohr, step_idx, "A")
                eb, grad_b = self._run_state(positions_bohr, step_idx, "B")
                energy_diff = ea - eb
                eff_grad = self._effective_gradient(
                    energy_diff=energy_diff, grad_a=grad_a, grad_b=grad_b
                )
                displacement = -self.settings.step_size * eff_grad
                displacement = self._apply_trust_radius(displacement)
                self._log_step(
                    report, step_idx, ea, eb, eff_grad, displacement
                )
                if self._is_converged(
                    energy_diff=energy_diff,
                    eff_grad=eff_grad,
                    displacement=displacement,
                ):
                    report.write(f"Converged at step {step_idx}.\n")
                    break
                positions_bohr = positions_bohr + displacement
            else:
                raise RuntimeError(
                    "MECP optimization did not converge within max_steps."
                )

        self.molecule.positions = positions_bohr * self.BOHR_TO_ANGSTROM
