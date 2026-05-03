"""
Gaussian Minimum Energy Cross Point (MECP) job implementation.

This module provides the GaussianMECPJob class for performing
Minimum Energy Cross Point calculations using Gaussian.
"""

import logging
import os
from typing import Type

import numpy as np
from ase import units

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
    MIN_DIFF_GRAD_NORM_SQ = 1.0e-20

    # MECP-specific attribute names that must be stripped when building
    # a GaussianLinkJobSettings for each sub-job (broken-symmetry mode).
    _MECP_ONLY_KEYS = frozenset(
        {
            "multiplicity_a",
            "multiplicity_b",
            "charge_a",
            "charge_b",
            "title_a",
            "title_b",
            "max_steps",
            "step_size",
            "trust_radius",
            "energy_diff_tol",
            "force_max_tol",
            "force_rms_tol",
            "disp_max_tol",
            "disp_rms_tol",
            "adaptive_step_size",
            "step_size_method",
            "step_size_grow",
            "step_size_shrink",
            "step_size_min",
            "step_size_max",
            "use_link",
            "num_alpha_a",
            "num_beta_a",
            "num_alpha_b",
            "num_beta_b",
            # 'stable' and 'guess' are kept: they are valid GaussianLinkJobSettings
            # params and will be overridden with state-specific values anyway.
        }
    )

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
            return any(line.startswith("Converged at step") for line in f)

    def _build_guess_string(self, state):
        """
        Build the Gaussian ``guess=`` parameter for a given state.

        For plain MECP sub-jobs (``use_link=False``) this is not called.
        For broken-symmetry link sub-jobs the base guess (e.g. ``"mix"``)
        is extended with ``nalpha``/``nbeta`` specifiers when the
        corresponding ``num_alpha_*`` / ``num_beta_*`` settings are set:

        * No nalpha/nbeta → ``"mix"``
        * With nalpha/nbeta → ``"(mix,nalpha=N,nbeta=M)"``

        Args:
            state (str): ``"A"`` or ``"B"``.

        Returns:
            str: The guess string for the Gaussian route section.
        """
        base = self.settings.guess  # e.g. "mix"
        if state == "A":
            nalpha = self.settings.num_alpha_a
            nbeta = self.settings.num_beta_a
        else:
            nalpha = self.settings.num_alpha_b
            nbeta = self.settings.num_beta_b

        if nalpha is None and nbeta is None:
            return base

        parts = [base] if base else []
        if nalpha is not None:
            parts.append(f"nalpha={nalpha}")
        if nbeta is not None:
            parts.append(f"nbeta={nbeta}")
        return "(" + ",".join(parts) + ")"

    def _state_settings(self, charge, multiplicity, title, state="A"):
        """
        Build per-state calculation settings.

        In plain mode (``use_link=False``) returns a ``GaussianJobSettings``
        copy configured for a SP+forces run.

        In broken-symmetry link mode (``use_link=True``) returns a
        ``GaussianLinkJobSettings`` that writes a two-section Gaussian input:

        * **Section 1** – ``stable=opt`` with ``guess=(mix[,nalpha=N,nbeta=M])``
          to converge to the broken-symmetry wavefunction.
        * **Section 2** – SP + ``force`` with ``geom=check guess=read`` to
          compute the energy and gradients on the stable BS solution.

        Args:
            charge (int): Formal charge for this state.
            multiplicity (int): Spin multiplicity for this state.
            title (str): Title string for the calculation.
            state (str): ``"A"`` or ``"B"`` (used for nalpha/nbeta lookup).

        Returns:
            GaussianJobSettings | GaussianLinkJobSettings: State settings.
        """
        if self.settings.use_link:
            from chemsmart.jobs.gaussian.settings import GaussianLinkJobSettings

            # Start from the MECP settings dict, strip private attrs and
            # MECP-only keys that GaussianLinkJobSettings does not expect.
            link_kwargs = {
                k: v
                for k, v in self.settings.__dict__.items()
                if not k.startswith("_") and k not in self._MECP_ONLY_KEYS
            }
            # Override with state-specific and link-specific values.
            link_kwargs.update(
                {
                    "jobtype": "sp",
                    "freq": False,
                    "numfreq": False,
                    "forces": True,
                    "charge": charge,
                    "multiplicity": multiplicity,
                    "title": title,
                    "stable": self.settings.stable,
                    "guess": self._build_guess_string(state),
                    "link": True,
                }
            )
            return GaussianLinkJobSettings(**link_kwargs)

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
        mol.positions = positions_bohr * units.Bohr
        settings = self._state_settings(
            charge=charge,
            multiplicity=multiplicity,
            title=f"{title} step {step_idx}",
            state=state,
        )
        state_label = f"{self.label}_step{step_idx:03d}_{state}"

        if self.settings.use_link:
            from chemsmart.jobs.gaussian.link import GaussianLinkJob

            job = GaussianLinkJob(
                molecule=mol,
                settings=settings,
                label=state_label,
                jobrunner=self.jobrunner,
                skip_completed=False,
            )
        else:
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
        forces = np.array(output.forces[-1], dtype=float)
        if forces.shape != np.array(mol.positions).shape:
            raise RuntimeError(
                f"Parsed force shape {forces.shape} does not match "
                f"coordinate shape {np.array(mol.positions).shape}."
            )
        gradient = -forces
        return energy, gradient

    def _write_trajectory_frame(self, positions_bohr, step_idx):
        positions = positions_bohr * units.Bohr
        mode = "w" if step_idx == 0 else "a"
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

    def _mecp_displacement(self, energy_diff, grad_a, grad_b, step_size):
        diff_grad = grad_a - grad_b
        diff_norm_sq = float(np.sum(diff_grad * diff_grad))

        if diff_norm_sq < self.MIN_DIFF_GRAD_NORM_SQ:
            raise RuntimeError(
                "Difference gradient is too small; cannot continue MECP step."
            )

        # Move toward the linearized crossing seam.
        seam_correction = -(energy_diff / diff_norm_sq) * diff_grad

        # Minimize state A projected onto the crossing seam.
        proj = np.sum(grad_a * diff_grad) / diff_norm_sq
        projected_grad = grad_a - proj * diff_grad

        # Downhill step along the seam.
        downhill_step = -step_size * projected_grad

        displacement = seam_correction + downhill_step

        return displacement, projected_grad

    def _adapt_step_size(self, current_step_size, prev_merit, current_merit):
        """
        Return an updated step size based on the merit function progress.

        The merit is dimensionless: ``|ΔE|/energy_diff_tol + RMS(g_perp)/force_rms_tol``.
        If merit decreased (progress), the step size is grown by ``step_size_grow``.
        If merit increased (overshoot/oscillation), it is shrunk by ``step_size_shrink``.
        The result is clamped to ``[step_size_min, step_size_max]``.
        """
        if current_merit < prev_merit:
            new_step = current_step_size * self.settings.step_size_grow
        else:
            new_step = current_step_size * self.settings.step_size_shrink
        return float(
            np.clip(
                new_step,
                self.settings.step_size_min,
                self.settings.step_size_max,
            )
        )

    def _bb_step_size(
        self, prev_positions, curr_positions, prev_proj_grad, curr_proj_grad
    ):
        """
        Compute a Barzilai-Borwein (BB2) step size from the secant condition.

        Uses the formula ``α = ||Δr||² / (Δr · Δg_⊥)`` where
        ``Δr = r_n − r_{n−1}`` and ``Δg_⊥ = g_⊥,n − g_⊥,n−1``.

        Falls back to the initial ``step_size`` when the denominator is
        non-positive (negative curvature or numerically zero step).
        The result is clamped to ``[step_size_min, step_size_max]``.
        """
        delta_r = (curr_positions - prev_positions).ravel()
        delta_g = (curr_proj_grad - prev_proj_grad).ravel()

        r_dot_r = float(np.dot(delta_r, delta_r))
        r_dot_g = float(np.dot(delta_r, delta_g))

        if r_dot_r < 1e-30 or r_dot_g <= 0.0:
            return float(
                np.clip(
                    self.settings.step_size,
                    self.settings.step_size_min,
                    self.settings.step_size_max,
                )
            )

        bb_step = r_dot_r / r_dot_g
        return float(
            np.clip(
                bb_step,
                self.settings.step_size_min,
                self.settings.step_size_max,
            )
        )

    def _apply_trust_radius(self, displacement):
        """
        Apply a per-atom Cartesian trust radius.

        The trust radius is interpreted as the maximum allowed
        Cartesian displacement norm for each atom, in Bohr, per step.
        """
        step = np.array(displacement, dtype=float)

        atom_step_norms = np.linalg.norm(step, axis=1)
        exceed = atom_step_norms > self.settings.trust_radius

        if np.any(exceed):
            scale = self.settings.trust_radius / atom_step_norms[exceed]
            step[exceed] *= scale[:, None]

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

    def _log_step(
        self, f, step_idx, ea, eb, eff_grad, displacement, step_size
    ):
        grad_max = float(np.max(np.abs(eff_grad)))
        grad_rms = self._rms(eff_grad)
        disp_max = float(np.max(np.abs(displacement)))
        disp_rms = self._rms(displacement)
        f.write(
            f"step={step_idx:03d} E_A={ea:.10f} E_B={eb:.10f} "
            f"dE={ea - eb:+.10e} grad_max={grad_max:.3e} "
            f"grad_rms={grad_rms:.3e} disp_max={disp_max:.3e} "
            f"disp_rms={disp_rms:.3e} step_size={step_size:.4e}\n"
        )

    def _run(self, **kwargs):
        positions_bohr = (
            np.array(self.molecule.positions, dtype=float) / units.Bohr
        )
        displacement = np.zeros_like(positions_bohr)
        logger.info(
            f"Starting MECP optimization for {self.label} at position: "
            f"{positions_bohr} Bohr and displacement: {displacement}\n"
        )

        current_step_size = self.settings.step_size
        prev_merit = None
        prev_positions = None
        prev_proj_grad = None

        with open(self.report_file, "w") as report:
            report.write("CHEMSMART self-contained MECP optimization\n")
            report.write(
                f"max_steps={self.settings.max_steps} "
                f"step_size={self.settings.step_size} "
                f"trust_radius={self.settings.trust_radius} "
                f"adaptive_step_size={self.settings.adaptive_step_size} "
                f"step_size_method={self.settings.step_size_method}\n"
            )
            for step_idx in range(self.settings.max_steps):
                self._write_trajectory_frame(positions_bohr, step_idx)
                ea, grad_a = self._run_state(positions_bohr, step_idx, "A")
                eb, grad_b = self._run_state(positions_bohr, step_idx, "B")

                energy_diff = ea - eb

                displacement, projected_grad = self._mecp_displacement(
                    energy_diff=energy_diff,
                    grad_a=grad_a,
                    grad_b=grad_b,
                    step_size=current_step_size,
                )

                displacement = self._apply_trust_radius(displacement)

                self._log_step(
                    report,
                    step_idx,
                    ea,
                    eb,
                    projected_grad,
                    displacement,
                    current_step_size,
                )

                if self._is_converged(
                    energy_diff=energy_diff,
                    eff_grad=projected_grad,
                    displacement=displacement,
                ):
                    report.write(f"Converged at step {step_idx}.\n")
                    break

                if self.settings.adaptive_step_size:
                    if self.settings.step_size_method == "bb":
                        if prev_positions is not None:
                            current_step_size = self._bb_step_size(
                                prev_positions,
                                positions_bohr,
                                prev_proj_grad,
                                projected_grad,
                            )
                        prev_positions = positions_bohr.copy()
                        prev_proj_grad = projected_grad.copy()
                    else:  # "grow_shrink"
                        current_merit = (
                            abs(energy_diff) / self.settings.energy_diff_tol
                            + self._rms(projected_grad)
                            / self.settings.force_rms_tol
                        )
                        if prev_merit is not None:
                            current_step_size = self._adapt_step_size(
                                current_step_size, prev_merit, current_merit
                            )
                        prev_merit = current_merit

                positions_bohr = positions_bohr + displacement
            else:
                raise RuntimeError(
                    "MECP optimization did not converge within max_steps."
                )

        self.molecule.positions = positions_bohr * units.Bohr
