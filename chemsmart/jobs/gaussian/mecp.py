"""
Gaussian Minimum Energy Cross Point (MECP) job implementation.

This module provides the GaussianMECPJob class for performing
Minimum Energy Cross Point calculations using Gaussian.
"""

import logging
import os
import re
from functools import partial
from typing import Type

import numpy as np
from ase import units

from chemsmart.jobs.gaussian.job import GaussianGeneralJob, GaussianJob
from chemsmart.jobs.gaussian.settings import GaussianMECPJobSettings
from chemsmart.utils.constants import (
    MECP_FREQUENCY_HEADER,
    MECP_FREQUENCY_TERMINATION_MARKER,
    amu_to_kg,
    bohr_to_meter,
    hartree_to_joules,
)
from chemsmart.utils.utils import remove_word_from_parenthesized_option

logger = logging.getLogger(__name__)

_MIN_DIFF_GRAD_NORM_SQ = 1.0e-20

_GROW_SHRINK_THRESHOLDS = {
    "improvement": 0.10,
    "regression": 0.02,
}

_BB_SAFEGUARDS = {
    "curvature_cosine_min": 1.0e-4,
    "relative_step_min": 0.5,
    "relative_step_max": 2.0,
}

# Settings consumed by the MECP driver rather than Gaussian link sub-jobs.
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
        "harvey_initial_hessian",
        "harvey_max_component_step",
        "harvey_max_condition",
        "use_link",
        "convergence_preset",
        "mecp_numfreq",
        "hess_step_size",
        "follow_seam_imaginary_mode",
        "seam_mode_displacement",
        "seam_mode_max_steps",
        "restart",
    }
)


class GaussianMECPJob(GaussianJob):
    """
    Gaussian job class for Minimum Energy Cross Point (MECP) calculations.

    Performs iterative molecular geometry optimization to find
    the minimum energy cross point between two spin multiplicities
    for a molecular species.

    At each iteration (step 1, 2, …, max_steps) two Gaussian sub-jobs are
    executed, labelled ``<label>_step<N>_A`` and
    ``<label>_step<N>_B`` where ``<N>`` is the **1-indexed** step number,
    supporting up to 999 999 steps.

    Attributes:
        TYPE (str): Job type identifier ('g16mecp').
        molecule (Molecule): Molecular structure to optimize.
        settings (GaussianJobSettings): Calculation configuration options.
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "g16mecp"

    @staticmethod
    def _route_without_guess_read(route):
        """Remove ``read`` from a Gaussian guess option for the first step."""
        route = re.sub(
            r"\b(guess)\s*=\s*\(([^)]*)\)",
            partial(remove_word_from_parenthesized_option, word="read"),
            route,
            flags=re.IGNORECASE,
        )
        route = re.sub(r"\bguess\s*=\s*read\b", "", route, flags=re.IGNORECASE)
        return " ".join(route.split())

    @staticmethod
    def _with_required_nosymm(route):
        """Ensure Gaussian forces remain in the MECP Cartesian frame."""
        route = route or ""
        if re.search(r"\bnosymm(?:etry)?\b", route, flags=re.IGNORECASE):
            return route
        if re.search(r"\bsymm(?:etry)?\b", route, flags=re.IGNORECASE):
            raise ValueError(
                "MECP requires nosymm so Cartesian forces remain aligned with "
                "the optimization coordinates; remove the explicit symmetry option."
            )
        return f"{route} nosymm".strip()

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
        self._last_spin_squared = {"A": None, "B": None}

    @classmethod
    def settings_class(cls) -> Type[GaussianMECPJobSettings]:
        return GaussianMECPJobSettings

    @property
    def report_file(self):
        return os.path.join(
            self.optimization_folder, f"{self.label}_report.log"
        )

    @property
    def optimization_folder(self):
        return os.path.join(self.folder, f"{self.label}_optimization")

    @property
    def numfreq_folder(self):
        return os.path.join(self.folder, f"{self.label}_numfreq")

    @property
    def final_report_file(self):
        return os.path.join(self.folder, f"{self.label}_final_report.log")

    @property
    def trajectory_file(self):
        return os.path.join(
            self.optimization_folder, f"{self.label}_traj.xyz"
        )

    @property
    def state_file(self):
        return os.path.join(self.folder, f"{self.label}_state.npz")

    @staticmethod
    def _optional_array(value):
        return (
            np.array([], dtype=float) if value is None else np.asarray(value)
        )

    @staticmethod
    def _restore_optional_array(value):
        return None if value.size == 0 else value

    def _save_optimizer_state(self, **state):
        """Atomically persist enough optimizer state to resume the next step."""
        temporary_file = f"{self.state_file}.tmp.npz"
        np.savez(
            temporary_file,
            version=np.array(1, dtype=int),
            symbols=np.asarray(self.molecule.symbols, dtype="U4"),
            step_size_method=np.array(self.settings.step_size_method),
            next_step=np.array(state["next_step"], dtype=int),
            positions_bohr=np.asarray(state["positions_bohr"], dtype=float),
            current_step_size=np.array(
                state["current_step_size"], dtype=float
            ),
            prev_merit=np.array(
                np.nan if state["prev_merit"] is None else state["prev_merit"],
                dtype=float,
            ),
            prev_positions=self._optional_array(state["prev_positions"]),
            prev_proj_grad=self._optional_array(state["prev_proj_grad"]),
            inv_hessian=self._optional_array(state["inv_hessian"]),
            prev_eff_grad=self._optional_array(state["prev_eff_grad"]),
            prev_positions_bfgs=self._optional_array(
                state["prev_positions_bfgs"]
            ),
        )
        os.replace(temporary_file, self.state_file)

    def _load_optimizer_state(self):
        """Load and validate a previously persisted optimizer state."""
        with np.load(self.state_file, allow_pickle=False) as saved:
            symbols = saved["symbols"].tolist()
            method = str(saved["step_size_method"])
            positions = np.asarray(saved["positions_bohr"], dtype=float)
            if symbols != list(self.molecule.symbols):
                raise RuntimeError(
                    "Cannot restart MECP: atom symbols differ from saved state."
                )
            if method != self.settings.step_size_method:
                raise RuntimeError(
                    "Cannot restart MECP: optimizer method differs from saved state "
                    f"({method!r} != {self.settings.step_size_method!r})."
                )
            if positions.shape != np.asarray(self.molecule.positions).shape:
                raise RuntimeError(
                    "Cannot restart MECP: coordinate shape differs from saved state."
                )
            prev_merit = float(saved["prev_merit"])
            return {
                "next_step": int(saved["next_step"]),
                "positions_bohr": positions,
                "current_step_size": float(saved["current_step_size"]),
                "prev_merit": None if np.isnan(prev_merit) else prev_merit,
                "prev_positions": self._restore_optional_array(
                    saved["prev_positions"]
                ),
                "prev_proj_grad": self._restore_optional_array(
                    saved["prev_proj_grad"]
                ),
                "inv_hessian": self._restore_optional_array(
                    saved["inv_hessian"]
                ),
                "prev_eff_grad": self._restore_optional_array(
                    saved["prev_eff_grad"]
                ),
                "prev_positions_bfgs": self._restore_optional_array(
                    saved["prev_positions_bfgs"]
                ),
            }

    def _job_is_complete(self):
        """Check MECP completion by looking for a 'Converged' marker in the report file."""
        if not os.path.isfile(self.report_file):
            return False
        with open(self.report_file, encoding="utf-8") as f:
            converged = any(line.startswith("Converged at step") for line in f)
        if not converged:
            return False
        if self.settings.mecp_numfreq:
            if not os.path.isfile(
                os.path.join(
                    self.numfreq_folder, f"{self.label}_seam_check.log"
                )
            ):
                return False
        if self.settings.mecp_numfreq:
            return os.path.isfile(
                os.path.join(self.folder, f"{self.label}_mecp_freq.log")
            )
        return True

    def _state_settings(self, charge, multiplicity, title, state="A"):
        """
        Build per-state calculation settings.

        In plain mode (``use_link=False``) returns a ``GaussianJobSettings``
        copy configured for a SP+forces run.

        In broken-symmetry link mode (``use_link=True``) returns a
        ``GaussianLinkJobSettings`` that writes a two-section Gaussian input:

        * **Section 1** – ``stable=opt`` with ``guess=mix`` (or the value of
          ``settings.guess``) to converge to the broken-symmetry wavefunction.
          The number of α/β electrons is determined by the charge/multiplicity
          line, not by Guess options.
        * **Section 2** – SP + ``force`` with ``geom=check guess=read`` to
          compute the energy and gradients on the stable BS solution.

        Args:
            charge (int): Formal charge for this state.
            multiplicity (int): Spin multiplicity for this state.
            title (str): Title string for the calculation.
            state (str): ``"A"`` or ``"B"`` (unused in base implementation,
                retained for API compatibility with subclasses).

        Returns:
            GaussianJobSettings | GaussianLinkJobSettings: State settings.
        """
        if self.settings.use_link:
            from chemsmart.jobs.gaussian.settings import (
                GaussianLinkJobSettings,
            )

            # Start from the MECP settings dict, strip private attrs and
            # MECP-only keys that GaussianLinkJobSettings does not expect.
            link_kwargs = {
                k: v
                for k, v in self.settings.__dict__.items()
                if not k.startswith("_") and k not in _MECP_ONLY_KEYS
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
                    "guess": self.settings.guess,
                    "link": True,
                }
            )
            link_kwargs["additional_route_parameters"] = (
                self._with_required_nosymm(
                    link_kwargs.get("additional_route_parameters")
                )
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
        state_settings.additional_route_parameters = (
            self._with_required_nosymm(
                state_settings.additional_route_parameters
            )
        )
        return state_settings

    def _run_state(self, positions_bohr, step_idx, state, checkpoint_tag=None):
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
        if isinstance(step_idx, str):
            state_label = f"{self.label}_{step_idx}_{state}"
        else:
            job_tag = f"{checkpoint_tag}_" if checkpoint_tag else ""
            state_label = f"{self.label}_{job_tag}step{step_idx}_{state}"

        checkpoint_key = (checkpoint_tag, state)
        oldchkfile = self._state_checkpoint_files.get(checkpoint_key)
        if not oldchkfile and not self.settings.use_link:
            route = settings.additional_route_parameters or ""
            settings.additional_route_parameters = (
                self._route_without_guess_read(route)
            )

        checkpoint_part = f"{checkpoint_tag}_" if checkpoint_tag else ""
        job_kwargs = {
            "molecule": mol,
            "settings": settings,
            "label": state_label,
            "jobrunner": self.jobrunner,
            "skip_completed": False,
            "scratch_parent_folder": (
                f"{self.label}_numfreq"
                if checkpoint_tag == "check"
                else f"{self.label}_optimization"
            ),
            "checkpoint_filename": (
                f"{self.label}_{checkpoint_part}{state}.chk"
            ),
        }

        if self.settings.use_link:
            from chemsmart.jobs.gaussian.link import GaussianLinkJob

            job = GaussianLinkJob(**job_kwargs)
        else:
            job = GaussianGeneralJob(**job_kwargs)
        job_folder = (
            self.numfreq_folder
            if checkpoint_tag == "check"
            else self.steps_folder
        )
        os.makedirs(job_folder, exist_ok=True)
        job.set_folder(job_folder)
        job.oldchkfile = oldchkfile
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
        if not hasattr(self, "_last_spin_squared"):
            self._last_spin_squared = {"A": None, "B": None}
        self._last_spin_squared[state] = output.spin_squared_after_annihilation
        if os.path.isfile(job.chkfile):
            self._state_checkpoint_files[checkpoint_key] = job.chkfile
        return energy, gradient

    def _write_trajectory_frame(self, positions_bohr, step_idx):
        positions = positions_bohr * units.Bohr
        mode = "w" if step_idx == 1 else "a"
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

        if diff_norm_sq < _MIN_DIFF_GRAD_NORM_SQ:
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

        return displacement, projected_grad, seam_correction

    @staticmethod
    def _update_inverse_hessian(
        inv_hessian, delta_x, delta_g, return_diagnostics=False
    ):
        """Return Harvey/easyMECP's inverse-BFGS update.

        easyMECP applies the rank-three update even for negative curvature.
        This implementation does the same and only skips a numerically
        singular or non-finite update.
        """
        h_del_g = inv_hessian @ delta_g
        fac = float(np.dot(delta_g, delta_x))
        fae = float(np.dot(delta_g, h_del_g))
        scale = max(
            float(np.linalg.norm(delta_x) * np.linalg.norm(delta_g)), 1.0
        )
        singular_tol = np.finfo(float).eps * scale

        if abs(fac) <= singular_tol or abs(fae) <= singular_tol:
            result = (inv_hessian, "SKIP_SINGULAR", fac, fae)
            return result if return_diagnostics else result[0]

        w = delta_x / fac - h_del_g / fae
        updated = (
            inv_hessian
            + np.outer(delta_x, delta_x) / fac
            - np.outer(h_del_g, h_del_g) / fae
            + fae * np.outer(w, w)
        )
        if not np.all(np.isfinite(updated)):
            result = (inv_hessian, "SKIP_NONFINITE", fac, fae)
            return result if return_diagnostics else result[0]

        updated = 0.5 * (updated + updated.T)
        result = (updated, "UPDATE", fac, fae)
        return result if return_diagnostics else result[0]

    def _bfgs_displacement(
        self,
        ea,
        eb,
        grad_a,
        grad_b,
        prev_positions,
        curr_positions,
        prev_eff_grad,
        inv_hessian,
    ):
        """
        Compute the displacement using Harvey's original BFGS
        quasi-Newton method.

        This implements the MECP optimization algorithm of J. N. Harvey
        (2003):

        1. Compute the effective gradient
           ``G_eff = (Ea-Eb)*facPP*PerpG + facP*ParG``.
        2. Maintain and update the inverse Hessian using the BFGS formula.
        3. Compute the displacement as ``-HI @ G_eff``.
        4. Apply the STPMX displacement limit.

        The implementation follows the ``UpdateX`` subroutine in
        ``easymecp``'s ``MECP_FORTRAN`` code.

        Args:
            ea, eb: Energies of the two states.
            grad_a, grad_b: Gradients of the two states (Hartree/Bohr).
            prev_positions: Previous geometry (Bohr), or ``None`` on the
                first step.
            curr_positions: Current geometry (Bohr).
            prev_eff_grad: Previous effective gradient as a one-dimensional
                array, or ``None`` on the first step.
            inv_hessian: Current inverse Hessian (N x N), or ``None`` on the
                first step.

        Returns:
            displacement: Cartesian displacement (Bohr).
            par_grad: Gradient tangent to the seam, used for convergence.
            seam_correction: Seam-correction term used for logging.
            inv_hessian: Updated inverse Hessian.
            update_status, fac, fae: BFGS update diagnostics.
            eff_grad: Harvey effective gradient used for the next BFGS
                history update.
        """
        n = grad_a.size

        # 1. Compute the effective gradient following Harvey's original
        # Effective_Gradient subroutine.
        diff_grad = grad_a - grad_b  # PerpG
        diff_norm_sq = float(np.sum(diff_grad * diff_grad))
        if diff_norm_sq < _MIN_DIFF_GRAD_NORM_SQ:
            raise RuntimeError(
                "Difference gradient is too small; cannot continue MECP step."
            )
        diff_norm = np.sqrt(diff_norm_sq)
        pp = float(np.sum(grad_a * diff_grad)) / diff_norm
        par_grad = grad_a - diff_grad / diff_norm * pp  # ParG

        # facPP=140 is an empirical value that gives an inverse-Hessian
        # component of approximately 1/140 along PerpG. facP=1 uses the
        # BFGS-maintained inverse Hessian along ParG.
        fac_pp = 140.0
        fac_p = 1.0
        eff_grad = (ea - eb) * fac_pp * diff_grad + fac_p * par_grad
        eff_grad_flat = eff_grad.ravel().copy()

        # 2. Update the inverse Hessian with BFGS and compute the displacement.
        # Harvey's original code uses Angstrom, whereas ChemSmart uses Bohr.
        # Initial inverse Hessian: 0.7 Angstrom^2/Hartree converted to
        # Bohr^2/Hartree.
        bohr_per_ang = 1.0 / units.Bohr
        initial_hi_val = self.settings.harvey_initial_hessian * (
            bohr_per_ang**2
        )
        initial_inv_hessian = initial_hi_val * np.eye(n)

        if prev_positions is None or prev_eff_grad is None:
            # First step: use the diagonal inverse Hessian from Harvey's
            # Initialize subroutine.
            inv_hess = initial_inv_hessian
            displacement_flat = -inv_hess @ eff_grad_flat
            update_status = "INITIAL"
            fac = float("nan")
            fae = float("nan")
        else:
            # Apply the BFGS update from Harvey's UpdateX subroutine.
            delta_x = (curr_positions - prev_positions).ravel()  # DelX
            delta_g = eff_grad_flat - prev_eff_grad  # DelG

            inv_hess, update_status, fac, fae = self._update_inverse_hessian(
                inv_hessian,
                delta_x,
                delta_g,
                return_diagnostics=True,
            )

            condition_number = float(np.linalg.cond(inv_hess))
            if (
                not np.isfinite(condition_number)
                or condition_number > self.settings.harvey_max_condition
            ):
                inv_hess = initial_inv_hessian
                update_status = "RESET_ILL_CONDITIONED"

            displacement_flat = -inv_hess @ eff_grad_flat
            if not np.all(np.isfinite(displacement_flat)):
                raise RuntimeError("Harvey BFGS produced a non-finite step.")
            if float(np.dot(displacement_flat, eff_grad_flat)) >= 0.0:
                inv_hess = initial_inv_hessian
                displacement_flat = -inv_hess @ eff_grad_flat
                update_status = "RESET_NON_DESCENT"

        # 3. Limit the displacement using the STPMX logic from Harvey's
        # UpdateX subroutine.
        # STPMX = 0.1 Å -> 0.1 * (1/Bohr) Bohr
        stpmx = self.settings.harvey_max_component_step * bohr_per_ang

        lgstst = float(np.max(np.abs(displacement_flat)))
        if lgstst > stpmx:
            displacement_flat = displacement_flat / lgstst * stpmx

        displacement = displacement_flat.reshape(grad_a.shape)

        # Exact linear seam correction, retained for log comparisons.
        seam_correction = -(ea - eb) / diff_norm_sq * diff_grad

        return (
            displacement,
            par_grad,
            seam_correction,
            inv_hess,
            update_status,
            fac,
            fae,
            eff_grad,
        )

    def _adapt_step_size(self, current_step_size, prev_merit, current_merit):
        """
        Return an updated step size based on the merit function progress.

        The merit is dimensionless: ``|ΔE|/energy_diff_tol + RMS(g_perp)/force_rms_tol``.
        A relative dead band prevents numerical noise from repeatedly growing and
        shrinking the step. Only an improvement above 10% grows the step; a
        regression above 2% shrinks it; otherwise the step is retained.
        """
        merit_scale = max(abs(prev_merit), np.finfo(float).tiny)
        relative_progress = (prev_merit - current_merit) / merit_scale
        if relative_progress > _GROW_SHRINK_THRESHOLDS["improvement"]:
            new_step = current_step_size * self.settings.step_size_grow
        elif relative_progress < -_GROW_SHRINK_THRESHOLDS["regression"]:
            new_step = current_step_size * self.settings.step_size_shrink
        else:
            new_step = current_step_size
        return float(
            np.clip(
                new_step,
                self.settings.step_size_min,
                self.settings.step_size_max,
            )
        )

    def _bb_step_size(
        self,
        current_step_size,
        prev_positions,
        curr_positions,
        prev_proj_grad,
        curr_proj_grad,
        constraint_gradient=None,
    ):
        """
        Compute a safeguarded Barzilai-Borwein step from the secant condition.

        Uses the formula ``α = ||Δr||² / (Δr · Δg_⊥)`` where
        ``Δr = r_n − r_{n−1}`` and ``Δg_⊥ = g_⊥,n − g_⊥,n−1``.

        When supplied, ``constraint_gradient`` projects the displacement onto
        the current seam tangent before pairing it with the projected-gradient
        change. Unreliable curvature damps the current step instead of resetting
        it. A relative safeguard also prevents a valid but noisy secant pair
        from changing the step by more than a factor of two in one iteration.
        """
        delta_r = (curr_positions - prev_positions).ravel()
        delta_g = (curr_proj_grad - prev_proj_grad).ravel()

        if constraint_gradient is not None:
            normal = np.asarray(constraint_gradient, dtype=float).ravel()
            normal_norm_sq = float(np.dot(normal, normal))
            if normal_norm_sq > np.finfo(float).tiny:
                delta_r = (
                    delta_r
                    - float(np.dot(delta_r, normal)) / normal_norm_sq * normal
                )

        r_dot_r = float(np.dot(delta_r, delta_r))
        r_dot_g = float(np.dot(delta_r, delta_g))
        g_dot_g = float(np.dot(delta_g, delta_g))
        curvature_scale = np.sqrt(max(r_dot_r * g_dot_g, 0.0))
        reliable_curvature = (
            r_dot_r >= 1.0e-30
            and g_dot_g >= 1.0e-30
            and r_dot_g
            > _BB_SAFEGUARDS["curvature_cosine_min"] * curvature_scale
        )

        if reliable_curvature:
            candidate = r_dot_r / r_dot_g
        else:
            candidate = current_step_size * self.settings.step_size_shrink

        relative_min = current_step_size * _BB_SAFEGUARDS["relative_step_min"]
        relative_max = current_step_size * _BB_SAFEGUARDS["relative_step_max"]
        safeguarded_step = np.clip(candidate, relative_min, relative_max)
        return float(
            np.clip(
                safeguarded_step,
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

    def log_step(
        self,
        f,
        step_idx,
        ea,
        eb,
        projected_grad,
        displacement,
        seam_correction,
        step_size,
    ):
        energy_diff = ea - eb
        pgrad_max = float(np.max(np.abs(projected_grad)))
        pgrad_rms = self._rms(projected_grad)
        disp_max = float(np.max(np.abs(displacement)))
        disp_rms = self._rms(displacement)
        seam_max = float(np.max(np.abs(seam_correction)))
        seam_rms = self._rms(seam_correction)
        f.write(
            f"step={step_idx} E_A={ea:.10f} E_B={eb:.10f} "
            f"dE={energy_diff:+.6e} "
            f"pgrad_max={pgrad_max:.3e} pgrad_rms={pgrad_rms:.3e} "
            f"disp_max={disp_max:.3e} disp_rms={disp_rms:.3e} "
            f"seam_max={seam_max:.3e} seam_rms={seam_rms:.3e} "
            f"step_size={step_size:.4e}\n"
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
        # BFGS state variables for the Harvey method
        inv_hessian = None
        prev_eff_grad = None
        prev_positions_bfgs = None
        self.steps_folder = self.optimization_folder
        os.makedirs(self.steps_folder, exist_ok=True)
        self._state_checkpoint_files = {}
        self._last_spin_squared = {"A": None, "B": None}

        start_step = 1
        restarting = self.settings.restart and os.path.isfile(self.state_file)
        if restarting:
            saved = self._load_optimizer_state()
            start_step = saved["next_step"]
            positions_bohr = saved["positions_bohr"]
            current_step_size = saved["current_step_size"]
            prev_merit = saved["prev_merit"]
            prev_positions = saved["prev_positions"]
            prev_proj_grad = saved["prev_proj_grad"]
            inv_hessian = saved["inv_hessian"]
            prev_eff_grad = saved["prev_eff_grad"]
            prev_positions_bfgs = saved["prev_positions_bfgs"]
            for state in ("A", "B"):
                checkpoint_file = os.path.join(
                    self.steps_folder, f"{self.label}_{state}.chk"
                )
                if os.path.isfile(checkpoint_file):
                    self._state_checkpoint_files[(None, state)] = (
                        checkpoint_file
                    )

        report_mode = "a" if restarting else "w"
        converged_step = None
        with open(self.report_file, report_mode) as report:
            if restarting:
                report.write(
                    f"Restarting from saved state at step {start_step}.\n"
                )
            else:
                report.write("CHEMSMART self-contained MECP optimization\n")
            if self.settings.step_size_method == "harvey":
                report.write(
                    f"max_steps={self.settings.max_steps} "
                    "step_size_method=harvey "
                    "harvey_initial_hessian="
                    f"{self.settings.harvey_initial_hessian} "
                    "harvey_max_component_step="
                    f"{self.settings.harvey_max_component_step} "
                    f"harvey_max_condition={self.settings.harvey_max_condition}\n"
                )
            else:
                report.write(
                    f"max_steps={self.settings.max_steps} "
                    f"step_size={self.settings.step_size} "
                    f"trust_radius={self.settings.trust_radius} "
                    f"adaptive_step_size={self.settings.adaptive_step_size} "
                    f"step_size_method={self.settings.step_size_method}\n"
                )
            for step_idx in range(start_step, self.settings.max_steps + 1):
                self._write_trajectory_frame(positions_bohr, step_idx)
                ea, grad_a = self._run_state(positions_bohr, step_idx, "A")
                eb, grad_b = self._run_state(positions_bohr, step_idx, "B")

                energy_diff = ea - eb

                if self.settings.step_size_method == "harvey":
                    (
                        displacement,
                        projected_grad,
                        seam_correction,
                        inv_hessian,
                        bfgs_status,
                        bfgs_fac,
                        bfgs_fae,
                        optimizer_gradient,
                    ) = self._bfgs_displacement(
                        ea=ea,
                        eb=eb,
                        grad_a=grad_a,
                        grad_b=grad_b,
                        prev_positions=prev_positions_bfgs,
                        curr_positions=positions_bohr,
                        prev_eff_grad=prev_eff_grad,
                        inv_hessian=inv_hessian,
                    )
                else:
                    displacement, projected_grad, seam_correction = (
                        self._mecp_displacement(
                            energy_diff=energy_diff,
                            grad_a=grad_a,
                            grad_b=grad_b,
                            step_size=current_step_size,
                        )
                    )
                    optimizer_gradient = projected_grad

                if self.settings.step_size_method != "harvey":
                    displacement = self._apply_trust_radius(displacement)

                self.log_step(
                    report,
                    step_idx,
                    ea,
                    eb,
                    projected_grad,
                    displacement,
                    seam_correction,
                    current_step_size,
                )
                if self.settings.step_size_method == "harvey":
                    report.write(
                        f"bfgs_status={bfgs_status} "
                        f"delta_g_dot_delta_x={bfgs_fac:+.8e} "
                        f"delta_g_dot_h_delta_g={bfgs_fae:+.8e}\n"
                    )
                spin_a = self._last_spin_squared["A"]
                spin_b = self._last_spin_squared["B"]
                if spin_a is not None or spin_b is not None:
                    value_a = "NA" if spin_a is None else f"{spin_a:.6f}"
                    value_b = "NA" if spin_b is None else f"{spin_b:.6f}"
                    report.write(
                        f"spin_squared_A={value_a} spin_squared_B={value_b}\n"
                    )

                if self._is_converged(
                    energy_diff=energy_diff,
                    eff_grad=projected_grad,
                    displacement=displacement,
                ):
                    self._final_optimization_steps = step_idx
                    self._final_convergence_metrics = {
                        "energy_diff": energy_diff,
                        "pgrad_max": float(np.max(np.abs(projected_grad))),
                        "pgrad_rms": self._rms(projected_grad),
                        "disp_max": float(np.max(np.abs(displacement))),
                        "disp_rms": self._rms(displacement),
                    }
                    report.write(
                        f"Optimization converged at step {step_idx}.\n"
                    )
                    converged_step = step_idx
                    break

                if self.settings.adaptive_step_size:
                    if self.settings.step_size_method == "bb":
                        if prev_positions is not None:
                            current_step_size = self._bb_step_size(
                                current_step_size,
                                prev_positions,
                                positions_bohr,
                                prev_proj_grad,
                                projected_grad,
                                grad_a - grad_b,
                            )
                        prev_positions = positions_bohr.copy()
                        prev_proj_grad = projected_grad.copy()
                    elif self.settings.step_size_method == "grow_shrink":
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

                if self.settings.step_size_method == "harvey":
                    prev_positions_bfgs = positions_bohr.copy()
                    prev_eff_grad = optimizer_gradient.ravel().copy()

                positions_bohr = positions_bohr + displacement
                self._save_optimizer_state(
                    next_step=step_idx + 1,
                    positions_bohr=positions_bohr,
                    current_step_size=current_step_size,
                    prev_merit=prev_merit,
                    prev_positions=prev_positions,
                    prev_proj_grad=prev_proj_grad,
                    inv_hessian=inv_hessian,
                    prev_eff_grad=prev_eff_grad,
                    prev_positions_bfgs=prev_positions_bfgs,
                )
            else:
                raise RuntimeError(
                    "MECP optimization did not converge within max_steps."
                )

        self.molecule.positions = positions_bohr * units.Bohr

        seam_result = None
        if self.settings.mecp_numfreq:
            try:
                seam_result = self._run_seam_minimum_check(positions_bohr)
            except Exception as error:
                with open(self.report_file, "a", encoding="utf-8") as report:
                    report.write(
                        f"Initial MECP optimization converged at step "
                        f"{converged_step}.\n"
                        f"Final status: FAILED SEAM VERIFICATION "
                        f"({type(error).__name__}: {error})\n"
                    )
                raise

        with open(self.report_file, "a", encoding="utf-8") as report:
            report.write(
                f"Initial MECP optimization converged at step "
                f"{converged_step}.\n"
            )
            if seam_result is not None:
                initial = self._initial_seam_result
                report.write(
                    "Initial seam status: "
                    f"{'MINIMUM' if initial['is_minimum'] else 'SADDLE'}; "
                    f"significant imaginary modes={initial['n_negative']}\n"
                )
                if not initial["is_minimum"] and "frequencies" in initial:
                    report.write(
                        "Initial lowest projected frequency="
                        f"{np.min(initial['frequencies']):+.6f} cm^-1\n"
                    )
                if not initial["is_minimum"]:
                    report.write(
                        "Seam following: completed; selected branch="
                        f"{self._selected_seam_follow_branch}\n"
                        "Seam-following macro steps="
                        f"{self._selected_seam_follow_macro_steps}\n"
                    )
                report.write(
                    f"Final MECP energy={seam_result['mecp_energy']:+.12f} "
                    "Hartree\n"
                    f"Final energy gap={seam_result['energy_diff']:+.6e} "
                    "Hartree\n"
                    "Final significant imaginary modes="
                    f"{seam_result['n_negative']}\n"
                    "Final status: VERIFIED MECP MINIMUM\n"
                )
                if self.settings.mecp_numfreq:
                    report.write(
                        f"Final geometry and frequencies: "
                        f"{self.label}_mecp_freq.log\n"
                    )
            report.write(f"Converged at step {converged_step}.\n")
        if not getattr(self, "_is_seam_follow_branch", False):
            self._write_final_report(
                initial_steps=converged_step,
                initial_energy_a=ea,
                initial_energy_b=eb,
                seam_result=seam_result,
            )
        if os.path.isfile(self.state_file):
            os.remove(self.state_file)

    def _write_final_report(
        self, initial_steps, initial_energy_a, initial_energy_b, seam_result
    ):
        """Summarize only the selected final structure and its convergence."""
        result = seam_result or {}
        metrics = result.get(
            "convergence_metrics", self._final_convergence_metrics
        )
        energy_a = result.get("energy_a", initial_energy_a)
        energy_b = result.get("energy_b", initial_energy_b)
        positions = np.asarray(
            result.get("positions_angstrom", self.molecule.positions),
            dtype=float,
        )
        thresholds = result.get("convergence_thresholds") or {
            "energy_diff": self.settings.energy_diff_tol,
            "pgrad_max": self.settings.force_max_tol,
            "pgrad_rms": self.settings.force_rms_tol,
            "disp_max": self.settings.disp_max_tol,
            "disp_rms": self.settings.disp_rms_tol,
        }
        with open(self.final_report_file, "w", encoding="utf-8") as report:
            report.write("CHEMSMART final MECP result\n")
            report.write(f"label={self.label}\n")
            report.write(f"initial_optimization_steps={initial_steps}\n")
            macro_steps = getattr(
                self, "_selected_seam_follow_macro_steps", 0
            )
            report.write(f"seam_follow_macro_steps={macro_steps}\n")
            if macro_steps:
                report.write(
                    "seam_follow_selected_branch="
                    f"{self._selected_seam_follow_branch}\n"
                )
                report.write(
                    "final_branch_optimization_steps="
                    f"{result.get('optimization_steps', 'NA')}\n"
                )
            report.write(f"energy_A={energy_a:+.12f} Hartree\n")
            report.write(f"energy_B={energy_b:+.12f} Hartree\n")
            report.write(
                f"mecp_energy={0.5 * (energy_a + energy_b):+.12f} "
                "Hartree\n"
            )
            report.write("\nFinal convergence criteria:\n")
            for name, threshold in thresholds.items():
                value = metrics[name]
                if name == "energy_diff":
                    value = abs(value)
                unit = "Hartree" if name == "energy_diff" else (
                    "Hartree/Bohr" if name.startswith("pgrad") else "Bohr"
                )
                report.write(
                    f"{name}: value={value:.6e} "
                    f"threshold={threshold:.6e} {unit} "
                    f"status={'PASS' if value <= threshold else 'FAIL'}\n"
                )
            if seam_result is not None:
                report.write(
                    "\nseam_minimum="
                    f"{'PASS' if result.get('is_minimum', True) else 'FAIL'}\n"
                    f"significant_imaginary_modes={result['n_negative']}\n"
                    f"frequency_file={self.label}_mecp_freq.log\n"
                )
            else:
                report.write("\nseam_minimum=NOT_CHECKED\n")
            report.write("\nFinal geometry (Angstrom):\n")
            for symbol, position in zip(self.molecule.symbols, positions):
                report.write(
                    f"{symbol:>3s} {position[0]:+16.8f} "
                    f"{position[1]:+16.8f} {position[2]:+16.8f}\n"
                )
            if "frequencies" in result:
                report.write("\nProjected MECP frequencies (cm^-1):\n")
                for index, frequency in enumerate(result["frequencies"], 1):
                    report.write(f"mode {index:4d}: {frequency:+14.6f}\n")

    # ------------------------------------------------------------------
    # Seam-minimum verification via effective Hessian analysis
    # ------------------------------------------------------------------

    def _build_projection_vectors(self, positions_bohr, diff_grad):
        """
        Return an orthonormal basis for the subspace to be removed from
        the Hessian: translations (3), rotations (up to 3), and the
        unit gradient-difference direction.

        The vectors are returned as a list of 1-D numpy arrays of
        length ``3*n_atoms``.  Linear dependencies (e.g. linear
        molecules have only 2 rotational modes) are handled by a
        Gram-Schmidt sweep that discards near-zero vectors.

        Args:
            positions_bohr (np.ndarray): Current geometry, shape (n_atoms, 3),
                in Bohr.
            diff_grad (np.ndarray): Gradient difference :math:`\\nabla E_A -
                \\nabla E_B`, shape (n_atoms, 3) or (3*n_atoms,).

        Returns:
            list[np.ndarray]: Orthonormal projection vectors.
        """
        n_atoms = len(self.molecule.symbols)
        n = 3 * n_atoms
        r = np.array(positions_bohr, dtype=float).reshape(n_atoms, 3)

        raw = []

        # Translational modes: unit displacement of all atoms along x/y/z.
        for k in range(3):
            v = np.zeros(n)
            v[k::3] = 1.0
            raw.append(v)

        # Rotational modes: r_i × e_k for each Cartesian axis.
        rot_axes = [
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([0.0, 0.0, 1.0]),
        ]
        for axis in rot_axes:
            v = np.zeros((n_atoms, 3))
            for i in range(n_atoms):
                v[i] = np.cross(r[i], axis)
            raw.append(v.ravel())

        # Gradient-difference direction.
        gd = np.array(diff_grad, dtype=float).ravel()
        raw.append(gd)

        # Gram-Schmidt orthonormalization; skip vectors that are linearly
        # dependent (norm < 1e-10 after projection).
        orthonormal = []
        for v in raw:
            v = v.astype(float)
            for u in orthonormal:
                v = v - np.dot(v, u) * u
            nv = np.linalg.norm(v)
            if nv > 1e-10:
                orthonormal.append(v / nv)

        return orthonormal

    @staticmethod
    def _reduced_hessian(hessian, projection_vectors):
        """Return the Hessian represented in the unprojected subspace.

        Diagonalising ``P @ H @ P`` leaves one numerical zero eigenvalue for
        every projected vector.  Removing a fixed number of eigenvalues after
        sorting is unsafe because genuine negative eigenvalues sort before the
        projected zeros.  Building an explicit orthonormal complement avoids
        that ambiguity.
        """
        if not projection_vectors:
            return np.array(hessian, dtype=float, copy=True)
        projected_basis = np.column_stack(projection_vectors)
        q_full, _ = np.linalg.qr(projected_basis, mode="complete")
        seam_basis = q_full[:, projected_basis.shape[1] :]
        return seam_basis.T @ hessian @ seam_basis

    @staticmethod
    def _frequency_from_mass_weighted_eigenvalue(eigenvalue):
        """Convert Hartree/(Bohr² amu) to a signed wavenumber in cm⁻¹."""
        angular_frequency = np.sqrt(
            abs(eigenvalue)
            * hartree_to_joules
            / (bohr_to_meter**2 * amu_to_kg)
        )
        wavenumber = angular_frequency / (2.0 * np.pi * units._c * 100.0)
        return float(np.copysign(wavenumber, eigenvalue))

    def _projected_frequencies_and_modes(
        self, hessian, positions_bohr, diff_grad
    ):
        """Return mass-weighted MECP frequencies and Cartesian normal modes.

        Translation, rotation, and the gradient-difference direction are
        removed in mass-weighted coordinates. Frequencies are returned as
        signed wavenumbers in cm⁻¹; a negative value denotes an imaginary
        mode on the crossing seam.
        """
        positions = np.asarray(positions_bohr, dtype=float)
        masses = np.asarray(self.molecule.most_abundant_masses, dtype=float)
        sqrt_masses = np.sqrt(masses)
        mass_vector = np.repeat(masses, 3)
        sqrt_mass_vector = np.sqrt(mass_vector)

        mass_weighted_hessian = hessian / np.sqrt(
            np.outer(mass_vector, mass_vector)
        )

        center_of_mass = np.average(positions, axis=0, weights=masses)
        centered = positions - center_of_mass
        raw_projection_vectors = []

        for axis_index in range(3):
            translation = np.zeros_like(positions)
            translation[:, axis_index] = sqrt_masses
            raw_projection_vectors.append(translation.ravel())

        for axis in np.eye(3):
            rotation = np.cross(centered, axis) * sqrt_masses[:, None]
            raw_projection_vectors.append(rotation.ravel())

        mass_weighted_diff_grad = (
            np.asarray(diff_grad, dtype=float).ravel() / sqrt_mass_vector
        )
        raw_projection_vectors.append(mass_weighted_diff_grad)

        projection_vectors = []
        for vector in raw_projection_vectors:
            vector = np.asarray(vector, dtype=float)
            for basis_vector in projection_vectors:
                vector -= np.dot(vector, basis_vector) * basis_vector
            norm = np.linalg.norm(vector)
            if norm > 1.0e-10:
                projection_vectors.append(vector / norm)

        if projection_vectors:
            projected_basis = np.column_stack(projection_vectors)
            full_basis, _ = np.linalg.qr(projected_basis, mode="complete")
            seam_basis = full_basis[:, projected_basis.shape[1] :]
        else:
            seam_basis = np.eye(hessian.shape[0])

        reduced_hessian = (
            seam_basis.T @ mass_weighted_hessian @ seam_basis
        )
        eigenvalues, reduced_modes = np.linalg.eigh(reduced_hessian)
        mass_weighted_modes = seam_basis @ reduced_modes
        cartesian_modes = mass_weighted_modes / sqrt_mass_vector[:, None]
        mode_norms = np.linalg.norm(cartesian_modes, axis=0)
        cartesian_modes /= mode_norms

        frequencies = np.array(
            [
                self._frequency_from_mass_weighted_eigenvalue(value)
                for value in eigenvalues
            ]
        )
        return frequencies, cartesian_modes.T, eigenvalues, len(
            projection_vectors
        )

    @classmethod
    def _lagrangian_hessian(cls, hessian_a, hessian_b, grad_a, grad_b):
        """Return the constrained MECP Lagrangian Hessian and multiplier."""
        diff_grad = np.asarray(grad_a).ravel() - np.asarray(grad_b).ravel()
        diff_norm_sq = float(np.dot(diff_grad, diff_grad))
        if diff_norm_sq < _MIN_DIFF_GRAD_NORM_SQ:
            raise RuntimeError(
                "Difference gradient is too small for seam verification."
            )
        multiplier = float(np.dot(np.asarray(grad_a).ravel(), diff_grad))
        multiplier /= diff_norm_sq
        hessian = (1.0 - multiplier) * hessian_a + multiplier * hessian_b
        return hessian, multiplier

    def _compute_numerical_hessian(
        self, positions_bohr, h, step_prefix, macro_step=None
    ):
        """
        Compute both state Hessians numerically via central finite differences
        of the Cartesian forces.

        For each coordinate ``j`` (0 … 3N−1) the geometry is displaced by
        ``±h`` Bohr and both states are evaluated:

        .. math::

            H_{ij} \\approx \\frac{g_i(+h_j) - g_i(-h_j)}{2h}

        where :math:`g_i` denotes the gradient component (force negated).
        Each Hessian is symmetrised as :math:`(H + H^T)/2`.

        .. warning::

            This method runs **4 × 3N** Gaussian sub-jobs (2 displacements
            × 2 spin states × 3N coordinates), which is expensive for large
            molecules (120 jobs for a 10-atom molecule).  Use the
            ``hess_step_size`` setting to control accuracy vs. cost.

        Args:
            positions_bohr (np.ndarray): Current geometry, shape (n_atoms, 3).
            h (float): Finite-difference step size in Bohr.
            step_prefix (int): First check-specific sub-job step index.

        Returns:
            tuple[np.ndarray, np.ndarray]: Symmetrised ``(H_A, H_B)`` matrices,
            each with shape (3N, 3N), in Hartree/Bohr².
        """
        n_atoms = len(self.molecule.symbols)
        n = 3 * n_atoms
        H_A = np.zeros((n, n))
        H_B = np.zeros((n, n))
        pos = np.array(positions_bohr, dtype=float).reshape(n_atoms, 3)

        for j in range(n):
            atom_idx = j // 3
            coord_idx = j % 3

            pos_plus = pos.copy()
            pos_plus[atom_idx, coord_idx] += h

            pos_minus = pos.copy()
            pos_minus[atom_idx, coord_idx] -= h

            if macro_step is None:
                step_p = step_prefix + 2 * j
                step_m = step_prefix + 2 * j + 1
            else:
                coord = f"macro{macro_step:02d}_check_coord{j + 1:02d}"
                step_p = f"{coord}_plus"
                step_m = f"{coord}_minus"

            _, g_A_plus = self._run_state(
                pos_plus, step_p, "A", checkpoint_tag="check"
            )
            _, g_B_plus = self._run_state(
                pos_plus, step_p, "B", checkpoint_tag="check"
            )
            _, g_A_minus = self._run_state(
                pos_minus, step_m, "A", checkpoint_tag="check"
            )
            _, g_B_minus = self._run_state(
                pos_minus, step_m, "B", checkpoint_tag="check"
            )

            H_A[:, j] = (g_A_plus.ravel() - g_A_minus.ravel()) / (2 * h)
            H_B[:, j] = (g_B_plus.ravel() - g_B_minus.ravel()) / (2 * h)

        H_A = (H_A + H_A.T) / 2
        H_B = (H_B + H_B.T) / 2
        return H_A, H_B

    def verify_seam_minimum(
        self,
        h=None,
        step_prefix=1,
        write_frequencies=False,
        macro_step=None,
    ):
        """
        Verify that the current MECP geometry is a **minimum on the crossing
        seam**, not merely a crossing point.

        The method computes both state Hessians numerically and forms the
        constrained Lagrangian Hessian. It then represents that Hessian in the
        subspace orthogonal to translations, rotations, and the
        gradient-difference direction:

        .. math::

            H_\\text{seam} = Q^T [(1-\\lambda)H_A + \\lambda H_B] Q

        The columns of the reduced-space basis span the seam tangent space.
        Any negative eigenvalue indicates a saddle point on the seam.

        Results are written to ``<label>_seam_check.log``.

        .. note::

            This analysis requires **4 × 3N** additional Gaussian sub-jobs
            (see :meth:`_compute_numerical_hessian`).  Call it only on the
            converged geometry.  Trigger automatically via the
            ``--mecp-numfreq`` CLI flag.

        Args:
            h (float, optional): Finite-difference step size in Bohr.
                Defaults to ``settings.hess_step_size`` (1 × 10⁻³ Bohr).
            step_prefix (int, optional): Starting check-specific sub-job step
                index (default: 1). The ``check`` name component prevents
                clashes with normal MECP steps.
            write_frequencies (bool, optional): Write mass-weighted projected
                frequencies and normal modes in addition to the seam check.

        Returns:
            dict: Keys ``"eigenvalues"`` (1-D array, non-projected modes),
            ``"n_negative"`` (int), ``"is_minimum"`` (bool),
            ``"energy_diff"`` (float, Hartree),
            ``"n_projected"`` (int, modes removed), and
            ``"lagrange_multiplier"`` (float).
        """
        if h is None:
            h = self.settings.hess_step_size

        positions_bohr = (
            np.array(self.molecule.positions, dtype=float) / units.Bohr
        )

        logger.info(
            f"verify_seam_minimum: computing gradient difference at {self.label}"
        )
        reference_step = (
            step_prefix
            if macro_step is None
            else f"macro{macro_step:02d}_check_reference"
        )
        ea, grad_a = self._run_state(
            positions_bohr, reference_step, "A", checkpoint_tag="check"
        )
        eb, grad_b = self._run_state(
            positions_bohr, reference_step, "B", checkpoint_tag="check"
        )
        diff_grad = grad_a - grad_b

        logger.info(
            f"verify_seam_minimum: computing numerical Hessian for {self.label} "
            f"(h={h} Bohr, {4 * 3 * len(self.molecule.symbols)} sub-jobs)"
        )
        H_A, H_B = self._compute_numerical_hessian(
            positions_bohr,
            h=h,
            step_prefix=step_prefix + 1,
            macro_step=macro_step,
        )

        H_lagrangian, lagrange_multiplier = self._lagrangian_hessian(
            H_A, H_B, grad_a, grad_b
        )
        proj_vecs = self._build_projection_vectors(positions_bohr, diff_grad)
        H_seam = self._reduced_hessian(H_lagrangian, proj_vecs)

        non_zero_evals = np.sort(np.linalg.eigvalsh(H_seam))
        n_proj = len(proj_vecs)
        # Use the same mass-weighted curvature tolerance for both CLI modes.
        # Raw Cartesian and mass-weighted eigenvalues have different units.
        frequencies, modes, frequency_eigenvalues, frequency_n_proj = (
            self._projected_frequencies_and_modes(
                H_lagrangian, positions_bohr, diff_grad
            )
        )
        n_negative = int(np.sum(frequency_eigenvalues < -1.0e-6))

        result = {
            "eigenvalues": non_zero_evals,
            "n_negative": n_negative,
            "is_minimum": n_negative == 0,
            "energy_a": ea,
            "energy_b": eb,
            "mecp_energy": 0.5 * (ea + eb),
            "energy_diff": ea - eb,
            "n_projected": n_proj,
            "lagrange_multiplier": lagrange_multiplier,
            "positions_angstrom": positions_bohr * units.Bohr,
            "atomic_masses": np.asarray(
                self.molecule.most_abundant_masses, dtype=float
            ),
            "multiplicity_a": self.settings.multiplicity_a,
            "multiplicity_b": self.settings.multiplicity_b,
            # Molecular point-group detection is not currently available on
            # Molecule. Use the conservative C1 symmetry number by default.
            "rotational_symmetry_number": 1,
        }

        if write_frequencies:
            result.update(
                frequency_eigenvalues=frequency_eigenvalues,
                frequencies=frequencies,
                modes=modes,
                n_projected=frequency_n_proj,
            )

        self._write_seam_check_log(result, h, step_prefix)
        if write_frequencies:
            self._write_mecp_frequency_log(result, h)
        self._remove_checkpoint_set("check")
        return result

    def _remove_checkpoint_set(self, checkpoint_tag):
        """Remove successful verification checkpoints, preserving MECP ones."""
        for state in ("A", "B"):
            checkpoint_key = (checkpoint_tag, state)
            checkpoint_file = self._state_checkpoint_files.pop(
                checkpoint_key, None
            )
            if checkpoint_file and os.path.isfile(checkpoint_file):
                os.remove(checkpoint_file)

    def _run_seam_minimum_check(self, positions_bohr):
        """Run the requested post-convergence MECP Hessian analysis."""
        logger.info(f"Starting seam-minimum verification for {self.label}")
        result = self.verify_seam_minimum(
            write_frequencies=self.settings.mecp_numfreq
        )
        self._initial_seam_result = result
        self._last_seam_result = result
        status = (
            "MINIMUM"
            if result["is_minimum"]
            else "NOT A MINIMUM (saddle point on seam)"
        )
        logger.info(
            f"Seam-minimum check for {self.label}: {status} "
            f"(n_negative={result['n_negative']})"
        )
        if not result["is_minimum"]:
            if getattr(self.settings, "follow_seam_imaginary_mode", False):
                result = self._follow_seam_imaginary_mode(result)
                self._last_seam_result = result
                return result
            raise RuntimeError(
                f"Converged crossing is not a minimum on the seam "
                f"({result['n_negative']} negative eigenvalue(s))."
            )
        return result

    def _write_mode_displacement_xyz(self, path, positions, comment):
        """Write a projected-mode displacement as an XYZ structure."""
        with open(path, "w", encoding="utf-8") as output:
            output.write(f"{len(self.molecule.symbols)}\n{comment}\n")
            for symbol, position in zip(self.molecule.symbols, positions):
                output.write(
                    f"{symbol:<3s} {position[0]:+16.10f} "
                    f"{position[1]:+16.10f} {position[2]:+16.10f}\n"
                )

    @staticmethod
    def _constrained_seam_displacement(
        energy_diff,
        progress_error,
        grad_a,
        grad_b,
        progress_mode,
        step_size,
    ):
        """Return a step tangent to the crossing seam and progress plane.

        The two linearized constraints are ``E_A - E_B = 0`` and
        ``mode . (x - x_plane) = 0``.  A minimum-norm correction restores
        both constraints, while the state-A gradient is projected into their
        common null space before taking the downhill step.
        """
        diff_grad = np.asarray(grad_a, dtype=float).ravel() - np.asarray(
            grad_b, dtype=float
        ).ravel()
        mode = np.asarray(progress_mode, dtype=float).ravel()
        mode_norm = float(np.linalg.norm(mode))
        if mode_norm <= np.finfo(float).tiny:
            raise RuntimeError("Seam-following mode has zero norm.")
        mode /= mode_norm

        constraints = np.vstack((diff_grad, mode))
        gram = constraints @ constraints.T
        if np.linalg.cond(gram) > 1.0e12:
            raise RuntimeError(
                "Energy-gap and mode-following constraints are linearly "
                "dependent."
            )

        correction_rhs = np.array([-energy_diff, -progress_error])
        correction = constraints.T @ np.linalg.solve(gram, correction_rhs)
        gradient = np.asarray(grad_a, dtype=float).ravel()
        gradient_multipliers = np.linalg.solve(gram, constraints @ gradient)
        tangent_gradient = gradient - constraints.T @ gradient_multipliers
        displacement = correction - step_size * tangent_gradient
        return (
            displacement.reshape(np.asarray(grad_a).shape),
            tangent_gradient.reshape(np.asarray(grad_a).shape),
            correction.reshape(np.asarray(grad_a).shape),
        )

    def _optimize_on_seam_progress_plane(
        self,
        positions_bohr,
        plane_point_bohr,
        progress_mode,
        trace,
        macro_step,
    ):
        """Optimize all seam coordinates except one fixed progress mode."""
        positions = np.asarray(positions_bohr, dtype=float).copy()
        plane_point = np.asarray(plane_point_bohr, dtype=float)
        mode = np.asarray(progress_mode, dtype=float).ravel()
        mode /= np.linalg.norm(mode)
        displacement = np.zeros_like(positions)
        current_step_size = self.settings.step_size

        self.steps_folder = self.optimization_folder
        os.makedirs(self.steps_folder, exist_ok=True)
        self._state_checkpoint_files = {}
        self._last_spin_squared = {"A": None, "B": None}

        for inner_step in range(1, self.settings.max_steps + 1):
            step_index = f"macro{macro_step:02d}_inner{inner_step:03d}"
            trace.write(
                f"macro_step={macro_step} inner_step={inner_step} "
                "status=RUNNING_STATE_A\n"
            )
            trace.flush()
            ea, grad_a = self._run_state(positions, step_index, "A")
            trace.write(
                f"macro_step={macro_step} inner_step={inner_step} "
                "status=RUNNING_STATE_B\n"
            )
            trace.flush()
            eb, grad_b = self._run_state(positions, step_index, "B")
            progress_error = float(
                np.dot((positions - plane_point).ravel(), mode)
            )
            displacement, tangent_gradient, _correction = (
                self._constrained_seam_displacement(
                    energy_diff=ea - eb,
                    progress_error=progress_error,
                    grad_a=grad_a,
                    grad_b=grad_b,
                    progress_mode=mode,
                    step_size=current_step_size,
                )
            )
            displacement = self._apply_trust_radius(displacement)
            trace.write(
                f"macro_step={macro_step} inner_step={inner_step} "
                f"E_A={ea:.10f} E_B={eb:.10f} dE={ea-eb:+.6e} "
                f"progress_error={progress_error:+.6e} "
                f"tangent_grad_max={np.max(np.abs(tangent_gradient)):.3e} "
                f"displacement_max={np.max(np.abs(displacement)):.3e}\n"
            )

            constraints_converged = (
                abs(ea - eb) <= self.settings.energy_diff_tol
                and abs(progress_error) <= self.settings.disp_max_tol
            )
            if constraints_converged and self._is_converged(
                energy_diff=ea - eb,
                eff_grad=tangent_gradient,
                displacement=displacement,
            ):
                self.molecule.positions = positions * units.Bohr
                return positions

            positions = positions + displacement

        raise RuntimeError(
            "Constrained seam optimization did not converge within "
            f"{self.settings.max_steps} steps."
        )

    @staticmethod
    def _overlap_tracked_negative_mode(result, previous_mode):
        """Select and orient the negative mode overlapping the prior mode."""
        modes = np.asarray(result["modes"], dtype=float)
        eigenvalues = np.asarray(result["frequency_eigenvalues"], dtype=float)
        negative = np.flatnonzero(eigenvalues < -1.0e-6)
        if not len(negative):
            return None, None, None

        previous = np.asarray(previous_mode, dtype=float).ravel()
        previous /= np.linalg.norm(previous)
        overlaps = np.array(
            [np.dot(modes[index].ravel(), previous) for index in negative]
        )
        selected = int(negative[np.argmax(np.abs(overlaps))])
        signed_overlap = float(overlaps[np.argmax(np.abs(overlaps))])
        mode = modes[selected].ravel()
        if signed_overlap < 0.0:
            mode = -mode
            signed_overlap = -signed_overlap
        mode /= np.linalg.norm(mode)
        return selected, mode, signed_overlap

    @staticmethod
    def _displace_along_mode(positions, mode, distance):
        """Displace an ``(N, 3)`` geometry along a flattened normal mode."""
        positions = np.asarray(positions, dtype=float)
        mode = np.asarray(mode, dtype=float)
        if mode.size != positions.size:
            raise ValueError(
                "Normal-mode size does not match the Cartesian geometry."
            )
        return positions + distance * mode.reshape(positions.shape)

    def _follow_seam_imaginary_mode(self, result):
        """Follow a negative seam mode with a constrained progress plane."""
        frequencies = np.asarray(result.get("frequencies", []), dtype=float)
        modes = np.asarray(result.get("modes", []), dtype=float)
        if not len(frequencies) or modes.shape[0] != len(frequencies):
            raise RuntimeError(
                "Seam-mode following requires projected frequency modes."
            )

        mode_index = int(np.argmin(frequencies))
        if frequencies[mode_index] >= 0.0:
            return result

        base_positions = np.asarray(result["positions_angstrom"], dtype=float)
        mode = modes[mode_index].reshape(base_positions.shape)
        candidates = []
        displacement_norm = self.settings.seam_mode_displacement
        displacement_bohr = displacement_norm / units.Bohr
        follow_folder = os.path.join(
            self.folder, f"{self.label}_seam_follow"
        )
        os.makedirs(follow_folder, exist_ok=True)

        summary_file = os.path.join(
            follow_folder, f"{self.label}_seam_follow.log"
        )
        with open(summary_file, "w", encoding="utf-8") as summary:
            summary.write("CHEMSMART constrained MECP seam-mode following\n")
            summary.write(
                f"source_frequency={frequencies[mode_index]:+.6f} cm^-1\n"
                f"progress_step={displacement_norm:.6f} Angstrom\n"
            )

        for suffix, sign in (("plus", 1.0), ("minus", -1.0)):
            tracked_mode = sign * mode.ravel()
            current_bohr = base_positions / units.Bohr
            branch_label = f"{self.label}_seam_follow_{suffix}"
            xyz_file = os.path.join(follow_folder, f"{branch_label}.xyz")
            molecule = self.molecule.copy()
            molecule.positions = base_positions
            settings = self.settings.copy()
            settings.follow_seam_imaginary_mode = False
            settings.mecp_numfreq = True
            settings.restart = False
            # A loose optimization can declare convergence before the
            # displaced structure has escaped the seam saddle. Preserve
            # stricter custom values, otherwise require the tight preset.
            tight = settings.CONVERGENCE_PRESETS["tight"]
            for name in (
                "energy_diff_tol",
                "force_max_tol",
                "force_rms_tol",
                "disp_max_tol",
                "disp_rms_tol",
                "trust_radius",
            ):
                setattr(
                    settings,
                    name,
                    min(getattr(settings, name), tight[name]),
                )
            settings.convergence_preset = "tight"

            branch = self.__class__(
                molecule=molecule,
                settings=settings,
                label=branch_label,
                jobrunner=self.jobrunner,
                skip_completed=False,
            )
            branch._is_seam_follow_branch = True
            branch.set_folder(follow_folder)
            trace_file = os.path.join(
                follow_folder, f"{branch_label}_mode_follow.log"
            )
            try:
                with open(trace_file, "w", encoding="utf-8") as trace:
                    trace.write(
                        "CHEMSMART constrained MECP seam-mode branch\n"
                        f"branch={suffix}\n"
                    )
                    for macro_step in range(
                        1, self.settings.seam_mode_max_steps + 1
                    ):
                        plane_point = self._displace_along_mode(
                            current_bohr,
                            tracked_mode,
                            displacement_bohr,
                        )
                        trace.write(
                            f"macro_step={macro_step} "
                            "status=CONSTRAINED_OPTIMIZATION\n"
                        )
                        trace.flush()
                        self._write_mode_displacement_xyz(
                            xyz_file,
                            plane_point * units.Bohr,
                            (
                                f"{self.label}: constrained {suffix} "
                                f"macro step {macro_step}"
                            ),
                        )
                        current_bohr = branch._optimize_on_seam_progress_plane(
                            positions_bohr=plane_point,
                            plane_point_bohr=plane_point,
                            progress_mode=tracked_mode,
                            trace=trace,
                            macro_step=macro_step,
                        )
                        branch.molecule.positions = current_bohr * units.Bohr
                        trace.write(
                            f"macro_step={macro_step} "
                            "status=COMPUTING_PROJECTED_HESSIAN\n"
                        )
                        trace.flush()
                        step_result = branch.verify_seam_minimum(
                            write_frequencies=True,
                            macro_step=macro_step,
                        )
                        selected, next_mode, overlap = (
                            self._overlap_tracked_negative_mode(
                                step_result, tracked_mode
                            )
                        )
                        trace.write(
                            f"macro_step={macro_step} "
                            f"mecp_energy={step_result['mecp_energy']:+.12f} "
                            f"n_negative={step_result['n_negative']}"
                        )
                        if next_mode is None:
                            trace.write(" status=NEGATIVE_MODE_CROSSED\n")
                            break
                        trace.write(
                            f" tracked_mode={selected + 1} "
                            f"frequency={step_result['frequencies'][selected]:+.6f} "
                            f"overlap={overlap:.6f}\n"
                        )
                        tracked_mode = next_mode
                    else:
                        raise RuntimeError(
                            "Negative seam mode remains after "
                            f"{self.settings.seam_mode_max_steps} constrained "
                            "following steps."
                        )

                # Release the progress coordinate only after the constrained
                # Hessian no longer contains a significant negative mode.
                branch.molecule.positions = current_bohr * units.Bohr
                branch.run()
            except RuntimeError as error:
                logger.warning(
                    "Seam-mode %s branch did not yield a minimum: %s",
                    suffix,
                    error,
                )
            branch_result = getattr(branch, "_last_seam_result", None)
            if branch_result is not None and branch_result["is_minimum"]:
                branch_result["convergence_metrics"] = (
                    branch._final_convergence_metrics
                )
                branch_result["convergence_thresholds"] = {
                    "energy_diff": branch.settings.energy_diff_tol,
                    "pgrad_max": branch.settings.force_max_tol,
                    "pgrad_rms": branch.settings.force_rms_tol,
                    "disp_max": branch.settings.disp_max_tol,
                    "disp_rms": branch.settings.disp_rms_tol,
                }
                branch_result["optimization_steps"] = (
                    branch._final_optimization_steps
                )
                candidates.append(
                    (
                        branch_result["mecp_energy"],
                        branch_result,
                        suffix,
                        macro_step,
                    )
                )

        if not candidates:
            raise RuntimeError(
                "Neither projected-mode displacement converged to a seam "
                "minimum. The +/- XYZ structures and branch reports were "
                "kept for inspection."
            )

        _, selected, selected_suffix, selected_macro_steps = min(
            candidates, key=lambda item: item[0]
        )
        self._selected_seam_follow_branch = selected_suffix
        self._selected_seam_follow_macro_steps = selected_macro_steps
        self.molecule.positions = np.asarray(
            selected["positions_angstrom"], dtype=float
        )
        self._write_seam_check_log(selected, self.settings.hess_step_size, 1)
        self._write_mecp_frequency_log(selected, self.settings.hess_step_size)
        with open(summary_file, "a", encoding="utf-8") as output:
            output.write(f"selected_branch={selected_suffix}\n")
            output.write(f"selected_macro_steps={selected_macro_steps}\n")
            output.write(
                f"selected_mecp_energy={selected['mecp_energy']:+.12f} "
                "Hartree\nstatus=MECP MINIMUM\n"
            )
        logger.info(
            "Seam-mode following selected the %s verified MECP minimum.",
            selected_suffix,
        )
        return selected

    def _write_seam_check_log(self, result, h, step_prefix):
        """Write the seam-minimum verification results to a log file."""
        os.makedirs(self.numfreq_folder, exist_ok=True)
        seam_check_file = os.path.join(
            self.numfreq_folder, f"{self.label}_seam_check.log"
        )
        with open(seam_check_file, "w", encoding="utf-8") as f:
            f.write("CHEMSMART MECP seam-minimum verification\n")
            f.write(
                f"label={self.label} hess_step={h:.2e} Bohr "
                f"check_step_start={step_prefix}\n"
            )
            f.write(
                f"energy_diff={result['energy_diff']:+.6e} Hartree "
                f"n_projected={result['n_projected']} "
                f"lagrange_multiplier={result['lagrange_multiplier']:+.8e}\n"
            )
            n_neg = result["n_negative"]
            is_min = result["is_minimum"]
            f.write(
                f"n_significant_negative_mass_weighted_eigenvalues={n_neg}  "
                f"{'MECP MINIMUM' if is_min else 'SADDLE POINT ON SEAM'}\n"
            )
            f.write("Minimum test tolerance: -1e-6 Hartree/(Bohr^2 amu).\n")
            f.write("\nEigenvalues of H_eff (Hartree/Bohr^2):\n")
            for i, ev in enumerate(result["eigenvalues"]):
                flag = "  ** NEGATIVE **" if ev < 0.0 else ""
                f.write(f"  mode {i + 1:4d}: {ev:+.6e}{flag}\n")
        logger.info(f"Seam-minimum check results written to {seam_check_file}")

    def _write_mecp_frequency_log(self, result, h):
        """Write projected MECP frequencies and Cartesian normal modes."""
        frequency_file = os.path.join(
            self.folder, f"{self.label}_mecp_freq.log"
        )
        symbols = list(self.molecule.symbols)
        positions = np.asarray(result["positions_angstrom"], dtype=float)
        masses = np.asarray(result["atomic_masses"], dtype=float)

        with open(frequency_file, "w", encoding="utf-8") as output:
            output.write(f"{MECP_FREQUENCY_HEADER}\n")
            output.write(f"label={self.label}\n")
            output.write(f"hess_step={h:.6e} Bohr\n")
            output.write(f"energy_A={result['energy_a']:+.12f} Hartree\n")
            output.write(f"energy_B={result['energy_b']:+.12f} Hartree\n")
            output.write(
                f"mecp_energy={result['mecp_energy']:+.12f} Hartree\n"
            )
            output.write(
                f"energy_difference={result['energy_diff']:+.6e} Hartree\n"
            )
            output.write(
                f"lagrange_multiplier="
                f"{result['lagrange_multiplier']:+.8e}\n"
            )
            output.write(f"n_atoms={len(symbols)}\n")
            output.write(f"n_projected={result['n_projected']}\n")
            output.write(f"n_modes={len(result['frequencies'])}\n")
            output.write(
                f"n_imaginary={int(np.sum(result['frequencies'] < 0.0))}\n"
                f"n_significant_imaginary={result['n_negative']}\n\n"
            )
            output.write(
                f"multiplicity_A={result['multiplicity_a']}\n"
                f"multiplicity_B={result['multiplicity_b']}\n"
                "rotational_symmetry_number="
                f"{result['rotational_symmetry_number']}\n\n"
            )

            output.write("Geometry (Angstrom) and masses (amu):\n")
            for symbol, position, mass in zip(symbols, positions, masses):
                output.write(
                    f"{symbol:>3s} {position[0]:+16.8f} "
                    f"{position[1]:+16.8f} {position[2]:+16.8f} "
                    f"{mass:12.8f}\n"
                )

            output.write("\nProjected MECP frequencies (cm^-1):\n")
            for index, frequency in enumerate(result["frequencies"], 1):
                output.write(f"mode {index:4d}: {frequency:+14.6f}\n")

            output.write(
                "\nCartesian normal modes "
                "(mass-unweighted, unit-normalized):\n"
            )
            for mode_index, mode in enumerate(result["modes"], 1):
                frequency = result["frequencies"][mode_index - 1]
                output.write(
                    f"mode {mode_index:4d} "
                    f"frequency={frequency:+.6f} cm^-1\n"
                )
                mode_vectors = np.asarray(mode).reshape(len(symbols), 3)
                for atom_index, (symbol, vector) in enumerate(
                    zip(symbols, mode_vectors), 1
                ):
                    output.write(
                        f"{atom_index:5d} {symbol:>3s} "
                        f"{vector[0]:+14.8f} {vector[1]:+14.8f} "
                        f"{vector[2]:+14.8f}\n"
                    )
                output.write("\n")

            output.write(f"{MECP_FREQUENCY_TERMINATION_MARKER}\n")

        logger.info(
            f"MECP projected frequencies written to {frequency_file}"
        )
