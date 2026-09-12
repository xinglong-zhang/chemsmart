"""
Gaussian Minimum Energy Cross Point (MECP) job implementation.

This module provides the GaussianMECPJob class for performing
Minimum Energy Cross Point calculations using Gaussian.
"""

import logging
import os
import re
from typing import Type

import numpy as np
from ase import units

from chemsmart.jobs.gaussian.job import GaussianGeneralJob, GaussianJob
from chemsmart.jobs.gaussian.settings import GaussianMECPJobSettings

logger = logging.getLogger(__name__)

MIN_DIFF_GRAD_NORM_SQ = 1.0e-20

# Settings consumed by the MECP driver rather than Gaussian link sub-jobs.
MECP_ONLY_KEYS = frozenset(
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
        "convergence_preset",
        "verify_seam_minimum",
        "hess_step_size",
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
<<<<<<< Updated upstream
    MIN_DIFF_GRAD_NORM_SQ = 1.0e-20

    @staticmethod
    def _without_unavailable_guess_read(route):
        """Remove ``read`` from a Gaussian guess option for the first step."""

        def strip_parenthesized_read(match):
            options = [
                option.strip()
                for option in match.group(1).split(",")
                if option.strip().lower() != "read"
            ]
            return f"guess=({','.join(options)})" if options else ""

        route = re.sub(
            r"\bguess\s*=\s*\(([^)]*)\)",
            strip_parenthesized_read,
            route,
            flags=re.IGNORECASE,
        )
        route = re.sub(r"\bguess\s*=\s*read\b", "", route, flags=re.IGNORECASE)
        return " ".join(route.split())

    @staticmethod
    def _with_required_nosymm(route):
        """Ensure Gaussian forces remain in the MECP Cartesian frame."""
        route = route or ""
        if re.search(r"\bnosymm\b", route, flags=re.IGNORECASE):
            return route
        if re.search(r"\bsymm(?:etry)?\b", route, flags=re.IGNORECASE):
            raise ValueError(
                "MECP requires nosymm so Cartesian forces remain aligned with "
                "the optimization coordinates; remove the explicit symmetry option."
            )
        return f"{route} nosymm".strip()

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
            "harvey_initial_hessian",
            "harvey_max_component_step",
            "harvey_max_condition",
            "use_link",
            "convergence_preset",
            "verify_seam_minimum",
            "hess_step_size",
            "restart",
            # 'stable' and 'guess' are kept: they are valid GaussianLinkJobSettings
            # params and will be overridden with state-specific values anyway.
        }
    )
=======
>>>>>>> Stashed changes

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
            return any(line.startswith("Converged at step") for line in f)

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
                if not k.startswith("_") and k not in MECP_ONLY_KEYS
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
        job_tag = f"{checkpoint_tag}_" if checkpoint_tag else ""
        state_label = f"{self.label}_{job_tag}step{step_idx}_{state}"

        checkpoint_key = (checkpoint_tag, state)
        oldchkfile = self._state_checkpoint_files.get(checkpoint_key)
        if not oldchkfile and not self.settings.use_link:
            route = settings.additional_route_parameters or ""
            settings.additional_route_parameters = (
                self._without_unavailable_guess_read(route)
            )

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
        job.set_folder(self.steps_folder)
        job.scratch_parent_folder = f"{self.label}_steps"
        checkpoint_part = f"{checkpoint_tag}_" if checkpoint_tag else ""
        job.checkpoint_filename = f"{self.label}_{checkpoint_part}{state}.chk"
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
        self._last_spin_squared[state] = getattr(output, "spin_squared", None)
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

        if diff_norm_sq < MIN_DIFF_GRAD_NORM_SQ:
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
        if diff_norm_sq < self.MIN_DIFF_GRAD_NORM_SQ:
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

    _GROW_SHRINK_IMPROVEMENT_THRESHOLD = 0.10
    _GROW_SHRINK_REGRESSION_THRESHOLD = 0.02
    _BB_CURVATURE_COSINE_MIN = 1.0e-4
    _BB_RELATIVE_STEP_MIN = 0.5
    _BB_RELATIVE_STEP_MAX = 2.0

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
        if relative_progress > self._GROW_SHRINK_IMPROVEMENT_THRESHOLD:
            new_step = current_step_size * self.settings.step_size_grow
        elif relative_progress < -self._GROW_SHRINK_REGRESSION_THRESHOLD:
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
            and r_dot_g > self._BB_CURVATURE_COSINE_MIN * curvature_scale
        )

        if reliable_curvature:
            candidate = r_dot_r / r_dot_g
        else:
            candidate = current_step_size * self.settings.step_size_shrink

        relative_min = current_step_size * self._BB_RELATIVE_STEP_MIN
        relative_max = current_step_size * self._BB_RELATIVE_STEP_MAX
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

<<<<<<< Updated upstream
    def _log_step(
=======
    def log_step(
>>>>>>> Stashed changes
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
        self.steps_folder = os.path.join(self.folder, f"{self.label}_steps")
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

        if self.settings.verify_seam_minimum:
            self._run_seam_minimum_check(positions_bohr)

        with open(self.report_file, "a", encoding="utf-8") as report:
            report.write(f"Converged at step {converged_step}.\n")
        if os.path.isfile(self.state_file):
            os.remove(self.state_file)

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

    @classmethod
    def _lagrangian_hessian(cls, hessian_a, hessian_b, grad_a, grad_b):
        """Return the constrained MECP Lagrangian Hessian and multiplier."""
        diff_grad = np.asarray(grad_a).ravel() - np.asarray(grad_b).ravel()
        diff_norm_sq = float(np.dot(diff_grad, diff_grad))
        if diff_norm_sq < cls.MIN_DIFF_GRAD_NORM_SQ:
            raise RuntimeError(
                "Difference gradient is too small for seam verification."
            )
        multiplier = float(np.dot(np.asarray(grad_a).ravel(), diff_grad))
        multiplier /= diff_norm_sq
        hessian = (1.0 - multiplier) * hessian_a + multiplier * hessian_b
        return hessian, multiplier

    def _compute_numerical_hessian(self, positions_bohr, h, step_prefix):
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

            step_p = step_prefix + 2 * j
            step_m = step_prefix + 2 * j + 1

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

    def verify_seam_minimum(self, h=None, step_prefix=1):
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
            ``--verify-seam-minimum`` CLI flag or by setting
            ``settings.verify_seam_minimum = True``.

        Args:
            h (float, optional): Finite-difference step size in Bohr.
                Defaults to ``settings.hess_step_size`` (1 × 10⁻³ Bohr).
            step_prefix (int, optional): Starting check-specific sub-job step
                index (default: 1). The ``check`` name component prevents
                clashes with normal MECP steps.

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
        ea, grad_a = self._run_state(
            positions_bohr, step_prefix, "A", checkpoint_tag="check"
        )
        eb, grad_b = self._run_state(
            positions_bohr, step_prefix, "B", checkpoint_tag="check"
        )
        diff_grad = grad_a - grad_b

        logger.info(
            f"verify_seam_minimum: computing numerical Hessian for {self.label} "
            f"(h={h} Bohr, {4 * 3 * len(self.molecule.symbols)} sub-jobs)"
        )
        H_A, H_B = self._compute_numerical_hessian(
            positions_bohr, h=h, step_prefix=step_prefix + 1
        )

        proj_vecs = self._build_projection_vectors(positions_bohr, diff_grad)
        H_lagrangian, lagrange_multiplier = self._lagrangian_hessian(
            H_A, H_B, grad_a, grad_b
        )
        H_seam = self._reduced_hessian(H_lagrangian, proj_vecs)

        non_zero_evals = np.sort(np.linalg.eigvalsh(H_seam))
        n_proj = len(proj_vecs)
        n_negative = int(np.sum(non_zero_evals < -1.0e-6))

        result = {
            "eigenvalues": non_zero_evals,
            "n_negative": n_negative,
            "is_minimum": n_negative == 0,
            "energy_diff": ea - eb,
            "n_projected": n_proj,
            "lagrange_multiplier": lagrange_multiplier,
        }

        self._write_seam_check_log(result, h, step_prefix)
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
        """Called at the end of ``_run()`` when ``verify_seam_minimum`` is set."""
        logger.info(f"Starting seam-minimum verification for {self.label}")
        result = self.verify_seam_minimum()
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
            raise RuntimeError(
                f"Converged crossing is not a minimum on the seam "
                f"({result['n_negative']} negative eigenvalue(s))."
            )

    def _write_seam_check_log(self, result, h, step_prefix):
        """Write the seam-minimum verification results to a log file."""
        seam_check_file = os.path.join(
            self.folder, f"{self.label}_seam_check.log"
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
                f"n_negative_eigenvalues={n_neg}  "
                f"{'MECP MINIMUM' if is_min else 'SADDLE POINT ON SEAM'}\n"
            )
            f.write("\nEigenvalues of H_eff (Hartree/Bohr^2):\n")
            for i, ev in enumerate(result["eigenvalues"]):
                flag = "  ** NEGATIVE **" if ev < -1.0e-6 else ""
                f.write(f"  mode {i + 1:4d}: {ev:+.6e}{flag}\n")
        logger.info(f"Seam-minimum check results written to {seam_check_file}")
