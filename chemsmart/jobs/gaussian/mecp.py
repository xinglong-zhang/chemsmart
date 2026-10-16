"""
Gaussian Minimum Energy Cross Point (MECP) job implementation.

This module provides the GaussianMECPJob class for performing
Minimum Energy Cross Point calculations using Gaussian.
"""

import logging
import os
import re
from functools import partial
from shutil import copyfile
from tempfile import TemporaryDirectory
from typing import Type

import numpy as np
from ase import units

from chemsmart.io.mecp import read_report_fields, state_job_label
from chemsmart.jobs.gaussian.job import GaussianGeneralJob, GaussianJob
from chemsmart.jobs.gaussian.settings import (
    GaussianJobSettings,
    GaussianLinkJobSettings,
    GaussianMECPJobSettings,
)
from chemsmart.jobs.mecp import MECPSeamMixin
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
        "mecp_numfreq_only",
        "hess_step_size",
        "hessian_method",
        "follow_seam_imaginary_mode",
        "seam_mode_displacement",
        "seam_mode_max_steps",
        "restart",
    }
)


class GaussianMECPJob(MECPSeamMixin, GaussianJob):
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
    STATE_OUTPUT_EXTENSION = ".log"

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
        self._reset_execution_state()

    @classmethod
    def settings_class(cls) -> Type[GaussianMECPJobSettings]:
        return GaussianMECPJobSettings

    @property
    def report_file(self):
        return os.path.join(
            self.optimization_folder, f"{self.label}_report.log"
        )

    @property
    def final_gaussian_logs(self):
        """Convenient single-state outputs at the selected final geometry."""
        label = (
            self.label
            if self.label.lower().endswith("_mecp")
            else f"{self.label}_mecp"
        )
        return {
            "A": os.path.join(self.folder, f"{label}.log"),
            "B": os.path.join(self.folder, f"{label}_B.log"),
        }

    def _export_final_gaussian_logs(self, sources):
        """Copy the selected geometry's original logs without rewriting them."""
        for state in ("A", "B"):
            if not os.path.isfile(sources[state]):
                raise RuntimeError(
                    f"Missing final state {state} Gaussian output: {sources[state]}"
                )
        with TemporaryDirectory(
            prefix="mecp-export-", dir=self.folder
        ) as temporary:
            for state in ("A", "B"):
                copyfile(
                    sources[state], os.path.join(temporary, f"{state}.log")
                )
            for state, destination in self.final_gaussian_logs.items():
                os.replace(
                    os.path.join(temporary, f"{state}.log"), destination
                )

    @property
    def trajectory_file(self):
        return os.path.join(self.optimization_folder, f"{self.label}_traj.xyz")

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
        if self.settings.mecp_numfreq:
            if os.path.isfile(self.frequency_file):
                saved_method = read_report_fields(self.frequency_file).get(
                    "hessian_method", "numerical"
                )
                if saved_method != self.hessian_method:
                    return False
        if getattr(self.settings, "mecp_numfreq_only", False):
            if not os.path.isfile(self.final_report_file):
                return False
            with open(self.final_report_file, encoding="utf-8") as report:
                complete = "NUMFREQ ONLY COMPLETED" in report.read()
            return complete and all(
                os.path.isfile(path)
                for path in (
                    os.path.join(self.folder, f"{self.label}_mecp_freq.log"),
                    os.path.join(
                        self.numfreq_folder, f"{self.label}_seam_check.log"
                    ),
                )
            )
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

    def _analytic_state_settings(self, settings):
        """Build the fixed-geometry analytic input from the reference settings."""
        if settings.input_string:
            raise ValueError(
                "anafreq cannot use a literal input_string; use structured settings."
            )
        # The reference force job has already stabilized the wavefunction
        # in link mode. Read that checkpoint in a fixed-geometry freq job.
        settings = GaussianJobSettings(
            **{
                k: v
                for k, v in settings.__dict__.items()
                if not k.startswith("_")
                and k not in _MECP_ONLY_KEYS
                and k not in {"link", "link_route"}
            }
        )
        settings.jobtype = "sp"
        settings.freq = False
        settings.numfreq = False
        settings.dieze_tag = "p"
        settings.forces = False
        settings.stable = None
        if settings.route_to_be_written:
            route = GaussianLinkJobSettings._without_route_options(
                settings.route_to_be_written,
                {"freq", "force", "opt", "irc", "stable"},
            )
            route = re.sub(r"^\s*#(?:[pnt])?\s*", "", route, flags=re.I)
            settings.route_to_be_written = (
                "#p " + self._with_required_nosymm(route) + " freq=analytic"
            )
        settings.additional_route_parameters = (
            GaussianLinkJobSettings._without_route_options(
                settings.additional_route_parameters or "",
                {"freq", "force", "opt", "irc", "stable"},
            )
            + " freq=analytic"
        )
        return settings

    def _run_state(
        self,
        positions_bohr,
        step_idx,
        state,
        checkpoint_tag=None,
        analytic=False,
    ):
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
        if analytic:
            settings = self._analytic_state_settings(settings)
        state_label = state_job_label(
            self.label, step_idx, state, checkpoint_tag
        )

        checkpoint_key = (checkpoint_tag, state)
        oldchkfile = self._state_checkpoint_files.get(checkpoint_key)
        if analytic and (not oldchkfile or not os.path.isfile(oldchkfile)):
            raise RuntimeError(
                f"Analytic Hessian for state {state} requires the reference checkpoint. "
                "Enable Gaussian checkpoint writing and rerun the reference calculation."
            )
        if analytic and oldchkfile:
            settings.additional_route_parameters = (
                GaussianLinkJobSettings._without_route_options(
                    settings.additional_route_parameters or "", {"guess"}
                )
                + " guess=read"
            )
            if settings.route_to_be_written:
                settings.route_to_be_written = (
                    GaussianLinkJobSettings._without_route_options(
                        settings.route_to_be_written, {"guess"}
                    )
                    + " guess=read"
                )
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
                os.path.basename(self.numfreq_folder)
                if checkpoint_tag == "check"
                else f"{self.label}_optimization"
            ),
            "checkpoint_filename": (
                f"{self.label}_{checkpoint_part}{state}.chk"
            ),
        }

        if self.settings.use_link and not analytic:
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
        if analytic:
            if not hasattr(self, "_analytic_hessians"):
                self._analytic_hessians = {}
            self._analytic_hessians[state] = (
                self.jobrunner.format_hessian_checkpoint(job)
            )
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
        self._reset_execution_state()
        if getattr(self.settings, "mecp_numfreq_only", False):
            return self._run_numfreq_only()
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

        if not getattr(self, "_is_seam_follow_branch", False):
            sources = (
                seam_result["gaussian_output_files"]
                if seam_result is not None
                else {
                    state: os.path.join(
                        self.optimization_folder,
                        f"{self.label}_step{converged_step}_{state}.log",
                    )
                    for state in ("A", "B")
                }
            )
            self._export_final_gaussian_logs(sources)

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

    def _run_numfreq_only(self):
        """Analyze the supplied stationary crossing without changing geometry."""
        return self._run_frequency_only()

    # ------------------------------------------------------------------
    # Seam-minimum verification via effective Hessian analysis
    # ------------------------------------------------------------------

    def _remove_checkpoint_set(self, checkpoint_tag):
        """Remove successful verification checkpoints, preserving MECP ones."""
        for state in ("A", "B"):
            checkpoint_key = (checkpoint_tag, state)
            checkpoint_file = self._state_checkpoint_files.pop(
                checkpoint_key, None
            )
            if checkpoint_file and os.path.isfile(checkpoint_file):
                os.remove(checkpoint_file)

    def verify_seam_minimum(self, *args, **kwargs):
        result = super().verify_seam_minimum(*args, **kwargs)
        # Preserve the Gaussian result API used by existing export consumers.
        result["gaussian_output_files"] = result["state_output_files"]
        return result

    def _final_report_outputs(self, result):
        return self.final_gaussian_logs

    def _final_report_legacy_fields(self, result, frequency_only):
        if frequency_only:
            return {}
        return {
            f"gaussian_output_{state}": os.path.basename(path)
            for state, path in self.final_gaussian_logs.items()
        }
