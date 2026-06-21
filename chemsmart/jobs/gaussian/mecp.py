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

    At each iteration (step 1, 2, …, max_steps) two Gaussian sub-jobs are
    executed, labelled ``<label>_step<NNNNNN>_A`` and
    ``<label>_step<NNNNNN>_B`` where ``<NNNNNN>`` is the **1-indexed**
    six-digit zero-padded step number (e.g. ``000001`` for the first step),
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
            "convergence_preset",
            "verify_seam_minimum",
            "hess_step_size",
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
                    "guess": self.settings.guess,
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
            title=f"{title} step {step_idx:06d}",
            state=state,
        )
        state_label = f"{self.label}_step{step_idx:06d}_{state}"

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

        return displacement, projected_grad, seam_correction

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
        self, f, step_idx, ea, eb, projected_grad, displacement, seam_correction, step_size
    ):
        energy_diff = ea - eb
        pgrad_max = float(np.max(np.abs(projected_grad)))
        pgrad_rms = self._rms(projected_grad)
        disp_max = float(np.max(np.abs(displacement)))
        disp_rms = self._rms(displacement)
        seam_max = float(np.max(np.abs(seam_correction)))
        seam_rms = self._rms(seam_correction)
        f.write(
            f"step={step_idx:06d} E_A={ea:.10f} E_B={eb:.10f} "
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

        with open(self.report_file, "w") as report:
            report.write("CHEMSMART self-contained MECP optimization\n")
            report.write(
                f"max_steps={self.settings.max_steps} "
                f"step_size={self.settings.step_size} "
                f"trust_radius={self.settings.trust_radius} "
                f"adaptive_step_size={self.settings.adaptive_step_size} "
                f"step_size_method={self.settings.step_size_method}\n"
            )
            for step_idx in range(1, self.settings.max_steps + 1):
                self._write_trajectory_frame(positions_bohr, step_idx)
                ea, grad_a = self._run_state(positions_bohr, step_idx, "A")
                eb, grad_b = self._run_state(positions_bohr, step_idx, "B")

                energy_diff = ea - eb

                displacement, projected_grad, seam_correction = self._mecp_displacement(
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
                    seam_correction,
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

        if self.settings.verify_seam_minimum:
            self._run_seam_minimum_check(positions_bohr)

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

    def _project_hessian(self, hessian, projection_vectors):
        """
        Apply the projector :math:`P = I - \\sum_i |v_i\\rangle\\langle v_i|`
        to the Hessian from both sides: :math:`H_\\text{eff} = P H P`.

        Args:
            hessian (np.ndarray): Square Hessian matrix (3N × 3N),
                Hartree/Bohr².
            projection_vectors (list[np.ndarray]): Orthonormal vectors to
                project out (each of length 3N).

        Returns:
            np.ndarray: Projected Hessian, same shape as input.
        """
        n = hessian.shape[0]
        P = np.eye(n)
        for v in projection_vectors:
            P -= np.outer(v, v)
        return P @ hessian @ P

    def _compute_numerical_hessian(self, positions_bohr, h, step_prefix):
        """
        Compute the average Hessian :math:`H = (H_A + H_B)/2` numerically
        via central finite differences of the Cartesian forces.

        For each coordinate ``j`` (0 … 3N−1) the geometry is displaced by
        ``±h`` Bohr and both states are evaluated:

        .. math::

            H_{ij} \\approx \\frac{g_i(+h_j) - g_i(-h_j)}{2h}

        where :math:`g_i` denotes the gradient component (force negated).
        The Hessian is symmetrised as :math:`(H + H^T)/2` before averaging.

        .. warning::

            This method runs **4 × 3N** Gaussian sub-jobs (2 displacements
            × 2 spin states × 3N coordinates), which is expensive for large
            molecules (120 jobs for a 10-atom molecule).  Use the
            ``hess_step_size`` setting to control accuracy vs. cost.

        Args:
            positions_bohr (np.ndarray): Current geometry, shape (n_atoms, 3).
            h (float): Finite-difference step size in Bohr.
            step_prefix (int): Starting sub-job step index (to avoid
                clashing with MECP optimisation step labels).

        Returns:
            np.ndarray: Average symmetrised Hessian, shape (3N, 3N),
            Hartree/Bohr².
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

            _, g_A_plus = self._run_state(pos_plus, step_p, "A")
            _, g_B_plus = self._run_state(pos_plus, step_p, "B")
            _, g_A_minus = self._run_state(pos_minus, step_m, "A")
            _, g_B_minus = self._run_state(pos_minus, step_m, "B")

            H_A[:, j] = (g_A_plus.ravel() - g_A_minus.ravel()) / (2 * h)
            H_B[:, j] = (g_B_plus.ravel() - g_B_minus.ravel()) / (2 * h)

        H_A = (H_A + H_A.T) / 2
        H_B = (H_B + H_B.T) / 2
        return (H_A + H_B) / 2

    def verify_seam_minimum(self, h=None, step_prefix=900000):
        """
        Verify that the current MECP geometry is a **minimum on the crossing
        seam**, not merely a crossing point.

        The method computes the average Hessian :math:`H = (H_A + H_B)/2`
        numerically via central finite differences, then projects out
        translational, rotational, and gradient-difference degrees of
        freedom:

        .. math::

            H_\\text{eff} = P H P, \\quad
            P = I - \\textstyle\\sum_i |v_i\\rangle\\langle v_i|

        where :math:`\\{v_i\\}` spans translations (3), rotations (up to 3),
        and the gradient-difference direction
        :math:`\\mathbf{g}_\\Delta / \\|\\mathbf{g}_\\Delta\\|`.

        The eigenvalues of :math:`H_\\text{eff}` are diagonalised after
        discarding the (approximately) zero eigenvalues corresponding to
        projected-out modes.  All remaining eigenvalues positive → confirmed
        MECP minimum; any negative eigenvalue indicates a saddle point on the
        seam.

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
            step_prefix (int, optional): Starting sub-job step index
                (default: 900000, well above any normal MECP step).

        Returns:
            dict: Keys ``"eigenvalues"`` (1-D array, non-projected modes),
            ``"n_negative"`` (int), ``"is_minimum"`` (bool),
            ``"energy_diff"`` (float, Hartree),
            ``"n_projected"`` (int, modes removed).
        """
        if h is None:
            h = self.settings.hess_step_size

        positions_bohr = (
            np.array(self.molecule.positions, dtype=float) / units.Bohr
        )

        logger.info(
            f"verify_seam_minimum: computing gradient difference at {self.label}"
        )
        ea, grad_a = self._run_state(positions_bohr, step_prefix - 1, "A")
        eb, grad_b = self._run_state(positions_bohr, step_prefix - 1, "B")
        diff_grad = grad_a - grad_b

        logger.info(
            f"verify_seam_minimum: computing numerical Hessian for {self.label} "
            f"(h={h} Bohr, {4 * 3 * len(self.molecule.symbols)} sub-jobs)"
        )
        H_avg = self._compute_numerical_hessian(
            positions_bohr, h=h, step_prefix=step_prefix
        )

        proj_vecs = self._build_projection_vectors(positions_bohr, diff_grad)
        H_eff = self._project_hessian(H_avg, proj_vecs)

        eigenvalues = np.sort(np.linalg.eigvalsh(H_eff))
        n_proj = len(proj_vecs)
        non_zero_evals = eigenvalues[n_proj:]
        n_negative = int(np.sum(non_zero_evals < -1.0e-6))

        result = {
            "eigenvalues": non_zero_evals,
            "n_negative": n_negative,
            "is_minimum": n_negative == 0,
            "energy_diff": ea - eb,
            "n_projected": n_proj,
        }

        self._write_seam_check_log(result, h, step_prefix)
        return result

    def _run_seam_minimum_check(self, positions_bohr):
        """Called at the end of ``_run()`` when ``verify_seam_minimum`` is set."""
        logger.info(
            f"Starting seam-minimum verification for {self.label}"
        )
        result = self.verify_seam_minimum()
        status = "MINIMUM" if result["is_minimum"] else "NOT A MINIMUM (saddle point on seam)"
        logger.info(
            f"Seam-minimum check for {self.label}: {status} "
            f"(n_negative={result['n_negative']})"
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
                f"step_prefix={step_prefix}\n"
            )
            f.write(
                f"energy_diff={result['energy_diff']:+.6e} Hartree "
                f"n_projected={result['n_projected']}\n"
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
