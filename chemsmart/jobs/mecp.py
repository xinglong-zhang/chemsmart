"""Backend-independent MECP seam analysis and constrained mode following.

Backends provide _run_state (Bohr, Hartree, Hartree/Bohr), analytic
Hessians (Hartree/Bohr squared), STATE_OUTPUT_EXTENSION, checkpoint
cleanup, and a run() implementation for final unconstrained relaxation.
Settings supply convergence thresholds and follow controls; branch jobs
are constructed with the same backend class. No backend modules are imported.
"""

import logging
import os

import numpy as np
from ase import units

from chemsmart.io.mecp import (
    MECPFiles,
    convergence_thresholds,
    write_final_report,
    write_frequency_log,
    write_seam_check_log,
)
from chemsmart.utils.constants import (
    amu_to_kg,
    bohr_to_meter,
    hartree_to_joules,
)

logger = logging.getLogger(__name__)
_MIN_DIFF_GRAD_NORM_SQ = 1.0e-20


MECP_CONVERGENCE_PRESETS = {
    "standard": {
        "energy_diff_tol": 5.0e-5,  # Hartree
        "force_max_tol": 7.0e-4,  # Hartree/Bohr
        "force_rms_tol": 5.0e-4,  # Hartree/Bohr
        "disp_max_tol": 4.0e-3,  # Bohr
        "disp_rms_tol": 2.5e-3,  # Bohr
        "trust_radius": 0.3,  # Bohr
    },
    "tight": {
        "energy_diff_tol": 1.0e-5,  # Hartree
        "force_max_tol": 3.0e-4,  # Hartree/Bohr
        "force_rms_tol": 1.0e-4,  # Hartree/Bohr
        "disp_max_tol": 2.0e-3,  # Bohr
        "disp_rms_tol": 1.0e-3,  # Bohr
        "trust_radius": 0.1,  # Bohr
    },
}


def configure_convergence(settings, preset_name, overrides):
    """Apply a preset and explicit overrides without changing their values.
    This should be used in constructor before any call that sets attributes
    after object instantiation."""
    if preset_name not in settings.CONVERGENCE_PRESETS:
        raise ValueError(
            f"Unknown convergence_preset={preset_name!r}. "
            f"Choose from: {list(settings.CONVERGENCE_PRESETS)}"
        )
    settings.convergence_preset = preset_name
    for name, value in settings.CONVERGENCE_PRESETS[preset_name].items():
        setattr(settings, name, overrides.get(name, value))


class MECPSeamMixin:
    """Shared algorithms for backends implementing the state-job contract."""

    @property
    def mecp_files(self):
        return MECPFiles(self.folder, self.label)

    @property
    def hessian_method(self):
        return self.settings.hessian_method

    def _reset_execution_state(self):
        """Initialize shared state for optimization and frequency workflows."""
        self._state_checkpoint_files = {}
        self._last_spin_squared = {"A": None, "B": None}
        for name in (
            "_final_report_context",
            "_initial_seam_result",
            "_last_seam_result",
            "_selected_seam_follow_branch",
            "_selected_seam_follow_macro_steps",
        ):
            self.__dict__.pop(name, None)

    @property
    def optimization_folder(self):
        return self.mecp_files.optimization_folder

    @property
    def frequency_type_folder(self):
        return self.mecp_files.hessian_folder(self.hessian_method)

    @property
    def final_report_file(self):
        return self.mecp_files.final_report

    @property
    def frequency_file(self):
        return self.mecp_files.frequencies

    def _reference_output_files(self, reference_step):
        return self.mecp_files.reference_outputs(
            reference_step, self.hessian_method, self.STATE_OUTPUT_EXTENSION
        )

    def _write_seam_check_log(self, result, h, step_prefix):
        filename = self.mecp_files.seam_check(self.hessian_method)
        os.makedirs(self.frequency_type_folder, exist_ok=True)
        write_seam_check_log(
            filename, self.label, result, h, step_prefix, self.hessian_method
        )
        logger.info("Seam-minimum check results written to %s", filename)

    def _write_mecp_frequency_log(self, result, h):
        write_frequency_log(
            self.frequency_file,
            self.label,
            list(self.molecule.symbols),
            result,
            h,
            self.hessian_method,
        )
        logger.info(
            "MECP projected frequencies written to %s", self.frequency_file
        )

    def _final_report_outputs(self, result):
        return result.get("state_output_files", {})

    def _final_report_legacy_fields(self, result, frequency_only):
        return {}

    def _run_frequency_only(self):
        """Analyze the supplied crossing and use the common result report."""
        os.makedirs(self.folder, exist_ok=True)
        with open(self.final_report_file, "w", encoding="utf-8") as report:
            report.write(
                "CHEMSMART MECP frequency-only analysis\nstatus=RUNNING\n"
            )
        result = self.verify_seam_minimum(write_frequencies=True)
        self._last_seam_result = result
        self._final_optimization_steps = 0
        self._final_convergence_metrics = {
            "energy_diff": result["energy_diff"]
        }
        self.log_results(
            0,
            result["energy_a"],
            result["energy_b"],
            result,
            frequency_only=True,
        )
        return result

    def log_results(
        self,
        initial_steps=None,
        initial_energy_a=None,
        initial_energy_b=None,
        seam_result=None,
        frequency_only=False,
        *,
        energy_gap_tolerance=None,
    ):
        """Write a final MECP report for either backend.

        Supply optimization energies and optional seam-analysis results to
        record a driver-based result. With no arguments, regenerate the
        recorded report or read native output. Return its path, or None
        when no reportable result is available.
        """
        if initial_steps is not None:
            if initial_energy_a is None or initial_energy_b is None:
                raise ValueError("Both initial state energies are required.")
            if energy_gap_tolerance is not None:
                raise ValueError(
                    "Driver-based reports use the settings' convergence thresholds."
                )
            context = (
                initial_steps,
                initial_energy_a,
                initial_energy_b,
                seam_result,
                frequency_only,
            )
            path = self._write_final_report(*context)
            self._final_report_context = context
            return path
        if (
            initial_energy_a is not None
            or initial_energy_b is not None
            or seam_result is not None
            or frequency_only
        ):
            raise ValueError("Driver-based reports require initial_steps.")
        context = getattr(self, "_final_report_context", None)
        if context is not None:
            if energy_gap_tolerance is not None:
                raise ValueError(
                    "Driver-based reports use the settings' convergence thresholds."
                )
            return self._write_final_report(*context)
        return self._log_native_results(
            energy_gap_tolerance=energy_gap_tolerance
        )

    def _log_native_results(self, energy_gap_tolerance=None):
        """Backends with native MECP output implement result reporting."""
        return None

    def _write_final_report(
        self,
        initial_steps,
        initial_energy_a,
        initial_energy_b,
        seam_result,
        frequency_only=False,
    ):
        """Report the selected result using one schema for either backend."""
        result = seam_result or {}
        metrics = result.get(
            "convergence_metrics",
            getattr(self, "_final_convergence_metrics", {}),
        )
        thresholds = result.get(
            "convergence_thresholds"
        ) or convergence_thresholds(self.settings)
        if frequency_only:
            outputs = result.get("state_output_files", {})
        else:
            outputs = self._final_report_outputs(result)
        mode = getattr(self.settings, "mode", None)
        if mode is None:
            if getattr(self.settings, "follow_seam_imaginary_mode", False):
                mode = "followseam"
            elif frequency_only or getattr(
                self.settings, "mecp_numfreq", False
            ):
                mode = (
                    "anafreq"
                    if self.hessian_method == "analytic"
                    else "numfreq"
                )
            else:
                mode = "opt"
        write_final_report(
            self.final_report_file,
            label=self.label,
            backend=self.PROGRAM.lower(),
            hessian_method=self.hessian_method,
            symbols=self.molecule.symbols,
            positions=result.get(
                "positions_angstrom", self.molecule.positions
            ),
            energy_a=result.get("energy_a", initial_energy_a),
            energy_b=result.get("energy_b", initial_energy_b),
            initial_steps=initial_steps,
            metrics=metrics,
            thresholds=thresholds,
            result=result,
            frequency_only=frequency_only,
            macro_steps=getattr(self, "_selected_seam_follow_macro_steps", 0),
            selected_branch=getattr(
                self, "_selected_seam_follow_branch", None
            ),
            state_outputs=outputs,
            legacy_fields=self._final_report_legacy_fields(
                result, frequency_only
            ),
            mode=mode,
        )
        return self.final_report_file

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

        reduced_hessian = seam_basis.T @ mass_weighted_hessian @ seam_basis
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
        return (
            frequencies,
            cartesian_modes.T,
            eigenvalues,
            len(projection_vectors),
        )

    @staticmethod
    def lagrangian_hessian(hessian_a, hessian_b, grad_a, grad_b):
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

            This method runs **4 × 3N** electronic-structure sub-jobs (2 displacements
            × 2 spin states × 3N coordinates), which is expensive for large
            molecules (120 jobs for a 10-atom molecule).  Use the
            ``hess_step_size`` setting to control accuracy vs. cost.

        Args:
            positions_bohr (np.ndarray): Current geometry, shape (n_atoms, 3).
            h (float): Finite-difference step size in Bohr.
            step_prefix (int): Legacy argument retained for API compatibility.
                File labels identify atoms, Cartesian axes, and displacement signs.

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

            prefix = "" if macro_step is None else f"macro{macro_step:02d}_"
            axis = "xyz"[coord_idx]
            coord = f"{prefix}hessian_a{atom_idx + 1}_{axis}"
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

    def _compute_analytic_hessians(self, positions_bohr, macro_step=None):
        """Obtain both state Hessians using their reference checkpoints."""
        prefix = "" if macro_step is None else f"macro{macro_step:02d}_"
        self._analytic_hessians = {}
        for state in ("A", "B"):
            self._run_state(
                positions_bohr,
                f"{prefix}hessian_analytic",
                state,
                checkpoint_tag="check",
                analytic=True,
            )
        return self._analytic_hessians["A"], self._analytic_hessians["B"]

    def _get_state_hessians(
        self, positions_bohr, h, step_prefix, macro_step=None
    ):
        """Single Hessian-source dispatch for all seam verification calls."""
        if self.hessian_method == "analytic":
            return self._compute_analytic_hessians(positions_bohr, macro_step)
        return self._compute_numerical_hessian(
            positions_bohr,
            h=h,
            step_prefix=step_prefix,
            macro_step=macro_step,
        )

    def verify_seam_minimum(
        self,
        h=None,
        step_prefix=1,
        write_frequencies=False,
        macro_step=None,
        require_stationary=False,
    ):
        """
        Verify that the current MECP geometry is a **minimum on the crossing
        seam**, not merely a crossing point.

        The method obtains both state Hessians and forms the
        constrained Lagrangian Hessian. It then represents that Hessian in the
        subspace orthogonal to translations, rotations, and the
        gradient-difference direction:

        .. math::

            H_\\text{seam} = Q^T [(1-\\lambda)H_A + \\lambda H_B] Q

        The columns of the reduced-space basis span the seam tangent space.
        Any negative eigenvalue indicates a saddle point on the seam.

        Results are written to ``<label>_seam_check.log``.

        .. note::

            Numerical analysis requires 12N displaced force jobs plus two
            reference jobs. Analytic analysis uses two reference force jobs
            and two fixed-geometry analytic Hessian jobs. Both routes share
            the same constrained projection and mass weighting.

        Args:
            h (float, optional): Finite-difference step size in Bohr.
                Defaults to ``settings.hess_step_size`` (1 × 10⁻³ Bohr).
            step_prefix (int, optional): Legacy report metadata (default: 1).
                Sub-job labels use atom and Cartesian-axis indices.
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
            "hessian_reference"
            if macro_step is None
            else f"macro{macro_step:02d}_hessian_reference"
        )
        ea, grad_a = self._run_state(
            positions_bohr, reference_step, "A", checkpoint_tag="check"
        )
        eb, grad_b = self._run_state(
            positions_bohr, reference_step, "B", checkpoint_tag="check"
        )
        diff_grad = grad_a - grad_b

        if require_stationary or getattr(
            self.settings, "mecp_numfreq_only", False
        ):
            _, projected_grad, _ = self._mecp_displacement(
                ea - eb, grad_a, grad_b, self.settings.step_size
            )
            gap = abs(ea - eb)
            grad_max = float(np.max(np.abs(projected_grad)))
            grad_rms = self._rms(projected_grad)
            if (
                gap > self.settings.energy_diff_tol
                or grad_max > self.settings.force_max_tol
                or grad_rms > self.settings.force_rms_tol
            ):
                raise RuntimeError(
                    "Seam frequency analysis requires a stationary crossing: "
                    f"|energy gap|={gap:.6e} Eh "
                    f"(limit {self.settings.energy_diff_tol:.6e}); "
                    f"projected gradient max={grad_max:.6e} Eh/Bohr "
                    f"(limit {self.settings.force_max_tol:.6e}), "
                    f"RMS={grad_rms:.6e} Eh/Bohr "
                    f"(limit {self.settings.force_rms_tol:.6e}). "
                    "Input energy gap or projected gradient exceeds the convergence "
                    "thresholds; Hessians have not been computed. Optimize the "
                    "crossing to these thresholds or explicitly choose tolerances "
                    "consistent with the intended optimization accuracy."
                )

        logger.info(
            "verify_seam_minimum: computing %s Hessians for %s",
            self.hessian_method,
            self.label,
        )
        H_A, H_B = self._get_state_hessians(
            positions_bohr,
            h=h,
            step_prefix=step_prefix + 1,
            macro_step=macro_step,
        )

        H_lagrangian, lagrange_multiplier = self.lagrangian_hessian(
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
            "state_output_files": self._reference_output_files(reference_step),
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
        diff_grad = (
            np.asarray(grad_a, dtype=float).ravel()
            - np.asarray(grad_b, dtype=float).ravel()
        )
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
            step_index = f"macro{macro_step:02d}_opt_step{inner_step:03d}"
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

    def _new_seam_follow_branch(self, base_positions, label, folder):
        """Create an independent backend branch with the existing tight limits."""
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
            label=label,
            jobrunner=self.jobrunner,
            skip_completed=False,
            is_seam_follow_branch=True,
        )
        branch.set_folder(folder)
        return branch

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
        follow_folder = self.mecp_files.follow_folder
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
            branch = self._new_seam_follow_branch(
                base_positions, branch_label, follow_folder
            )
            trace_file = os.path.join(
                follow_folder, f"{branch_label}_record.log"
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
                branch_result["convergence_thresholds"] = (
                    convergence_thresholds(branch.settings)
                )
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
