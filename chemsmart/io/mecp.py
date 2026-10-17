"""ChemSmart MECP file conventions and backend-independent reports.

State A/B means state 1/2 for both engines. Energies are in Hartree,
Cartesian coordinates in Angstrom, gradients in Hartree/Bohr, and
projected frequencies in cm^-1. Native engine files retain their formats.
"""

import os
from dataclasses import dataclass

import numpy as np

from chemsmart.utils.constants import (
    MECP_FREQUENCY_HEADER,
    MECP_FREQUENCY_TERMINATION_MARKER,
)

CONVERGENCE_CRITERIA = {
    "energy_diff": ("energy_diff_tol", "Hartree"),
    "pgrad_max": ("force_max_tol", "Hartree/Bohr"),
    "pgrad_rms": ("force_rms_tol", "Hartree/Bohr"),
    "disp_max": ("disp_max_tol", "Bohr"),
    "disp_rms": ("disp_rms_tol", "Bohr"),
}


def state_job_label(label, step, state, checkpoint_tag=None):
    if not isinstance(step, str):
        prefix = f"{checkpoint_tag}_" if checkpoint_tag else ""
        step = f"{prefix}step{step}"
    return f"{label}_{step}_{state}"


@dataclass(frozen=True)
class MECPFiles:
    """One naming convention for ChemSmart files and engine sub-jobs."""

    folder: str
    label: str

    def path(self, suffix):
        return os.path.join(self.folder, f"{self.label}_{suffix}")

    @property
    def optimization_folder(self):
        return self.path("optimization")

    @property
    def final_report(self):
        return self.path("final_report.log")

    @property
    def frequencies(self):
        return self.path("mecp_freq.log")

    @property
    def follow_folder(self):
        return self.path("seam_follow")

    def hessian_folder(self, method):
        return self.path(f"{method}_frequency")

    def seam_check(self, method):
        return os.path.join(
            self.hessian_folder(method), f"{self.label}_seam_check.log"
        )

    def state_label(self, step, state, checkpoint_tag=None):
        return state_job_label(self.label, step, state, checkpoint_tag)

    def reference_outputs(self, reference_step, method, extension):
        return {
            state: os.path.join(
                self.hessian_folder(method),
                self.state_label(reference_step, state) + extension,
            )
            for state in ("A", "B")
        }


def convergence_thresholds(settings):
    return {
        name: getattr(settings, attribute)
        for name, (attribute, _) in CONVERGENCE_CRITERIA.items()
    }


def read_report_fields(filename):
    """Read key=value metadata; later fields override earlier ones."""
    if not os.path.isfile(filename):
        return {}
    with open(filename, encoding="utf-8") as stream:
        return dict(
            line.strip().split("=", 1) for line in stream if "=" in line
        )


def state_energy_lines(energy_a, energy_b, energy_difference=None):
    """Canonical two-state energy fields, including explicit units."""
    both_known = energy_a is not None and energy_b is not None
    if energy_difference is None and both_known:
        energy_difference = energy_a - energy_b
    values = {
        "energy_A": energy_a,
        "energy_B": energy_b,
        "energy_difference": energy_difference,
        "mecp_energy": 0.5 * (energy_a + energy_b) if both_known else None,
    }
    lines = []
    for name, value in values.items():
        if value is None:
            lines.append(f"{name}=NA")
        elif name == "energy_difference":
            lines.append(f"{name}={value:+.6e} Hartree")
        else:
            lines.append(f"{name}={value:+.12f} Hartree")
    return lines


def write_final_report(
    filename,
    *,
    label,
    backend,
    hessian_method,
    symbols,
    positions,
    energy_a,
    energy_b,
    initial_steps,
    metrics,
    thresholds,
    result=None,
    frequency_only=False,
    macro_steps=0,
    selected_branch=None,
    state_outputs=None,
    legacy_fields=None,
    mode=None,
):
    """Write one final-result schema for optimization and frequency-only jobs.

    Missing native convergence metrics are recorded as NA. Legacy fields
    preserve existing consumers without giving the canonical fields different
    meanings between backends.
    """
    result = result or {}
    checked = "is_minimum" in result
    minimum = result.get("is_minimum")
    seam_status = ("PASS" if minimum else "FAIL") if checked else "NOT_CHECKED"
    with open(filename, "w", encoding="utf-8") as output:
        output.write("CHEMSMART final MECP result\n")
        output.write(
            f"label={label}\nbackend={backend}\n"
            f"hessian_method={hessian_method}\n"
            f"optimization={'SKIPPED' if frequency_only else 'CONVERGED'}\n"
            f"initial_optimization_steps={initial_steps}\n"
            f"seam_follow_macro_steps={macro_steps}\n"
            f"seam_follow_selected_branch={selected_branch or 'none'}\n"
        )
        if mode is not None:
            output.write(f"mode={mode}\n")
        for state, path in (state_outputs or {}).items():
            output.write(
                f"state_output_{state}={os.path.relpath(path, os.path.dirname(filename))}\n"
            )
        for name, value in (legacy_fields or {}).items():
            output.write(f"{name}={value}\n")
        if macro_steps:
            output.write(
                f"final_branch_optimization_steps={result.get('optimization_steps', 'NA')}\n"
            )
        output.write("\n".join(state_energy_lines(energy_a, energy_b)) + "\n")
        output.write("\nFinal convergence criteria:\n")
        for name, threshold in thresholds.items():
            value = metrics.get(name)
            unit = CONVERGENCE_CRITERIA[name][1]
            if value is None:
                output.write(
                    f"{name}: value=NA threshold={threshold:.6e} {unit} status=NA\n"
                )
                continue
            value = abs(value) if name == "energy_diff" else value
            output.write(
                f"{name}: value={value:.6e} threshold={threshold:.6e} {unit} "
                f"status={'PASS' if value <= threshold else 'FAIL'}\n"
            )
        output.write(f"\nseam_minimum={seam_status}\n")
        if checked:
            output.write(
                f"is_minimum={minimum}\n"
                f"significant_imaginary_modes={result['n_negative']}\n"
                f"frequency_file={label}_mecp_freq.log\n"
            )
        output.write("\nFinal geometry (Angstrom):\n")
        for symbol, position in zip(symbols, np.asarray(positions)):
            output.write(
                f"{symbol:>3s} {position[0]:+16.8f} "
                f"{position[1]:+16.8f} {position[2]:+16.8f}\n"
            )
        if "frequencies" in result:
            output.write("\nProjected MECP frequencies (cm^-1):\n")
            for index, frequency in enumerate(result["frequencies"], 1):
                output.write(f"mode {index:4d}: {frequency:+14.6f}\n")
        # Retain the existing Gaussian frequency-only completion marker.
        if frequency_only:
            output.write("NUMFREQ ONLY COMPLETED\n")
        output.write("status=COMPLETED\n")


def write_seam_check_log(filename, label, result, h, step_prefix, method):
    """Write the seam-minimum verification results to a log file."""
    with open(filename, "w", encoding="utf-8") as f:
        f.write("CHEMSMART MECP seam-minimum verification\n")
        f.write(f"hessian_method={method}\n")
        f.write(
            f"label={label} "
            + (f"hess_step={h:.2e} Bohr " if method == "numerical" else "")
            + f"check_step_start={step_prefix}\n"
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


def write_frequency_log(filename, label, symbols, result, h, method):
    """Write projected MECP frequencies and Cartesian normal modes."""
    positions = np.asarray(result["positions_angstrom"], dtype=float)
    masses = np.asarray(result["atomic_masses"], dtype=float)

    with open(filename, "w", encoding="utf-8") as output:
        output.write(f"{MECP_FREQUENCY_HEADER}\n")
        output.write(f"label={label}\n")
        output.write(f"hessian_method={method}\n")
        if method == "numerical":
            output.write(f"hess_step={h:.6e} Bohr\n")
        output.write(f"energy_A={result['energy_a']:+.12f} Hartree\n")
        output.write(f"energy_B={result['energy_b']:+.12f} Hartree\n")
        output.write(f"mecp_energy={result['mecp_energy']:+.12f} Hartree\n")
        output.write(
            f"energy_difference={result['energy_diff']:+.6e} Hartree\n"
        )
        output.write(
            f"lagrange_multiplier=" f"{result['lagrange_multiplier']:+.8e}\n"
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
            "\nCartesian normal modes " "(mass-unweighted, unit-normalized):\n"
        )
        for mode_index, mode in enumerate(result["modes"], 1):
            frequency = result["frequencies"][mode_index - 1]
            output.write(
                f"mode {mode_index:4d} " f"frequency={frequency:+.6f} cm^-1\n"
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
