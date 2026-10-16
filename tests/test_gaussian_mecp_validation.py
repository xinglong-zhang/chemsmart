"""Independent analytical checks for MECP curvature and thermochemistry.

These synthetic benchmarks validate implementation, not electronic-structure
accuracy or the physical accuracy of crossing-seam thermal corrections.
"""

from types import SimpleNamespace

import numpy as np
import pytest
from scipy.constants import Avogadro, Boltzmann, c, h, physical_constants

from chemsmart.analysis.thermochemistry import thermochemistry_from_file
from chemsmart.jobs.gaussian.mecp import GaussianMECPJob


@pytest.mark.parametrize(
    "a,b,crossing_scale,alpha",
    [
        (1.0, 4.0, 2.0, 0.5),
        (-1.0, 4.0, 2.0, 0.5),
        (-1.0, 0.5, 2.0, 0.5),
        (1.0, 4.0, 2.0, 0.0),
    ],
)
def test_curved_seam_matches_exact_energy(a, b, crossing_scale, alpha):
    # EA = a*x + b*y**2/2; EB = EA - c*(x-alpha*y**2).
    # On x=alpha*y**2 the exact curvature is b+2*a*alpha.
    grad_a = np.array([a, 0.0])
    grad_b = np.array([a - crossing_scale, 0.0])
    hessian_a = np.diag([0.0, b])
    hessian_b = np.diag([0.0, b + 2 * crossing_scale * alpha])
    actual, multiplier = GaussianMECPJob.lagrangian_hessian(
        hessian_a,
        hessian_b,
        grad_a,
        grad_b,
    )
    assert multiplier == pytest.approx(a / crossing_scale)
    assert actual[1, 1] == pytest.approx(b + 2 * a * alpha)


def _known_spectrum_model(linear=False):
    """Build a Hessian with known seam eigenvalues using independent SVD."""
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0] if linear else [0.0, 1.0, 0.0],
        ]
    )
    masses = np.array([1.0, 12.0, 16.0])
    sqrt_mass = np.repeat(np.sqrt(masses), 3)
    centered = positions - np.average(positions, axis=0, weights=masses)
    external = []
    for axis in np.eye(3):
        external.append(
            (np.tile(axis, (3, 1)) * np.sqrt(masses)[:, None]).ravel()
        )
        external.append(
            (np.cross(centered, axis) * np.sqrt(masses)[:, None]).ravel()
        )
    u, singular, _ = np.linalg.svd(
        np.column_stack(external), full_matrices=True
    )
    rank = np.count_nonzero(singular > 1.0e-10)
    normal = u[:, rank]
    basis = u[:, rank + 1 :]
    spectrum = np.linspace(0.01, 0.03, basis.shape[1])
    weighted = (basis * spectrum) @ basis.T
    hessian = weighted * np.outer(sqrt_mass, sqrt_mass)
    diff_grad = (normal * sqrt_mass).reshape(3, 3)
    return positions, masses, hessian, diff_grad, spectrum, rank + 1


def _frequencies(positions, masses, hessian, diff_grad):
    job = object.__new__(GaussianMECPJob)
    job.molecule = SimpleNamespace(most_abundant_masses=masses)
    return job._projected_frequencies_and_modes(hessian, positions, diff_grad)


@pytest.mark.parametrize("linear", [False, True])
def test_mass_weighted_projection_has_known_spectrum(linear):
    positions, masses, hessian, diff_grad, expected, excluded = (
        _known_spectrum_model(linear)
    )
    frequencies, modes, eigenvalues, n_projected = _frequencies(
        positions, masses, hessian, diff_grad
    )
    np.testing.assert_allclose(
        eigenvalues, expected, rtol=1.0e-10, atol=1.0e-12
    )
    assert n_projected == excluded
    assert len(frequencies) == (3 if linear else 2)
    hartree = physical_constants["Hartree energy"][0]
    bohr = physical_constants["Bohr radius"][0]
    amu = physical_constants["atomic mass constant"][0]
    reference = np.sqrt(expected * hartree / (bohr**2 * amu)) / (
        2 * np.pi * c * 100
    )
    # ASE and scipy may use different editions of fundamental constants.
    np.testing.assert_allclose(frequencies, reference, rtol=2.0e-6)
    np.testing.assert_allclose(np.linalg.norm(modes, axis=1), 1.0)
    np.testing.assert_allclose(modes @ diff_grad.ravel(), 0.0, atol=1.0e-10)


def test_frequency_spectrum_is_translation_and_rotation_invariant():
    positions, masses, hessian, diff_grad, _, _ = _known_spectrum_model()
    original = _frequencies(positions, masses, hessian, diff_grad)
    rotation, _ = np.linalg.qr(
        np.array([[1.0, 2.0, 3.0], [3.0, 1.0, 2.0], [2.0, 3.0, 1.0]])
    )
    if np.linalg.det(rotation) < 0:
        rotation[:, 0] *= -1
    transform = np.kron(np.eye(3), rotation)
    changed = _frequencies(
        positions @ rotation.T + np.array([2.0, -3.0, 4.0]),
        masses,
        transform @ hessian @ transform.T,
        diff_grad @ rotation.T,
    )
    assert original[3] == changed[3]
    np.testing.assert_allclose(original[0], changed[0], rtol=1.0e-10)


def test_state_exchange_preserves_lagrangian_hessian_and_frequencies():
    positions, masses, target, diff_grad, _, _ = _known_spectrum_model()
    multiplier = -0.5
    grad_a = multiplier * diff_grad
    grad_b = grad_a - diff_grad
    perturbation = np.diag(np.linspace(0.01, 0.09, 9))
    hessian_a = target - multiplier * perturbation
    hessian_b = target + (1 - multiplier) * perturbation
    first, lam = GaussianMECPJob.lagrangian_hessian(
        hessian_a,
        hessian_b,
        grad_a,
        grad_b,
    )
    second, swapped = GaussianMECPJob.lagrangian_hessian(
        hessian_b,
        hessian_a,
        grad_b,
        grad_a,
    )
    assert swapped == pytest.approx(1 - lam)
    np.testing.assert_allclose(first, target, atol=1.0e-12)
    np.testing.assert_allclose(second, first, atol=1.0e-12)
    np.testing.assert_allclose(
        _frequencies(positions, masses, first, diff_grad)[0],
        _frequencies(positions, masses, second, -diff_grad)[0],
        rtol=1.0e-10,
    )


@pytest.mark.parametrize("linear", [False, True])
def test_mecp_file_thermal_contributions_match_independent_formulas(
    tmp_path, linear
):
    positions, masses, hessian, diff_grad, _, _ = _known_spectrum_model(linear)
    frequencies, modes, eigenvalues, excluded = _frequencies(
        positions, masses, hessian, diff_grad
    )
    job = object.__new__(GaussianMECPJob)
    job.folder = str(tmp_path)
    job.label = "reference"
    job.molecule = SimpleNamespace(symbols=["H", "C", "O"])
    result = dict(
        positions_angstrom=positions,
        atomic_masses=masses,
        energy_a=-1.0,
        energy_b=-1.0,
        mecp_energy=-1.0,
        energy_diff=0.0,
        lagrange_multiplier=0.5,
        n_projected=excluded,
        n_negative=0,
        frequencies=frequencies,
        modes=modes,
        multiplicity_a=1,
        multiplicity_b=3,
        rotational_symmetry_number=1,
    )
    job._write_mecp_frequency_log(result, 1.0e-3)
    analysis = thermochemistry_from_file(
        str(tmp_path / "reference_mecp_freq.log"), temperature=298.15
    )
    # Use rounded file frequencies, as the real reader does.
    stored = np.asarray(analysis.file_object.vibrational_frequencies)
    assert len(analysis.cleaned_frequencies) == len(frequencies)
    theta = h * c * 100 * stored / Boltzmann
    gas_constant = Avogadro * Boltzmann
    x = theta / 298.15
    zpe = gas_constant * np.sum(theta / 2)
    internal = gas_constant * np.sum(theta * (0.5 + 1 / np.expm1(x)))
    entropy = gas_constant * np.sum(x / np.expm1(x) - np.log(-np.expm1(-x)))
    assert analysis.zero_point_energy == pytest.approx(zpe, rel=2.0e-6)
    assert analysis.vibrational_internal_energy == pytest.approx(
        internal, rel=2.0e-6
    )
    assert analysis.vibrational_entropy == pytest.approx(entropy, rel=2.0e-6)


def test_frequency_only_analyzes_known_stationary_crossing(tmp_path):
    from ase import units

    from chemsmart.jobs.gaussian.settings import GaussianMECPJobSettings

    positions, masses, hessian, diff_grad, expected, _ = (
        _known_spectrum_model()
    )
    job = object.__new__(GaussianMECPJob)
    job.folder = str(tmp_path)
    job.label = "known"
    job.settings = GaussianMECPJobSettings(
        mecp_numfreq_only=True, multiplicity_a=3, multiplicity_b=1
    )
    job.molecule = SimpleNamespace(
        positions=positions * units.Bohr,
        most_abundant_masses=masses,
        symbols=["H", "C", "O"],
    )
    initial = job.molecule.positions.copy()
    job._state_checkpoint_files = {}
    job._run_state = lambda positions, step, state, **kwargs: (
        0.0,
        0.5 * diff_grad if state == "A" else -0.5 * diff_grad,
    )
    job._compute_numerical_hessian = lambda *args, **kwargs: (hessian, hessian)
    result = job._run()
    np.testing.assert_allclose(result["frequency_eigenvalues"], expected)
    np.testing.assert_array_equal(job.molecule.positions, initial)
    assert job._job_is_complete()
    assert not (tmp_path / "known_optimization").exists()
