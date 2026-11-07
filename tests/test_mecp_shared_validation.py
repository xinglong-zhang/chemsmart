"""Independent checks for shared MECP curvature and thermochemistry algorithms.

These synthetic benchmarks validate implementation, not electronic-structure
accuracy or the physical accuracy of crossing-seam thermal corrections.
"""

from types import SimpleNamespace

import numpy as np
import pytest
from scipy.constants import Avogadro, Boltzmann, c, h, physical_constants

from chemsmart.analysis.thermochemistry import thermochemistry_from_file
from chemsmart.jobs.gaussian.mecp import GaussianMECPJob
from chemsmart.jobs.gaussian.settings import GaussianMECPJobSettings


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
        hessian_a, hessian_b, grad_a, grad_b
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
        hessian_a, hessian_b, grad_a, grad_b
    )
    second, swapped = GaussianMECPJob.lagrangian_hessian(
        hessian_b, hessian_a, grad_b, grad_a
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
    job.settings = GaussianMECPJobSettings(hessian_method="numerical")
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


@pytest.mark.parametrize(
    "difference,mode,message",
    [
        ([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], "zero norm"),
        ([1.0, 0.0, 0.0], [1.0, 0.0, 0.0], "linearly dependent"),
    ],
)
def test_seam_progress_rejects_degenerate_constraints(
    difference, mode, message
):
    from chemsmart.jobs.mecp import MECPSeamMixin

    with pytest.raises(RuntimeError, match=message):
        MECPSeamMixin._constrained_seam_displacement(
            energy_diff=0.1,
            grad_a=np.asarray(difference).reshape(1, 3),
            grad_b=np.zeros((1, 3)),
            progress_error=0.0,
            progress_mode=np.asarray(mode),
            step_size=0.1,
        )


def test_seam_mode_displacement_rejects_incompatible_geometry():
    from chemsmart.jobs.mecp import MECPSeamMixin

    with pytest.raises(ValueError, match="does not match"):
        MECPSeamMixin._displace_along_mode(np.zeros((2, 3)), np.ones(3), 0.1)


def test_new_execution_invalidates_previous_final_report(make_driver_job):
    from pathlib import Path

    job = make_driver_job()
    Path(job.final_report_file).write_text("status=COMPLETED\n")
    job._last_seam_result = {"energy_a": -1.0}
    job._begin_execution()
    assert not Path(job.final_report_file).exists()
    assert job._last_seam_result is None
    with pytest.raises(RuntimeError, match="Difference gradient is too small"):
        job._mecp_displacement(0.0, np.zeros((1, 3)), np.zeros((1, 3)), 0.1)
    with pytest.raises(ValueError, match="convergence thresholds"):
        job.log_results(1, -1.0, -1.0, energy_gap_tolerance=0.1)


@pytest.mark.parametrize("result", [{}, {"frequencies": [-20.0], "modes": []}])
def test_seam_follow_requires_complete_projected_modes(
    make_driver_job, result
):
    job = make_driver_job()
    with pytest.raises(
        RuntimeError, match="requires projected frequency modes"
    ):
        job._follow_seam_imaginary_mode(result)


def test_constrained_seam_relaxation_converges_and_reports_exhaustion(
    make_driver_job,
):
    job = make_driver_job(max_steps=1)
    job._run_state = lambda *args, **kwargs: (
        (-1.0, np.array([[1.0, 0.0, 0.0]]))
        if args[2] == "A"
        else (-1.0, np.array([[-1.0, 0.0, 0.0]]))
    )
    coordinates = np.zeros((1, 3))
    mode = np.array([0.0, 1.0, 0.0])
    from io import StringIO

    from ase import units

    relaxed = job._optimize_on_seam_progress_plane(
        coordinates, coordinates, mode, StringIO(), 1
    )
    np.testing.assert_allclose(relaxed, coordinates)
    np.testing.assert_allclose(
        job.molecule.positions, coordinates * units.Bohr
    )
    job._is_converged = lambda **kwargs: False
    with pytest.raises(RuntimeError, match="did not converge"):
        job._optimize_on_seam_progress_plane(
            coordinates, coordinates, mode, StringIO(), 1
        )


def test_seam_follow_stops_when_negative_mode_persists(
    make_driver_job, monkeypatch
):
    from pathlib import Path

    job = make_driver_job(seam_mode_max_steps=2)
    result = {
        "positions_angstrom": np.zeros((1, 3)),
        "frequencies": np.array([-20.0]),
        "frequency_eigenvalues": np.array([-0.01]),
        "modes": np.array([[1.0, 0.0, 0.0]]),
        "mecp_energy": -1.0,
        "n_negative": 1,
    }
    calls = []

    def branch_factory(*args):
        branch = SimpleNamespace(
            molecule=SimpleNamespace(positions=np.zeros((1, 3))),
            _last_seam_result=None,
        )
        branch._optimize_on_seam_progress_plane = lambda **kwargs: kwargs[
            "positions_bohr"
        ]

        def verify(**kwargs):
            calls.append(kwargs["macro_step"])
            return result

        branch.verify_seam_minimum = verify
        return branch

    monkeypatch.setattr(job, "_new_seam_follow_branch", branch_factory)
    with pytest.raises(
        RuntimeError, match="Neither projected-mode displacement"
    ):
        job._follow_seam_imaginary_mode(result)
    assert calls == [1, 2, 1, 2]
    records = list(Path(job.mecp_files.follow_folder).glob("*_record.log"))
    assert len(records) == 2
    assert all("tracked_mode=1" in path.read_text() for path in records)


def test_seam_tracking_orients_an_opposite_signed_mode():
    from chemsmart.jobs.mecp import MECPSeamMixin

    result = {
        "modes": np.array([[-1.0, 0.0, 0.0]]),
        "frequency_eigenvalues": np.array([-0.1]),
    }
    index, mode, overlap = MECPSeamMixin._overlap_tracked_negative_mode(
        result, np.array([1.0, 0.0, 0.0])
    )
    assert index == 0
    np.testing.assert_allclose(mode, [1.0, 0.0, 0.0])
    assert overlap == pytest.approx(1.0)


def test_final_report_without_optional_metadata_is_readable(tmp_path):
    from chemsmart.io.mecp import read_final_result, write_final_report

    path = tmp_path / "report.log"
    write_final_report(
        str(path),
        label="crossing",
        backend="gaussian",
        hessian_method="numerical",
        symbols=["H"],
        positions=np.zeros((1, 3)),
        energy_a=-1.0,
        energy_b=-1.0,
        initial_steps=1,
        metrics={},
        thresholds={},
    )
    text = path.read_text()
    assert "status=COMPLETED" in text
    assert "mode=" not in text
    result = read_final_result(str(path))
    assert result.energy_a == -1.0
