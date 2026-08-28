import os
from io import StringIO
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from chemsmart.cli.gaussian.mecp_options import add_mecp_method_suffix

from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.jobs.gaussian.mecp import GaussianMECPJob
from chemsmart.jobs.gaussian.runner import GaussianJobRunner
from chemsmart.jobs.gaussian.settings import GaussianMECPJobSettings
from chemsmart.jobs.gaussian.writer import GaussianInputWriter


def test_first_mecp_step_removes_unavailable_guess_read():
    remove_read = GaussianMECPJob._without_unavailable_guess_read

    assert remove_read("scf=xqc guess=read nosymm") == "scf=xqc nosymm"
    assert remove_read("guess=(mix,read) nosymm") == "guess=(mix) nosymm"
    assert remove_read("scf=xqc nosymm") == "scf=xqc nosymm"


def test_mecp_requires_nosymm_for_force_alignment():
    ensure_nosymm = GaussianMECPJob._with_required_nosymm

    assert ensure_nosymm("scf=xqc") == "scf=xqc nosymm"
    assert ensure_nosymm("NoSymm scf=xqc") == "NoSymm scf=xqc"
    with pytest.raises(ValueError, match="requires nosymm"):
        ensure_nosymm("symmetry=loose")


def test_gaussian_spin_squared_parser_returns_final_annihilated_value():
    output = object.__new__(Gaussian16Output)
    output.__dict__["contents"] = [
        "S**2 before annihilation     2.1040, after     2.0060",
        "S**2 before annihilation     2.0550, after     2.0015",
    ]

    assert output.spin_squared == pytest.approx(2.0015)


def test_gaussian_header_reuses_rolling_checkpoint():
    job = SimpleNamespace(
        label="mecp_step2_A",
        chkfile="steps/mecp_A.chk",
        oldchkfile="steps/mecp_A.chk",
        settings=SimpleNamespace(chk=True),
        jobrunner=SimpleNamespace(num_cores=4, mem_gb=8),
    )
    output = StringIO()

    GaussianInputWriter(job)._write_gaussian_header(output)

    assert output.getvalue().splitlines()[0] == "%chk=mecp_A.chk"
    assert "%oldchk" not in output.getvalue().lower()


def test_seam_checkpoint_cleanup_preserves_final_mecp_checkpoints(
    monkeypatch,
):
    job = object.__new__(GaussianMECPJob)
    job._state_checkpoint_files = {
        (None, "A"): "mecp_A.chk",
        (None, "B"): "mecp_B.chk",
        ("check", "A"): "mecp_check_A.chk",
        ("check", "B"): "mecp_check_B.chk",
    }
    removed = []
    monkeypatch.setattr("os.path.isfile", lambda path: True)
    monkeypatch.setattr("os.remove", removed.append)

    job._remove_checkpoint_set("check")

    assert removed == ["mecp_check_A.chk", "mecp_check_B.chk"]
    assert job._state_checkpoint_files == {
        (None, "A"): "mecp_A.chk",
        (None, "B"): "mecp_B.chk",
    }


def test_mecp_scratch_jobs_are_grouped_in_steps_folder():
    runner = object.__new__(GaussianJobRunner)
    runner._scratch_dir = "/scratch/project"
    job = SimpleNamespace(
        label="mecp_step2_A", scratch_parent_folder="mecp_steps"
    )

    scratch_directory = runner._scratch_job_directory(job)

    assert scratch_directory == os.path.join(
        "/scratch/project", "mecp_steps", "mecp_step2_A"
    )


def test_inverse_bfgs_update_satisfies_secant_condition():
    inv_hessian = np.diag([0.7, 1.1, 1.6])
    delta_x = np.array([0.12, -0.04, 0.08])
    delta_g = np.array([0.20, -0.03, 0.10])

    updated = GaussianMECPJob._update_inverse_hessian(
        inv_hessian, delta_x, delta_g
    )

    np.testing.assert_allclose(updated @ delta_g, delta_x, atol=1.0e-12)
    np.testing.assert_allclose(updated, updated.T, atol=1.0e-14)
    assert np.all(np.linalg.eigvalsh(updated) > 0.0)


def test_inverse_bfgs_update_accepts_negative_curvature_like_easymecp():
    inv_hessian = np.eye(3)
    delta_x = np.array([1.0, 0.0, 0.0])
    delta_g = np.array([-1.0, 0.0, 0.0])

    updated = GaussianMECPJob._update_inverse_hessian(
        inv_hessian, delta_x, delta_g
    )

    np.testing.assert_allclose(updated @ delta_g, delta_x)
    assert np.linalg.eigvalsh(updated).min() < 0.0


def test_inverse_bfgs_update_skips_singular_secant_pair():
    inv_hessian = np.eye(3)
    delta_x = np.zeros(3)
    delta_g = np.ones(3)

    updated = GaussianMECPJob._update_inverse_hessian(
        inv_hessian, delta_x, delta_g
    )

    np.testing.assert_array_equal(updated, inv_hessian)


def test_removed_harvey_bfgs_setting_is_rejected():
    with pytest.raises(ValueError, match="Unknown MECP step_size_method"):
        GaussianMECPJobSettings(step_size_method="harvey_bfgs")


@pytest.mark.parametrize(
    "kwargs",
    [
        {"max_steps": 0},
        {"hess_step_size": 0.0},
        {"step_size_min": 2.0, "step_size_max": 1.0},
        {"multiplicity_a": 0},
        {"multiplicity_a": 1, "multiplicity_b": 1},
    ],
)
def test_invalid_mecp_settings_are_rejected(kwargs):
    with pytest.raises(ValueError):
        GaussianMECPJobSettings(**kwargs)


def test_reduced_hessian_preserves_negative_seam_curvature():
    hessian = np.diag([-2.0, 0.0, 3.0])
    projected = [np.array([0.0, 1.0, 0.0])]

    reduced = GaussianMECPJob._reduced_hessian(hessian, projected)
    eigenvalues = np.linalg.eigvalsh(reduced)

    np.testing.assert_allclose(eigenvalues, [-2.0, 3.0], atol=1.0e-12)


def test_lagrangian_hessian_weight_is_not_forced_to_average():
    grad_a = np.array([2.0, 0.0])
    grad_b = np.array([-1.0, 0.0])
    hessian_a = np.diag([1.0, 2.0])
    hessian_b = np.diag([4.0, 8.0])
    hessian, multiplier = GaussianMECPJob._lagrangian_hessian(
        hessian_a, hessian_b, grad_a, grad_b
    )

    assert multiplier == pytest.approx(2.0 / 3.0)
    np.testing.assert_allclose(
        hessian, (1.0 / 3.0) * hessian_a + (2.0 / 3.0) * hessian_b
    )


def test_optimizer_state_round_trip(tmp_path):
    job = object.__new__(GaussianMECPJob)
    job.folder = str(tmp_path)
    job.label = "restartable"
    job.molecule = SimpleNamespace(
        symbols=["H", "H"], positions=np.zeros((2, 3))
    )
    job.settings = SimpleNamespace(step_size_method="harvey")
    positions = np.arange(6, dtype=float).reshape(2, 3)
    inverse_hessian = np.eye(6)

    job._save_optimizer_state(
        next_step=12,
        positions_bohr=positions,
        current_step_size=0.1,
        prev_merit=None,
        prev_positions=None,
        prev_proj_grad=None,
        inv_hessian=inverse_hessian,
        prev_eff_grad=np.ones(6),
        prev_positions_bfgs=positions - 0.1,
    )
    restored = job._load_optimizer_state()

    assert restored["next_step"] == 12
    assert restored["prev_merit"] is None
    assert restored["prev_positions"] is None
    np.testing.assert_array_equal(restored["positions_bohr"], positions)
    np.testing.assert_array_equal(restored["inv_hessian"], inverse_hessian)


def test_optimizer_state_rejects_method_change(tmp_path):
    job = object.__new__(GaussianMECPJob)
    job.folder = str(tmp_path)
    job.label = "restartable"
    job.molecule = SimpleNamespace(symbols=["H"], positions=np.zeros((1, 3)))
    job.settings = SimpleNamespace(step_size_method="harvey")
    job._save_optimizer_state(
        next_step=2,
        positions_bohr=np.zeros((1, 3)),
        current_step_size=0.1,
        prev_merit=None,
        prev_positions=None,
        prev_proj_grad=None,
        inv_hessian=np.eye(3),
        prev_eff_grad=np.zeros(3),
        prev_positions_bfgs=np.zeros((1, 3)),
    )
    job.settings.step_size_method = "bb"

    with pytest.raises(RuntimeError, match="optimizer method differs"):
        job._load_optimizer_state()


def test_job_is_complete_requires_post_verification_marker(tmp_path):
    job = object.__new__(GaussianMECPJob)
    job.folder = str(tmp_path)
    job.label = "verified"
    report = Path(job.report_file)
    report.write_text("Optimization converged at step 8.\n", encoding="utf-8")

    assert job._job_is_complete() is False

    report.write_text("Converged at step 8.\n", encoding="utf-8")
    assert job._job_is_complete() is True


@pytest.mark.parametrize("swap_states", [False, True])
def test_harvey_optimizer_converges_on_analytic_crossing(swap_states):
    job = object.__new__(GaussianMECPJob)
    job.settings = GaussianMECPJobSettings()
    positions = np.array([[0.7, 1.0, 0.0]])
    prev_positions = None
    prev_gradient = None
    inverse_hessian = None

    for _ in range(100):
        x, y, _ = positions[0]
        energy_a = 0.5 * ((x - 1.0) ** 2 + y**2)
        energy_b = 0.5 * ((x + 1.0) ** 2 + y**2)
        gradient_a = np.array([[x - 1.0, y, 0.0]])
        gradient_b = np.array([[x + 1.0, y, 0.0]])
        if swap_states:
            energy_a, energy_b = energy_b, energy_a
            gradient_a, gradient_b = gradient_b, gradient_a

        (
            displacement,
            effective_gradient,
            _,
            inverse_hessian,
            _,
            _,
            _,
            optimizer_gradient,
        ) = job._bfgs_displacement(
            energy_a,
            energy_b,
            gradient_a,
            gradient_b,
            prev_positions,
            positions,
            prev_gradient,
            inverse_hessian,
        )
        if job._is_converged(
            energy_a - energy_b, effective_gradient, displacement
        ):
            break
        prev_positions = positions.copy()
        prev_gradient = optimizer_gradient.ravel().copy()
        positions = positions + displacement
    else:
        pytest.fail("Harvey optimizer did not converge on analytic surfaces")

    np.testing.assert_allclose(positions, np.zeros((1, 3)), atol=1.0e-5)


@pytest.mark.parametrize(
    ("label", "method", "expected"),
    [
        ("sample", "bb", "sample_bb_mecp"),
        ("sample_mecp", "grow_shrink", "sample_grow_shrink_mecp"),
        ("sample_harvey_mecp", "harvey", "sample_harvey_mecp"),
    ],
)
def test_add_mecp_method_suffix(label, method, expected):
    assert add_mecp_method_suffix(label, method) == expected
