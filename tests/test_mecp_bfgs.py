from io import StringIO
from types import SimpleNamespace

import numpy as np

from chemsmart.jobs.gaussian.mecp import GaussianMECPJob
from chemsmart.jobs.gaussian.writer import GaussianInputWriter


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
