import numpy as np

from chemsmart.jobs.gaussian.mecp import GaussianMECPJob


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


def test_inverse_bfgs_update_skips_negative_curvature():
    inv_hessian = np.eye(3)
    delta_x = np.array([1.0, 0.0, 0.0])
    delta_g = np.array([-1.0, 0.0, 0.0])

    updated = GaussianMECPJob._update_inverse_hessian(
        inv_hessian, delta_x, delta_g
    )

    np.testing.assert_array_equal(updated, inv_hessian)
