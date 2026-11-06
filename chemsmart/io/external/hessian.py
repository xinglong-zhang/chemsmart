"""
Hessians for hosts that ask the external program for second derivatives
(e.g. Gaussian with ``opt=calcfc``, ``freq`` or ``irc``).

Analytic engine Hessians are used when available and consistent;
otherwise the Hessian is built from central finite differences of
gradients.
"""

import logging

import numpy as np

from chemsmart.utils.constants import bohr_to_angstrom

logger = logging.getLogger(__name__)

HESSIAN_MODES = ("auto", "analytic", "fd")
DEFAULT_FD_STEP = 0.005  # Angstrom


def _fd_column(engine, request, k, step):
    """Central-difference column k of the Hessian, in Hartree/Bohr^2."""
    flat = request.positions.reshape(-1)
    grads = []
    for sign in (1.0, -1.0):
        displaced = flat.copy()
        displaced[k] += sign * step
        _, grad = engine.energy_and_gradient(
            request.numbers,
            displaced.reshape(-1, 3),
            request.charge,
            request.multiplicity,
        )
        grads.append(np.asarray(grad).reshape(-1))
    step_bohr = step / bohr_to_angstrom
    return (grads[0] - grads[1]) / (2.0 * step_bohr)


def finite_difference_hessian(engine, request, step=DEFAULT_FD_STEP):
    """Symmetrised central-difference Hessian (6N gradient calls)."""
    n = 3 * request.natoms
    hessian = np.empty((n, n))
    for k in range(n):
        hessian[:, k] = _fd_column(engine, request, k, step)
    return 0.5 * (hessian + hessian.T)


def compute_hessian(engine, request, mode="auto", step=DEFAULT_FD_STEP):
    """
    Return (hessian, source) for the request.

    ``mode``:
        - ``"fd"``: finite differences only.
        - ``"analytic"``: engine Hessian; error if unavailable.
        - ``"auto"``: engine Hessian when available and consistent with a
          finite-difference spot check of one column, otherwise finite
          differences.
    """
    if mode not in HESSIAN_MODES:
        raise ValueError(f"hessian mode must be one of {HESSIAN_MODES}.")
    if mode == "fd":
        return finite_difference_hessian(engine, request, step), "fd"

    analytic = engine.analytic_hessian(
        request.numbers,
        request.positions,
        request.charge,
        request.multiplicity,
    )
    if analytic is None:
        if mode == "analytic":
            raise RuntimeError(
                "hessian='analytic' requested but the engine provides none."
            )
        return finite_difference_hessian(engine, request, step), "fd"
    analytic = 0.5 * (analytic + analytic.T)
    if mode == "analytic":
        return analytic, "analytic"

    # Spot check the stiffest column against finite differences. This
    # catches sign, unit or mass-weighting mistakes in engine Hessians.
    k = int(np.argmax(np.abs(np.diag(analytic))))
    reference = _fd_column(engine, request, k, step)
    tolerance = 0.05 * np.max(np.abs(reference)) + 1e-4
    if np.max(np.abs(analytic[:, k] - reference)) <= tolerance:
        return analytic, "analytic"
    logger.warning(
        "Analytic Hessian failed the finite-difference spot check; "
        "using finite differences."
    )
    return finite_difference_hessian(engine, request, step), "fd"
