"""
Run Gaussian optimisers (opt, TS, modred, scan, IRC, ...) on an external
level of theory: machine-learning potentials, xTB, or any ASE calculator.

:class:`GaussianExternalCalculator` is set as the ``external`` attribute
of :class:`~chemsmart.jobs.gaussian.settings.GaussianJobSettings`. When
present, it replaces the QM method in the route with
``External="<python> <launcher>"``, and the Gaussian job runner writes the
launcher script next to the ``.com`` file. The launcher hands each
request to :func:`chemsmart.io.gaussian.external.main`.
"""

from chemsmart.io.external import DEFAULT_FD_STEP, HESSIAN_MODES
from chemsmart.jobs.external import ExternalCalculator

__all__ = ["GaussianExternalCalculator"]


class GaussianExternalCalculator(ExternalCalculator):
    """
    Backend configuration for Gaussian's ``External=`` keyword.

    Takes the arguments of
    :class:`~chemsmart.jobs.external.ExternalCalculator`, plus:

    Args:
        hessian (str): ``auto`` (engine Hessian if available and
            consistent, else finite differences), ``analytic`` or ``fd``.
            Gaussian asks the External program for the Hessian
            (``opt=calcfc``, TS searches, ``freq``, IRC).
        fd_step (float): Finite-difference displacement in Angstrom.
    """

    host = "Gaussian"
    runtime_module = "chemsmart.io.gaussian.external"
    script_prefix = "chemsmart_external"

    def __init__(
        self,
        backend,
        calculator=None,
        options=None,
        hessian="auto",
        fd_step=DEFAULT_FD_STEP,
        python_executable=None,
        keep_backend_output=False,
    ):
        super().__init__(
            backend,
            calculator=calculator,
            options=options,
            python_executable=python_executable,
            keep_backend_output=keep_backend_output,
        )
        if hessian not in HESSIAN_MODES:
            raise ValueError(
                f"hessian must be one of {', '.join(HESSIAN_MODES)}."
            )
        if fd_step <= 0:
            raise ValueError("fd_step must be positive.")
        self.hessian = hessian
        self.fd_step = float(fd_step)

    def _extra_fields(self):
        return {"hessian": self.hessian, "fd_step": self.fd_step}

    @property
    def route_value(self):
        """Value for ``External=`` in the route, always quoted."""
        return f'"{self.python} {self.script_filename}"'
