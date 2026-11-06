"""
Host protocol base: how a host program (Gaussian, ORCA, ...) exchanges
one energy/gradient request with a chemsmart launcher.

A host subclasses :class:`ExternalProtocol` and implements three methods:

- :meth:`~ExternalProtocol.parse_arguments`: the launcher's command line
  into an invocation object (no file I/O, so failures can still be
  reported to the host);
- :meth:`~ExternalProtocol.read_request`: the host's input files into an
  :class:`ExternalRequest`;
- :meth:`~ExternalProtocol.write_result`: an :class:`ExternalResult`
  into the files the host reads back.

Optionally it overrides :meth:`~ExternalProtocol.report_success` and
:meth:`~ExternalProtocol.report_failure` (default: stdout/stderr). The
shared driver (:meth:`~ExternalProtocol.run`,
:meth:`~ExternalProtocol.main`) builds the engine, evaluates energy and
gradient, and supplies a Hessian when the request asks for one.
"""

import sys
import time
import traceback
from dataclasses import dataclass
from typing import Optional

import numpy as np

from chemsmart.io.external.backends import build_engine
from chemsmart.io.external.hessian import DEFAULT_FD_STEP, compute_hessian
from chemsmart.io.external.runlog import (
    ExternalRunLog,
    capture_warnings,
    safely,
)


@dataclass
class ExternalRequest:
    """One energy/gradient request, with positions in Angstrom.

    ``deriv`` is the derivative order wanted: 0 energy, 1 gradient,
    2 Hessian.
    """

    deriv: int
    charge: int
    multiplicity: int
    numbers: np.ndarray
    positions: np.ndarray

    @property
    def natoms(self):
        return len(self.numbers)


@dataclass
class ExternalResult:
    """Answer to an :class:`ExternalRequest`, in atomic units.

    ``gradient`` (N, 3) is ``None`` for energy-only requests; ``hessian``
    (3N, 3N) is set only for ``deriv == 2``, with ``hessian_source``
    ``"analytic"`` or ``"fd"``.
    """

    energy: float
    gradient: Optional[np.ndarray] = None
    hessian: Optional[np.ndarray] = None
    hessian_source: Optional[str] = None

    def summary(self, backend):
        message = f"{backend}: E = {self.energy:.10f} Eh"
        if self.hessian_source is not None:
            message += f"; Hessian from {self.hessian_source}"
        return message


def evaluate(engine, request, hessian="auto", fd_step=DEFAULT_FD_STEP):
    """Evaluate ``request`` with ``engine``; return an ExternalResult."""
    energy, gradient = engine.energy_and_gradient(
        request.numbers,
        request.positions,
        request.charge,
        request.multiplicity,
    )
    result = ExternalResult(
        energy=energy, gradient=gradient if request.deriv >= 1 else None
    )
    if request.deriv >= 2:
        result.hessian, result.hessian_source = compute_hessian(
            engine, request, mode=hessian, step=fd_step
        )
    return result


class ExternalProtocol:
    """Exchange of one request between a host program and the launcher.

    Subclass once per host; see the module docstring. ``usage`` is shown
    when the launcher is called with unusable arguments.
    """

    host = None
    usage = None

    def parse_arguments(self, argv):
        """Return an invocation object; raise ``ValueError`` if unusable."""
        raise NotImplementedError

    def read_request(self, invocation):
        """Return the :class:`ExternalRequest` the host wrote."""
        raise NotImplementedError

    def write_result(self, invocation, request, result):
        """Write ``result`` where the host expects it."""
        raise NotImplementedError

    def report_success(self, invocation, message):
        sys.stdout.write(message + "\n")

    def report_failure(self, invocation, error):
        sys.stderr.write("chemsmart external calculation failed:\n")
        sys.stderr.write(error)

    def answer(self, config, invocation, engine=None):
        """Read, evaluate and write one request.

        Returns ``(request, result, engine)``. When
        ``config["log_backend_output"]`` is set, the engine keeps its
        programs' output for the per-job log.
        """
        request = self.read_request(invocation)
        engine = engine or build_engine(config)
        engine.keep_output = bool(config.get("log_backend_output"))
        result = evaluate(
            engine,
            request,
            hessian=config.get("hessian", "auto"),
            fd_step=config.get("fd_step", DEFAULT_FD_STEP),
        )
        self.write_result(invocation, request, result)
        return request, result, engine

    def run(self, config, invocation, engine=None):
        """Answer one request; return a status message."""
        _, result, _ = self.answer(config, invocation, engine)
        return result.summary(config["backend"])

    def main(self, config, argv=None):
        """Launcher entry point; returns the process exit code.

        Each request is recorded in the per-job external log when the
        host job was started with one (see
        :mod:`chemsmart.io.external.runlog`).
        """
        argv = sys.argv[1:] if argv is None else list(argv)
        try:
            invocation = self.parse_arguments(argv)
        except ValueError as exc:
            sys.stderr.write(f"{exc}\n")
            return 1
        run_log = ExternalRunLog.from_environment(self.host, config)
        start = time.perf_counter()
        with capture_warnings() as warnings:
            try:
                request, result, engine = self.answer(config, invocation)
            except Exception:
                error = traceback.format_exc()
                if run_log is not None:
                    safely(
                        run_log.record_failure,
                        error,
                        time.perf_counter() - start,
                        warnings,
                    )
                self.report_failure(invocation, error)
                return 1
        if run_log is not None:
            safely(
                run_log.record,
                request,
                result,
                time.perf_counter() - start,
                warnings,
                engine.drain_output(),
            )
        self.report_success(invocation, result.summary(config["backend"]))
        return 0
