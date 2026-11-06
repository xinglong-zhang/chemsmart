"""
Per-job log of an external level of theory.

The host program (Gaussian, ORCA, ...) starts a new launcher process for
every energy/gradient request, so each request appends to a log file
named by the :data:`EXTERNAL_LOG_ENV` environment variable. Chemsmart's
job runners set it to ``<running directory>/<label>_external.log`` for
external jobs; the file is copied back from scratch with the job's other
``<label>*`` files. Without the variable nothing is logged.

Each request adds one line::

    2026-11-06 22:10:03  step    3  deriv 1  E = -40.123456789012 Eh  max|g| = 1.234e-03  rms|g| = 4.567e-04 Eh/Bohr  hessian -  0.84 s

followed, when present, by indented backend warnings and (with
``log_backend_output``) the backend program's output. Failed requests
are logged with their traceback. Records are appended under a file lock,
so concurrent requests (e.g. ORCA computing a numerical Hessian in
parallel) neither interleave nor reuse step numbers.
"""

import json
import logging
import os
import re
import sys
from contextlib import contextmanager
from datetime import datetime

import numpy as np

from chemsmart.utils.repattern import external_log_step_pattern

try:  # POSIX; on other platforms records are appended without a lock
    import fcntl
except ImportError:  # pragma: no cover - Windows
    fcntl = None

EXTERNAL_LOG_ENV = "CHEMSMART_EXTERNAL_LOG"

_INDENT = "    "


@contextmanager
def capture_warnings(logger_name="chemsmart"):
    """Collect warnings logged under ``logger_name`` while active.

    Yields a list that fills with the formatted messages. The messages
    are still handled as usual (e.g. printed to stderr).
    """
    messages = []

    class _Collector(logging.Handler):
        def emit(self, record):
            messages.append(
                f"{record.levelname.lower()}: {record.getMessage()}"
            )

    handler = _Collector(level=logging.WARNING)
    target = logging.getLogger(logger_name)
    target.addHandler(handler)
    try:
        yield messages
    finally:
        target.removeHandler(handler)


def gradient_statistics(gradient):
    """Return (max |component|, RMS) of a gradient, or (None, None)."""
    if gradient is None:
        return None, None
    values = np.asarray(gradient, dtype=float).reshape(-1)
    if values.size == 0:
        return None, None
    return float(np.max(np.abs(values))), float(np.sqrt(np.mean(values**2)))


class ExternalRunLog:
    """Append-only log of the requests answered for one host job.

    Args:
        path (str): Log file.
        host (str): Host program name, for the header.
        config (dict): Launcher configuration, for the header.
    """

    def __init__(self, path, host, config):
        self.path = path
        self.host = host
        self.config = config

    @classmethod
    def from_environment(cls, host, config, environ=None):
        """The log named by :data:`EXTERNAL_LOG_ENV`, or ``None``."""
        environ = os.environ if environ is None else environ
        path = environ.get(EXTERNAL_LOG_ENV)
        return cls(path, host, config) if path else None

    @property
    def include_backend_output(self):
        return bool(self.config.get("log_backend_output"))

    # ── records ──────────────────────────────────────────────────────────

    def record(self, request, result, elapsed, warnings=(), outputs=()):
        """Log one answered request."""
        max_g, rms_g = gradient_statistics(result.gradient)
        gradient = (
            "max|g| = -          rms|g| = -         "
            if max_g is None
            else f"max|g| = {max_g:.3e}  rms|g| = {rms_g:.3e}"
        )
        line = (
            f"deriv {request.deriv}  E = {result.energy:.12f} Eh  "
            f"{gradient} Eh/Bohr  hessian {result.hessian_source or '-'}  "
            f"{elapsed:.2f} s"
        )
        self._append(line, self._details(warnings, outputs))

    def record_failure(self, error, elapsed, warnings=(), outputs=()):
        """Log a request that raised; ``error`` is the traceback text."""
        last = error.strip().splitlines()[-1] if error.strip() else "error"
        details = self._details(warnings, outputs)
        details += [f"{_INDENT}{line}" for line in error.rstrip().splitlines()]
        self._append(f"FAILED after {elapsed:.2f} s: {last}", details)

    def _details(self, warnings, outputs):
        details = [f"{_INDENT}{message}" for message in warnings]
        if self.include_backend_output:
            for label, text in outputs:
                details.append(f"{_INDENT}--- output of: {label} ---")
                details += [
                    f"{_INDENT}{line}" for line in text.rstrip().splitlines()
                ]
                details.append(f"{_INDENT}--- end of output ---")
        return details

    # ── file handling ────────────────────────────────────────────────────

    def _header(self):
        config = json.dumps(self.config, sort_keys=True, default=str)
        return [
            f"# chemsmart external log: host {self.host}, "
            f"backend {self.config.get('backend')}",
            f"# config: {config}",
            "# energies in Hartree, gradients in Hartree/Bohr; one line per "
            "request from the host program",
        ]

    def _append(self, line, details):
        """Append a record, numbering it under an exclusive lock."""
        directory = os.path.dirname(os.path.abspath(self.path))
        os.makedirs(directory, exist_ok=True)
        with open(self.path, "a+") as fh:
            if fcntl is not None:
                fcntl.flock(fh, fcntl.LOCK_EX)
            try:
                fh.seek(0)
                existing = fh.read()
                step = 1 + len(
                    re.findall(external_log_step_pattern, existing, re.M)
                )
                lines = [] if existing else self._header()
                timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                lines.append(f"{timestamp}  step {step:4d}  {line}")
                lines += details
                fh.seek(0, os.SEEK_END)
                fh.write("\n".join(lines) + "\n")
                fh.flush()
            finally:
                if fcntl is not None:
                    fcntl.flock(fh, fcntl.LOCK_UN)


def safely(write, *args, **kwargs):
    """Run a log write; a logging problem must never fail the request."""
    try:
        write(*args, **kwargs)
    except Exception as exc:
        sys.stderr.write(f"chemsmart external log not written: {exc}\n")
