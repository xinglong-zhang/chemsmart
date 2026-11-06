"""Fixtures for the program-neutral external backends
(:mod:`chemsmart.io.external`), shared by the Gaussian and ORCA external
interfaces."""

import sys

import numpy as np
import pytest

from chemsmart.io.external import ExternalRequest


@pytest.fixture()
def xtb_output_directory(tmp_path):
    """Directory holding xtb ``gradient`` and ``hessian`` files for H2."""
    (tmp_path / "gradient").write_text(
        "$grad\n"
        "  cycle =      1    SCF energy =    -0.9822404321   "
        "|dE/xyz| =  0.012345\n"
        "    0.00000000000000      0.00000000000000     -0.69926098000000  h\n"
        "    0.00000000000000      0.00000000000000      0.69926098000000  h\n"
        "   0.0000000000000E+00   0.0000000000000E+00  -1.2345678901234E-02\n"
        "   0.0000000000000E+00   0.0000000000000E+00   1.2345678901234E-02\n"
        "$end\n"
    )
    hessian = [[float(6 * i + j) for j in range(6)] for i in range(6)]
    body = "\n".join(" ".join(f"{v:.6f}" for v in row) for row in hessian)
    (tmp_path / "hessian").write_text(f"$hessian\n{body}\n$end\n")
    return str(tmp_path)


@pytest.fixture()
def fake_xtb_executable(tmp_path, xtb_output_directory, monkeypatch):
    """Stand-in ``xtb`` binary for :class:`XTBEngine` subprocess tests.

    It copies the canned ``gradient``/``hessian`` files from
    :func:`xtb_output_directory` into its working directory, records its
    arguments and ``coord.xyz`` in ``calls.log``, and exits with status 3
    when given ``--fail``. The test runs from ``tmp_path`` so the engine's
    temporary ``xtb_*`` folders stay out of the repository.
    """
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    log = tmp_path / "calls.log"
    script = bin_dir / "xtb"
    script.write_text(
        f"#!{sys.executable}\n"
        "import shutil, sys\n"
        f"log = {str(log)!r}\n"
        f"source = {xtb_output_directory!r}\n"
        "with open(log, 'a') as fh:\n"
        "    fh.write(' '.join(sys.argv[1:]) + '\\n')\n"
        "    fh.write(open('coord.xyz').read())\n"
        "if '--fail' in sys.argv:\n"
        "    print('xtb: boom')\n"
        "    sys.exit(3)\n"
        "for name in ('gradient', 'hessian'):\n"
        "    shutil.copy(f'{source}/{name}', name)\n"
    )
    script.chmod(0o755)
    monkeypatch.chdir(tmp_path)
    return str(script), log


@pytest.fixture()
def harmonic_hessian():
    """Stiff, coupled 2-atom Hessian in Hartree/Bohr^2."""
    hessian = np.diag([0.5, 0.4, 0.3, 0.6, 0.2, 0.7])
    hessian[0, 3] = hessian[3, 0] = -0.1
    return hessian


@pytest.fixture()
def harmonic_external_request():
    """A Hessian (deriv=2) request for two H atoms 1 Angstrom apart."""
    return ExternalRequest(
        deriv=2,
        charge=0,
        multiplicity=1,
        numbers=np.array([1, 1]),
        positions=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    )


@pytest.fixture()
def isolated_backend_registry(monkeypatch):
    """The external backend registry, restored after the test.

    Backends registered during the test disappear afterwards, and entry
    point plugins are treated as already loaded.
    """
    from chemsmart.io.external import backends

    monkeypatch.setattr(backends, "_REGISTRY", dict(backends._REGISTRY))
    monkeypatch.setattr(backends, "_plugins_loaded", True)
    return backends


@pytest.fixture()
def harmonic_backend(isolated_backend_registry, harmonic_hessian):
    """A registered ``harmonic`` backend returning a QuadraticEngine with
    an exact analytic Hessian (needs no external software)."""
    from tests.helpers import QuadraticEngine

    @isolated_backend_registry.register_backend(
        "harmonic", description="Harmonic test potential"
    )
    def build(options, calculator):
        return QuadraticEngine(harmonic_hessian, analytic=harmonic_hessian)

    return "harmonic"
