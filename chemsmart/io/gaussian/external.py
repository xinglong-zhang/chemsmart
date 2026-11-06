"""
Runtime side of Gaussian's ``External=`` interface.

Gaussian runs a generated launcher script (see
:class:`chemsmart.jobs.gaussian.external.GaussianExternalCalculator`)
once per energy/gradient/Hessian request. The launcher calls
:func:`main` with its backend configuration, which:

1. reads Gaussian's External input file (geometry, charge, multiplicity
   and the derivative order requested),
2. evaluates the energy and gradient with an engine (any ASE calculator,
   a machine-learning-potential preset, or the ``xtb`` binary),
3. supplies a Hessian when Gaussian asks for second derivatives
   (``opt=calcfc``, ``freq``, ``irc``), analytically when the engine
   provides one and by central finite differences of gradients otherwise,
4. writes the results back in Gaussian's External output format.

Gaussian invokes the launcher as::

    <launcher> [extra_args] layer InFile OutFile MsgFile FChkFile MatElFile

The six standard arguments are always the last six.

Input file (atomic units)::

    natoms  deriv  charge  multiplicity              (4I10)
    Z  x  y  z  mm_charge                            (I10,4F20.12) x natoms

Output file (atomic units)::

    energy  dipole_x  dipole_y  dipole_z             (4D20.12)
    dE/dx  dE/dy  dE/dz                              (3D20.12) x natoms
    # only when deriv == 2:
    polarizability (6 values, lower triangle)        (3D20.12)
    dipole derivatives (9 * natoms values)           (3D20.12)
    force constants (lower triangle, 3N(3N+1)/2)     (3D20.12)
"""

import logging
import sys
import traceback

import numpy as np

from chemsmart.io.external import (
    BOHR_TO_ANG,
    DEFAULT_FD_STEP,
    ExternalRequest,
    build_engine,
    compute_hessian,
)

logger = logging.getLogger(__name__)


# ── File formats ──────────────────────────────────────────────────────────


def read_external_input(filename):
    """Parse a Gaussian External input file into an ExternalRequest."""
    with open(filename) as fh:
        lines = fh.readlines()
    natoms, deriv, charge, multiplicity = (
        int(v) for v in lines[0].split()[:4]
    )
    numbers = []
    positions = []
    for line in lines[1 : natoms + 1]:
        parts = line.replace("D", "E").split()
        numbers.append(int(float(parts[0])))
        positions.append([float(v) * BOHR_TO_ANG for v in parts[1:4]])
    return ExternalRequest(
        deriv=deriv,
        charge=charge,
        multiplicity=multiplicity,
        numbers=np.array(numbers, dtype=int),
        positions=np.array(positions, dtype=float),
    )


def _write_values(fh, values, per_line=3):
    values = list(values)
    for i in range(0, len(values), per_line):
        fh.write("".join(f"{v:20.12E}" for v in values[i : i + per_line]))
        fh.write("\n")


def write_external_output(
    filename, energy, gradient=None, hessian=None, dipole=None
):
    """
    Write a Gaussian External output file.

    Args:
        filename (str): Output path Gaussian will read.
        energy (float): Energy in Hartree.
        gradient (array, optional): dE/dx in Hartree/Bohr, shape (N, 3).
        hessian (array, optional): Cartesian Hessian in Hartree/Bohr^2,
            shape (3N, 3N). Its presence triggers the second-derivative
            block (zero polarizability and dipole derivatives, then the
            lower-triangle force constants).
        dipole (array, optional): Dipole moment; zeros when omitted.
    """
    dipole = np.zeros(3) if dipole is None else np.asarray(dipole)
    with open(filename, "w") as fh:
        _write_values(fh, [energy, *dipole], per_line=4)
        if gradient is None:
            return
        gradient = np.asarray(gradient).reshape(-1, 3)
        for row in gradient:
            _write_values(fh, row)
        if hessian is None:
            return
        natoms = gradient.shape[0]
        _write_values(fh, np.zeros(6))  # polarizability
        _write_values(fh, np.zeros(9 * natoms))  # dipole derivatives
        hessian = np.asarray(hessian)
        lower = hessian[np.tril_indices(3 * natoms)]
        _write_values(fh, lower)


# ── Driver ────────────────────────────────────────────────────────────────


def run(config, input_file, output_file, engine=None):
    """Answer one Gaussian External request. Returns a status message."""
    request = read_external_input(input_file)
    engine = engine or build_engine(config)
    energy, gradient = engine.energy_and_gradient(
        request.numbers,
        request.positions,
        request.charge,
        request.multiplicity,
    )
    hessian = None
    message = f"{config['backend']}: E = {energy:.10f} Ha"
    if request.deriv >= 2:
        hessian, source = compute_hessian(
            engine,
            request,
            mode=config.get("hessian", "auto"),
            step=config.get("fd_step", DEFAULT_FD_STEP),
        )
        message += f"; Hessian from {source}"
    write_external_output(
        output_file,
        energy,
        gradient if request.deriv >= 1 else None,
        hessian,
    )
    return message


def main(config, argv=None):
    """Entry point for generated launcher scripts. Returns an exit code."""
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) < 6:
        sys.stderr.write(
            "Usage: <script> [extra_args] layer InputFile OutputFile "
            "MsgFile FChkFile MatElFile\n"
        )
        return 1
    _layer, input_file, output_file, msg_file, _fchk, _matel = argv[-6:]
    try:
        message = run(config, input_file, output_file)
    except Exception:
        with open(msg_file, "w") as fh:
            fh.write("chemsmart External calculation failed:\n")
            fh.write(traceback.format_exc())
        traceback.print_exc()
        return 1
    with open(msg_file, "w") as fh:
        fh.write(message + "\n")
    return 0
