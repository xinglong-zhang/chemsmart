"""
Runtime side of ORCA's ``ExtOpt`` interface.

With ``! ExtOpt`` in the route and ``%method ProgExt "<launcher>" end``,
ORCA keeps its own optimiser (Opt, OptTS, ScanTS, IRC, constrained
optimisations and scans) and asks an external program for every energy
and gradient. Numerical Hessians (TS searches, ``NumFreq``, IRC initial
Hessians) are assembled by ORCA itself from those gradients.

For each request ORCA writes ``<base>.xyz`` (Angstrom) and an input file
``<base>.extinp.tmp``, then runs::

    <ProgExt> <base>.extinp.tmp [Ext_Params]

Input file, one value per line (anything after the value is a comment)::

    <base>.xyz          # geometry file, Angstrom
    0                   # charge
    1                   # multiplicity
    4                   # number of cores
    1                   # 1: gradient wanted, 0: energy only
    <base>.pc           # point charges (optional; not supported here)

The launcher must write ``<base>.engrad`` (atomic units)::

    #
    # Number of atoms
    #
    natoms
    #
    # The current total energy in Eh
    #
    energy
    #
    # The current gradient in Eh/bohr
    #
    dE/dx1
    dE/dy1
    ...                 # 3 * natoms values, one per line

The launcher configuration is baked into the generated script (see
:class:`chemsmart.jobs.orca.external.ORCAExternalCalculator`), so
``Ext_Params`` is not used.
"""

import os
from dataclasses import dataclass
from typing import Optional

import numpy as np

from chemsmart.io.external import ExternalProtocol, ExternalRequest

# ── File formats ──────────────────────────────────────────────────────────


@dataclass
class ORCAExternalInput:
    """Contents of ORCA's ``<base>.extinp.tmp`` file."""

    xyz_file: str
    charge: int
    multiplicity: int
    ncores: int
    gradient: bool
    pointcharge_file: Optional[str] = None

    @property
    def basename(self):
        """``<base>`` shared by the ``.xyz`` and ``.engrad`` files."""
        root, ext = os.path.splitext(self.xyz_file)
        return root if ext.lower() == ".xyz" else self.xyz_file

    @property
    def engrad_file(self):
        return f"{self.basename}.engrad"


def _first_token(line):
    tokens = line.split()
    return tokens[0] if tokens else ""


def read_orca_external_input(filename):
    """Parse ORCA's ``.extinp.tmp`` file.

    Relative file names inside it are resolved against the folder of
    ``filename``, which is ORCA's working directory.
    """
    with open(filename) as fh:
        values = [_first_token(line) for line in fh]
    if len(values) < 5 or not all(values[:5]):
        raise ValueError(
            f"{filename} is not an ORCA ExtOpt input file: expected xyz "
            "file, charge, multiplicity, cores and gradient flag."
        )
    folder = os.path.dirname(os.path.abspath(filename))

    def resolve(name):
        return name if os.path.isabs(name) else os.path.join(folder, name)

    pointcharges = values[5] if len(values) > 5 and values[5] else None
    return ORCAExternalInput(
        xyz_file=resolve(values[0]),
        charge=int(values[1]),
        multiplicity=int(values[2]),
        ncores=int(values[3]),
        gradient=bool(int(values[4])),
        pointcharge_file=resolve(pointcharges) if pointcharges else None,
    )


def read_xyz_geometry(filename):
    """Return (atomic numbers, positions in Angstrom) from an XYZ file."""
    from ase.data import atomic_numbers

    with open(filename) as fh:
        lines = fh.readlines()
    natoms = int(lines[0].split()[0])
    numbers = []
    positions = []
    for line in lines[2 : 2 + natoms]:
        symbol, x, y, z = line.split()[:4]
        # ORCA may write labelled symbols such as "C1"; keep the element
        element = "".join(c for c in symbol if c.isalpha()).capitalize()
        numbers.append(atomic_numbers[element])
        positions.append([float(x), float(y), float(z)])
    if len(numbers) != natoms:
        raise ValueError(
            f"{filename} lists {natoms} atoms but has {len(numbers)} "
            "coordinate lines."
        )
    return np.array(numbers, dtype=int), np.array(positions, dtype=float)


def write_engrad(filename, natoms, energy, gradient=None):
    """Write an ORCA ``.engrad`` file (Hartree, Hartree/Bohr).

    ``gradient`` has shape (natoms, 3); when omitted only the energy
    section is written (ORCA asked for the energy only).
    """
    if gradient is not None:
        gradient = np.asarray(gradient, dtype=float)
        if gradient.size != 3 * natoms:
            raise ValueError(
                f"Gradient has {gradient.size} values; expected "
                f"{3 * natoms}."
            )
    with open(filename, "w") as fh:
        fh.write("#\n# Number of atoms\n#\n")
        fh.write(f"{natoms:d}\n")
        fh.write("#\n# The current total energy in Eh\n#\n")
        fh.write(f"{energy:22.12f}\n")
        if gradient is None:
            return
        fh.write("#\n# The current gradient in Eh/bohr\n#\n")
        for value in gradient.reshape(-1):
            fh.write(f"{value:22.12f}\n")


# ── Protocol ──────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class ORCAExternalInvocation:
    """ORCA's ``ProgExt`` arguments: the input file, then ``Ext_Params``."""

    input_file: str
    ext_params: tuple = ()


class ORCAExternalProtocol(ExternalProtocol):
    """ORCA ``ExtOpt`` file exchange.

    Errors go to stderr, which the ORCA job runner keeps in the ``.err``
    file; the non-zero exit code makes ORCA stop.
    """

    host = "ORCA"
    usage = "Usage: <script> <base>.extinp.tmp [Ext_Params]"

    def parse_arguments(self, argv):
        if not argv:
            raise ValueError(self.usage)
        return ORCAExternalInvocation(argv[0], tuple(argv[1:]))

    def read_request(self, invocation):
        request_input = read_orca_external_input(invocation.input_file)
        if request_input.pointcharge_file and os.path.getsize(
            request_input.pointcharge_file
        ):
            raise NotImplementedError(
                "Point charges (ORCA QM/MM embedding) are not supported by "
                "the chemsmart external backends."
            )
        numbers, positions = read_xyz_geometry(request_input.xyz_file)
        return ExternalRequest(
            deriv=int(request_input.gradient),
            charge=request_input.charge,
            multiplicity=request_input.multiplicity,
            numbers=numbers,
            positions=positions,
        )

    def write_result(self, invocation, request, result):
        request_input = read_orca_external_input(invocation.input_file)
        write_engrad(
            request_input.engrad_file,
            request.natoms,
            result.energy,
            result.gradient,
        )


def main(config, argv=None):
    """Entry point for generated launcher scripts. Returns an exit code."""
    return ORCAExternalProtocol().main(config, argv)
