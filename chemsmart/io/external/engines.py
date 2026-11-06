"""
Engines: compute energies, gradients and (optionally) Hessians for an
external level of theory.

All public methods take positions in Angstrom and return atomic units:
Hartree, Hartree/Bohr and Hartree/Bohr^2. Engines know nothing about the
host program (Gaussian, ORCA, ...) that asked for the numbers.
"""

import logging
import os
import shutil
import subprocess
import tempfile

import numpy as np

from chemsmart.utils.constants import (
    ev_per_angstrom2_to_hartree_per_bohr2,
    ev_per_angstrom_to_hartree_per_bohr,
    hartree_to_ev,
)

logger = logging.getLogger(__name__)


class ExternalEngine:
    """Computes energies, gradients and (optionally) Hessians.

    All public methods take positions in Angstrom and return atomic
    units: Hartree, Hartree/Bohr and Hartree/Bohr^2.
    """

    def energy_and_gradient(self, numbers, positions, charge, multiplicity):
        raise NotImplementedError

    def analytic_hessian(self, numbers, positions, charge, multiplicity):
        """Return the analytic Hessian, or None if unavailable."""
        return None


class ASEEngine(ExternalEngine):
    """Engine backed by an ASE calculator.

    Args:
        factory (callable): ``factory(charge, multiplicity)`` returning an
            ASE calculator. Called once per (charge, multiplicity).

    The charge and multiplicity are also stored in ``atoms.info["charge"]``
    and ``atoms.info["spin"]``, the convention used by charge/spin-aware
    models such as FAIRChem UMA.
    """

    def __init__(self, factory):
        self.factory = factory
        self._calculators = {}

    def _atoms(self, numbers, positions, charge, multiplicity):
        from ase import Atoms

        key = (charge, multiplicity)
        if key not in self._calculators:
            self._calculators[key] = self.factory(charge, multiplicity)
        atoms = Atoms(numbers=numbers, positions=positions)
        atoms.info["charge"] = charge
        atoms.info["spin"] = multiplicity
        atoms.calc = self._calculators[key]
        return atoms

    def energy_and_gradient(self, numbers, positions, charge, multiplicity):
        atoms = self._atoms(numbers, positions, charge, multiplicity)
        energy = atoms.get_potential_energy() / hartree_to_ev
        gradient = -atoms.get_forces() * ev_per_angstrom_to_hartree_per_bohr
        return energy, gradient

    def analytic_hessian(self, numbers, positions, charge, multiplicity):
        atoms = self._atoms(numbers, positions, charge, multiplicity)
        get_hessian = getattr(atoms.calc, "get_hessian", None)
        if get_hessian is None:
            return None
        n = 3 * len(numbers)
        try:
            hessian = np.asarray(get_hessian(atoms=atoms), dtype=float)
            hessian = hessian.reshape(n, n)
        except Exception as exc:  # backend-specific failures
            logger.warning(f"Analytic Hessian unavailable: {exc}")
            return None
        return hessian * ev_per_angstrom2_to_hartree_per_bohr2


class XTBEngine(ExternalEngine):
    """Engine that runs the ``xtb`` executable.

    Args:
        executable (str): Path to (or name of) the xtb binary.
        method (str): One of ``gfn0``, ``gfn1``, ``gfn2``, ``gfnff``.
        alpb (str, optional): ALPB implicit solvent name.
        gbsa (str, optional): GBSA implicit solvent name.
        args (str or list, optional): Extra xtb command-line arguments.
    """

    METHOD_ARGS = {
        "gfn0": ["--gfn", "0"],
        "gfn1": ["--gfn", "1"],
        "gfn2": ["--gfn", "2"],
        "gfnff": ["--gfnff"],
    }

    def __init__(
        self, executable="xtb", method="gfn2", alpb=None, gbsa=None, args=None
    ):
        method = str(method).lower().replace("-xtb", "").replace("_", "")
        if method not in self.METHOD_ARGS:
            raise ValueError(
                f"Unknown xtb method {method!r}; choose from "
                f"{sorted(self.METHOD_ARGS)}."
            )
        if alpb and gbsa:
            raise ValueError("Specify only one of alpb and gbsa.")
        self.executable = executable
        self.method = method
        self.alpb = alpb
        self.gbsa = gbsa
        if isinstance(args, str):
            args = args.split()
        self.args = list(args or [])

    def _command(self, runtype, charge, multiplicity):
        command = [
            self.executable,
            "coord.xyz",
            runtype,
            "--chrg",
            str(charge),
            "--uhf",
            str(multiplicity - 1),
            *self.METHOD_ARGS[self.method],
        ]
        if self.alpb:
            command += ["--alpb", str(self.alpb)]
        if self.gbsa:
            command += ["--gbsa", str(self.gbsa)]
        return command + self.args

    def _run(self, runtype, numbers, positions, charge, multiplicity, parse):
        from ase.data import chemical_symbols

        workdir = tempfile.mkdtemp(prefix="xtb_", dir=os.getcwd())
        try:
            with open(os.path.join(workdir, "coord.xyz"), "w") as fh:
                fh.write(f"{len(numbers)}\n\n")
                for z, (x, y, zc) in zip(numbers, positions):
                    fh.write(
                        f"{chemical_symbols[z]:2s} "
                        f"{x:18.10f} {y:18.10f} {zc:18.10f}\n"
                    )
            command = self._command(runtype, charge, multiplicity)
            with open(os.path.join(workdir, "xtb.out"), "w") as out:
                result = subprocess.run(
                    command,
                    cwd=workdir,
                    stdout=out,
                    stderr=subprocess.STDOUT,
                )
            if result.returncode != 0:
                with open(os.path.join(workdir, "xtb.out")) as fh:
                    tail = fh.read()[-2000:]
                raise RuntimeError(
                    f"xtb failed (exit {result.returncode}): "
                    f"{' '.join(command)}\n{tail}"
                )
            return parse(workdir, len(numbers))
        finally:
            shutil.rmtree(workdir, ignore_errors=True)

    def energy_and_gradient(self, numbers, positions, charge, multiplicity):
        return self._run(
            "--grad",
            numbers,
            positions,
            charge,
            multiplicity,
            read_xtb_gradient,
        )

    def analytic_hessian(self, numbers, positions, charge, multiplicity):
        try:
            return self._run(
                "--hess",
                numbers,
                positions,
                charge,
                multiplicity,
                read_xtb_hessian,
            )
        except (RuntimeError, OSError, ValueError) as exc:
            logger.warning(f"xtb Hessian unavailable: {exc}")
            return None


def read_xtb_gradient(workdir, natoms):
    """Return (energy [Ha], gradient [Ha/Bohr]) from xtb's ``gradient``
    file in ``workdir``, checked against ``natoms``."""
    from chemsmart.io.xtb.file import XTBGradientFile

    gradient_file = XTBGradientFile(os.path.join(workdir, "gradient"))
    gradients = gradient_file.gradients
    gradient = None if gradients is None else gradients[-1]
    if gradient_file.energy is None or gradient is None:
        raise ValueError(f"No energy/gradient in xtb output in {workdir}.")
    if gradient.shape != (natoms, 3):
        raise ValueError(f"Unexpected xtb gradient shape {gradient.shape}.")
    return gradient_file.energy, gradient


def read_xtb_hessian(workdir, natoms):
    """Return the Cartesian Hessian [Ha/Bohr^2] from xtb's ``hessian``
    file in ``workdir``, checked against ``natoms``."""
    from chemsmart.io.xtb.file import XTBHessianFile

    hessian = XTBHessianFile(os.path.join(workdir, "hessian")).hessian
    n = 3 * natoms
    if hessian is None or hessian.shape != (n, n):
        shape = None if hessian is None else hessian.shape
        raise ValueError(f"xtb hessian has shape {shape}; expected {(n, n)}.")
    return hessian
