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

import importlib
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import traceback
from dataclasses import dataclass

import numpy as np

logger = logging.getLogger(__name__)

# CODATA 2018
BOHR_TO_ANG = 0.529177210903
EV_PER_HA = 27.211386245988
# eV/Angstrom -> Hartree/Bohr
EV_ANG_TO_HA_BOHR = BOHR_TO_ANG / EV_PER_HA
# eV/Angstrom^2 -> Hartree/Bohr^2
EV_ANG2_TO_HA_BOHR2 = BOHR_TO_ANG**2 / EV_PER_HA

HESSIAN_MODES = ("auto", "analytic", "fd")
DEFAULT_FD_STEP = 0.005  # Angstrom


# ── File formats ──────────────────────────────────────────────────────────


@dataclass
class ExternalRequest:
    """One Gaussian External request, with positions in Angstrom."""

    deriv: int
    charge: int
    multiplicity: int
    numbers: np.ndarray
    positions: np.ndarray

    @property
    def natoms(self):
        return len(self.numbers)


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


# ── Engines ───────────────────────────────────────────────────────────────


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
        energy = atoms.get_potential_energy() / EV_PER_HA
        gradient = -atoms.get_forces() * EV_ANG_TO_HA_BOHR
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
        return hessian * EV_ANG2_TO_HA_BOHR2


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
            parse_xtb_gradient_file,
        )

    def analytic_hessian(self, numbers, positions, charge, multiplicity):
        try:
            return self._run(
                "--hess",
                numbers,
                positions,
                charge,
                multiplicity,
                parse_xtb_hessian_file,
            )
        except (RuntimeError, OSError, ValueError) as exc:
            logger.warning(f"xtb Hessian unavailable: {exc}")
            return None


def parse_xtb_gradient_file(workdir, natoms):
    """Return (energy [Ha], gradient [Ha/Bohr]) from xtb's ``gradient``."""
    with open(os.path.join(workdir, "gradient")) as fh:
        lines = [line for line in fh if line.strip()]
    start = next(i for i, line in enumerate(lines) if "$grad" in line)
    energy = float(lines[start + 1].split("SCF energy =")[1].split()[0])
    grad_lines = lines[start + 2 + natoms : start + 2 + 2 * natoms]
    gradient = np.array(
        [
            [float(v) for v in line.replace("D", "E").split()[:3]]
            for line in grad_lines
        ]
    )
    if gradient.shape != (natoms, 3):
        raise ValueError(f"Unexpected xtb gradient shape {gradient.shape}.")
    return energy, gradient


def parse_xtb_hessian_file(workdir, natoms):
    """Return the Cartesian Hessian [Ha/Bohr^2] from xtb's ``hessian``."""
    values = []
    with open(os.path.join(workdir, "hessian")) as fh:
        in_block = False
        for line in fh:
            if line.strip().startswith("$hessian"):
                in_block = True
                continue
            if in_block and line.strip().startswith("$"):
                break
            if in_block:
                values.extend(float(v) for v in line.replace("D", "E").split())
    n = 3 * natoms
    if len(values) != n * n:
        raise ValueError(
            f"xtb hessian has {len(values)} values; expected {n * n}."
        )
    return np.array(values).reshape(n, n)


# ── Hessians ──────────────────────────────────────────────────────────────


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
    step_bohr = step / BOHR_TO_ANG
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
    args = (
        request.numbers,
        request.positions,
        request.charge,
        request.multiplicity,
    )
    if mode == "fd":
        return finite_difference_hessian(engine, request, step), "fd"

    analytic = engine.analytic_hessian(*args)
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


# ── Backends ──────────────────────────────────────────────────────────────


def _import_object(dotted_path):
    module_name, _, attr = dotted_path.rpartition(".")
    if not module_name:
        raise ValueError(
            f"{dotted_path!r} is not a dotted path like "
            "'ase.calculators.emt.EMT'."
        )
    return getattr(importlib.import_module(module_name), attr)


def _mace_factory(function_name):
    def factory(options):
        defaults = {
            "model": "medium",
            "device": "cpu",
            "default_dtype": "float64",
        }
        kwargs = {**defaults, **options}

        def make(charge, multiplicity):
            from mace import calculators

            return getattr(calculators, function_name)(**kwargs)

        return make

    return factory


def _aimnet2_factory(options):
    options = dict(options)
    model = options.pop("model", "aimnet2")

    def make(charge, multiplicity):
        from aimnet.calculators import AIMNet2ASE

        return AIMNet2ASE(model, charge=charge, mult=multiplicity, **options)

    return make


def _uma_factory(options):
    options = dict(options)
    model = options.pop("model", "uma-s-1")
    device = options.pop("device", "cpu")
    task_name = options.pop("task_name", "omol")

    def make(charge, multiplicity):
        from fairchem.core import FAIRChemCalculator, pretrained_mlip

        predictor = pretrained_mlip.get_predict_unit(model, device=device)
        return FAIRChemCalculator(predictor, task_name=task_name, **options)

    return make


ASE_PRESETS = {
    "mace-off": _mace_factory("mace_off"),
    "mace-mp": _mace_factory("mace_mp"),
    "aimnet2": _aimnet2_factory,
    "uma": _uma_factory,
}


def build_engine(config):
    """Create the engine described by a launcher configuration dict."""
    backend = config["backend"]
    options = dict(config.get("options") or {})
    if backend == "xtb":
        return XTBEngine(**options)
    if backend == "ase":
        target = _import_object(config["calculator"])

        def make(charge, multiplicity):
            return target(**options)

        return ASEEngine(make)
    if backend in ASE_PRESETS:
        return ASEEngine(ASE_PRESETS[backend](options))
    raise ValueError(f"Unknown External backend {backend!r}.")


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
