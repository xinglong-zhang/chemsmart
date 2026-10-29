"""Read molecular structures from NumPy archives, including revQM9."""

import pickle
from functools import cached_property
from zipfile import ZipFile

import numpy as np
from ase.symbols import Symbols

from chemsmart.utils.mixins import FileMixin
from chemsmart.utils.utils import parse_index_specification


class NPZFile(FileMixin):
    """Molecular NPZ parser; coordinates must be in Angstrom.

    Supports DGL QM9 (R/Z/N) and revQM9 (coords/charges).
    """

    def __init__(self, filename):
        self.filename = filename

    @cached_property
    def molecule(self):
        """Return the last structure, like XYZFile.molecule."""
        return self.get_molecules(index="-1")

    def get_molecules(self, index=":", return_list=False):
        """Read selected structures using the existing 1-based index syntax."""
        from chemsmart.io.molecules.structure import Molecule

        with ZipFile(self.filepath) as archive:
            keys = {
                name[:-4]
                for name in archive.namelist()
                if name.endswith(".npy")
            }
            if {"R", "Z", "N"} <= keys:
                position_key, atom_key = "R", "Z"
                energy_key = None
                property_keys = (
                    "id",
                    "A",
                    "B",
                    "C",
                    "mu",
                    "alpha",
                    "homo",
                    "lumo",
                    "gap",
                    "r2",
                    "zpve",
                    "U0",
                    "U",
                    "H",
                    "G",
                    "Cv",
                )
            elif {"coords", "charges"} <= keys:
                position_key, atom_key = "coords", "charges"
                # Assume revQM9 total energies are in Hartree.
                energy_key = "energies"
                property_keys = (
                    "atomization",
                    "homo",
                    "lumo",
                    "dipole",
                    "inds",
                )
            else:
                raise ValueError(
                    "NPZ requires DGL QM9 fields R/Z/N or revQM9 fields coords/charges."
                )
            positions = self._read_array(archive, position_key)
            atoms = self._read_array(archive, atom_key)
            # Check explicit identity fields when present instead of silently
            # selecting one contradictory set of elements.
            other_atoms = {
                key: self._read_array(archive, key)
                for key in ("elements", "atomic_numbers")
                if key in keys
            }
            # The DimeNet/DGL QM9 archive concatenates all atoms in R and Z;
            # N contains the number of atoms for each molecule.
            offsets = None
            if positions.ndim == 2 and atoms.ndim == 1 and "N" in keys:
                counts = self._read_array(archive, "N")
                if (
                    counts.ndim != 1
                    or counts.dtype.kind not in "iu"
                    or np.any(counts <= 0)
                ):
                    raise ValueError(
                        "NPZ N must contain positive integer atom counts."
                    )
                if positions.shape[1:] != (3,) or len(positions) != len(atoms):
                    raise ValueError(
                        "NPZ concatenated coordinates and atoms must have matching atom counts."
                    )
                if sum(map(int, counts)) != len(atoms):
                    raise ValueError(
                        "NPZ N atom counts must sum to the number of atoms."
                    )
                offsets = np.concatenate(
                    ([0], np.cumsum(counts, dtype=np.int64))
                )
                num_molecules = len(counts)
            else:
                single = positions.ndim == 2 and positions.shape[-1] == 3
                if single:
                    positions = positions[None, ...]
                    atoms = atoms[None, ...]
                    other_atoms = {
                        key: value[None, ...]
                        for key, value in other_atoms.items()
                    }
                if positions.ndim not in (1, 3) or len(positions) != len(
                    atoms
                ):
                    raise ValueError(
                        "NPZ coordinate and atom arrays must have matching molecule counts."
                    )
                num_molecules = len(positions)
            if not num_molecules:
                raise ValueError("NPZ contains no molecular structures.")
            properties = {}
            for key in property_keys:
                if key in keys:
                    values = self._read_array(archive, key)
                    if values.ndim != 1 or len(values) != num_molecules:
                        raise ValueError(
                            f"NPZ {key} must contain one value per molecule."
                        )
                    if (
                        values.dtype.kind not in "iuf"
                        or not np.isfinite(values).all()
                    ):
                        raise ValueError(
                            f"NPZ {key} must contain finite numeric values."
                        )
                    properties[key] = values
            energies = None
            if energy_key is not None and energy_key in keys:
                energies = self._read_array(archive, energy_key)
                if (
                    energies.shape != (num_molecules,)
                    or energies.dtype.kind not in "iuf"
                ):
                    raise ValueError(
                        f"NPZ {energy_key} must contain one numeric energy per molecule."
                    )
            selector = parse_index_specification(str(index))
            if isinstance(selector, slice):
                indices = list(range(num_molecules))[selector]
            elif isinstance(selector, int):
                indices = [range(num_molecules)[selector]]
            else:
                indices = [range(num_molecules)[i] for i in selector]
            molecules = []
            for i in indices:
                if offsets is None:
                    xyz = np.asarray(positions[i], dtype=float)
                    atom_values = np.asarray(atoms[i])
                else:
                    atom_slice = slice(offsets[i], offsets[i + 1])
                    xyz = np.asarray(positions[atom_slice], dtype=float)
                    atom_values = np.asarray(atoms[atom_slice])
                if xyz.ndim != 2 or xyz.shape[1] != 3:
                    raise ValueError(
                        f"NPZ structure {i + 1}: coordinates must have shape (N, 3)."
                    )
                if atom_values.ndim != 1 or len(atom_values) != len(xyz):
                    raise ValueError(
                        f"NPZ structure {i + 1}: atom and coordinate counts differ."
                    )
                numbers = np.asarray(atom_values, dtype=float)
                if (
                    not np.isfinite(numbers).all()
                    or (numbers != np.floor(numbers)).any()
                    or (numbers < 0).any()
                    or (numbers > 118).any()
                ):
                    raise ValueError(
                        f"NPZ structure {i + 1}: invalid atomic numbers."
                    )
                if (numbers == 0).any():
                    raise ValueError(
                        f"NPZ structure {i + 1}: zero-padded atomic numbers are not supported."
                    )
                symbols = Symbols(numbers.astype(int))
                for key, other in other_atoms.items():
                    values = other[i] if offsets is None else other[atom_slice]
                    if key == "elements":
                        try:
                            other_numbers = np.asarray(
                                Symbols.fromsymbols(
                                    np.asarray(values).tolist()
                                ).numbers
                            )
                        except (TypeError, ValueError, KeyError) as exc:
                            raise ValueError(
                                f"NPZ structure {i + 1}: invalid elements."
                            ) from exc
                    else:
                        other_numbers = np.asarray(values)
                    if (
                        other_numbers.shape != numbers.shape
                        or not np.array_equal(other_numbers, numbers)
                    ):
                        raise ValueError(
                            f"NPZ structure {i + 1}: conflicting atomic identities in {atom_key} and {key}."
                        )
                if not np.isfinite(xyz).all():
                    raise ValueError(
                        f"NPZ structure {i + 1}: coordinates must be finite."
                    )
                energy = None
                if energies is not None:
                    energy = float(energies[i])
                    if not np.isfinite(energy):
                        raise ValueError(
                            f"NPZ structure {i + 1}: {energy_key} must be finite."
                        )
                molecule = Molecule(
                    symbols=symbols,
                    positions=xyz,
                    energy=energy,
                    structure_index_in_file=int(i) + 1,
                )
                molecule.dataset_source = (
                    "qm9_eV" if energy_key is None else "revQM9"
                )
                if "id" in properties:
                    molecule.id = properties["id"][i].item()
                if "A" in properties:
                    molecule.A = properties["A"][i].item()
                if "B" in properties:
                    molecule.B = properties["B"][i].item()
                if "C" in properties:
                    molecule.C = properties["C"][i].item()
                if "mu" in properties:
                    molecule.mu = properties["mu"][i].item()
                if "alpha" in properties:
                    molecule.alpha = properties["alpha"][i].item()
                if "homo" in properties:
                    molecule.homo = properties["homo"][i].item()
                if "lumo" in properties:
                    molecule.lumo = properties["lumo"][i].item()
                if "gap" in properties:
                    molecule.gap = properties["gap"][i].item()
                if "r2" in properties:
                    molecule.r2 = properties["r2"][i].item()
                if "zpve" in properties:
                    molecule.zpve = properties["zpve"][i].item()
                if "U0" in properties:
                    molecule.U0 = properties["U0"][i].item()
                if "U" in properties:
                    molecule.U = properties["U"][i].item()
                if "H" in properties:
                    molecule.H = properties["H"][i].item()
                if "G" in properties:
                    molecule.G = properties["G"][i].item()
                if "Cv" in properties:
                    molecule.Cv = properties["Cv"][i].item()
                if "atomization" in properties:
                    molecule.atomization = properties["atomization"][i].item()
                if "dipole" in properties:
                    molecule.dipole = properties["dipole"][i].item()
                if "inds" in properties:
                    molecule.inds = properties["inds"][i].item()
                molecules.append(molecule)
        if return_list or not isinstance(selector, int):
            return molecules
        return molecules[0]

    @staticmethod
    def _read_array(archive, key):
        """Read NPY data, restricting object arrays to NumPy reconstruction."""
        with archive.open(f"{key}.npy") as stream:
            version = np.lib.format.read_magic(stream)
            readers = {
                (1, 0): np.lib.format.read_array_header_1_0,
                (2, 0): np.lib.format.read_array_header_2_0,
            }
            if version not in readers:
                raise ValueError(f"Unsupported NPY version: {version}")
            shape, _, dtype = readers[version](stream)
            if not dtype.hasobject:
                stream.seek(0)
                return np.lib.format.read_array(stream, allow_pickle=False)
            try:
                array = _NumpyArrayUnpickler(stream).load()
            except pickle.UnpicklingError as exc:
                raise ValueError(
                    f"Unsupported object in NPZ array {key}: {exc}"
                ) from exc
            if (
                not isinstance(array, np.ndarray)
                or array.shape != shape
                or array.dtype != dtype
            ):
                raise ValueError(
                    f"NPZ array {key} does not match its NPY header."
                )
            return array


class _NumpyArrayUnpickler(pickle.Unpickler):
    """Reject arbitrary globals while restoring NumPy's ragged arrays."""

    def find_class(self, module, name):
        if module == "numpy" and name in ("ndarray", "dtype"):
            return getattr(np, name)
        if module in (
            "numpy.core.multiarray",
            "numpy._core.multiarray",
        ) and name in (
            "_reconstruct",
            "scalar",
        ):
            return getattr(np.core.multiarray, name)
        raise pickle.UnpicklingError(
            f"Disallowed pickle global: {module}.{name}"
        )
