"""Read molecular structures from NumPy archives, including revQM9."""

import pickle
from dataclasses import dataclass, field
from functools import cached_property
from typing import Optional
from zipfile import ZipFile

import numpy as np
from ase.symbols import Symbols
from ase.units import Hartree

from chemsmart.utils.mixins import FileMixin
from chemsmart.utils.utils import parse_index_specification


@dataclass(frozen=True)
class NPZSchema:
    """Array layout and property metadata of a supported NPZ dataset.

    ``property_units`` describes values when the dataset is explicitly declared.
    """

    name: str
    dataset: str
    position_key: str
    atom_key: str
    property_units: dict
    dipole_key: str
    orbital_energy_to_eV: float = 1.0
    zpve_to_hartree: float = 1.0
    rotational_constant_keys: tuple = ()
    rotational_constant_to_Hz: float = 1.0
    count_key: Optional[str] = None
    energy_key: Optional[str] = None
    meta: dict = field(default_factory=dict)

    @property
    def required_keys(self):
        """Array names identifying this dataset layout."""
        keys = {self.position_key, self.atom_key, self.count_key}
        return keys - {None}


QM9_SCHEMA = NPZSchema(
    name="qm9_eV",
    dataset="qm9",
    position_key="R",
    atom_key="Z",
    count_key="N",
    dipole_key="mu",
    zpve_to_hartree=1.0 / Hartree,
    rotational_constant_keys=("A", "B", "C"),
    rotational_constant_to_Hz=1.0e9,
    meta={"method": "B3LYP", "basis": "6-31G(2df,p)"},
    property_units={
        "id": None,
        "A": "GHz",
        "B": "GHz",
        "C": "GHz",
        "mu": "Debye",
        "alpha": "bohr^3",
        "homo": "eV",
        "lumo": "eV",
        "gap": "eV",
        "r2": "bohr^2",
        "zpve": "eV",
        "U0": "eV",
        "U": "eV",
        "H": "eV",
        "G": "eV",
        "Cv": "cal/(mol K)",
    },
)

REVQM9_SCHEMA = NPZSchema(
    name="revQM9",
    dataset="revqm9",
    position_key="coords",
    atom_key="charges",
    # Assume revQM9 total energies are in Hartree.
    energy_key="energies",
    dipole_key="dipole",
    orbital_energy_to_eV=Hartree,
    property_units={
        "inds": None,
        "atomization": "Hartree",
        "homo": "Hartree",
        "lumo": "Hartree",
        "dipole": "Debye",
    },
    meta={"method": "aPBE0", "basis": "cc-pVTZ"},
)

NPZ_SCHEMAS = (QM9_SCHEMA, REVQM9_SCHEMA)
NPZ_DATASETS = tuple(schema.dataset for schema in NPZ_SCHEMAS)
DATASET_PROPERTY_KEYS = frozenset(
    key for schema in NPZ_SCHEMAS for key in schema.property_units
)


@dataclass(frozen=True)
class _NPZData:
    """Validated arrays loaded from one NPZ archive."""

    positions: np.ndarray
    atoms: np.ndarray
    other_atoms: dict
    offsets: Optional[np.ndarray]
    properties: dict
    energies: Optional[np.ndarray]

    @property
    def num_molecules(self):
        """Number of structures in the archive."""
        if self.offsets is None:
            return len(self.positions)
        return len(self.offsets) - 1

    def structure_arrays(self, i):
        """Return coordinates, atoms, and other atom identities for one index."""
        if self.offsets is None:
            selection = i
        else:
            selection = slice(self.offsets[i], self.offsets[i + 1])
        return (
            np.asarray(self.positions[selection], dtype=float),
            np.asarray(self.atoms[selection]),
            {
                key: values[selection]
                for key, values in self.other_atoms.items()
            },
        )


class NPZFile(FileMixin):
    """Molecular NPZ parser; coordinates must be in Angstrom.

    Supports DGL QM9 (R/Z/N) and revQM9 (coords/charges).
    """

    def __init__(self, filename, dataset=None):
        self.filename = filename
        if dataset is not None and dataset not in NPZ_DATASETS:
            raise ValueError("Dataset must be 'revqm9' or 'qm9'.")
        self.dataset = dataset

    @cached_property
    def molecule(self):
        """Return the last structure, like XYZFile.molecule."""
        return self.get_molecules(index="-1")

    @cached_property
    def schema(self):
        """Return the supported dataset schema matching this archive."""
        for schema in NPZ_SCHEMAS:
            if schema.required_keys <= self._array_names:
                if self.dataset not in (None, schema.dataset):
                    raise ValueError(
                        f"NPZ layout does not match dataset {self.dataset}."
                    )
                return schema
        raise ValueError(
            "NPZ requires DGL QM9 fields R/Z/N or revQM9 fields coords/charges."
        )

    @cached_property
    def _array_names(self):
        with ZipFile(self.filepath) as archive:
            return {
                name[:-4]
                for name in archive.namelist()
                if name.endswith(".npy")
            }

    def get_molecules(self, index=":", return_list=False):
        """Read selected structures using the existing 1-based index syntax."""
        data = self._read_data()
        selector = parse_index_specification(str(index))
        indices = self._select_indices(selector, data.num_molecules, index)
        molecules = [self._build_molecule(data, i) for i in indices]
        if return_list or not isinstance(selector, int):
            return molecules
        return molecules[0]

    def iter_molecule_results(self, index=":"):
        """Yield one-based indices and either a molecule or its validation error.

        Archive layout and selection errors are raised before iteration.
        """
        data = self._read_data()
        selector = parse_index_specification(str(index))
        indices = self._select_indices(selector, data.num_molecules, index)
        for i in indices:
            try:
                yield i + 1, self._build_molecule(data, i), None
            except (ValueError, TypeError, IndexError) as error:
                yield i + 1, None, error

    def _select_indices(self, selector, num_molecules, index):
        """Return the zero-based structures chosen by a parsed index."""
        available = range(num_molecules)
        try:
            if isinstance(selector, slice):
                indices = list(available[selector])
            elif isinstance(selector, int):
                indices = [available[selector]]
            else:
                indices = [available[i] for i in selector]
        except IndexError:
            indices = []
        if not indices:
            raise IndexError(
                f"Index {index} selects no structures in {self.filename}, "
                f"which has {num_molecules}."
            )
        return indices

    def _read_data(self):
        """Load and validate the arrays shared by all selected structures."""
        schema = self.schema
        with ZipFile(self.filepath) as archive:
            positions = self._read_array(archive, schema.position_key)
            atoms = self._read_array(archive, schema.atom_key)
            # Check explicit identity fields when present instead of silently
            # selecting one contradictory set of elements.
            other_atoms = {
                key: self._read_array(archive, key)
                for key in ("elements", "atomic_numbers")
                if key in self._array_names
            }
            # The DimeNet/DGL QM9 archive concatenates all atoms in R and Z;
            # N contains the number of atoms for each molecule.
            offsets = None
            if schema.count_key:
                if positions.ndim != 2 or atoms.ndim != 1:
                    raise ValueError(
                        f"NPZ {schema.name} must concatenate all atoms in "
                        f"{schema.position_key} and {schema.atom_key}; "
                        f"padded arrays are not supported."
                    )
                counts = self._read_array(archive, schema.count_key)
                offsets = self._get_atom_offsets(counts, positions, atoms)
                num_molecules = len(counts)
            else:
                if positions.ndim == 2 and positions.shape[-1] == 3:
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
            return _NPZData(
                positions=positions,
                atoms=atoms,
                other_atoms=other_atoms,
                offsets=offsets,
                properties=self._read_properties(archive, num_molecules),
                energies=self._read_energies(archive, num_molecules),
            )

    @staticmethod
    def _get_atom_offsets(counts, positions, atoms):
        """Return per-molecule atom boundaries of concatenated arrays."""
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
        return np.concatenate(([0], np.cumsum(counts, dtype=np.int64)))

    def _read_properties(self, archive, num_molecules):
        """Read the schema's per-molecule properties present in the archive."""
        properties = {}
        for key in self.schema.property_units:
            if key not in self._array_names:
                continue
            values = self._read_array(archive, key)
            if values.ndim != 1 or len(values) != num_molecules:
                raise ValueError(
                    f"NPZ {key} must contain one value per molecule."
                )
            if values.dtype.kind not in "iuf" or not np.isfinite(values).all():
                raise ValueError(
                    f"NPZ {key} must contain finite numeric values."
                )
            properties[key] = values
        return properties

    def _read_energies(self, archive, num_molecules):
        """Read total energies when the schema and archive provide them."""
        energy_key = self.schema.energy_key
        if energy_key is None or energy_key not in self._array_names:
            return None
        energies = self._read_array(archive, energy_key)
        if (
            energies.shape != (num_molecules,)
            or energies.dtype.kind not in "iuf"
        ):
            raise ValueError(
                f"NPZ {energy_key} must contain one numeric energy per molecule."
            )
        return energies

    def _build_molecule(self, data, i):
        """Validate one structure and build its Molecule."""
        from chemsmart.io.molecules.structure import Molecule

        xyz, atom_values, other_atoms = data.structure_arrays(i)
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
            raise ValueError(f"NPZ structure {i + 1}: invalid atomic numbers.")
        if (numbers == 0).any():
            raise ValueError(
                f"NPZ structure {i + 1}: zero-padded atomic numbers are not supported."
            )
        self._check_atom_identities(numbers, other_atoms, i)
        if not np.isfinite(xyz).all():
            raise ValueError(
                f"NPZ structure {i + 1}: coordinates must be finite."
            )
        energy = None
        if data.energies is not None:
            energy = float(data.energies[i])
            if not np.isfinite(energy):
                raise ValueError(
                    f"NPZ structure {i + 1}: {self.schema.energy_key} must be finite."
                )
        return Molecule(
            symbols=Symbols(numbers.astype(int)),
            positions=xyz,
            energy=energy,
            charge=0,
            multiplicity=1,
            structure_index_in_file=int(i) + 1,
            dataset_source=self.schema.name,
            **{
                key: values[i].item()
                for key, values in data.properties.items()
            },
        )

    def _check_atom_identities(self, numbers, other_atoms, i):
        """Reject explicit identity fields contradicting the atom array."""
        for key, values in other_atoms.items():
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
            if other_numbers.shape != numbers.shape or not np.array_equal(
                other_numbers, numbers
            ):
                raise ValueError(
                    f"NPZ structure {i + 1}: conflicting atomic identities "
                    f"in {self.schema.atom_key} and {key}."
                )

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
