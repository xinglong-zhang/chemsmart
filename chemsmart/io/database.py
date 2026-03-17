from functools import cached_property

import numpy as np

from chemsmart.assembler.database import Database
from chemsmart.utils.mixins import FileMixin
from chemsmart.utils.utils import string2index_1based


class DatabaseFile(FileMixin):
    """Chemsmart database file object."""

    def __init__(self, filename):
        self.filename = filename

    @cached_property
    def last_structure(self):
        """Return the last molecular structure from the record(s)."""
        return self.get_molecules(index="-1")

    @property
    def molecule(self):
        """Alias for the last molecular structure."""
        return self.last_structure

    def get_molecules(
        self,
        index=":",
        return_list=False,
        record_index=None,
        record_id=None,
    ):
        from chemsmart.io.molecules.structure import Molecule

        def build_molecule_from_database(mol_dict):
            """Convert a molecule dictionary to a Molecule object."""
            vibrational_modes = mol_dict.get("vibrational_modes")
            if vibrational_modes is not None:
                vibrational_modes = [
                    np.array(mode) for mode in vibrational_modes
                ]
            return Molecule(
                symbols=mol_dict.get("chemical_symbols"),
                positions=np.array(mol_dict.get("positions")),
                charge=mol_dict.get("charge"),
                multiplicity=mol_dict.get("multiplicity"),
                energy=mol_dict.get("energy"),
                forces=(
                    np.array(mol_dict["forces"])
                    if mol_dict.get("forces")
                    else None
                ),
                frozen_atoms=mol_dict.get("frozen_atoms"),
                vibrational_frequencies=mol_dict.get(
                    "vibrational_frequencies"
                ),
                vibrational_modes=vibrational_modes,
                mulliken_atomic_charges=mol_dict.get(
                    "mulliken_atomic_charges"
                ),
                rotational_symmetry_number=mol_dict.get(
                    "rotational_symmetry_number"
                ),
                is_optimized_structure=mol_dict.get("is_optimized_structure"),
                structure_index_in_file=mol_dict.get(
                    "structure_index_in_file"
                ),
            )

        db = Database(self.filename)

        # Resolve record selection
        if record_index is not None:
            record = db.get_record(record_index=record_index)
            if record is None:
                raise ValueError(f"No record found with index {record_index}.")
            records = [record]
        elif record_id is not None:
            full_id = db.get_record_by_partial_id(record_id)
            record = db.get_record(record_id=full_id)
            if record is None:
                raise ValueError(f"No record found with ID '{record_id}'.")
            records = [record]
        else:
            records = db.get_all_records()

        molecules = []
        index = string2index_1based(index)
        for record in records:
            mol_dicts = record.get("molecules", [])
            if not mol_dicts:
                continue
            if isinstance(index, int):
                try:
                    molecule = build_molecule_from_database(mol_dicts[index])
                    molecules.append(molecule)
                except IndexError:
                    raise ValueError(
                        f"Molecule index '{index}' out of range for record {record.get('record_index')}"
                    )
            elif isinstance(index, slice):
                molecule = [
                    build_molecule_from_database(mol_dict)
                    for mol_dict in mol_dicts[index]
                ]
                molecules.extend(molecule)
        if return_list:
            return molecules
        else:
            if len(molecules) == 1:
                return molecules[0]
            return molecules
