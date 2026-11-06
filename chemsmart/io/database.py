import numpy as np

from chemsmart.database.database import Database
from chemsmart.database.utils import (
    resolve_record,
    sort_structure_dicts_by_energy,
)
from chemsmart.utils.mixins import FileMixin
from chemsmart.utils.utils import string2index_1based


class DatabaseFile(FileMixin):
    """CHEMSMART database file object."""

    def __init__(self, filename):
        self.filename = filename

    def build_molecule_from_database(self, struct_dict):
        """Convert a structure dictionary from the database to a Molecule object."""
        from chemsmart.io.molecules.structure import Molecule

        vibrational_modes = struct_dict.get("vibrational_modes")
        if vibrational_modes is not None:
            vibrational_modes = [np.array(mode) for mode in vibrational_modes]
        positions = struct_dict.get("positions")
        return Molecule(
            symbols=struct_dict.get("chemical_symbols"),
            positions=(np.array(positions) if positions is not None else None),
            charge=struct_dict.get("charge"),
            multiplicity=struct_dict.get("multiplicity"),
            frozen_atoms=struct_dict.get("frozen_atoms"),
            energy=struct_dict.get("energy"),
            forces=(
                np.array(struct_dict["forces"])
                if struct_dict.get("forces")
                else None
            ),
            vibrational_frequencies=struct_dict.get("vibrational_frequencies"),
            vibrational_modes=vibrational_modes,
            structure_index_in_file=struct_dict.get("structure_index_in_file"),
            rotational_symmetry_number=struct_dict.get(
                "rotational_symmetry_number"
            ),
            mulliken_atomic_charges=struct_dict.get("mulliken_atomic_charges"),
            mulliken_spin_densities=struct_dict.get("mulliken_spin_densities"),
            is_optimized_structure=struct_dict.get("is_optimized_structure"),
            dipole_moment=(
                np.array(struct_dict["dipole_moment"])
                if struct_dict.get("dipole_moment") is not None
                else None
            ),
            dipole_moment_magnitude=struct_dict.get("dipole_moment_magnitude"),
            rotational_constants=(
                np.array(struct_dict["rotational_constants"])
                if struct_dict.get("rotational_constants") is not None
                else None
            ),
            point_group=struct_dict.get("point_group"),
        )

    def get_molecule_by_structure_id(self, structure_id, return_list=False):
        """Get a molecule by its global structure ID or prefix."""
        db = Database(self.filename)
        full_sid = db.get_structure_by_partial_id(structure_id)
        struct = db.get_structure(full_sid)
        if struct is None:
            raise ValueError(f"No structure found with ID '{structure_id}'.")
        molecule = self.build_molecule_from_database(struct)
        if return_list:
            return [molecule]
        return molecule

    def get_molecules_by_molecule_id(self, molecule_id, return_list=False):
        """Get molecules by their molecule ID or prefix."""
        db = Database(self.filename)
        full_mid = db.get_molecule_by_partial_id(molecule_id)
        struct_dicts = db.get_structures_for_molecule(full_mid)
        if not struct_dicts:
            raise ValueError(
                f"No structures found for molecule '{molecule_id}'."
            )
        struct_dicts = sort_structure_dicts_by_energy(
            self.filename, struct_dicts
        )
        molecules = [
            self.build_molecule_from_database(struct)
            for struct in struct_dicts
        ]
        if return_list:
            return molecules
        return molecules if len(molecules) != 1 else molecules[0]

    def get_molecules_by_record(
        self,
        record_index=None,
        record_id=None,
        structure_index=":",
        return_list=False,
    ):
        """Get molecules by record index or ID and structure index."""
        db = Database(self.filename)
        records = resolve_record(
            db,
            record_index=record_index,
            record_id=record_id,
            return_list=True,
        )
        molecules = []
        selected_index = string2index_1based(str(structure_index))
        for record in records:
            mol_dicts = record.get("molecules", [])
            if not mol_dicts:
                continue
            if isinstance(selected_index, int):
                try:
                    molecule = self.build_molecule_from_database(
                        mol_dicts[selected_index]
                    )
                    molecules.append(molecule)
                except IndexError as exc:
                    raise ValueError(
                        f"Structure index {structure_index} out of "
                        "range for selected record."
                    ) from exc
            elif isinstance(selected_index, slice):
                molecule = [
                    self.build_molecule_from_database(mol_dict)
                    for mol_dict in mol_dicts[selected_index]
                ]
                molecules.extend(molecule)
        if return_list:
            return molecules
        return molecules if len(molecules) != 1 else molecules[0]

    def get_all_molecules(self, return_list=False):
        """Get all molecules from the chemsmart database."""
        db = Database(self.filename)
        all_structs = db.get_all_structures()
        molecules = [self.build_molecule_from_database(s) for s in all_structs]
        if return_list:
            return molecules
        return molecules if len(molecules) != 1 else molecules[0]
