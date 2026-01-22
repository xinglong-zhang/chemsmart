import hashlib
import logging

from chemsmart.assembler.provenance import build_provenance
from chemsmart.io.molecules.structure import Molecule

logger = logging.getLogger(__name__)


def get_record_id(filename):
    """Generate a stable record ID from filename."""
    return hashlib.sha256(str(filename).encode()).hexdigest()


class BaseAssembler:
    output_class = None  # To be defined in subclasses
    program_name = "unknown"  # To be defined in subclasses

    def __init__(self, filename, index="-1"):
        self.filename = filename
        self.index = index
        self.output = self.output_class(filename)

    @property
    def molecules_list(self):
        return Molecule.from_filepath(
            self.filename, index=self.index, return_list=True
        )

    def assemble(self):
        if not self.output.normal_termination:
            logger.warning(
                f"Calculation in {self.filename} did not terminate normally, skip assembling..."
            )
            return None
        if not self.molecules_list:
            logger.error(f"No molecules parsed from {self.filename}.")
            return None

        record_id = get_record_id(self.filename)

        data = {
            "record_id": record_id,
            "program": self.program_name,
            "meta": self.get_meta_data(),
            "results": self.get_calculation_results(),
            "molecules": [],
            "provenance": build_provenance(self.filename, self.output),
        }

        for i, mol in enumerate(self.molecules_list):
            mol_entry = {"index": i + 1, **self.get_molecule_info(mol)}
            data["molecules"].append(mol_entry)

        return data

    def get_meta_data(self):
        meta_data = {
            "filename": self.filename,
            "version": self.output.version,
            "date": self.output.date,
            "functional": self.output.functional,
            "basis_set": self.output.basis,
            "num_basis_functions": self.output.num_basis_functions,
            "spin": self.output.spin,
            "job_type": self.output.job_type,
            "solvent_on": self.output.solvent_on,
            "route_string": self.output.route_string,
        }
        if self.output.solvent_on:
            meta_data.update(
                {
                    "solvent_model": self.output.solvent_model,
                    "solvent_id": self.output.solvent_id,
                }
            )
        if self.output.freq:
            meta_data.update(
                {
                    "temperature_in_K": self.output.temperature_in_K,
                    "pressure_in_atm": self.output.pressure_in_atm,
                }
            )
        return meta_data

    def get_molecule_info(self, mol):
        molecule_info = {
            "structure_index_in_file": mol.structure_index_in_file,
            "is_optimized_structure": mol.is_optimized_structure,
            "charge": mol.charge,
            "multiplicity": mol.multiplicity,
            "coordinates": list(zip(mol.chemical_symbols, mol.positions)),
            "chemical_formula": mol.chemical_formula,
            "number_of_atoms": mol.num_atoms,
            "mass": mol.mass,
            "elements": mol.elements,
            "element_counts": mol.element_counts,
            "center_of_mass": mol.center_of_mass,
            "is_chiral": mol.is_chiral,
            "is_ring": mol.is_ring,
            "is_monoatomic": mol.is_monoatomic,
            "is_diatomic": mol.is_diatomic,
            "is_linear": mol.is_linear,
            "smiles": mol.to_smiles(),
            "moments_of_inertia": mol.moments_of_inertia,
            "frozen_atoms": mol.frozen_atoms,
            "energy": mol.energy,
            "forces": mol.forces,
        }
        if mol.mulliken_atomic_charges is not None:
            molecule_info["mulliken_atomic_charges"] = (
                mol.mulliken_atomic_charges
            )
        if mol.rotational_symmetry_number is not None:
            molecule_info["rotational_symmetry_number"] = (
                mol.rotational_symmetry_number
            )
        if mol.has_vibrations:
            molecule_info.update(
                {
                    "num_vibrational_modes": mol.num_vib_modes,
                    "vibrational_frequencies": mol.vibrational_frequencies,
                    "vibrational_modes": mol.vibrational_modes,
                }
            )
        return molecule_info

    def get_calculation_results(self):
        calculation_results = {
            "total_energy": self.output.energies[-1],
            "homo_energy": self.output.homo_energy,
            "lumo_energy": self.output.lumo_energy,
            "fmo_gap": self.output.fmo_gap,
            "total_core_hours": self.output.total_core_hours,
            "total_elapsed_walltime": self.output.total_elapsed_walltime,
        }
        if self.output.freq:
            calculation_results.update(
                {
                    "rotational_temperatures_in_K": self.output.rotational_temperatures,
                    "rotational_constants_in_Hz": self.output.rotational_constants_in_Hz,
                    "zero_point_energy": self.output.zero_point_energy,
                    "thermal_vibration_correction": self.output.thermal_vibration_correction,
                    "thermal_rotation_correction": self.output.thermal_rotation_correction,
                    "thermal_translation_correction": self.output.thermal_translation_correction,
                    "thermal_energy_correction": self.output.thermal_energy_correction,
                    "thermal_enthalpy_correction": self.output.thermal_enthalpy_correction,
                    "thermal_free_energy_correction": self.output.thermal_gibbs_free_energy_correction,
                    "internal_energy": self.output.internal_energy,
                    "enthalpy": self.output.enthalpy,
                    "electronic_entropy": self.output.electronic_entropy,
                    "vibrational_entropy": self.output.vibrational_entropy,
                    "rotational_entropy": self.output.rotational_entropy,
                    "translational_entropy": self.output.translational_entropy,
                    "entropy": self.output.entropy,
                    "entropy_times_temperature": self.output.entropy_times_temperature,
                    "gibbs_free_energy": self.output.gibbs_free_energy,
                }
            )
        return calculation_results
