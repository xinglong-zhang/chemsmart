import logging
from functools import cached_property

from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.orca.output import ORCAOutput

logger = logging.getLogger(__name__)


class BaseAssembler:
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
        if not self.molecules_list:
            logger.warning(f"No molecules parsed from {self.filename}.")
            return []

        assemble_data = []

        if self.index == "-1" and self.output.normal_termination:
            mol = self.molecules_list[0]
            entry = {
                **self.get_meta_data(),
                **self.get_molecule_info(mol),
                **self.get_calculation_results(),
            }
            assemble_data.append(entry)
        elif self.index == ":":
            for i, mol in enumerate(self.molecules_list):
                entry = {
                    "structure_index": i,
                    **self.get_meta_data(),
                    **self.get_molecule_info(mol),
                }
                assemble_data.append(entry)
        else:
            mol = self.molecules_list[0]
            entry = {
                "structure_index": self.index,
                **self.get_meta_data(),
                **self.get_molecule_info(mol),
            }
            assemble_data.append(entry)
        return assemble_data

    def get_meta_data(self):
        meta_data = {
            "filename": self.filename,
            "program": self.program_name,
            "version": self.output.version,
            "date": self.output.date,
            "functional": self.output.functional,
            "basis_set": self.output.basis,
            #           'spin': self.output.spin,
            "job_type": self.output.job_type,
            "solvent_on": self.output.solvent_on,
            "solvent_model": self.output.solvent_model,
            "solvent_id": self.output.solvent_id,
        }
        return meta_data

    def get_molecule_info(self, mol):
        molecule_info = {
            "charge": mol.charge,
            "multiplicity": mol.multiplicity,
            "chemical_symbols": mol.chemical_symbols,
            "positions": mol.positions,
            "chemical_formula": mol.chemical_formula,
            "number_of_atoms": mol.num_atoms,
            "mass": mol.mass,
            "elements": mol.elements,
            "element_counts": mol.element_counts,
            "center_of_mass": mol.center_of_mass,
            "is_chiral": mol.is_chiral,
            #           'is_aromatic': mol.is_aromatic,
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
        return molecule_info

    def get_calculation_results(self):
        calculation_results = {
            "zero_point_energy": self.output.zero_point_energy,
            "homo_energy": self.output.homo_energy,
            "lumo_energy": self.output.lumo_energy,
            "fmo_energy": self.output.fmo_gap,
            "mulliken_atomic_charges": self.output.mulliken_atomic_charges,
            "total_core_hours": self.output.total_core_hours,
            "total_elapsed_walltime": self.output.total_elapsed_walltime,
        }
        if self.output.freq:
            calculation_results["vibrational_frequencies"] = (
                self.output.vibrational_frequencies
            )
            calculation_results["num_vib_frequencies"] = (
                self.output.num_vib_frequencies
            )
            calculation_results["rotational_symmetry_number"] = (
                self.output.rotational_symmetry_number
            )
        return calculation_results


class GaussianAssembler(BaseAssembler):
    output_class = Gaussian16Output
    program_name = "Gaussian"


class ORCAAssembler(BaseAssembler):
    output_class = ORCAOutput
    program_name = "ORCA"


class SingleFileAssembler:
    def __init__(self, filename, index="-1", database_file="database"):
        self.filename = filename
        self.index = index
        self.database_file = database_file

    @cached_property
    def assemble_data(self):
        assembler = self._get_assembler(self.filename)
        try:
            data = assembler.assemble()
        except Exception as e:
            logger.error(f"Error assembling {self.filename}: {e}")
            return []
        return data

    def _get_assembler(self, file):
        if str(file).endswith(".log"):
            assembler = GaussianAssembler(file, index=self.index)
        elif str(file).endswith(".out"):
            assembler = ORCAAssembler(file, index=self.index)
        else:
            raise ValueError(
                "Unsupported file format. Use .log or .out files."
            )
        return assembler

    def query(self, key, value):
        results = []
        for entry in self.assemble_data:
            if key in entry and entry[key] == value:
                results.append(entry)
        return results
