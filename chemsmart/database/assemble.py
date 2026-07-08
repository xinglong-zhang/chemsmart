"""
Assembler module for parsing quantum chemistry output files and
assembling structured records for database storage.

This module consolidates the assembler pipeline:
- BaseAssembler: core assembly logic shared across programs
- GaussianAssembler: Gaussian-specific assembly
- ORCAAssembler: ORCA-specific assembly
- XTBAssembler: xTB folder-based assembly
- SingleFileAssembler: dispatcher for file-based program outputs
- SingleFolderAssembler: dispatcher for folder-based program outputs
- build_provenance: provenance metadata builder
"""

import logging
from functools import cached_property

from chemsmart import __version__ as chemsmart_version
from chemsmart.database.records import AssembledRecord
from chemsmart.database.utils import (
    canonical_json_hash,
    canonicalize_route_string,
    compute_trajectory_id,
    directory_size,
    file_size,
    get_record_id,
    is_custom_basis,
    is_custom_solvent,
    sha256_content,
    standardize_basis_set,
    utcnow_iso,
)
from chemsmart.io.folder import BaseFolder
from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.orca.output import ORCAOutput
from chemsmart.io.xtb.output import XTBOutput
from chemsmart.utils.io import get_program_type_from_file

logger = logging.getLogger(__name__)


class BaseAssembler:
    OUTPUT_CLASS = None  # To be defined in subclasses
    PROGRAM = "unknown"  # To be defined in subclasses
    FOLDER_BASED = False

    def __init__(
        self, filename=None, folder=None, index=":", include_failed=False
    ):
        self.filename = filename
        self.folder = folder
        self.index = index
        self.include_failed = include_failed
        self.target = self.folder if self.FOLDER_BASED else self.filename
        self.output = self.OUTPUT_CLASS(self.target)

    @property
    def molecules_list(self):
        if self.FOLDER_BASED:
            molecules = Molecule.from_directorypath(
                self.folder,
                program=self.PROGRAM.lower(),
                index=self.index,
            )
            return molecules if isinstance(molecules, list) else [molecules]
        return Molecule.from_filepath(
            self.filename, index=self.index, return_list=True
        )

    def build_provenance(self):
        """Build provenance metadata for an assembled record."""
        target = self.target
        output = self.output

        if self.FOLDER_BASED:
            hash_source = output.main_out
            source_size = directory_size(target)
        else:
            hash_source = output
            source_size = file_size(target)

        return {
            "source": target,
            "source_file_hash": (
                sha256_content(hash_source) if hash_source else None
            ),
            "source_size": source_size,
            "source_date": output.file_date,
            "program": self.PROGRAM,
            "program_version": output.version,
            "parser": output.__class__.__name__,
            "chemsmart_version": chemsmart_version,
            "assembled_at": utcnow_iso(),
            "normal_termination": output.normal_termination,
        }

    def assemble(self):
        if not self.output.normal_termination:
            if not self.include_failed:
                logger.warning(
                    f"Calculation in {self.target} did not terminate normally, "
                    f"skip assembling..."
                )
                return None
            else:
                logger.warning(
                    f"Calculation in {self.target} did not terminate normally, "
                    f"assembling partial data..."
                )
        if not self.molecules_list:
            logger.error(f"No molecules parsed from {self.target}.")
            return None

        meta = self.get_meta_data()
        results = self.get_calculation_results()
        provenance = self.build_provenance()
        molecules = []
        for i, mol in enumerate(self.molecules_list):
            mol_entry = {"index": i + 1, **self.get_molecule_info(mol)}
            molecules.append(mol_entry)
        trajectory_id = compute_trajectory_id(
            [m.structure_id for m in self.molecules_list]
        )
        meta = {**meta, "trajectory_id": trajectory_id}
        custom_basis_hash = canonical_json_hash(meta.get("custom_basis"))
        custom_solvent_hash = canonical_json_hash(meta.get("custom_solvent"))
        route_tokens = canonicalize_route_string(
            self.output.route_string,
            drop_terms=[
                self.output.method,
                self.output.basis,
                self.output.jobtype,
                self.output.solvent_model,
                self.output.solvent_id,
            ],
        )
        route_hash = canonical_json_hash(route_tokens)
        record_id = get_record_id(
            program=provenance.get("program") or None,
            method=meta.get("method") or None,
            basis=meta.get("basis") or None,
            jobtype=meta.get("jobtype") or None,
            trajectory_id=trajectory_id,
            custom_basis_hash=custom_basis_hash,
            solvent_model=meta.get("solvent_model"),
            solvent_id=meta.get("solvent_id"),
            custom_solvent_hash=custom_solvent_hash,
            route_hash=route_hash,
        )

        return AssembledRecord(
            record_id=record_id,
            meta=meta,
            results=results,
            molecules=molecules,
            provenance=provenance,
        )

    def get_meta_data(self):
        basis = self.output.basis
        if is_custom_basis(basis):
            basis = "customized_basis"
        else:
            basis = standardize_basis_set(basis)
        meta_data = {
            "method": self.output.method,
            "basis": basis,
            "num_basis_functions": self.output.num_basis_functions,
            "spin": self.output.spin,
            "jobtype": self.output.jobtype,
            "solvent_on": self.output.solvent_on,
            "route_string": self.output.route_string,
        }
        if self.output.solvent_on:
            solvent_id = self.output.solvent_id
            if is_custom_solvent(solvent_id):
                solvent_id = "customized_solvent"
            meta_data.update(
                {
                    "solvent_model": self.output.solvent_model,
                    "solvent_id": solvent_id,
                    "custom_solvent": self.output.custom_solvent,
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
            "molecule_label": mol.molecule_label,
            "molecule_id": mol.molecule_id,
            "structure_index_in_file": mol.structure_index_in_file,
            "is_optimized_structure": mol.is_optimized_structure,
            "charge": mol.charge,
            "multiplicity": mol.multiplicity,
            "structure_id": mol.structure_id,
            "structure_label": mol.structure_label,
            "chemical_symbols": mol.chemical_symbols,
            "positions": mol.positions,
            "empirical_formula": mol.empirical_formula,
            "chemical_formula": mol.chemical_formula,
            "number_of_atoms": mol.num_atoms,
            "mass": mol.mass,
            "elements": mol.elements,
            "element_counts": mol.element_counts,
            "center_of_mass": mol.center_of_mass,
            "is_chiral": mol.is_chiral,
            "is_ring": mol.is_ring,
            "is_aromatic": mol.is_aromatic,
            "is_monoatomic": mol.is_monoatomic,
            "is_diatomic": mol.is_diatomic,
            "is_linear": mol.is_linear,
            "is_multicomponent": mol.is_multicomponent,
            "num_components": mol.num_components,
            "smiles": mol.smiles,
            "inchi": mol.inchi,
            "chiral_centers": mol.chiral_centers,
            "moments_of_inertia": mol.moments_of_inertia,
            "frozen_atoms": mol.frozen_atoms,
            "energy": mol.energy,
            "forces": mol.forces,
        }
        if mol.mulliken_atomic_charges is not None:
            molecule_info["mulliken_atomic_charges"] = (
                mol.mulliken_atomic_charges
            )
        if mol.mulliken_spin_densities is not None:
            molecule_info["mulliken_spin_densities"] = (
                mol.mulliken_spin_densities
            )
        if mol.rotational_symmetry_number is not None:
            molecule_info["rotational_symmetry_number"] = (
                mol.rotational_symmetry_number
            )
        if mol.rotational_constants is not None:
            molecule_info["rotational_constants"] = mol.rotational_constants
        if mol.point_group is not None:
            molecule_info["point_group"] = mol.point_group
        if mol.dipole_moment is not None:
            molecule_info["dipole_moment"] = mol.dipole_moment
        if mol.dipole_moment_magnitude is not None:
            molecule_info["dipole_moment_magnitude"] = (
                mol.dipole_moment_magnitude
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
            "total_energy": (
                self.output.energies[-1] if self.output.energies else None
            ),
            "num_unpaired_electrons": self.output.num_unpaired_electrons,
            "homo_energy": self.output.homo_energy,
            "lumo_energy": self.output.lumo_energy,
            "alpha_homo_energy": self.output.alpha_homo_energy,
            "beta_homo_energy": self.output.beta_homo_energy,
            "alpha_lumo_energy": self.output.alpha_lumo_energy,
            "beta_lumo_energy": self.output.beta_lumo_energy,
            "somo_energies": self.output.somo_energies,
            "lowest_somo_energy": self.output.lowest_somo_energy,
            "highest_somo_energy": self.output.highest_somo_energy,
            "fmo_gap": self.output.fmo_gap,
            "alpha_fmo_gap": self.output.alpha_fmo_gap,
            "beta_fmo_gap": self.output.beta_fmo_gap,
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


class GaussianAssembler(BaseAssembler):
    OUTPUT_CLASS = Gaussian16Output
    PROGRAM = "Gaussian"

    def get_meta_data(self):
        meta_data = super().get_meta_data()
        meta_data.update(
            {
                "num_primitive_gaussians": self.output.num_primitive_gaussians,
                "num_cartesian_basis_functions": self.output.num_cartesian_basis_functions,
            }
        )
        if self.output.modredundant_group is not None:
            meta_data["modredundant_group"] = self.output.modredundant_group
        if is_custom_basis(self.output.basis):
            meta_data["custom_basis"] = {
                "heavy_elements": self.output.heavy_elements,
                "heavy_elements_basis": self.output.heavy_elements_basis,
                "heavy_elements_ecp": self.output.heavy_elements_ecp,
                "light_elements": self.output.light_elements,
                "light_elements_basis": self.output.light_elements_basis,
            }
        return meta_data

    def get_calculation_results(self):
        calculation_results = super().get_calculation_results()
        calculation_results.update(
            {
                "optimized_steps": self.output.optimized_steps,
                "reduced_masses": self.output.reduced_masses,
                "force_constants": self.output.force_constants,
                "ir_intensities": self.output.ir_intensities,
                "vibrational_mode_symmetries": self.output.vibrational_mode_symmetries,
            }
        )
        return calculation_results


class ORCAAssembler(BaseAssembler):
    OUTPUT_CLASS = ORCAOutput
    PROGRAM = "ORCA"

    def get_meta_data(self):
        meta_data = super().get_meta_data()
        meta_data.update(
            {
                "num_shells": self.output.num_shells,
            }
        )
        return meta_data

    def get_calculation_results(self):
        calculation_results = super().get_calculation_results()
        calculation_results.update(
            {
                "molar_absorption_coefficients": self.output.molar_absorption_coefficients,
            }
        )
        return calculation_results


class XTBAssembler(BaseAssembler):
    OUTPUT_CLASS = XTBOutput
    PROGRAM = "xTB"
    FOLDER_BASED = True

    def get_meta_data(self):
        meta_data = super().get_meta_data()
        meta_data.update(
            {
                "num_atomic_orbital": self.output.num_atomic_orbital,
                "num_shells": self.output.num_shells,
                "num_electrons": self.output.num_electrons,
                "max_iteration": self.output.max_iter,
                "pc_potential": self.output.pc_potential,
                "accuracy": self.output.accuracy,
            }
        )
        return meta_data

    def get_calculation_results(self):
        calculation_results = super().get_calculation_results()
        calculation_results.update(
            {
                "c6_coefficient": self.output.c6_coefficient,
                "c8_coefficient": self.output.c8_coefficient,
                "alpha_coefficient": self.output.alpha_coefficient,
            }
        )
        return calculation_results


class SingleFileAssembler:
    """Auto-detect program type and delegate to the appropriate file assembler."""

    def __init__(self, filename, index=":", include_failed=False):
        self.filename = filename
        self.index = index
        self.include_failed = include_failed

    @cached_property
    def assemble_data(self):
        assembler = self._get_assembler(self.filename)
        try:
            data = assembler.assemble()
        except Exception as e:
            logger.error(f"Error assembling {self.filename}: {e}")
            return None
        return data

    def _get_assembler(self, path):
        program = get_program_type_from_file(path)
        if program == "gaussian":
            return GaussianAssembler(
                filename=path,
                index=self.index,
                include_failed=self.include_failed,
            )
        if program == "orca":
            return ORCAAssembler(
                filename=path,
                index=self.index,
                include_failed=self.include_failed,
            )
        raise ValueError(
            f"Unsupported file '{path}'. "
            "Only 'gaussian' and 'orca' output files are supported."
        )


class SingleFolderAssembler:
    """Auto-detect program type and delegate to the appropriate folder assembler."""

    def __init__(self, folder, index=":", include_failed=False):
        self.folder = folder
        self.index = index
        self.include_failed = include_failed

    @cached_property
    def assemble_data(self):
        assembler = self._get_assembler(self.folder)
        try:
            data = assembler.assemble()
        except Exception as e:
            logger.error(f"Error assembling {self.folder}: {e}")
            return None
        return data

    def _get_assembler(self, path):
        program = BaseFolder(folder=path).get_program_type_from_folder()
        if program == "xtb":
            return XTBAssembler(
                folder=path,
                index=self.index,
                include_failed=self.include_failed,
            )
        if program == "mixed":
            raise ValueError(
                f"Folder '{path}' contains outputs from multiple programs."
            )
        raise ValueError(
            f"Unsupported folder '{path}'. "
            f"Only 'xtb' output folders are supported."
        )
