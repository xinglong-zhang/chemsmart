"""
Assembler module for parsing quantum chemistry output files and
assembling structured records for database storage.

This module consolidates the assembler pipeline:
- BaseAssembler: core assembly logic shared across programs
- GaussianAssembler: Gaussian-specific assembly
- ORCAAssembler: ORCA-specific assembly
- SingleFileAssembler: dispatcher that auto-detects program type
- build_provenance: provenance metadata builder
"""

import logging
from functools import cached_property

from chemsmart import __version__ as chemsmart_version
from chemsmart.assembler.records import AssembledRecord
from chemsmart.assembler.utils import (
    canonical_geometry_string,
    file_size,
    get_record_id,
    is_custom_basis,
    is_custom_solvent,
    sha256_content,
    utcnow_iso,
)
from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.orca.output import ORCAOutput
from chemsmart.utils.io import get_program_type_from_file

logger = logging.getLogger(__name__)


def build_provenance(filename, output):
    """Build provenance metadata for an assembled record."""
    return {
        "source_file": filename,
        "source_file_hash": sha256_content(output),
        "source_file_size": file_size(filename),
        "source_file_date": output.date,
        "program": get_program_type_from_file(filename),
        "program_version": output.version,
        "parser": output.__class__.__name__,
        "chemsmart_version": chemsmart_version,
        "assembled_at": utcnow_iso(),
    }


class BaseAssembler:
    OUTPUT_CLASS = None  # To be defined in subclasses
    PROGRAM = "unknown"  # To be defined in subclasses

    def __init__(self, filename, index=":"):
        self.filename = filename
        self.index = index
        self.output = self.OUTPUT_CLASS(filename)

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

        meta = self.get_meta_data()
        results = self.get_calculation_results()
        provenance = build_provenance(self.filename, self.output)
        molecules = []
        for i, mol in enumerate(self.molecules_list):
            mol_entry = {"index": i + 1, **self.get_molecule_info(mol)}
            molecules.append(mol_entry)

        # Use the last molecule (typically the optimized structure) for the ID
        ref_mol = self.molecules_list[-1]
        canon_geom = canonical_geometry_string(
            ref_mol.chemical_symbols, ref_mol.positions
        )
        record_id = get_record_id(
            canonical_geometry=canon_geom,
            charge=ref_mol.charge,
            multiplicity=ref_mol.multiplicity,
            program=provenance.get("program", "unknown"),
            functional=meta.get("functional", ""),
            basis=meta.get("basis", ""),
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
        meta_data = {
            "functional": self.output.functional,
            "basis": "customized_basis" if is_custom_basis(basis) else basis,
            "num_basis_functions": self.output.num_basis_functions,
            "spin": self.output.spin,
            "jobtype": self.output.jobtype,
            "solvent_on": self.output.solvent_on,
            "route_string": self.output.route_string,
        }
        if self.output.solvent_on:
            solvent_id = self.output.solvent_id
            meta_data.update(
                {
                    "solvent_model": self.output.solvent_model,
                    "solvent_id": (
                        "customized_solvent"
                        if is_custom_solvent(solvent_id)
                        else solvent_id
                    ),
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
            "structure_index_in_file": mol.structure_index_in_file,
            "is_optimized_structure": mol.is_optimized_structure,
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
            "is_ring": mol.is_ring,
            "is_aromatic": mol.is_aromatic,
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


class SingleFileAssembler:
    """Auto-detect program type and delegate to the appropriate assembler."""

    def __init__(self, filename, index=":"):
        self.filename = filename
        self.index = index

    @cached_property
    def assemble_data(self):
        assembler = self._get_assembler(self.filename)
        try:
            data = assembler.assemble()
        except Exception as e:
            logger.error(f"Error assembling {self.filename}: {e}")
            return None
        return data

    def _get_assembler(self, file):
        program = get_program_type_from_file(self.filename)
        if program == "gaussian":
            assembler = GaussianAssembler(file, index=self.index)
        elif program == "orca":
            assembler = ORCAAssembler(file, index=self.index)
        else:
            raise ValueError(
                "Unsupported format. Only 'gaussian' and 'orca' are supported."
            )
        return assembler
