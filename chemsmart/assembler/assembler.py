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
            "num_basis_functions": self.output.num_basis_functions,
            "job_type": self.output.job_type,
            "solvent_on": self.output.solvent_on,
            "solvent_model": self.output.solvent_model,
            "solvent_id": self.output.solvent_id,
            "route_string": self.output.route_string,
        }
        return meta_data

    def get_molecule_info(self, mol):
        molecule_info = {
            "charge": mol.charge,
            "multiplicity": mol.multiplicity,
            "chemical_symbols": mol.chemical_symbols,
            "positions": mol.positions,
            "chemical_formula": mol.chemical_formula,
            "empirical_formula": mol.empirical_formula,
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
            "input_coordinates_block": self.output.input_coordinates_block,
            "all_structures": self.output.all_structures,
            "optimized_structure": self.output.optimized_structure,
            "last_structure": self.output.last_structure,
            "molecule": self.output.molecule,
            "optimized_steps_indices": self.output.optimized_steps_indices,
            "zero_point_energy": self.output.zero_point_energy,
            "homo_energy": self.output.homo_energy,
            "lumo_energy": self.output.lumo_energy,
            "fmo_gap": self.output.fmo_gap,
            "cpu_runtime_by_jobs_core_hours": self.output.cpu_runtime_by_jobs_core_hours,
            "service_units_by_jobs": self.output.service_units_by_jobs,
            "total_core_hours": self.output.total_core_hours,
            "total_service_unit": self.output.total_service_unit,
            "elapsed_walltime_by_jobs": self.output.elapsed_walltime_by_jobs,
            "total_elapsed_walltime": self.output.total_elapsed_walltime,
            "vibrational_frequencies": self.output.vibrational_frequencies,
            "num_vib_frequencies": self.output.num_vib_frequencies,
            "rotational_symmetry_number": self.output.rotational_symmetry_number,
            "energies": self.output.energies,
            "has_forces": self.output.has_forces,
            "num_forces": self.output.num_forces,
            "mulliken_atomic_charges": self.output.mulliken_atomic_charges,
            "hirshfeld_charges": self.output.hirshfeld_charges,
            "hirshfeld_spin_densities": self.output.hirshfeld_spin_densities,
        }
        return calculation_results


class GaussianAssembler(BaseAssembler):
    output_class = Gaussian16Output
    program_name = "Gaussian"

    def get_meta_data(self):
        meta_data = super().get_meta_data()
        meta_data.update(
            {
                "spin": self.output.spin,
                "gen_genecp": self.output.gen_genecp,
                "modredundant_group": self.output.modredundant_group,
                "num_primitive_gaussians": self.output.num_primitive_gaussians,
                "num_cartesian_basis_functions": self.output.num_cartesian_basis_functions,
            }
        )
        return meta_data

    def get_calculation_results(self):
        calculation_results = super().get_calculation_results()
        calculation_results.update(
            {
                "num_steps": self.output.num_steps,
                "intermediate_steps": self.output.intermediate_steps,
                "optimized_steps": self.output.optimized_steps,
                "reduced_masses": self.output.reduced_masses,
                "force_constants": self.output.force_constants,
                "ir_intensities": self.output.ir_intensities,
                "vibrational_mode_symmetries": self.output.vibrational_mode_symmetries,
                "vibrational_modes": self.output.vibrational_modes,
                "num_vib_modes": self.output.num_vib_modes,
                "has_frozen_coordinates": self.output.has_frozen_coordinates,
                "frozen_coordinate_indices": self.output.frozen_coordinate_indices,
                "free_coordinate_indices": self.output.free_coordinate_indices,
                "frozen_elements": self.output.frozen_elements,
                "free_elements": self.output.free_elements,
                "frozen_atoms_masks": self.output.frozen_atoms_masks,
                "scf_energies": self.output.scf_energies,
                "mp2_energies": self.output.mp2_energies,
                "oniom_energies": self.output.onionom_energies,
                "convergence_criterion_not_met": self.output.convergence_criterion_not_met,
                "input_orientations": self.output.input_orientations,
                "input_orientations_pbc": self.output.input_orientations_pbc,
                "standard_orientations": self.output.standard_orientations,
                "standard_orientations_pbc": self.output.standard_orientations_pbc,
                "tddft_transitions": self.output.tddft_transitions,
                "excitation_energies_eV": self.output.excitation_energies_eV,
                "absorptions_in_nm": self.output.absorptions_in_nm,
                "oscillatory_strengths": self.output.oscillatory_strengths,
                "transitions": self.output.transitions,
                "contribution_coefficients": self.output.contribution_coefficients,
                "alpha_occ_eigenvalues": self.output.alpha_occ_eigenvalues,
                "alpha_virtual_eigenvalues": self.output.alpha_virtual_eigenvalues,
                "beta_occ_eigenvalues": self.output.beta_occ_eigenvalues,
                "beta_virtual_eigenvalues": self.output.beta_virtual_eigenvalues,
                "num_unpaired_electrons": self.output.num_unpaired_electrons,
                "somo_energy": self.output.somo_energy,
                "mulliken_spin_densities": self.output.mulliken_spin_densities,
                "mulliken_atomic_charges_heavy_atoms": self.output.mulliken_atomic_charges_heavy_atoms,
                "mulliken_spin_densities_heavy_atoms": self.output.mulliken_spin_densities_heavy_atoms,
                "hirshfeld_dipoles": self.output.hirshfeld_dipoles,
                "hirshfeld_cm5_charges": self.output.hirshfeld_cm5_charges,
                "hirshfeld_charges_heavy_atoms": self.output.hirshfeld_charges_heavy_atoms,
                "hirshfeld_spin_densities_heavy_atoms": self.output.hirshfeld_spin_densities_heavy_atoms,
                "hirshfeld_cm5_charges_heavy_atoms": self.output.hirshfeld_cm5_charges_heavy_atoms,
                "moments_of_inertia_principal_axes": self.output.moments_of_inertia_principal_axes,
                "rotational_temperatures": self.output.rotational_temperatures,
                "rotational_constants_in_Hz": self.output.rotational_constants_in_Hz,
            }
        )
        return calculation_results


class ORCAAssembler(BaseAssembler):
    output_class = ORCAOutput
    program_name = "ORCA"

    def get_meta_data(self):
        meta_data = super().get_meta_data()
        meta_data.update(
            {
                "num_shells": self.output.num_shells,
                "max_ang_mom": self.output.max_ang_mom,
                "contraction_scheme": self.output.contraction_scheme,
                "coulomb_range_seperation": self.output.coulomb_range_seperation,
                "exchange_range_seperation": self.output.exchange_range_seperation,
                "finite_nucleus_model": self.output.finite_nucleus_model,
                "aux_j_fitting_basis": self.output.aux_j_fitting_basis,
                "aux_j_num_basis_functions": self.output.aux_j_num_basis_functions,
                "aux_j_num_shells": self.output.aux_j_num_shells,
                "aux_j_max_ang_mom": self.output.aux_j_max_ang_mom,
                "aux_jk_fitting_basis": self.output.aux_jk_fitting_basis,
                "aux_k_fitting_basis": self.output.aux_k_fitting_basis,
                "aux_external_fitting_basis": self.output.aux_external_fitting_basis,
                "integral_threshold": self.output.integral_threshold,
                "primitive_cutoff": self.output.primitive_cutoff,
                "primitive_pair_threshold": self.output.primitive_pair_threshold,
                "ri_approx": self.output.ri_approx,
                "rij_cosx": self.output.rij_cosx,
                "num_electrons": self.output.num_electrons,
                "basis_dim": self.output.basis_dim,
                "diis_acceleration": self.output.diis_acceleration,
                "scf_maxiter": self.output.scf_maxiter,
                "scf_convergence": self.output.scf_convergence,
                "dipole": self.output.dipole,
                "quadrupole": self.output.quadrupole,
                "temperature_in_K": self.output.temperature_in_K,
                "pressure_in_atm": self.output.pressure_in_atm,
            }
        )
        return meta_data

    def get_calculation_results(self):
        results = super().get_calculation_results()
        results.update(
            {
                "constrained_bond_lengths": self.output.constrained_bond_lengths,
                "constrained_bond_angles": self.output.constrained_bond_angles,
                "constrained_dihedral_angles": self.output.constrained_dihedral_angles,
                "converged": self.output.converged,
                "final_scf_energy": self.output.final_scf_energy,
                "final_energy": self.output.final_energy,
                "single_point_energy": self.output.single_point_energy,
                "single_point_energy_eV": self.output.single_point_energy_eV,
                "final_nuclear_repulsion": self.output.final_nuclear_repulsion,
                "final_nuclear_repulsion_eV": self.output.final_nuclear_repulsion_eV,
                "final_electronic_energy": self.output.final_electronic_energy,
                "final_electronic_energy_eV": self.output.final_electronic_energy_eV,
                "one_electron_energy": self.output.one_electron_energy,
                "one_electron_energy_eV": self.output.one_electron_energy_eV,
                "two_electron_energy": self.output.two_electron_energy,
                "two_electron_energy_eV": self.output.two_electron_energy_eV,
                "max_cosx_asymmetry_energy": self.output.max_cosx_asymmetry_energy,
                "max_cosx_asymmetry_energy_eV": self.output.max_cosx_asymmetry_energy_eV,
                "potential_energy": self.output.potential_energy,
                "potential_energy_eV": self.output.potential_energy_eV,
                "kinetic_energy": self.output.kinetic_energy,
                "kinetic_energy_eV": self.output.kinetic_energy_eV,
                "virial_ratio": self.output.virial_ratio,
                "xc_energy": self.output.xc_energy,
                "xc_energy_eV": self.output.xc_energy_eV,
                "dfet_embed_energy": self.output.dfet_embed_energy,
                "dfet_embed_energy_eV": self.output.dfet_embed_energy_eV,
                "orbital_occupancy": self.output.orbital_occupancy,
                "orbital_energies": self.output.orbital_energies,
                "loewdin_atomic_charges": self.output.loewdin_atomic_charges,
                "mayer_mulliken_gross_atomic_population": self.output.mayer_mulliken_gross_atomic_population,
                "mayer_total_nuclear_charge": self.output.mayer_total_nuclear_charge,
                "mayer_mulliken_gross_atomic_charge": self.output.mayer_total_nuclear_charge,
                "mayer_total_valence": self.output.mayer_total_valence,
                "mayer_bonded_valence": self.output.mayer_bonded_valence,
                "mayer_free_valence": self.output.mayer_free_valence,
                "mayer_bond_orders_larger_than_zero_point_one": self.output.mayer_bond_orders_larger_than_zero_point_one,
                "total_integrated_alpha_density": self.output.total_integrated_alpha_density,
                "total_integrated_beta_density": self.output.total_integrated_beta_density,
                "dipole_moment_electric_contribution": self.output.dipole_moment_electric_contribution,
                "dipole_moment_nuclear_contribution": self.output.dipole_moment_nuclear_contribution,
                "total_dipole_moment": self.output.total_dipole_moment,
                "dipole_moment_in_au": self.output.dipole_moment_in_au,
                "dipole_moment_in_debye": self.output.dipole_moment_in_debye,
                "dipole_moment_along_axis_in_au": self.output.dipole_moment_along_axis_in_au,
                "dipole_moment_along_axis_in_debye": self.output.dipole_moment_along_axis_in_debye,
                "rotational_constants_in_wavenumbers": self.output.rotational_constants_in_wavenumbers,
                "rotational_constants_in_MHz": self.output.rotational_constants_in_MHz,
                "vib_freq_scale_factor": self.output.vib_freq_scale_factor,
                "molar_absorption_coefficients": self.output.molar_absorption_coefficients,
                "integrated_absorption_coefficients": self.output.integrated_absorption_coefficients,
                "transition_dipole_deriv_norm": self.output.transition_dipole_deriv_norm,
                "num_translation_and_rotation_modes": self.output.num_translation_and_rotation_modes,
                "num_vibration_modes": self.output.num_vibration_modes,
                "internal_energy": self.output.internal_energy,
                "internal_energy_in_eV": self.output.internal_energy_in_eV,
                "electronic_energy": self.output.electronic_energy,
                "electronic_energy_in_eV": self.output.electronic_energy_in_eV,
                "zero_point_energy_in_eV": self.output.zero_point_energy_in_eV,
                "thermal_vibration_correction": self.output.thermal_vibration_correction,
                "thermal_vibration_correction_in_eV": self.output.thermal_vibration_correction_in_eV,
                "thermal_rotation_correction": self.output.thermal_rotation_correction,
                "thermal_rotation_correction_in_eV": self.output.thermal_rotation_correction_in_eV,
                "thermal_translation_correction": self.output.thermal_translation_correction,
                "thermal_translation_correction_in_eV": self.output.thermal_translation_correction_in_eV,
                "total_thermal_correction_due_to_trans_rot_vib": self.output.total_thermal_correction_due_to_trans_rot_vib,
                "total_correction": self.output.total_correction,
                "enthalpy": self.output.enthalpy,
                "enthalpy_in_eV": self.output.enthalpy_in_eV,
                "thermal_enthalpy_correction": self.output.thermal_enthalpy_correction,
                "thermal_enthalpy_correction_in_eV": self.output.thermal_enthalpy_correction_in_eV,
                "electronic_entropy_no_temperature_in_SI": self.output.electronic_entropy_no_temperature_in_SI,
                "vibrational_entropy_no_temperature_in_SI": self.output.vibrational_entropy_no_temperature_in_SI,
                "rotational_entropy_no_temperature_in_SI": self.output.rotational_entropy_no_temperature_in_SI,
                "translational_entropy_no_temperature_in_SI": self.output.translational_entropy_no_temperature_in_SI,
                "entropy_in_J_per_mol_per_K": self.output.entropy_in_J_per_mol_per_K,
                "entropy_TS": self.output.entropy_TS,
                "rotational_entropy_symmetry_correction_J_per_mol_per_K": self.output.rotational_entropy_symmetry_correction_J_per_mol_per_K,
                "gibbs_free_energy": self.output.gibbs_free_energy,
                "gibbs_free_energy_in_eV": self.output.gibbs_free_energy_in_eV,
            }
        )
        return results


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
