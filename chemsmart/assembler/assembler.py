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

    @cached_property
    def assemble(self):
        if not self.output.normal_termination:
            logger.warning(
                f"Calculation in {self.filename} did not terminate normally, skip assembling..."
            )
            return {}
        if not self.molecules_list:
            logger.error(f"No molecules parsed from {self.filename}.")
            return {}

        data = {
            "meta": self.get_meta_data(),
            "results": self.get_calculation_results(),
            "molecules": [],
        }

        for i, mol in enumerate(self.molecules_list):
            mol_entry = {"index": i + 1, **self.get_molecule_info(mol)}
            data["molecules"].append(mol_entry)

        return data

    def get_meta_data(self):
        meta_data = {
            "filename": self.filename,
            "program": self.program_name,
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


class GaussianAssembler(BaseAssembler):
    output_class = Gaussian16Output
    program_name = "Gaussian"

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
                # "num_steps": self.output.num_steps,
                # "intermediate_steps": self.output.intermediate_steps,
                # "optimized_steps": self.output.optimized_steps,
                # "reduced_masses": self.output.reduced_masses,
                # "force_constants": self.output.force_constants,
                # "ir_intensities": self.output.ir_intensities,
                # "vibrational_mode_symmetries": self.output.vibrational_mode_symmetries,
                # "has_frozen_coordinates": self.output.has_frozen_coordinates,
                # "frozen_coordinate_indices": self.output.frozen_coordinate_indices,
                # "free_coordinate_indices": self.output.free_coordinate_indices,
                # "frozen_elements": self.output.frozen_elements,
                # "free_elements": self.output.free_elements,
                # "frozen_atoms_masks": self.output.frozen_atoms_masks,
                # "scf_energies": self.output.scf_energies,
                # "mp2_energies": self.output.mp2_energies,
                # "oniom_energies": self.output.oniom_energies,
                # "convergence_criterion_not_met": self.output.convergence_criterion_not_met,
                # "input_orientations": self.output.input_orientations,
                # "input_orientations_pbc": self.output.input_orientations_pbc,
                # "standard_orientations": self.output.standard_orientations,
                # "standard_orientations_pbc": self.output.standard_orientations_pbc,
                # "tddft_transitions": self.output.tddft_transitions,
                # "excitation_energies_eV": self.output.excitation_energies_eV,
                # "absorptions_in_nm": self.output.absorptions_in_nm,
                # "oscillatory_strengths": self.output.oscillatory_strengths,
                # "transitions": self.output.transitions,
                # "contribution_coefficients": self.output.contribution_coefficients,
                # "alpha_occ_eigenvalues": self.output.alpha_occ_eigenvalues,
                # "alpha_virtual_eigenvalues": self.output.alpha_virtual_eigenvalues,
                # "beta_occ_eigenvalues": self.output.beta_occ_eigenvalues,
                # "beta_virtual_eigenvalues": self.output.beta_virtual_eigenvalues,
                # "num_unpaired_electrons": self.output.num_unpaired_electrons,
                # "somo_energy": self.output.somo_energy,
                # "mulliken_spin_densities": self.output.mulliken_spin_densities,
                # "moments_of_inertia_principal_axes": self.output.moments_of_inertia_principal_axes,
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
                # "max_ang_mom": self.output.max_ang_mom,
                # "contraction_scheme": self.output.contraction_scheme,
                # "coulomb_range_seperation": self.output.coulomb_range_seperation,
                # "exchange_range_seperation": self.output.exchange_range_seperation,
                # "finite_nucleus_model": self.output.finite_nucleus_model,
                # "aux_j_fitting_basis": self.output.aux_j_fitting_basis,
                # "aux_j_num_basis_functions": self.output.aux_j_num_basis_functions,
                # "aux_j_num_shells": self.output.aux_j_num_shells,
                # "aux_j_max_ang_mom": self.output.aux_j_max_ang_mom,
                # "aux_jk_fitting_basis": self.output.aux_jk_fitting_basis,
                # "aux_k_fitting_basis": self.output.aux_k_fitting_basis,
                # "aux_external_fitting_basis": self.output.aux_external_fitting_basis,
                # "integral_threshold": self.output.integral_threshold,
                # "primitive_cutoff": self.output.primitive_cutoff,
                # "primitive_pair_threshold": self.output.primitive_pair_threshold,
                # "ri_approx": self.output.ri_approx,
                # "rij_cosx": self.output.rij_cosx,
                # "num_electrons": self.output.num_electrons,
                # "basis_dim": self.output.basis_dim,
                # "diis_acceleration": self.output.diis_acceleration,
                # "scf_maxiter": self.output.scf_maxiter,
                # "scf_convergence": self.output.scf_convergence,
                # "dipole": self.output.dipole,
                # "quadrupole": self.output.quadrupole,
            }
        )
        return meta_data

    def get_calculation_results(self):
        results = super().get_calculation_results()
        results.update(
            {
                # "constrained_bond_lengths": self.output.constrained_bond_lengths,
                # "constrained_bond_angles": self.output.constrained_bond_angles,
                # "constrained_dihedral_angles": self.output.constrained_dihedral_angles,
                # "converged": self.output.converged,
                # "final_scf_energy": self.output.final_scf_energy,
                # "final_energy": self.output.final_energy,
                # "single_point_energy": self.output.single_point_energy,
                # "single_point_energy_eV": self.output.single_point_energy_eV,
                # "final_nuclear_repulsion": self.output.final_nuclear_repulsion,
                # "final_nuclear_repulsion_eV": self.output.final_nuclear_repulsion_eV,
                # "final_electronic_energy": self.output.final_electronic_energy,
                # "final_electronic_energy_eV": self.output.final_electronic_energy_eV,
                # "one_electron_energy": self.output.one_electron_energy,
                # "one_electron_energy_eV": self.output.one_electron_energy_eV,
                # "two_electron_energy": self.output.two_electron_energy,
                # "two_electron_energy_eV": self.output.two_electron_energy_eV,
                # "max_cosx_asymmetry_energy": self.output.max_cosx_asymmetry_energy,
                # "max_cosx_asymmetry_energy_eV": self.output.max_cosx_asymmetry_energy_eV,
                # "potential_energy": self.output.potential_energy,
                # "potential_energy_eV": self.output.potential_energy_eV,
                # "kinetic_energy": self.output.kinetic_energy,
                # "kinetic_energy_eV": self.output.kinetic_energy_eV,
                # "virial_ratio": self.output.virial_ratio,
                # "xc_energy": self.output.xc_energy,
                # "xc_energy_eV": self.output.xc_energy_eV,
                # "dfet_embed_energy": self.output.dfet_embed_energy,
                # "dfet_embed_energy_eV": self.output.dfet_embed_energy_eV,
                # "orbital_occupancy": self.output.orbital_occupancy,
                # "orbital_energies": self.output.orbital_energies,
                # "loewdin_atomic_charges": self.output.loewdin_atomic_charges,
                # "mayer_mulliken_gross_atomic_population": self.output.mayer_mulliken_gross_atomic_population,
                # "mayer_total_nuclear_charge": self.output.mayer_total_nuclear_charge,
                # "mayer_mulliken_gross_atomic_charge": self.output.mayer_total_nuclear_charge,
                # "mayer_total_valence": self.output.mayer_total_valence,
                # "mayer_bonded_valence": self.output.mayer_bonded_valence,
                # "mayer_free_valence": self.output.mayer_free_valence,
                # "mayer_bond_orders_larger_than_zero_point_one": self.output.mayer_bond_orders_larger_than_zero_point_one,
                # "total_integrated_alpha_density": self.output.total_integrated_alpha_density,
                # "total_integrated_beta_density": self.output.total_integrated_beta_density,
                # "dipole_moment_electric_contribution": self.output.dipole_moment_electric_contribution,
                # "dipole_moment_nuclear_contribution": self.output.dipole_moment_nuclear_contribution,
                # "total_dipole_moment": self.output.total_dipole_moment,
                # "dipole_moment_in_au": self.output.dipole_moment_in_au,
                # "dipole_moment_in_debye": self.output.dipole_moment_in_debye,
                # "dipole_moment_along_axis_in_au": self.output.dipole_moment_along_axis_in_au,
                # "dipole_moment_along_axis_in_debye": self.output.dipole_moment_along_axis_in_debye,
                # "vib_freq_scale_factor": self.output.vib_freq_scale_factor,
                # "molar_absorption_coefficients": self.output.molar_absorption_coefficients,
                # "integrated_absorption_coefficients": self.output.integrated_absorption_coefficients,
                # "transition_dipole_deriv_norm": self.output.transition_dipole_deriv_norm,
            }
        )
        return results


class SingleFileAssembler:
    def __init__(self, filename, index="-1"):
        self.filename = filename
        self.index = index

    @cached_property
    def assemble_data(self):
        assembler = self._get_assembler(self.filename)
        try:
            data = assembler.assemble
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
