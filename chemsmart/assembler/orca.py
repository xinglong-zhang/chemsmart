from chemsmart.assembler.base import BaseAssembler
from chemsmart.io.orca.output import ORCAOutput


class ORCAAssembler(BaseAssembler):
    OUTPUT_CLASS = ORCAOutput
    PROGRAM = "ORCA"

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
