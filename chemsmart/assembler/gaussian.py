from chemsmart.assembler.base import BaseAssembler
from chemsmart.io.gaussian.output import Gaussian16Output


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
                "optimized_steps": self.output.optimized_steps,
                "reduced_masses": self.output.reduced_masses,
                "force_constants": self.output.force_constants,
                "ir_intensities": self.output.ir_intensities,
                "vibrational_mode_symmetries": self.output.vibrational_mode_symmetries,
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
