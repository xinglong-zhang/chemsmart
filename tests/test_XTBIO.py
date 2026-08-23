import os.path
from unittest.mock import MagicMock, PropertyMock, patch

import numpy as np
import pytest

from chemsmart.io.xtb.file import (
    XTBChargesFile,
    XTBEnergyFile,
    XTBEngradFile,
    XTBG98File,
    XTBGradientFile,
    XTBHessianFile,
    XTBMainOut,
    XTBVibSpectrumFile,
    XTBWibergBondOrderFile,
)
from chemsmart.io.xtb.folder import XTBFolder
from chemsmart.io.xtb.input import XTBInput
from chemsmart.io.xtb.output import XTBOutput
from chemsmart.io.xtb.route import XTBRoute
from chemsmart.utils.constants import kcal_per_mol_to_hartree


class TestXTBRoute:
    """Tests for XTBRoute class."""

    def test_read_route(self):
        s1 = "xtb mol.xyz --gfn 2"
        r1 = XTBRoute(route_string=s1)
        assert r1.method == "gfn2"
        assert r1.gfn_version == "gfn2"
        assert r1.basis == "default"
        assert r1.charge is None
        assert r1.uhf is None
        assert r1.jobtype == "sp"
        assert r1.optimization_level is None
        assert r1.solvent_model is None
        assert r1.solvent_id is None
        assert r1.freq is False
        assert r1.grad is False

        s2 = (
            "xtb mol.xyz "
            "--gfn2 "
            "--chrg -1 "
            "--uhf 2 "
            "--alpb water "
            "--opt tight "
            "--grad "
            "--acc 0.1 "
            "--etemp 400"
        )
        r2 = XTBRoute(route_string=s2)
        assert r2.method == "gfn2"
        assert r2.gfn_version == "gfn2"
        assert r2.charge == -1
        assert r2.uhf == 2
        assert r2.jobtype == "opt"
        assert r2.solvent_model == "alpb"
        assert r2.solvent_id == "water"
        assert r2.optimization_level == "tight"
        assert r2.grad is True
        assert r2.freq is False
        assert r2.accuracy == 0.1
        assert r2.electronic_temperature == 400.0

        # GBSA solvent and Hessian job
        s3 = "xtb mol.xyz --gfn 1 --gbsa toluene --hess"
        r3 = XTBRoute(route_string=s3)
        assert r3.method == "gfn1"
        assert r3.solvent_model == "gbsa"
        assert r3.solvent_id == "toluene"
        assert r3.jobtype == "hess"
        assert r3.freq is True

        # Short charge/uhf flags, COSMO, and property printouts
        s4 = (
            "xtb mol.xyz "
            "--gfn 0 "
            "-c 1 "
            "-u 1 "
            "--cosmo thf "
            "--pop "
            "--wbo "
            "--dipole"
        )
        r4 = XTBRoute(route_string=s4)
        assert r4.method == "gfn0"
        assert r4.charge == 1
        assert r4.uhf == 1
        assert r4.solvent_model == "cosmo"
        assert r4.solvent_id == "thf"
        assert r4.jobtype == "sp"
        assert r4.mulliken_population is True
        assert r4.wbo is True
        assert r4.dipole is True

        # verytight normalized to vtight; CPCMX solvent
        s5 = "xtb mol.xyz --gfn 2 --opt verytight --cpcmx water"
        r5 = XTBRoute(route_string=s5)
        assert r5.optimization_level == "vtight"
        assert r5.jobtype == "opt"
        assert r5.solvent_model == "cpcmx"
        assert r5.solvent_id == "water"

        # GFN-FF and combined optimize+Hessian flags
        s6 = "xtb mol.xyz --gfnff --ohess"
        r6 = XTBRoute(route_string=s6)
        assert r6.method == "gfnff"
        assert r6.gfn_version == "gfnff"
        assert r6.jobtype == "opt"
        assert r6.freq is True

        # --gff alias and MD job type
        s7 = "xtb mol.xyz --gff --md"
        r7 = XTBRoute(route_string=s7)
        assert r7.gfn_version == "gfnff"
        assert r7.method == "gfnff"
        assert r7.jobtype == "md"

        # Boolean property flags
        s8 = (
            "xtb mol.xyz "
            "--gfn 2 "
            "--ptb "
            "--spinpol "
            "--ceh "
            "--molden "
            "--lmo "
            "--fod "
            "--esp "
            "--stm "
            "--vip "
            "--vea "
            "--vipea "
            "--vfukui "
            "--vomega "
            "--alpha "
            "--cma"
        )
        r8 = XTBRoute(route_string=s8)
        assert r8.jobtype == "sp"
        assert r8.ptb is True
        assert r8.spin_polarization is True
        assert r8.charge_extended_hueckel is True
        assert r8.molden_file is True
        assert r8.localized_molecular_orbitals is True
        assert r8.fractional_occupation_density is True
        assert r8.electrostatic_potential is True
        assert r8.stm_image is True
        assert r8.vertical_ionization_potential is True
        assert r8.vertical_electron_affinity is True
        assert r8.vertical_ionization_and_affinity is True
        assert r8.fukui_indices is True
        assert r8.electrophilicity_index is True
        assert r8.polarizability is True
        assert r8.center_of_mass_transform is True


class TestXTBInput:
    def test_default_input(self, xtb_default_inputfile):
        assert os.path.exists(xtb_default_inputfile)
        xtb_input = XTBInput(filename=xtb_default_inputfile)
        assert xtb_input.charge == 0
        assert xtb_input.spin == 0
        assert xtb_input.method == "GFN2-xTB"
        assert xtb_input.scc
        assert xtb_input.periodic is False
        assert xtb_input.dispersion_energy_scale == 1.0
        assert xtb_input.max_iterations == 250
        assert xtb_input.electronic_temperature == 300.0
        assert xtb_input.broyden_damping == 0.4
        assert xtb_input.guess_charges is None
        assert xtb_input.engine == "rf"
        assert xtb_input.optimization_level == "normal"
        assert xtb_input.anc_microcycles == 20
        assert xtb_input.max_optcycles == 0
        assert xtb_input.max_displacement == 1.0
        assert xtb_input.low_frequency_cutoff == 0.01
        assert xtb_input.hessian_model == "old"
        assert xtb_input.s6_in_model_hessian == 20.0
        assert xtb_input.stretch_force_constant == 0.4
        assert xtb_input.bend_force_constant == 0.13
        assert xtb_input.torsion_force_constant == 0.75e-2
        assert xtb_input.out_of_plane_force_constant == 0.0
        assert xtb_input.additional_vdw_contribution == 0.0
        assert xtb_input.electrostatic_contribution == 0.0
        assert xtb_input.distance_cutoff == 8.366600265340756
        assert xtb_input.exact_rational_function is False
        assert xtb_input.average_convergence is False
        assert xtb_input.thermo_temperature == 298.15
        assert xtb_input.rotor_cutoff == 50.0
        assert xtb_input.imaginary_frequency_cutoff == -20.0
        assert xtb_input.scaling_factor == 1.0
        assert xtb_input.md_temperature == 298.15
        assert xtb_input.md_time == 50.0
        assert xtb_input.dump_structure == 50.0
        assert xtb_input.velocity_in_trj is False
        assert xtb_input.nvt_ensemble
        assert xtb_input.skip_interval == 500
        assert xtb_input.md_step == 4.0
        assert xtb_input.hydrogen_mass == 4
        assert xtb_input.shake_algorithm == 2
        assert xtb_input.md_scc_accuracy == 2.0
        assert xtb_input.force_writing_restart is False
        assert xtb_input.hess_scc_accuracy == 0.3
        assert xtb_input.hess_step == 0.005
        assert xtb_input.hess_scale == 1.0
        assert xtb_input.modef_n == 31
        assert xtb_input.modef_step == 1.0
        assert xtb_input.modef_update == 0.2
        assert xtb_input.modef_local == 0
        assert xtb_input.modef_threshold == 0.0
        assert xtb_input.projected_mode == 0
        assert xtb_input.mode_following == 7
        assert xtb_input.cube_step == 0.4
        assert xtb_input.density_matrix_threshold == 0.05
        assert xtb_input.boundary_offset == 3.0
        assert xtb_input.cube_output == 1
        assert xtb_input.symmetry_threshold == 0.1
        assert xtb_input.symmetry_max_atoms == 200
        assert xtb_input.atom_type == 7
        assert xtb_input.isotropic_electrostatic
        assert xtb_input.pathfinder_runs == 3
        assert xtb_input.path_points == 50
        assert xtb_input.path_optimization_steps == 3
        assert xtb_input.rmsd_push_factor == 0.05
        assert xtb_input.rmsd_pull_factor == -0.04
        assert xtb_input.rmsd_width == 0.7
        assert xtb_input.wall_potential == "polynomial"
        assert xtb_input.wall_potential_exponent == 30
        assert xtb_input.logfermi_bias_exponent == 6.0
        assert xtb_input.wall_temperature == 300.0
        assert xtb_input.auto_scale == 1.0
        assert xtb_input.axis_shift == 3.5

    def test_sp_alpb_input(self, xtb_sp_alpb_inputfile):
        assert os.path.exists(xtb_sp_alpb_inputfile)
        xtb_input = XTBInput(filename=xtb_sp_alpb_inputfile)
        assert xtb_input.charge == 0
        assert xtb_input.spin == 0


class TestXTBMainOut:
    """Tests for XTBMainOut class."""

    def test_main_out_co2(self, xtb_co2_outfolder):
        """Test parsing main output from CO2 ohess calculation."""
        xtb_main_out_file = os.path.join(xtb_co2_outfolder, "co2_ohess.out")
        assert os.path.exists(xtb_main_out_file)
        co2_main_out = XTBMainOut(xtb_main_out_file)
        assert co2_main_out.version == "6.7.1"
        assert co2_main_out.normal_termination
        assert (
            co2_main_out.route_string
            == "xtb co2.xyz --ohess vtight --grad --copy"
        )
        assert not co2_main_out.solvent_on
        # GFN2-xTB Setup
        assert co2_main_out.num_basis_functions == 12
        assert co2_main_out.num_atomic_orbital == 12
        assert co2_main_out.num_shells == 6
        assert co2_main_out.num_electrons == 16
        assert co2_main_out.max_iter == 250
        assert co2_main_out.hamiltonian == "GFN2-xTB"
        assert co2_main_out.method == "gfn2"
        assert co2_main_out.basis == "default"
        assert co2_main_out.custom_solvent is None
        assert not co2_main_out.restart
        assert not co2_main_out.solvent_on
        assert not co2_main_out.pc_potential
        assert co2_main_out.electronic_temperature == 300.0
        assert co2_main_out.temperature_in_K == 298.15
        assert co2_main_out.accuracy == 1.0
        assert co2_main_out.integral_cutoff == 25.0
        assert co2_main_out.integral_neglect == 1e-8
        assert co2_main_out.scf_convergence == 1.0e-6
        assert co2_main_out.wf_convergence == 1.0e-4
        assert co2_main_out.broyden_damping == 0.4
        assert co2_main_out.net_charge == 0
        assert co2_main_out.unpaired_electrons == 0
        # Geometry Optimization Setup
        assert co2_main_out.optimization_level == "vtight"
        assert co2_main_out.max_optcycles == 200
        assert co2_main_out.anc_microcycles == 20
        assert co2_main_out.degrees_of_freedom == 4
        assert co2_main_out.rf_solver == "davidson"
        assert co2_main_out.write_all_intermediate_geometries
        assert co2_main_out.is_linear
        assert co2_main_out.energy_convergence == 1.0e-7
        assert co2_main_out.gradient_convergence == 2.0e-4
        assert co2_main_out.max_rf_displacement == 1.0
        assert co2_main_out.low_frequency_cutoff == 0.01
        assert co2_main_out.max_frequency_cutoff == 5.0
        assert co2_main_out.s6_in_model_hessian == 20.0
        # Geometry Optimization Results
        assert co2_main_out.geometry_optimization_converged
        assert co2_main_out.optimized_structure_block == [
            "3",
            "xtb: 6.7.1 (edcfbbe)",
            "O           -1.14365140481883        0.00000000000000        0.00000000000000",
            "O            1.14365140481883       -0.00000000000000        0.00000000000000",
            "C            0.00000000000000       -0.00000000000000       -0.00000000000000",
            "",
        ]
        assert co2_main_out.scc_energy == -10.430605117263
        assert co2_main_out.isotropic_es == 0.032324567807
        assert co2_main_out.anisotropic_es == 0.003405663023
        assert co2_main_out.anisotropic_xc == 0.000432280404
        assert co2_main_out.dispersion_energy == -0.000687152300
        assert co2_main_out.solvation_energy_gsolv is None
        assert co2_main_out.electronic_solvation_energy_gelec is None
        assert co2_main_out.surface_area_solvation_energy_gsasa is None
        assert co2_main_out.hydrogen_bonding_solvation_energy_ghb is None
        assert co2_main_out.empirical_shift_correction_gshift is None
        assert co2_main_out.repulsion_energy == 0.122152828089
        assert co2_main_out.additional_restraining_energy == 0.0
        assert co2_main_out.total_charge == 0
        assert co2_main_out.energies == [
            -10.2973989,
            -10.3084470,
            -10.3084521,
            -10.3084522,
            -10.3084523,
        ]
        # Hessian Setup
        assert co2_main_out.numfreq
        assert co2_main_out.hessian_step_length == 0.00500
        assert co2_main_out.scc_accuracy == 0.30000
        assert co2_main_out.hessian_scale_factor == 1.00000
        assert co2_main_out.rms_gradient == 0.00000
        # Hessian Results
        assert co2_main_out.homo_energy == -14.5428
        assert co2_main_out.lumo_energy == -6.0942
        assert co2_main_out.c6_coefficient == 174.800200
        assert co2_main_out.c8_coefficient == 4029.884814
        assert co2_main_out.alpha_coefficient == 19.088396
        assert np.allclose(
            co2_main_out.molecular_dipole_qonly, [0.0, -0.0, -0.0]
        )
        assert np.allclose(
            co2_main_out.molecular_dipole_full, [0.0, 0.0, -0.0]
        )
        assert co2_main_out.total_molecular_dipole_moment == 0.0
        assert np.allclose(
            co2_main_out.molecular_quadrupole_qonly,
            [[-2.169, 0.0, 0.0], [0.0, 1.084, -0.0], [0.0, -0.0, 1.084]],
        )
        assert np.allclose(
            co2_main_out.molecular_quadrupole_q_dip,
            [[-3.107, 0.0, 0.0], [0.0, 1.553, -0.0], [0.0, -0.0, 1.553]],
        )
        assert np.allclose(
            co2_main_out.molecular_quadrupole_full,
            [[-4.360, 0.0, 0.0], [0.0, 2.180, -0.0], [0.0, -0.0, 2.180]],
        )
        assert co2_main_out.molecular_mass == 44.0095457
        assert co2_main_out.center_of_mass == [-0.0, 0.0, 0.0]
        assert co2_main_out.moments_of_inertia == [
            -0.3040259e-14,
            0.4185248e02,
            0.4185248e02,
        ]
        assert co2_main_out.rotational_constants_in_wavenumbers == [
            -0.5544801e16,
            0.4027869,
            0.4027869,
        ]
        assert co2_main_out.all_vibrational_frequencies == [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            600.31,
            600.31,
            1424.78,
            2593.07,
        ]
        assert co2_main_out.vibrational_frequencies == [
            600.31,
            600.31,
            1424.78,
            2593.07,
        ]
        assert co2_main_out.ir_intensities == [
            0.00,
            68.69,
            68.69,
            0.00,
        ]
        assert co2_main_out.raman_intensities == [
            0.00,
            0.00,
            0.00,
            0.00,
        ]
        # Thermodynamic Setup
        assert co2_main_out.num_vib_frequencies == 4
        assert co2_main_out.num_imaginary_frequencies == 0
        assert not co2_main_out.only_rot_calc
        assert co2_main_out.symmetry == "din"
        assert co2_main_out.rotational_symmetry_number == 2
        assert co2_main_out.scaling_factor == 1.0
        assert co2_main_out.rotor_cutoff == 50.0
        assert co2_main_out.imaginary_frequency_cutoff == -20.0
        # Thermodynamic Results
        assert co2_main_out.zero_point_energy == 0.011888572359
        assert co2_main_out.grrho_without_zpve == -0.020688920356
        assert co2_main_out.grrho_contribution == -0.008800347997
        assert co2_main_out.total_energy_without_gsasa_hb is None
        assert co2_main_out.total_energy == -10.308452289174
        assert co2_main_out.enthalpy == -10.292932643737
        assert co2_main_out.gibbs_free_energy == -10.317252637172
        assert co2_main_out.gradient_norm == 0.000000274582
        assert co2_main_out.fmo_gap == 8.448655866329
        assert np.isclose(
            co2_main_out.electronic_entropy,
            (51.1857 - 0.915 - 13.017 - 37.255)
            * 1e-3
            * kcal_per_mol_to_hartree,
        )
        assert np.isclose(
            co2_main_out.vibrational_entropy,
            0.915 * 1e-3 * kcal_per_mol_to_hartree,
        )
        assert np.isclose(
            co2_main_out.rotational_entropy,
            13.017 * 1e-3 * kcal_per_mol_to_hartree,
        )
        assert np.isclose(
            co2_main_out.translational_entropy,
            37.255 * 1e-3 * kcal_per_mol_to_hartree,
        )
        assert np.isclose(
            co2_main_out.entropy,
            51.1857 * 1e-3 * kcal_per_mol_to_hartree,
        )
        assert np.isclose(co2_main_out.entropy_times_temperature, 0.243200e-01)
        # Time Information
        assert np.isclose(co2_main_out.total_elapsed_walltime * 3600, 0.468)
        assert np.isclose(co2_main_out.total_core_hours * 3600, 0.088)
        assert np.isclose(co2_main_out.scf_wall_time * 3600, 0.094)
        assert np.isclose(co2_main_out.scf_cpu_time * 3600, 0.015)
        assert np.isclose(co2_main_out.optimizer_wall_time * 3600, 0.302)
        assert np.isclose(co2_main_out.optimizer_cpu_time * 3600, 0.047)
        assert np.isclose(co2_main_out.hessian_wall_time * 3600, 0.040)
        assert np.isclose(co2_main_out.hessian_cpu_time * 3600, 0.017)
        # Reduced masses (not exercised by any other test)
        assert co2_main_out.reduced_masses == [13.1, 13.1, 16.0, 13.1]
        # Solvent-off branches: co2 is a gas-phase (no solvent) calculation,
        # so every solvent-dependent property must short-circuit to None.
        assert not co2_main_out.solvent_on
        assert co2_main_out.solvent_model is None
        assert co2_main_out.solvent_id is None
        assert co2_main_out.dielectric_constant is None
        assert co2_main_out.free_energy_shift is None
        assert co2_main_out.solvent_temperature is None
        assert co2_main_out.density is None
        assert co2_main_out.solvent_mass is None
        assert co2_main_out.h_bond_correction is None
        assert co2_main_out.ion_screening is None
        assert co2_main_out.surface_tension is None
        # Closed-shell default multiplicity (unpaired_electrons == 0)
        assert co2_main_out.unpaired_electrons == 0
        assert co2_main_out.multiplicity == 1
        assert co2_main_out.spin == "restricted"
        # Special analysis properties (--vip/--vipea/--vomega/--vfukui) are
        # not requested for this run, so these must all be None.
        assert co2_main_out.vertical_ionization_potential is None
        assert co2_main_out.vertical_electron_affinity is None
        assert co2_main_out.global_electrophilicity_index is None
        assert co2_main_out.fukui_index is None
        assert co2_main_out.incomplete_optimized_geometry is False

    def test_main_out_acetaldehyde_hess(self, xtb_acetaldehyde_outfolder):
        """Test parsing main output from acetaldehyde hess-only calculation
        (no geometry optimization was performed)."""
        xtb_main_out_file = os.path.join(
            xtb_acetaldehyde_outfolder, "acetaldehyde_hess.out"
        )
        assert os.path.exists(xtb_main_out_file)
        main_out = XTBMainOut(xtb_main_out_file)
        assert main_out.normal_termination
        # No geometry optimization was requested for this hess-only job.
        assert main_out.geometry_optimization_converged is False
        assert main_out.optimized_structure_block is None
        assert main_out.molecular_mass is None
        assert main_out.center_of_mass is None
        assert main_out.moments_of_inertia is None
        assert main_out.rotational_constants_in_wavenumbers is None
        # No ANC optimizer timing since no optimization was run.
        assert main_out.optimizer_wall_time is None
        assert main_out.optimizer_cpu_time is None
        # Partition function table (VIB/ROT/INT/TR), not exercised elsewhere.
        partition_function = main_out.partition_function
        assert partition_function == {
            "vibrational": 2.23,
            "rotational": 0.117e05,
            "internal": 0.262e05,
            "translational": 0.283e27,
        }

    def test_main_out_p_benzyne_sp(self, xtb_p_benzyne_sp_outfolder):
        """Test parsing main output from a single-point calculation: no
        geometry optimization and no (numerical) Hessian were requested,
        so all opt-/hessian-dependent properties must be None/False."""
        xtb_main_out_file = os.path.join(
            xtb_p_benzyne_sp_outfolder, "p_benzyne_sp_alpb_toluene.out"
        )
        assert os.path.exists(xtb_main_out_file)
        main_out = XTBMainOut(xtb_main_out_file)
        assert main_out.normal_termination
        # No geometry optimization was performed.
        assert main_out.geometry_optimization_converged is False
        assert main_out.optimized_structure_block is None
        assert main_out.molecular_mass is None
        assert main_out.center_of_mass is None
        assert main_out.moments_of_inertia is None
        assert main_out.rotational_constants_in_wavenumbers is None
        assert main_out.optimizer_wall_time is None
        assert main_out.optimizer_cpu_time is None
        # No Hessian was computed for this sp job.
        assert main_out.numerical_hessian_block is None
        assert main_out.numfreq is False
        assert main_out.hessian_step_length is None
        assert main_out.scc_accuracy is None
        assert main_out.hessian_scale_factor is None
        assert main_out.rms_gradient is None
        assert main_out.hessian_wall_time is None
        assert main_out.hessian_cpu_time is None
        # No thermochemistry section (no frequency calculation performed).
        assert main_out.temperature_in_K is None
        assert main_out.entropy_no_temperature_in_SI is None
        assert main_out.electronic_entropy_no_temperature_in_SI is None
        assert main_out.electronic_entropy is None
        assert main_out.vibrational_entropy_no_temperature_in_SI is None
        assert main_out.vibrational_entropy is None
        assert main_out.rotational_entropy_no_temperature_in_SI is None
        assert main_out.rotational_entropy is None
        assert main_out.translational_entropy_no_temperature_in_SI is None
        assert main_out.translational_entropy is None
        assert main_out.entropy is None
        assert main_out.entropy_times_temperature is None
        # HOMO/LUMO/dispersion coefficients are still printed for sp jobs.
        assert main_out.homo_energy == -8.4001
        assert main_out.lumo_energy == -6.3910
        assert main_out.c6_coefficient == 1511.991497
        assert main_out.c8_coefficient == 38995.680442
        assert main_out.alpha_coefficient == 63.027243
        # No vibrational analysis at all for this sp-only job.
        assert main_out.reduced_masses is None
        assert main_out.ir_intensities is None
        assert main_out.raman_intensities is None
        # This sp job's SETUP block lacks the opt-related keys entirely.
        assert main_out.write_all_intermediate_geometries is None
        assert main_out.is_linear is None

    def test_main_out_p_benzyne_opt(self, xtb_p_benzyne_opt_outfolder):
        xtb_main_out_file = os.path.join(
            xtb_p_benzyne_opt_outfolder, "p_benzyne_opt_alpb_toluene.out"
        )
        assert os.path.exists(xtb_main_out_file)
        p_benzyne_opt_main_out = XTBMainOut(xtb_main_out_file)
        assert p_benzyne_opt_main_out.version == "6.7.1"
        assert p_benzyne_opt_main_out.normal_termination
        assert (
            p_benzyne_opt_main_out.route_string
            == "xtb p_benzyne.xyz --opt loose --alpb toluene --uhf 2 --grad --json"
        )
        # Open-shell (--uhf 2): unpaired electrons, multiplicity, and FMO levels
        assert p_benzyne_opt_main_out.net_charge == 0
        assert p_benzyne_opt_main_out.unpaired_electrons == 2
        assert p_benzyne_opt_main_out.multiplicity == 3
        assert p_benzyne_opt_main_out.spin == "unrestricted"
        assert p_benzyne_opt_main_out.homo_energy == -8.4
        assert p_benzyne_opt_main_out.lumo_energy == -6.391
        assert p_benzyne_opt_main_out.fmo_gap == 2.009059622602
        assert p_benzyne_opt_main_out.alpha_occ_eigenvalues is None
        assert p_benzyne_opt_main_out.beta_occ_eigenvalues is None
        assert p_benzyne_opt_main_out.solvent_on
        assert p_benzyne_opt_main_out.solvent_model == "ALPB"
        assert p_benzyne_opt_main_out.solvent_id == "toluene"
        assert p_benzyne_opt_main_out.dielectric_constant == 7.0
        assert p_benzyne_opt_main_out.free_energy_shift == 2.2081e-03
        assert p_benzyne_opt_main_out.solvent_temperature == 298.15
        assert p_benzyne_opt_main_out.density == 0.867
        assert p_benzyne_opt_main_out.solvent_mass == 78.11
        assert not p_benzyne_opt_main_out.h_bond_correction
        assert not p_benzyne_opt_main_out.ion_screening
        assert p_benzyne_opt_main_out.surface_tension == 1.000e-05
        assert p_benzyne_opt_main_out.solvation_energy_gsolv == -0.006398930400
        assert (
            p_benzyne_opt_main_out.electronic_solvation_energy_gelec
            == -0.000119108995
        )
        assert (
            p_benzyne_opt_main_out.surface_area_solvation_energy_gsasa
            == -0.008487964543
        )
        assert (
            p_benzyne_opt_main_out.hydrogen_bonding_solvation_energy_ghb
            == 0.000000000000
        )
        assert (
            p_benzyne_opt_main_out.empirical_shift_correction_gshift
            == 0.002208143139
        )


class TestXTBMainOutSyntheticEdgeCases:
    """Covers XTBMainOut branches that no real fixture output naturally
    exercises: abnormal terminations, missing route information, rarely
    requested analyses (--vip/--vipea/--vomega/--vfukui), and malformed
    SUMMARY/Numerical-Hessian blocks. Each synthetic file follows the
    real xTB output formatting conventions seen in the fixtures under
    tests/data/XTBTests/outputs/."""

    def test_error_terminated_output(self, tmp_path):
        """A run that dies mid-calculation: contents end in an [ERROR]
        line with no '* finished run' marker, and never reaches the
        HOMO/LUMO, dipole/quadrupole, or SUMMARY sections."""
        content = (
            "          program call               : xtb bad.xyz --opt\n"
            "\n"
            "          ...................................................\n"
            "          :                      SETUP                      :\n"
            "          :.................................................:\n"
            "          :  # basis functions                  12          :\n"
            "          :  net charge                          0          :\n"
            "          :  unpaired electrons                  0          :\n"
            "          ...................................................\n"
            "\n"
            "[ERROR] SCF did not converge\n"
        )
        out_file = tmp_path / "error_terminated.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))

        assert main_out.normal_termination is False
        assert main_out.homo_energy is None
        assert main_out.lumo_energy is None
        assert main_out.c6_coefficient is None
        assert main_out.c8_coefficient is None
        assert main_out.alpha_coefficient is None
        assert main_out.molecular_dipole_lines is None
        assert main_out.molecular_dipole_qonly is None
        assert main_out.molecular_dipole_full is None
        assert main_out.total_molecular_dipole_moment is None
        assert main_out.molecular_quadrupole_lines is None
        assert main_out.molecular_quadrupole_qonly is None
        assert main_out.molecular_quadrupole_q_dip is None
        assert main_out.molecular_quadrupole_full is None
        # No SUMMARY block at all in a run that errored out this early.
        assert main_out.get_all_summary_blocks() is None
        # No timing information was ever printed.
        assert main_out.total_elapsed_walltime is None
        assert main_out.total_core_hours is None
        assert main_out.scf_wall_time is None
        assert main_out.scf_cpu_time is None
        assert main_out.optimizer_wall_time is None
        assert main_out.optimizer_cpu_time is None
        assert main_out.hessian_wall_time is None
        assert main_out.hessian_cpu_time is None
        assert main_out.partition_function is None

    def test_no_termination_marker(self, tmp_path):
        """A run that was killed/truncated: no '[ERROR]' line and no
        '* finished run' marker anywhere -- normal_termination must
        fall back to False after scanning the whole file."""
        content = (
            "          program call               : xtb water.xyz --opt\n"
            "\n"
            "          ...................................................\n"
            "          :                      SETUP                      :\n"
            "          :.................................................:\n"
            "          :  net charge                          0          :\n"
            "          ...................................................\n"
        )
        out_file = tmp_path / "no_termination_marker.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))
        assert main_out.normal_termination is False

    def test_no_program_call_line(self, tmp_path):
        """A file lacking the 'program call' line entirely -- route
        parsing must fall back to None even though the run otherwise
        terminated normally."""
        content = (
            "          ...................................................\n"
            "          :                      SETUP                      :\n"
            "          :.................................................:\n"
            "          :  net charge                          0          :\n"
            "          ...................................................\n"
            "\n"
            "* finished run on 2026/01/01 at 00:00:00.000\n"
        )
        out_file = tmp_path / "no_program_call.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))
        assert main_out.route_string is None
        assert main_out.normal_termination is True

    def test_vip_ea_omega_and_fukui_analysis(self, tmp_path):
        """Synthetic --vip/--vipea/--vomega/--vfukui analysis printout,
        a feature none of the real fixtures were run with."""
        content = (
            "          program call               : xtb mol.xyz --vfukui --vomega\n"
            "\n"
            "delta SCC IP (eV): 12.345\n"
            "delta SCC EA (eV): -1.234\n"
            "Global electrophilicity index (eV): 3.456\n"
            "\n"
            "Fukui functions:\n"
            "    #        f(+)     f(-)     f(0)\n"
            "    1O      -0.113   -0.098   -0.106\n"
            "    2H       0.056    0.049    0.053\n"
            "------------------------------------------------------\n"
        )
        out_file = tmp_path / "vip_ea_fukui.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))

        assert main_out.vertical_ionization_potential == 12.345
        assert main_out.vertical_electron_affinity == -1.234
        assert main_out.global_electrophilicity_index == 3.456
        assert main_out.fukui_index == [
            "#        f(+)     f(-)     f(0)",
            "1O      -0.113   -0.098   -0.106",
            "2H       0.056    0.049    0.053",
        ]

    def test_partition_function_table_absent(self, tmp_path):
        """A file with no VIB/ROT/INT/TR partition-function table at
        all -- partition_function must return None."""
        content = (
            "          program call               : xtb mol.xyz\n"
            "\n"
            "          ...................................................\n"
            "          :                      SETUP                      :\n"
            "          :.................................................:\n"
            "          :  net charge                          0          :\n"
            "          ...................................................\n"
            "\n"
            "* finished run on 2026/01/01 at 00:00:00.000\n"
        )
        out_file = tmp_path / "no_partition_table.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))
        assert main_out.partition_function is None

    def test_partition_function_heading_without_vib_row(self, tmp_path):
        """The 'partition function'+'entropy' heading is present, but
        the row table breaks (hits the 'T/K'+'H(T)' marker) before any
        VIB/ROT/INT/TR row is found -- temperature_in_K and the
        entropy properties derived from the same scan must all fall
        back to None."""
        content = (
            "   temp. (K)  partition function   enthalpy   heat capacity  entropy\n"
            "       T/K    H(0)-H(T)+PV         H(T)/Eh          T*S/Eh         G(T)/Eh\n"
        )
        out_file = tmp_path / "partition_heading_no_rows.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))

        assert main_out.temperature_in_K is None
        assert main_out.entropy_no_temperature_in_SI is None
        assert main_out.electronic_entropy_no_temperature_in_SI is None
        assert main_out.electronic_entropy is None
        assert main_out.vibrational_entropy_no_temperature_in_SI is None
        assert main_out.vibrational_entropy is None
        assert main_out.rotational_entropy_no_temperature_in_SI is None
        assert main_out.rotational_entropy is None
        assert main_out.translational_entropy_no_temperature_in_SI is None
        assert main_out.translational_entropy is None
        assert main_out.entropy is None
        assert main_out.entropy_times_temperature is None

    def test_summary_block_with_no_trailing_blank_line(self, tmp_path):
        """A SUMMARY block that is the very last thing in the file (no
        trailing blank line) -- the inner block-collection loop must
        exhaust the file's contents rather than break on a blank line."""
        content = (
            "         ::                     SUMMARY                     ::\n"
            "         :::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
            "         :: SCC energy               -1.000000000000 Eh    ::\n"
            "         :::::::::::::::::::::::::::::::::::::::::::::::::::::"
        )
        out_file = tmp_path / "summary_no_trailing_blank.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))

        blocks = main_out.get_all_summary_blocks()
        assert blocks == [
            [":: SCC energy               -1.000000000000 Eh    ::"]
        ]
        assert main_out.scc_energy == -1.0

    def test_summary_block_immediately_empty(self, tmp_path):
        """A SUMMARY heading immediately followed by a blank line (no
        content) -- get_all_summary_blocks() must return a non-None
        list whose last block is an empty (falsy) list, and any
        property built on top of it must gracefully return None
        instead of indexing into that empty block."""
        content = (
            "         ::                     SUMMARY                     ::\n"
            "         :::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
            "\n"
        )
        out_file = tmp_path / "summary_empty_block.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))

        blocks = main_out.get_all_summary_blocks()
        assert blocks == [[]]
        assert main_out.scc_energy is None
        assert main_out.total_charge is None

    def test_numerical_hessian_block_runs_to_eof_and_step_length_not_first(
        self, tmp_path
    ):
        """The Numerical Hessian block is the last thing in the file
        (no trailing blank line), so its collection loop must exhaust
        the file rather than break. It also places an unrelated line
        before 'step length', so hessian_step_length's own search loop
        must iterate past a non-matching line before it finds a match."""
        content = (
            "          |                Numerical Hessian                |\n"
            "           ------------------------------------------------- \n"
            "some other diagnostic line\n"
            "step length          :   0.00500\n"
            "SCC accuracy         :   0.30000\n"
            "Hessian scale factor :   1.00000"
        )
        out_file = tmp_path / "hessian_block_no_trailing_blank.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))

        assert main_out.numerical_hessian_block == [
            "some other diagnostic line",
            "step length          :   0.00500",
            "SCC accuracy         :   0.30000",
            "Hessian scale factor :   1.00000",
        ]
        assert main_out.hessian_step_length == 0.00500

    def test_failed_convergence_and_incomplete_geometry(self, tmp_path):
        """Both the explicit 'FAILED TO CONVERGE' marker and the
        'INCOMPLETELY OPTIMIZED GEOMETRY' warning, neither of which
        appears in any successfully-converged real fixture."""
        content = (
            "some earlier output\n"
            "FAILED TO CONVERGE GEOMETRY OPTIMIZATION\n"
            "INCOMPLETELY OPTIMIZED GEOMETRY\n"
        )
        out_file = tmp_path / "failed_convergence.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))

        assert main_out.geometry_optimization_converged is False
        assert main_out.incomplete_optimized_geometry is True

    def test_converged_but_missing_derived_geometry_fields(self, tmp_path):
        """A run that reports GEOMETRY OPTIMIZATION CONVERGED and a
        'final structure:' block, but never prints the 'Bond Distances'
        marker, nor the molecular mass/center-of-mass/moments-of-inertia
        /rotational-constants lines that a real ohess fixture always
        includes -- every one of those derived properties must fall
        back to None (or the fully-collected block, for the structure
        itself) instead of crashing."""
        content = (
            "GEOMETRY OPTIMIZATION CONVERGED\n"
            "\n"
            "final structure:\n"
            "\n"
            "3\n"
            "xtb: 6.7.1 (edcfbbe)\n"
            "O   0.0 0.0 0.0\n"
            "O   1.0 0.0 0.0\n"
        )
        out_file = tmp_path / "converged_minimal.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))

        assert main_out.geometry_optimization_converged is True
        assert main_out.optimized_structure_block == [
            "3",
            "xtb: 6.7.1 (edcfbbe)",
            "O   0.0 0.0 0.0",
            "O   1.0 0.0 0.0",
        ]
        assert main_out.molecular_mass is None
        assert main_out.center_of_mass is None
        assert main_out.moments_of_inertia is None
        assert main_out.rotational_constants_in_wavenumbers is None

    def test_only_rot_calc_key_never_matches_real_xtb_output(
        self, xtb_co2_outfolder
    ):
        """`only_rot_calc` searches the SETUP block for the literal key
        'only rotational calc.', but real xTB output prints 'only rotor
        calc.' instead (see the Hessian SETUP block in co2_ohess.out).
        The two never match, so only_rot_calc always returns None for
        any real xTB output -- see BUGS_FOUND.md."""
        xtb_main_out_file = os.path.join(xtb_co2_outfolder, "co2_ohess.out")
        main_out = XTBMainOut(xtb_main_out_file)
        assert any("only rotor calc." in line for line in main_out.contents)
        assert not any(
            "only rotational calc." in line for line in main_out.contents
        )
        assert main_out.only_rot_calc is None

    def test_setup_information_scan_reaches_true_eof(self, tmp_path):
        """`_get_setup_information`'s inner block-scan must exhaust the
        file's contents (rather than break on a blank line) when the
        SETUP block is the very last thing in the file and the
        requested keyword isn't present in it."""
        content = (
            "          :                      SETUP                      :\n"
            "          :.................................................:\n"
            "          :  net charge                          0          :\n"
            "          :.................................................:"
        )
        out_file = tmp_path / "setup_truncated_no_blank.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))
        assert main_out.num_basis_functions is None

    def test_partition_function_heading_with_no_matching_rows_at_all(
        self, tmp_path
    ):
        """The 'partition function'+'entropy' heading is found, but the
        remainder of the file has neither a VIB/ROT/INT/TR row nor a
        'T/K'+'H(T)' closing marker -- the inner row-scan must exhaust
        the file's contents naturally instead of ever breaking."""
        content = (
            "   temp. (K)  partition function   enthalpy   heat capacity  entropy\n"
            "some unrelated trailing line with no special markers\n"
        )
        out_file = tmp_path / "partition_heading_exhausts.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))

        assert main_out.temperature_in_K is None
        assert main_out.vibrational_entropy_no_temperature_in_SI is None
        assert main_out.rotational_entropy_no_temperature_in_SI is None
        assert main_out.translational_entropy_no_temperature_in_SI is None
        assert main_out.entropy_no_temperature_in_SI is None

    def test_fukui_block_runs_to_eof_without_closing_marker(self, tmp_path):
        """A Fukui functions block that is the last thing in the file,
        with no closing '------' marker before EOF."""
        content = (
            "Fukui functions:\n"
            "    #        f(+)     f(-)     f(0)\n"
            "    1O      -0.113   -0.098   -0.106\n"
        )
        out_file = tmp_path / "fukui_no_closing.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))
        assert main_out.fukui_index == [
            "#        f(+)     f(-)     f(0)",
            "1O      -0.113   -0.098   -0.106",
        ]

    def test_numerical_hessian_block_present_but_all_keys_missing(
        self, tmp_path
    ):
        """The Numerical Hessian block exists (non-empty) but contains
        none of step length / SCC accuracy / Hessian scale factor / RMS
        gradient -- each of their own search loops must exhaust the
        block without ever matching."""
        content = (
            "          |                Numerical Hessian                |\n"
            "           ------------------------------------------------- \n"
            "some unrelated diagnostic line\n"
            "\n"
        )
        out_file = tmp_path / "hessian_block_no_keys.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))

        assert main_out.numerical_hessian_block == [
            "some unrelated diagnostic line"
        ]
        assert main_out.hessian_step_length is None
        assert main_out.scc_accuracy is None
        assert main_out.hessian_scale_factor is None
        assert main_out.rms_gradient is None

    def test_solvent_on_but_detail_lines_missing(self, tmp_path):
        """GBSA solvation is flagged on in the SETUP block, but none of
        the individual solvent detail lines (Dielectric constant, Free
        energy shift, Temperature, Density, Solvent mass, H-bond
        correction, Ion screening, Surface tension) are present -- each
        property's own content scan must exhaust without matching."""
        content = (
            "          :                      SETUP                      :\n"
            "          :.................................................:\n"
            "          :  GBSA solvation                   true          :\n"
            "          :.................................................:\n"
        )
        out_file = tmp_path / "solvent_on_no_details.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))

        assert main_out.solvent_on is True
        assert main_out.dielectric_constant is None
        assert main_out.free_energy_shift is None
        assert main_out.solvent_temperature is None
        assert main_out.density is None
        assert main_out.solvent_mass is None
        assert main_out.h_bond_correction is None
        assert main_out.ion_screening is None
        assert main_out.surface_tension is None

    def test_dipole_full_and_total_none_when_full_line_truncated(
        self, tmp_path
    ):
        """The dipole block is truncated right after the 'q only:' row
        (no 'full:' row follows) -- molecular_dipole_full and
        total_molecular_dipole_moment's own search loops must exhaust
        the (non-empty) dipole_lines without ever matching, while
        molecular_dipole_qonly still succeeds normally."""
        content = (
            "molecular dipole:\n"
            "                 x           y           z       tot (Debye)\n"
            " q only:        0.100       0.200       0.300\n"
        )
        out_file = tmp_path / "dipole_truncated.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))

        assert main_out.molecular_dipole_lines == [
            "q only:        0.100       0.200       0.300"
        ]
        assert np.allclose(
            main_out.molecular_dipole_qonly, [0.100, 0.200, 0.300]
        )
        assert main_out.molecular_dipole_full is None
        assert main_out.total_molecular_dipole_moment is None

    def test_quadrupole_all_none_when_block_truncated_to_labels_only(
        self, tmp_path
    ):
        """The quadrupole block is truncated right after the heading
        (only the column-labels row survives before EOF) -- none of
        qonly/q+dip/full's search loops can find their marker, so each
        must exhaust the non-empty quadrupole_lines without matching."""
        content = (
            "molecular quadrupole (traceless):\n"
            "                xx          xy          yy          xz          yz          zz\n"
        )
        out_file = tmp_path / "quadrupole_truncated.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))

        assert main_out.molecular_quadrupole_lines == [
            "xx          xy          yy          xz          yz          zz"
        ]
        assert main_out.molecular_quadrupole_qonly is None
        assert main_out.molecular_quadrupole_q_dip is None
        assert main_out.molecular_quadrupole_full is None

    def test_all_vibrational_frequencies_block_runs_to_eof(self, tmp_path):
        """The vibrational-frequencies row is the last thing in the
        file, with no 'reduced masses (amu)' marker ever following --
        the frequency-collection loop must exhaust the file naturally."""
        content = (
            "           |               Frequency Printout                |\n"
            "vibrational frequencies (cm⁻¹)\n"
            "eigval :        0.00   100.00\n"
        )
        out_file = tmp_path / "vib_freq_no_closing.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))
        assert main_out.all_vibrational_frequencies == [0.00, 100.00]
        assert main_out.vibrational_frequencies == [100.00]

    def test_vibrational_frequencies_none_source_returns_empty_list(
        self, tmp_path
    ):
        """When all_vibrational_frequencies is None (no Frequency
        Printout section at all), vibrational_frequencies must fall
        back to an empty list rather than None."""
        content = "no frequency section here at all\n"
        out_file = tmp_path / "no_freq_section.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))
        assert main_out.all_vibrational_frequencies is None
        assert main_out.vibrational_frequencies == []

    def test_reduced_masses_block_runs_to_eof(self, tmp_path):
        """The reduced-masses row is the last thing in the file, with
        no 'IR intensities (km·mol⁻¹)' marker ever following."""
        content = (
            "          :                      SETUP                      :\n"
            "          :.................................................:\n"
            "          :  # frequencies                       2          :\n"
            "          :.................................................:\n"
            "\n"
            "reduced masses (amu)\n"
            "   1: 12.345   2: 67.890\n"
        )
        out_file = tmp_path / "reduced_masses_no_closing.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))
        assert main_out.reduced_masses == [12.345, 67.890]

    def test_ir_intensities_block_runs_to_eof(self, tmp_path):
        """The IR-intensities row is the last thing in the file, with
        no 'Raman intensities' marker ever following."""
        content = (
            "          :                      SETUP                      :\n"
            "          :.................................................:\n"
            "          :  # frequencies                       2          :\n"
            "          :.................................................:\n"
            "\n"
            "IR intensities (km*mol-1)\n"
            "   1: 10.000   2: 20.000\n"
        )
        out_file = tmp_path / "ir_intensities_no_closing.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))
        assert main_out.ir_intensities == [10.000, 20.000]

    def test_raman_intensities_block_runs_to_eof_with_bad_token(
        self, tmp_path
    ):
        """The Raman-intensities row is the last thing in the file (no
        closing 'output can be read by thermo' marker), and it also
        contains a non-numeric token that must be silently skipped."""
        content = (
            "          :                      SETUP                      :\n"
            "          :.................................................:\n"
            "          :  # frequencies                       2          :\n"
            "          :.................................................:\n"
            "\n"
            "Raman intensities (A**4/amu)\n"
            "   1: notanumber   2: 30.000\n"
        )
        out_file = tmp_path / "raman_intensities_no_closing.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))
        assert main_out.raman_intensities == [30.000]

    def test_thermodynamics_block_runs_to_eof(self, tmp_path):
        """The THERMODYNAMIC block is the last thing in the file, with
        no blank line ever terminating it."""
        content = (
            "                   :::::::::::::::::::::::::::::::::::::::::::::::::\n"
            "                   ::                THERMODYNAMIC                ::\n"
            "                   :::::::::::::::::::::::::::::::::::::::::::::::::\n"
            "                   :: zero point energy         0.011888572359 Eh   ::"
        )
        out_file = tmp_path / "thermo_block_no_blank.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))
        assert main_out.zero_point_energy == 0.011888572359
        # Block is present (non-empty) but doesn't contain this keyword,
        # so the search loop must exhaust it without ever returning.
        assert main_out.grrho_without_zpve is None

    def test_thermodynamics_block_absent(self, tmp_path):
        """No THERMODYNAMIC heading at all -- the block getter must
        return None, and any property derived from it must gracefully
        fall back to None as well."""
        content = "nothing relevant in this file\n"
        out_file = tmp_path / "no_thermo_block.out"
        out_file.write_text(content)
        main_out = XTBMainOut(str(out_file))
        assert main_out._get_thermodynamics_block() is None
        assert main_out.zero_point_energy is None
        assert main_out.total_energy is None

    def test_always_none_trivial_properties(self, xtb_co2_outfolder):
        """A handful of XTBMainOut properties are unconditionally None
        (placeholders kept for interface parity with other QM-program
        parsers) -- exercise them all via any real fixture instance."""
        xtb_main_out_file = os.path.join(xtb_co2_outfolder, "co2_ohess.out")
        main_out = XTBMainOut(xtb_main_out_file)
        assert main_out.alpha_occ_eigenvalues is None
        assert main_out.beta_occ_eigenvalues is None
        assert main_out.alpha_virtual_eigenvalues is None
        assert main_out.beta_virtual_eigenvalues is None
        assert main_out.rotational_temperatures is None
        assert main_out.rotational_constants_in_Hz is None
        assert main_out.thermal_vibration_correction is None
        assert main_out.thermal_rotation_correction is None
        assert main_out.thermal_translation_correction is None
        assert main_out.thermal_energy_correction is None
        assert main_out.thermal_enthalpy_correction is None
        assert main_out.thermal_gibbs_free_energy_correction is None
        assert main_out.internal_energy is None
        assert main_out.pressure_in_atm is None


class TestXTBChargesFile:
    """Tests for XTBChargesFile class."""

    def test_charges_co2(self, xtb_co2_outfolder):
        """Test parsing charges from CO2 ohess calculation."""
        charges_file = os.path.join(xtb_co2_outfolder, "charges")
        assert os.path.exists(charges_file)
        co2_charges = XTBChargesFile(charges_file)
        assert co2_charges.partial_charges == [
            -0.23213972,
            -0.23213972,
            0.46427944,
        ]
        assert np.isclose(co2_charges.total_charge, 0, atol=1e-8)

    def test_charges_cyclopentadienyl_anion(
        self, xtb_cyclopentadienyl_anion_outfolder
    ):
        """Test parsing charges from cyclopentadienyl anion opt calculation."""
        charges_file = os.path.join(
            xtb_cyclopentadienyl_anion_outfolder, "charges"
        )
        assert os.path.exists(charges_file)
        cyclopentadienyl_anion_charges = XTBChargesFile(charges_file)
        assert cyclopentadienyl_anion_charges.partial_charges == [
            -0.10690174,
            -0.10686507,
            -0.10690237,
            -0.10686693,
            -0.10687190,
            -0.09313128,
            -0.09311500,
            -0.09310077,
            -0.09314340,
            -0.09310154,
        ]
        assert np.isclose(
            cyclopentadienyl_anion_charges.total_charge, -1, atol=1e-8
        )

    def test_charges_p_benzyne_sp(self, xtb_p_benzyne_sp_outfolder):
        """Test parsing charges from p-benzyne sp calculation."""
        charges_file = os.path.join(xtb_p_benzyne_sp_outfolder, "charges")
        assert os.path.exists(charges_file)
        p_benzyne_sp_charges = XTBChargesFile(charges_file)
        assert p_benzyne_sp_charges.partial_charges == [
            -0.06012030,
            -0.01778337,
            -0.01766299,
            -0.06009510,
            -0.01770971,
            -0.01769064,
            0.04789461,
            0.04763761,
            0.04779053,
            0.04773936,
        ]
        assert np.isclose(p_benzyne_sp_charges.total_charge, 0, atol=1e-8)

    def test_partial_charges_none_for_blank_and_invalid_lines(self, tmp_path):
        """A charges file that has no valid numeric lines at all -- just
        a blank line (skipped via `continue`) and a non-numeric line
        (skipped via the `except ValueError: continue` branch) -- must
        yield partial_charges=None and, in turn, total_charge=None."""
        charges_file = tmp_path / "charges"
        charges_file.write_text("\nnot_a_number\n")
        charges = XTBChargesFile(str(charges_file))
        assert charges.partial_charges is None
        assert charges.total_charge is None


class TestXTBEnergyFile:
    """Tests for XTBEnergyFile class."""

    def test_energy_co2(self, xtb_co2_outfolder):
        """Test parsing energy from CO2 ohess calculation."""
        energy_file = os.path.join(xtb_co2_outfolder, "energy")
        assert os.path.exists(energy_file)
        co2_energy = XTBEnergyFile(energy_file)
        assert co2_energy.last_energy == -10.30845228917

    def test_energy_water(self, xtb_water_outfolder):
        """Test parsing energy from water ohess calculation."""
        energy_file = os.path.join(xtb_water_outfolder, "energy")
        assert os.path.exists(energy_file)
        water_energy = XTBEnergyFile(energy_file)
        assert water_energy.last_energy == -5.07054444346

    def test_energy_p_benzyne_opt(self, xtb_p_benzyne_opt_outfolder):
        """Test parsing energy from p-benzyne opt calculation."""
        energy_file = os.path.join(xtb_p_benzyne_opt_outfolder, "energy")
        assert os.path.exists(energy_file)
        p_benzyne_opt_energy = XTBEnergyFile(energy_file)
        assert p_benzyne_opt_energy.last_energy == -14.66185695901

    def test_last_energy_none_when_no_data_line(self, tmp_path):
        """A file consisting only of $-marker and blank lines (both
        skipped) with no actual data line -- last_energy must fall
        back to None instead of raising."""
        energy_file = tmp_path / "energy"
        energy_file.write_text("$energy\n\n$end\n")
        energy = XTBEnergyFile(str(energy_file))
        assert energy.last_energy is None


class TestXTBEngradFile:
    """Tests for XTBEngradFile class."""

    def test_engrad_co2(self, xtb_co2_outfolder):
        """Test parsing energy gradient from CO2 ohess calculation."""
        engrad_file = os.path.join(xtb_co2_outfolder, "co2.engrad")
        assert os.path.exists(engrad_file)
        co2_engrad = XTBEngradFile(engrad_file)
        assert co2_engrad.num_atoms == 3
        assert co2_engrad.total_energy == -10.308452289174
        assert np.allclose(
            co2_engrad.forces[0][0],
            [0.000000194263, -0.000000000000, -0.000000000000],
        )
        assert np.allclose(
            co2_engrad.forces[0][1],
            [-0.000000194263, 0.000000000000, -0.000000000000],
        )
        assert np.allclose(
            co2_engrad.forces[0][2],
            [0.000000000000, 0.000000000000, 0.000000000000],
        )

    def test_engrad_p_benzyne_opt(self, xtb_p_benzyne_opt_outfolder):
        """Test parsing energy gradient from p-benzyne opt calculation."""
        engrad_file = os.path.join(
            xtb_p_benzyne_opt_outfolder, "p_benzyne.engrad"
        )
        assert os.path.exists(engrad_file)
        p_benzyne_opt_engrad = XTBEngradFile(engrad_file)
        assert p_benzyne_opt_engrad.num_atoms == 10
        assert p_benzyne_opt_engrad.total_energy == -14.661856959008
        assert np.allclose(
            p_benzyne_opt_engrad.forces[0][0],
            [-0.000151885154, 0.000008080973, -0.000083466044],
        )
        assert np.allclose(
            p_benzyne_opt_engrad.forces[0][1],
            [0.000550798511, -0.000200234055, -0.000078750880],
        )
        assert np.allclose(
            p_benzyne_opt_engrad.forces[0][2],
            [-0.000344033432, -0.000092282346, 0.000528937526],
        )
        assert np.allclose(
            p_benzyne_opt_engrad.forces[0][3],
            [0.000265720803, -0.000032115473, -0.000048689390],
        )
        assert np.allclose(
            p_benzyne_opt_engrad.forces[0][4],
            [-0.000041353408, 0.000168583076, 0.000075947816],
        )
        assert np.allclose(
            p_benzyne_opt_engrad.forces[0][5],
            [0.000006914570, 0.000063463873, -0.000186368586],
        )
        assert np.allclose(
            p_benzyne_opt_engrad.forces[0][6],
            [-0.000446235650, 0.000190941706, 0.000274234357],
        )
        assert np.allclose(
            p_benzyne_opt_engrad.forces[0][7],
            [-0.000021789087, -0.000005645447, -0.000313642685],
        )
        assert np.allclose(
            p_benzyne_opt_engrad.forces[0][8],
            [0.000151189558, -0.000003081686, -0.000132662651],
        )
        assert np.allclose(
            p_benzyne_opt_engrad.forces[0][9],
            [0.000030673290, -0.000097710620, -0.000035539462],
        )

    def test_no_headers_at_all(self, tmp_path):
        """A file with none of the '# Number of atoms' / 'current
        total energy' / 'current gradient' markers at all -- num_atoms,
        total_energy, and (via num_atoms being None) forces must all
        gracefully return None."""
        engrad_file = tmp_path / "empty.engrad"
        engrad_file.write_text("# nothing relevant here\n# just text\n")
        engrad = XTBEngradFile(str(engrad_file))
        assert engrad.num_atoms is None
        assert engrad.total_energy is None
        assert engrad.forces is None

    def test_headers_present_but_unparseable(self, tmp_path):
        """Both headers are present, but none of the lines following
        them parse as numbers -- the inner per-heading scan loop must
        exhaust its 3-line window without ever returning, so the outer
        loop keeps searching and eventually falls back to None."""
        engrad_file = tmp_path / "malformed.engrad"
        engrad_file.write_text(
            "#\n"
            "# Number of atoms\n"
            "#\n"
            "not_a_number\n"
            "also_not\n"
            "still_not\n"
            "#\n"
            "# The current total energy in Eh\n"
            "#\n"
            "not_a_number\n"
            "also_not\n"
            "still_not\n"
        )
        engrad = XTBEngradFile(str(engrad_file))
        assert engrad.num_atoms is None
        assert engrad.total_energy is None

    def test_gradient_component_count_mismatch(self, tmp_path):
        """num_atoms and total_energy parse fine, but the gradient
        section has too few numeric values for 3*num_atoms -- the
        gradient-count validation must fail, the search loop must
        continue past it (finding no other match) and exhaust to
        None, and forces must fall back to None as well."""
        engrad_file = tmp_path / "mismatch.engrad"
        engrad_file.write_text(
            "#\n"
            "# Number of atoms\n"
            "#\n"
            "         2\n"
            "#\n"
            "# The current total energy in Eh\n"
            "#\n"
            "     -1.234567890123\n"
            "#\n"
            "# The current gradient in Eh/bohr\n"
            "#\n"
            "       0.000000000100\n"
            "      -0.000000000050\n"
        )
        engrad = XTBEngradFile(str(engrad_file))
        assert engrad.num_atoms == 2
        assert engrad.total_energy == -1.234567890123
        assert engrad.forces is None


class TestXTBGradientFile:
    """Tests for XTBGradientFile class."""

    def test_gradient_co2(self, xtb_water_outfolder):
        gradient_file = os.path.join(xtb_water_outfolder, "gradient")
        grad = XTBGradientFile(gradient_file)
        assert np.isclose(grad.energy, -5.07054444346)
        assert grad.gradients[0].shape == (3, 3)
        assert np.allclose(
            grad.gradients[-1][0],
            [1.2982149851656e-10, -3.3293838441823e-18, 5.7137032470948e-05],
        )
        assert np.allclose(
            grad.gradients[-1][1],
            [-1.9065815509550e-05, -1.9010909790127e-17, -2.8568564065665e-05],
        )
        assert np.allclose(
            grad.gradients[-1][2],
            [1.9065685688053e-05, 2.2340293634309e-17, -2.8568468405277e-05],
        )
        assert np.allclose(grad.forces[-1], -grad.gradients[-1])

    def test_gradients_edge_cases_combined(self, tmp_path):
        """A single synthetic $grad file exercising several branches a
        well-formed real gradient file never takes:
        - content before the first '$grad' marker (skipped)
        - a junk line inside the block, before 'SCF energy', that
          matches neither the coordinate- nor gradient-line shape
          (loop continues without appending)
        - a malformed coordinate line (non-numeric tokens -- caught
          and skipped)
        - a malformed gradient-value line (non-numeric token -- caught
          and skipped)
        - no closing '$end' marker at all (the file just ends), so the
          block-collection loop must exhaust the file rather than
          break.
        """
        content = (
            "# leading comment before any $grad block\n"
            "$grad\n"
            " unrelated junk text here\n"
            "  cycle =      1    SCF energy =    -1.000000000000   "
            "|dE/dxyz| =  0.000100\n"
            "   BAD   BAD   BAD      O\n"
            "   -0.00000250190431     -0.00000125099553     "
            "-0.71677514611431      O\n"
            "   1.2982149851656E-10  -3.3293838441823E-18   BADVALUE\n"
            "   -1.9065815509550E-05  -1.9010909790127E-17  "
            "-2.8568564065665E-05\n"
        )
        gradient_file = tmp_path / "gradient"
        gradient_file.write_text(content)
        grad = XTBGradientFile(str(gradient_file))

        assert grad.gradients is not None
        assert grad.gradients[0].shape == (1, 3)
        assert np.allclose(
            grad.gradients[0][0],
            [-1.9065815509550e-05, -1.9010909790127e-17, -2.8568564065665e-05],
        )
        assert np.allclose(grad.forces[0], -grad.gradients[0])

    def test_gradients_none_when_no_grad_block(self, tmp_path):
        """No '$grad' marker anywhere -- gradients (and, in turn,
        forces) must fall back to None."""
        gradient_file = tmp_path / "gradient"
        gradient_file.write_text("nothing relevant here\n")
        grad = XTBGradientFile(str(gradient_file))
        assert grad.gradients is None
        assert grad.forces is None


class TestXTBHessianFile:
    """Tests for XTBHessianFile class."""

    def test_hessian_co2(self, xtb_co2_outfolder):
        hessian_file = os.path.join(xtb_co2_outfolder, "hessian")
        hess = XTBHessianFile(hessian_file)
        assert hess.hessian.shape == (9, 9)
        assert np.allclose(
            hess.hessian[0],
            np.array(
                [
                    1.1701657400,
                    -0.0000059805,
                    -0.0000039870,
                    -0.0589478867,
                    0.0000004551,
                    0.0000003034,
                    -1.1111279332,
                    0.0000069650,
                    0.0000046434,
                ]
            ),
        )
        assert np.allclose(
            hess.hessian[-1],
            np.array(
                [
                    0.0000046434,
                    -0.0000000000,
                    -0.0595361130,
                    0.0000046434,
                    -0.0000000000,
                    -0.0595361130,
                    -0.0000092863,
                    0.0000000001,
                    0.1191272426,
                ]
            ),
        )

    def test_hessian_none_for_empty_file(self, tmp_path):
        """A hessian file containing only the '$hessian' marker (and no
        numeric values) -- hessian must be None rather than an empty
        array or a crash."""
        hessian_file = tmp_path / "hessian"
        hessian_file.write_text("$hessian\n")
        hess = XTBHessianFile(str(hessian_file))
        assert hess.hessian is None

    def test_hessian_raises_for_non_square_value_count(self, tmp_path):
        """A hessian file whose flat value count isn't a perfect square
        can't be reshaped into an (N, N) matrix -- must raise
        ValueError rather than silently truncating/reshaping wrongly."""
        hessian_file = tmp_path / "hessian"
        hessian_file.write_text("$hessian\n1.0 2.0 3.0\n")
        hess = XTBHessianFile(str(hessian_file))
        with pytest.raises(ValueError, match="not a square matrix"):
            hess.hessian


class TestXTBVibSpectrumFile:
    """Tests for XTBVibSpectrumFile class."""

    def test_vibspectrum_co2(self, xtb_acetaldehyde_outfolder):
        vib_file = os.path.join(xtb_acetaldehyde_outfolder, "vibspectrum")
        vib = XTBVibSpectrumFile(vib_file)
        assert vib.vibrational_frequencies == [
            151.34,
            501.81,
            769.05,
            947.27,
            1045.68,
            1107.27,
            1355.34,
            1389.38,
            1446.60,
            1447.86,
            1798.58,
            2748.94,
            3018.34,
            3026.55,
            3059.76,
        ]
        assert vib.ir_intensities == [
            0.04462,
            11.45411,
            2.91256,
            15.16489,
            16.51912,
            60.12317,
            12.40125,
            46.09384,
            16.02810,
            10.90711,
            300.93898,
            142.63122,
            6.53630,
            3.41138,
            14.20614,
        ]
        assert vib.vibrational_mode_symmetries == [
            "a",
            "a",
            "a",
            "a",
            "a",
            "a",
            "a",
            "a",
            "a",
            "a",
            "a",
            "a",
            "a",
            "a",
            "a",
        ]
        # frequencies (all_modes=True) includes the zero-frequency
        # translational/rotational modes that vibrational_frequencies
        # filters out.
        assert vib.frequencies == [-0.00] * 3 + [0.00] * 3 + [
            151.34,
            501.81,
            769.05,
            947.27,
            1045.68,
            1107.27,
            1355.34,
            1389.38,
            1446.60,
            1447.86,
            1798.58,
            2748.94,
            3018.34,
            3026.55,
            3059.76,
        ]

    def test_frequencies_data_none_and_malformed_row_skipped(self, tmp_path):
        """A vibspectrum file with a comment/marker-only body (so
        _frequencies_data collects zero rows and returns None) must
        make every derived property return an empty list. A separate
        file with one malformed data row (first token not an integer)
        must have that row skipped via the except/continue branch."""
        empty_file = tmp_path / "vibspectrum_empty"
        empty_file.write_text(
            "$vibrational spectrum\n"
            "#  mode     symmetry     wave number   IR intensity\n"
            "$end\n"
        )
        empty_vib = XTBVibSpectrumFile(str(empty_file))
        assert empty_vib.vibrational_frequencies == []
        assert empty_vib.ir_intensities == []
        assert empty_vib.vibrational_mode_symmetries == []
        assert empty_vib.frequencies == []

        malformed_file = tmp_path / "vibspectrum_malformed"
        malformed_file.write_text(
            "$vibrational spectrum\n"
            "not_a_mode_index   garbage   row\n"
            "     1        a            151.34         0.04462         YES\n"
            "$end\n"
        )
        malformed_vib = XTBVibSpectrumFile(str(malformed_file))
        assert malformed_vib.vibrational_frequencies == [151.34]


class TestXTBWibergBondOrderFile:
    """Tests for XTBWibergBondOrderFile class."""

    def test_wbo_water(self, xtb_water_outfolder):
        wbo_file = os.path.join(xtb_water_outfolder, "wbo")
        wbo = XTBWibergBondOrderFile(wbo_file)
        assert len(wbo.bond_orders) == 2
        assert wbo.bond_orders[0] == (1, 2, 0.92021379026732564)
        assert wbo.bond_orders[1] == (1, 3, 0.92021379039282269)
        assert wbo.bond_order_matrix.shape == (3, 3)
        assert np.isclose(wbo.bond_order_matrix[0, 1], 0.92021379026732564)
        assert np.isclose(wbo.bond_order_matrix[1, 0], 0.92021379026732564)
        assert np.isclose(wbo.bond_order_matrix[0, 2], 0.92021379039282269)
        assert np.isclose(wbo.bond_order_matrix[2, 0], 0.92021379039282269)

    def test_bond_orders_skips_blank_and_malformed_lines(self, tmp_path):
        """Blank lines (too few tokens) and lines with non-numeric
        tokens must both be silently skipped rather than raising."""
        wbo_file = tmp_path / "wbo"
        wbo_file.write_text("\n" "a   b   not_a_number\n" "1   2   0.5\n")
        wbo = XTBWibergBondOrderFile(str(wbo_file))
        assert wbo.bond_orders == [(1, 2, 0.5)]

    def test_bond_order_matrix_none_when_no_bond_orders(self, tmp_path):
        """A wbo file with no valid bond-order lines at all -- both
        bond_orders and bond_order_matrix must be empty/None."""
        wbo_file = tmp_path / "wbo"
        wbo_file.write_text("\ntoo few\n")
        wbo = XTBWibergBondOrderFile(str(wbo_file))
        assert wbo.bond_orders == []
        assert wbo.bond_order_matrix is None


class TestXTBG98File:
    """Tests for XTBG98File class."""

    def test_g98_co2(self, xtb_co2_outfolder):
        """Test parsing G98 output from co2 ohess calculation."""
        g98_file = os.path.join(xtb_co2_outfolder, "g98.out")
        assert os.path.exists(g98_file)
        co2_g98 = XTBG98File(g98_file)
        assert co2_g98.standard_orientation == [
            [-1.143652, 0.000006, 0.000004],
            [1.143652, -0.000006, -0.000004],
            [0.000000, -0.000000, -0.000000],
        ]
        assert co2_g98.vibrational_frequencies == [
            600.3117,
            600.3117,
            1424.7819,
            2593.0748,
        ]
        assert co2_g98.reduced_masses == [13.0986, 13.0986, 15.9994, 13.0994]
        assert co2_g98.force_constants == [0.0, 0.0, 0.0, 0.0]
        assert co2_g98.ir_intensities == [68.6947, 68.6947, 0.0, 1046.6649]
        assert co2_g98.raman_activities == [0.0, 0.0, 0.0, 0.0]
        assert co2_g98.depolarization_ratios == [0.0, 0.0, 0.0, 0.0]
        assert co2_g98.vibrational_mode_symmetries == ["a", "a", "a", "a"]
        assert np.allclose(
            co2_g98.vibrational_modes[0],
            [
                [-0.00, -0.37, -0.00],
                [-0.00, -0.37, -0.00],
                [0.00, 0.85, 0.00],
            ],
        )
        assert np.allclose(
            co2_g98.vibrational_modes[1],
            [
                [-0.00, 0.00, -0.37],
                [-0.00, 0.00, -0.37],
                [0.00, -0.00, 0.85],
            ],
        )
        assert np.allclose(
            co2_g98.vibrational_modes[2],
            [
                [-0.71, 0.00, 0.00],
                [0.71, -0.00, -0.00],
                [-0.00, 0.00, -0.00],
            ],
        )
        assert np.allclose(
            co2_g98.vibrational_modes[3],
            [
                [0.37, -0.00, -0.00],
                [0.37, -0.00, -0.00],
                [-0.85, 0.00, 0.00],
            ],
        )
        assert co2_g98.num_vib_modes == 4
        assert co2_g98.num_vib_frequencies == 4

    def test_g98_acetaldehyde(self, xtb_acetaldehyde_outfolder):
        """Test parsing G98 output from acetaldehyde hess calculation."""
        g98_file = os.path.join(xtb_acetaldehyde_outfolder, "g98.out")
        assert os.path.exists(g98_file)
        acetaldehyde_g98 = XTBG98File(g98_file)
        assert acetaldehyde_g98.standard_orientation == [
            [4.238362, 0.515099, 0.187252],
            [1.871052, 0.670294, -0.172258],
            [3.153451, 0.008898, 0.220004],
            [2.046523, 1.693345, -0.493507],
            [1.187651, 0.658353, 0.675123],
            [1.407937, 0.102450, -0.977432],
            [3.022239, -1.038528, 0.560754],
        ]
        assert acetaldehyde_g98.symbols == ["O", "C", "C", "H", "H", "H", "H"]
        assert acetaldehyde_g98.vibrational_frequencies == [
            151.3396,
            501.8141,
            769.0542,
            947.2656,
            1045.6822,
            1107.2721,
            1355.3353,
            1389.3815,
            1446.6021,
            1447.8603,
            1798.5806,
            2748.9403,
            3018.3384,
            3026.5465,
            3059.7618,
        ]
        assert acetaldehyde_g98.reduced_masses == [
            3.5771,
            9.8091,
            2.4123,
            7.4606,
            6.1400,
            7.7131,
            3.0566,
            2.9306,
            1.7558,
            1.6047,
            13.3745,
            1.6770,
            1.9421,
            1.3886,
            1.8344,
        ]
        assert acetaldehyde_g98.force_constants == [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
        assert acetaldehyde_g98.ir_intensities == [
            0.0446,
            11.4541,
            2.9126,
            15.1649,
            16.5191,
            60.1232,
            12.4013,
            46.0938,
            16.0281,
            10.9071,
            300.9390,
            142.6312,
            6.5363,
            3.4114,
            14.2061,
        ]
        assert acetaldehyde_g98.raman_activities == [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
        assert acetaldehyde_g98.depolarization_ratios == [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
        assert acetaldehyde_g98.vibrational_mode_symmetries == [
            "a",
            "a",
            "a",
            "a",
            "a",
            "a",
            "a",
            "a",
            "a",
            "a",
            "a",
            "a",
            "a",
            "a",
            "a",
        ]
        assert np.allclose(
            acetaldehyde_g98.vibrational_modes[0],
            [
                [0.04, -0.09, -0.28],
                [0.00, -0.00, -0.01],
                [-0.04, 0.11, 0.32],
                [-0.05, 0.12, 0.36],
                [-0.27, -0.40, -0.23],
                [0.35, 0.20, -0.35],
                [-0.04, 0.09, 0.28],
            ],
        )
        assert np.allclose(
            acetaldehyde_g98.vibrational_modes[14],
            [
                [0.00, 0.00, -0.00],
                [-0.07, -0.25, 0.08],
                [-0.00, -0.01, 0.00],
                [0.16, 0.89, -0.28],
                [0.06, -0.03, -0.08],
                [0.04, 0.03, 0.09],
                [0.00, 0.02, -0.01],
            ],
        )
        assert acetaldehyde_g98.num_vib_modes == 15
        assert acetaldehyde_g98.num_vib_frequencies == 15

    def test_standard_orientation_block_runs_to_eof(self, tmp_path):
        """The Standard orientation table is the last thing in the
        file, with no closing dashed-line marker -- the row-collection
        loop must exhaust the file's contents naturally rather than
        break."""
        content = (
            "                     Standard orientation:\n"
            " --------------------------------------------------------------------\n"
            "  Center     Atomic     Atomic              Coordinates (Angstroms)\n"
            "  Number     Number      Type              X           Y           Z\n"
            " --------------------------------------------------------------------\n"
            "    1          8             0       -1.143652    0.000006    0.000004\n"
        )
        g98_file = tmp_path / "g98.out"
        g98_file.write_text(content)
        g98 = XTBG98File(str(g98_file))
        assert g98.standard_orientation == [[-1.143652, 0.000006, 0.000004]]
        assert g98.symbols == ["O"]

    def test_standard_orientation_none_when_heading_absent(self, tmp_path):
        """No 'Standard orientation:' heading anywhere in the file --
        both standard_orientation and symbols must be None."""
        g98_file = tmp_path / "g98.out"
        g98_file.write_text("nothing relevant here\n")
        g98 = XTBG98File(str(g98_file))
        assert g98.standard_orientation is None
        assert g98.symbols is None

    def test_vibrational_modes_empty_when_no_matching_rows(self, tmp_path):
        """A 'Frequencies --' line with nothing (or non-matching
        content) in its normal-mode data window -- first_col_vib_modes
        (and therefore the whole per-block collection) stays empty, so
        nothing is appended to the overall modes list for that block."""
        g98_file = tmp_path / "g98.out"
        g98_file.write_text("Frequencies --  1539.3017\n")
        g98 = XTBG98File(str(g98_file))
        assert g98.vibrational_modes == []
        assert g98.num_vib_modes == 0

    def test_attach_vib_metadata_sets_all_attributes(self, xtb_co2_outfolder):
        """Direct test of the private `_attach_vib_metadata` helper,
        which isn't exercised by any Molecule-building code path in
        the existing fixture-driven tests."""
        g98_file = os.path.join(xtb_co2_outfolder, "g98.out")
        co2_g98 = XTBG98File(g98_file)

        class _DummyMolecule:
            pass

        mol = _DummyMolecule()
        result = co2_g98._attach_vib_metadata(mol)

        assert result is mol
        assert mol.vibrational_frequencies == co2_g98.vibrational_frequencies
        assert mol.vibrational_reduced_masses == co2_g98.reduced_masses
        assert mol.vibrational_force_constants == co2_g98.force_constants
        assert mol.vibrational_ir_intensities == co2_g98.ir_intensities
        assert (
            mol.vibrational_mode_symmetries
            == co2_g98.vibrational_mode_symmetries
        )
        assert len(mol.vibrational_modes) == len(co2_g98.vibrational_modes)


class TestXTBFolder:
    """Tests for XTBFolder class."""

    def test_folder_co2(self, xtb_co2_outfolder):
        """Test XTBFolder with CO2 ohess calculation output."""
        assert os.path.exists(xtb_co2_outfolder)
        co2_folder = XTBFolder(xtb_co2_outfolder)
        assert co2_folder.is_xtb_calculation_directory

        assert co2_folder.xtb_out_filepath is not None
        assert os.path.basename(co2_folder.xtb_out_filepath) == "co2_ohess.out"

        assert co2_folder.xtbopt_log_filepath is not None
        assert os.path.basename(co2_folder.xtbopt_log_filepath) == "xtbopt.log"

        assert co2_folder.charges_filepath is not None
        assert os.path.basename(co2_folder.charges_filepath) == "charges"

        assert co2_folder.energy_filepath is not None
        assert os.path.basename(co2_folder.energy_filepath) == "energy"

        assert co2_folder.engrad_filepath is not None
        assert os.path.basename(co2_folder.engrad_filepath) == "co2.engrad"

        assert co2_folder.g98_out_filepath is not None
        assert os.path.basename(co2_folder.g98_out_filepath) == "g98.out"

        assert co2_folder.gradient_filepath is not None
        assert os.path.basename(co2_folder.gradient_filepath) == "gradient"

        assert co2_folder.hessian_filepath is not None
        assert os.path.basename(co2_folder.hessian_filepath) == "hessian"

        assert co2_folder.vibspectrum_filepath is not None
        assert (
            os.path.basename(co2_folder.vibspectrum_filepath) == "vibspectrum"
        )

        assert co2_folder.wbo_filepath is not None
        assert os.path.basename(co2_folder.wbo_filepath) == "wbo"

        assert co2_folder.input_geometry_filepath is not None
        assert (
            os.path.basename(co2_folder.input_geometry_filepath) == "co2.xyz"
        )

        assert co2_folder.xtbopt_geometry_filepath is not None
        assert (
            os.path.basename(co2_folder.xtbopt_geometry_filepath)
            == "xtbopt.xyz"
        )

        assert co2_folder.xtbtopo_mol_filepath is not None
        assert (
            os.path.basename(co2_folder.xtbtopo_mol_filepath) == "xtbtopo.mol"
        )

    def test_folder_cyclopentadienyl_anion(
        self, xtb_cyclopentadienyl_anion_outfolder
    ):
        """Test XTBFolder with cyclopentadienyl anion opt calculation output."""
        assert os.path.exists(xtb_cyclopentadienyl_anion_outfolder)
        cyclopentadienyl_anion_folder = XTBFolder(
            xtb_cyclopentadienyl_anion_outfolder
        )
        assert cyclopentadienyl_anion_folder.is_xtb_calculation_directory

        assert cyclopentadienyl_anion_folder.xtb_out_filepath is not None
        assert (
            os.path.basename(cyclopentadienyl_anion_folder.xtb_out_filepath)
            == "cyclopentadienyl_anion_opt.out"
        )

        assert cyclopentadienyl_anion_folder.xtbopt_log_filepath is not None
        assert cyclopentadienyl_anion_folder.charges_filepath is not None
        assert (
            cyclopentadienyl_anion_folder.energy_filepath is None
        )  # --grad calculation is not enabled
        assert (
            cyclopentadienyl_anion_folder.engrad_filepath is None
        )  # --grad calculation is not enabled
        assert (
            cyclopentadienyl_anion_folder.g98_out_filepath is None
        )  # --hess calculation is not enabled
        assert (
            cyclopentadienyl_anion_folder.gradient_filepath is None
        )  # --grad calculation is not enabled
        assert (
            cyclopentadienyl_anion_folder.hessian_filepath is None
        )  # --hess calculation is not enabled
        assert (
            cyclopentadienyl_anion_folder.vibspectrum_filepath is None
        )  # --hess calculation is not enabled
        assert cyclopentadienyl_anion_folder.wbo_filepath is not None
        assert (
            os.path.basename(
                cyclopentadienyl_anion_folder.input_geometry_filepath
            )
            == "cyclopentadienyl_anion.coord"
        )
        assert (
            os.path.basename(
                cyclopentadienyl_anion_folder.xtbopt_geometry_filepath
            )
            == "xtbopt.coord"
        )
        assert cyclopentadienyl_anion_folder.xtbtopo_mol_filepath is not None

    def test_folder_p_benzyne_sp(self, xtb_p_benzyne_sp_outfolder):
        """Test XTBFolder with p-benzyne sp calculation output."""
        assert os.path.exists(xtb_p_benzyne_sp_outfolder)
        p_benzyne_sp_folder = XTBFolder(xtb_p_benzyne_sp_outfolder)
        assert p_benzyne_sp_folder.is_xtb_calculation_directory

        assert p_benzyne_sp_folder.xtb_out_filepath is not None
        assert (
            os.path.basename(p_benzyne_sp_folder.xtb_out_filepath)
            == "p_benzyne_sp_alpb_toluene.out"
        )

        assert (
            p_benzyne_sp_folder.xtbopt_log_filepath is None
        )  # no optimization performed
        assert p_benzyne_sp_folder.charges_filepath is not None
        assert (
            p_benzyne_sp_folder.energy_filepath is None
        )  # --grad calculation is not enabled
        assert (
            p_benzyne_sp_folder.engrad_filepath is None
        )  # --grad calculation is not enabled
        assert (
            p_benzyne_sp_folder.g98_out_filepath is None
        )  # --hess calculation is not enabled
        assert (
            p_benzyne_sp_folder.gradient_filepath is None
        )  # --grad calculation is not enabled
        assert (
            p_benzyne_sp_folder.hessian_filepath is None
        )  # --hess calculation is not enabled
        assert (
            p_benzyne_sp_folder.vibspectrum_filepath is None
        )  # --hess calculation is not enabled
        assert p_benzyne_sp_folder.wbo_filepath is not None
        assert (
            p_benzyne_sp_folder.xtbopt_geometry_filepath is None
        )  # no optimization performed
        assert p_benzyne_sp_folder.input_geometry_filepath is not None
        assert p_benzyne_sp_folder.xtbtopo_mol_filepath is not None

    def test_is_xtb_calculation_directory_false_without_output(self, tmp_path):
        folder = XTBFolder(str(tmp_path))
        assert folder._xtb_out() is None
        assert folder.is_xtb_calculation_directory is False

    def test_is_xtb_calculation_directory_false_without_markers(
        self, tmp_path
    ):
        """An xTB *.out file is present, but none of the common
        auxiliary marker files (xtbrestart, charges, xtbtopo.mol,
        wbo) exist -- must not be considered a valid xTB directory."""
        (tmp_path / "job.out").write_text("x T B\nsome xtb output\n")
        folder = XTBFolder(str(tmp_path))
        assert folder._xtb_out() is not None
        assert folder.is_xtb_calculation_directory is False

    def test_input_geometry_skips_xtbopt_prefixed_parseable_files(
        self, tmp_path
    ):
        """A parseable-format file whose name starts with "xtbopt" is
        the optimized-geometry output, not the input geometry, and
        must be skipped by _input_geometry."""
        (tmp_path / "xtbopt.xyz").write_text("optimized geometry\n")
        real_input = tmp_path / "structure.xyz"
        real_input.write_text("input geometry\n")

        folder = XTBFolder(str(tmp_path))
        assert folder._input_geometry() == str(real_input)

    def test_input_geometry_skips_xtbopt_prefixed_unsupported_files(
        self, tmp_path
    ):
        (tmp_path / "xtbopt.coord").write_text("optimized geometry\n")
        real_input = tmp_path / "structure.coord"
        real_input.write_text("input geometry\n")

        folder = XTBFolder(str(tmp_path))
        assert folder._input_geometry() == str(real_input)

    def test_input_geometry_none_when_nothing_found(self, tmp_path):
        folder = XTBFolder(str(tmp_path))
        assert folder._input_geometry() is None

    def test_xtb_out_raises_for_multiple_output_files(self, tmp_path):
        import pytest

        (tmp_path / "job1.out").write_text("x T B\nsome xtb output\n")
        (tmp_path / "job2.out").write_text("x T B\nsome xtb output\n")

        folder = XTBFolder(str(tmp_path))
        with pytest.raises(ValueError, match="Multiple xTB main output"):
            folder._xtb_out()


class TestXTBOutput:

    def test_ohess_output(self, xtb_co2_outfolder):
        assert os.path.exists(xtb_co2_outfolder)
        xtb_co2_output = XTBOutput(folder=xtb_co2_outfolder)
        assert xtb_co2_output.normal_termination
        assert xtb_co2_output.geometry_optimization_converged
        assert xtb_co2_output.charge == 0
        assert xtb_co2_output.multiplicity == 1
        assert xtb_co2_output.mass == 44.0095457
        assert xtb_co2_output.final_energy == -10.308452289174
        assert np.allclose(
            xtb_co2_output.final_forces,
            [
                [0.000000194263, -0.000000000000, -0.000000000000],
                [-0.000000194263, 0.000000000000, -0.000000000000],
                [0.000000000000, 0.000000000000, 0.000000000000],
            ],
        )
        assert xtb_co2_output.symbols == ["O", "O", "C"]
        assert xtb_co2_output.partial_charges == {
            "O1": -0.23213972,
            "O2": -0.23213972,
            "C1": 0.46427944,
        }
        optimized_flags = [
            mol.is_optimized_structure for mol in xtb_co2_output.all_structures
        ]
        assert len(xtb_co2_output.all_structures) == 5
        assert optimized_flags == [False] * 4 + [True]


class TestXTBOutputGetattrDelegation:
    """Direct tests for XTBOutput.__getattr__'s parser-search loop."""

    def test_delegates_to_a_later_parser_after_skipping_none_and_no_attr(
        self, xtb_p_benzyne_sp_outfolder
    ):
        """p_benzyne_sp lacks energy/g98/gradient/hessian/vibspectrum
        files (all None, skipped), main_out doesn't have `bond_orders`
        (AttributeError, skipped), so delegation must reach wbo_file."""
        output = XTBOutput(folder=xtb_p_benzyne_sp_outfolder)
        assert output.energy_file is None
        assert output.g98_file is None
        assert output.gradient_file is None
        assert output.hessian_file is None
        assert output.vibspectrum_file is None
        assert output.bond_orders is not None

    def test_unknown_attribute_raises_attribute_error(self, xtb_co2_outfolder):
        output = XTBOutput(folder=xtb_co2_outfolder)
        with pytest.raises(AttributeError, match="no attribute"):
            output.totally_made_up_attribute_xyz

    def test_non_attribute_error_from_a_parser_is_logged_and_skipped(
        self, xtb_co2_outfolder
    ):
        """A parser that raises something other than AttributeError
        (e.g. a bug in a third-party file parser) must not abort the
        whole delegation search -- it's logged and the next parser in
        the priority list is tried instead."""

        class _BrokenParser:
            def __getattr__(self, name):
                raise RuntimeError("simulated parser bug")

        output = XTBOutput(folder=xtb_co2_outfolder)
        output.__dict__["main_out"] = _BrokenParser()
        assert output.bond_orders is not None


class TestXTBOutputFallbackProperties:
    """Direct tests for XTBOutput's properties that fall back across
    multiple underlying parsers, using __dict__ overrides on cached
    properties (functools.cached_property stores its value directly in
    the instance __dict__, so assigning there bypasses recomputation)
    to force otherwise-unreachable-via-real-fixtures branches."""

    def test_normal_termination_and_convergence_false_without_main_out(
        self, xtb_co2_outfolder
    ):
        output = XTBOutput(folder=xtb_co2_outfolder)
        output.__dict__["main_out"] = None
        assert output.normal_termination is False
        assert output.geometry_optimization_converged is False

    def test_charge_falls_back_to_charges_file_then_none(
        self, xtb_co2_outfolder
    ):
        output = XTBOutput(folder=xtb_co2_outfolder)
        output.__dict__["main_out"] = None
        assert output.charge == 0.0

        output2 = XTBOutput(folder=xtb_co2_outfolder)
        output2.__dict__["main_out"] = None
        output2.__dict__["charges_file"] = None
        assert output2.charge is None

    def test_multiplicity_none_without_main_out(self, xtb_co2_outfolder):
        output = XTBOutput(folder=xtb_co2_outfolder)
        output.__dict__["main_out"] = None
        assert output.multiplicity is None

    def test_mass_falls_back_to_molecule_then_none(self, xtb_co2_outfolder):
        output = XTBOutput(folder=xtb_co2_outfolder)
        output.__dict__["main_out"] = None
        assert output.mass == pytest.approx(44.009, abs=0.01)

        output2 = XTBOutput(folder=xtb_co2_outfolder)
        output2.__dict__["main_out"] = None
        output2.__dict__["all_structures"] = []
        assert output2.mass is None

    def test_num_atoms_falls_back_to_molecule_then_none(
        self, xtb_co2_outfolder
    ):
        output = XTBOutput(folder=xtb_co2_outfolder)
        output.__dict__["engrad_file"] = None
        assert output.num_atoms == 3

        output2 = XTBOutput(folder=xtb_co2_outfolder)
        output2.__dict__["engrad_file"] = None
        output2.__dict__["all_structures"] = []
        assert output2.num_atoms is None

    def test_final_energy_falls_back_through_energy_then_engrad_then_none(
        self, xtb_co2_outfolder
    ):
        output = XTBOutput(folder=xtb_co2_outfolder)
        output.__dict__["main_out"] = None
        assert output.final_energy == pytest.approx(-10.30845228917)

        output2 = XTBOutput(folder=xtb_co2_outfolder)
        output2.__dict__["main_out"] = None
        output2.__dict__["energy_file"] = None
        assert output2.final_energy == pytest.approx(
            output2.engrad_file.total_energy
        )

        output3 = XTBOutput(folder=xtb_co2_outfolder)
        output3.__dict__["main_out"] = None
        output3.__dict__["energy_file"] = None
        output3.__dict__["engrad_file"] = None
        assert output3.final_energy is None

    def test_final_forces_none_without_engrad_file(
        self, xtb_p_benzyne_sp_outfolder
    ):
        output = XTBOutput(folder=xtb_p_benzyne_sp_outfolder)
        assert output.engrad_file is None
        assert output.final_forces is None

    def test_symbols_falls_back_to_all_structures_then_none(
        self, xtb_co2_outfolder
    ):
        output = XTBOutput(folder=xtb_co2_outfolder)
        output.__dict__["xtbopt_geometry"] = None
        assert output.symbols == ["O", "O", "C"]

        output2 = XTBOutput(folder=xtb_co2_outfolder)
        output2.__dict__["xtbopt_geometry"] = None
        output2.__dict__["all_structures"] = []
        assert output2.symbols is None

    def test_partial_charges_none_without_charges_file_or_symbols(
        self, xtb_co2_outfolder
    ):
        output = XTBOutput(folder=xtb_co2_outfolder)
        output.__dict__["charges_file"] = None
        assert output.partial_charges is None

        output2 = XTBOutput(folder=xtb_co2_outfolder)
        output2.__dict__["symbols"] = None
        assert output2.partial_charges is None

    def test_vibrational_frequencies_monoatomic_is_empty(
        self, xtb_he_outfolder
    ):
        output = XTBOutput(folder=xtb_he_outfolder)
        assert output.molecule.is_monoatomic
        assert output.vibrational_frequencies == []

    def test_vibrational_frequencies_falls_back_to_g98_then_none(
        self, xtb_co2_outfolder
    ):
        output = XTBOutput(folder=xtb_co2_outfolder)
        output.__dict__["main_out"] = None
        assert (
            output.vibrational_frequencies
            == output.g98_file.vibrational_frequencies
        )

        output2 = XTBOutput(folder=xtb_co2_outfolder)
        output2.__dict__["main_out"] = None
        output2.__dict__["g98_file"] = None
        assert output2.vibrational_frequencies is None

    def test_optimized_structure_falls_back_to_last_structure(
        self, xtb_co2_outfolder
    ):
        output = XTBOutput(folder=xtb_co2_outfolder)
        output.__dict__["xtbopt_geometry"] = None
        assert output.optimized_structure is output.all_structures[-1]

    def test_xtbopt_geometry_enrichment_skipped_when_all_none(
        self, xtb_co2_outfolder
    ):
        """Covers the "if X is not None" guards' False arms in
        xtbopt_geometry's enrichment block (charge/multiplicity/
        final_energy/final_forces), which every real fixture always
        supplies (via main_out), making the skip-arms otherwise
        unreachable."""
        output = XTBOutput(folder=xtb_co2_outfolder)
        with (
            patch.object(
                XTBOutput,
                "charge",
                new_callable=PropertyMock,
                return_value=None,
            ),
            patch.object(
                XTBOutput,
                "multiplicity",
                new_callable=PropertyMock,
                return_value=None,
            ),
            patch.object(
                XTBOutput,
                "final_energy",
                new_callable=PropertyMock,
                return_value=None,
            ),
            patch.object(
                XTBOutput,
                "final_forces",
                new_callable=PropertyMock,
                return_value=None,
            ),
        ):
            mol = output.xtbopt_geometry
        assert mol is not None

    def test_input_geometry_enrichment_skipped_when_all_none(
        self, xtb_co2_outfolder
    ):
        output = XTBOutput(folder=xtb_co2_outfolder)
        with (
            patch.object(
                XTBOutput,
                "charge",
                new_callable=PropertyMock,
                return_value=None,
            ),
            patch.object(
                XTBOutput,
                "multiplicity",
                new_callable=PropertyMock,
                return_value=None,
            ),
        ):
            mol = output.input_geometry
        assert mol is not None


class TestXTBOutputUnsupportedGeometryFormats:
    """cyclopentadienyl_anion_opt only has xtbopt.coord/*.coord
    geometry files, an unsupported format for both the optimized and
    input geometry files -- these still get located at the folder
    level (with a warning) but XTBOutput's cached properties must
    return None for them since only .xyz/.sdf/.pdb are dispatched."""

    def test_xtbopt_geometry_file_none_for_unsupported_format(
        self, xtb_cyclopentadienyl_anion_outfolder
    ):
        output = XTBOutput(folder=xtb_cyclopentadienyl_anion_outfolder)
        assert output.folder._xtbopt_geometry() is not None
        assert output.xtbopt_geometry_file is None

    def test_input_geometry_file_none_for_unsupported_format(
        self, xtb_cyclopentadienyl_anion_outfolder
    ):
        output = XTBOutput(folder=xtb_cyclopentadienyl_anion_outfolder)
        assert output.folder._input_geometry() is not None
        assert output.input_geometry_file is None

    def test_xtbopt_geometry_file_none_when_no_xtbopt_file_at_all(
        self, xtb_p_benzyne_sp_outfolder
    ):
        output = XTBOutput(folder=xtb_p_benzyne_sp_outfolder)
        assert output.folder._xtbopt_geometry() is None
        assert output.xtbopt_geometry_file is None

    def test_input_geometry_file_none_when_no_geometry_file_at_all(
        self, xtb_co2_outfolder
    ):
        output = XTBOutput(folder=xtb_co2_outfolder)
        with patch.object(
            type(output.folder), "_input_geometry", return_value=None
        ):
            assert output.input_geometry_file is None


class TestXTBOutputReadGeometryFile:
    def test_missing_file_returns_none_and_logs_error(self, xtb_co2_outfolder):
        output = XTBOutput(folder=xtb_co2_outfolder)
        assert output._read_geometry_file("/no/such/path/geometry.xyz") is None

    def test_none_path_returns_none(self, xtb_co2_outfolder):
        output = XTBOutput(folder=xtb_co2_outfolder)
        assert output._read_geometry_file(None) is None

    def test_unsupported_extension_returns_none(
        self, xtb_cyclopentadienyl_anion_outfolder
    ):
        output = XTBOutput(folder=xtb_cyclopentadienyl_anion_outfolder)
        path = output.folder._input_geometry()
        assert path.endswith(".coord")
        assert output._read_geometry_file(path) is None

    def test_input_geometry_none_for_unsupported_format(
        self, xtb_cyclopentadienyl_anion_outfolder
    ):
        """Covers input_geometry's "if molecule:" False arm, reached
        when _read_geometry_file returns None for an unsupported
        format."""
        output = XTBOutput(folder=xtb_cyclopentadienyl_anion_outfolder)
        assert output.input_geometry is None


class TestXTBOutputChooseOrientationsTiers:
    """Covers _choose_orientations' fallback tiers: xtbopt.log
    (tier 1, already covered via the ohess fixtures elsewhere),
    g98 standard orientation (tier 3), input geometry for a completed
    single-point job (tier 4), and the final empty-tuple fallback."""

    def test_tier2_xtbopt_geometry_used_when_log_unavailable(
        self, xtb_co2_outfolder
    ):
        """co2_ohess has both xtbopt.log (tier 1) and xtbopt.xyz; force
        tier 1 to be skipped so tier 2's own branch body runs."""
        output = XTBOutput(folder=xtb_co2_outfolder)
        assert output.xtbopt_geometry is not None
        output.__dict__["xtbopt_log_file"] = None
        orientations, symbols, energies = output._choose_orientations()
        assert len(orientations) == 1
        assert symbols == list(output.xtbopt_geometry.symbols)

    def test_tier1_skipped_when_log_has_no_molecules(self, xtb_co2_outfolder):
        output = XTBOutput(folder=xtb_co2_outfolder)
        real_xtbopt_geometry = output.xtbopt_geometry
        fake_log_file = MagicMock()
        fake_log_file.get_molecules.return_value = []
        output.__dict__["xtbopt_log_file"] = fake_log_file
        orientations, symbols, energies = output._choose_orientations()
        # falls through to tier 2 (xtbopt_geometry) instead
        assert len(orientations) == 1
        assert symbols == list(real_xtbopt_geometry.symbols)

    def test_tier3_g98_standard_orientation_used_for_hess_job(
        self, xtb_acetaldehyde_outfolder
    ):
        output = XTBOutput(folder=xtb_acetaldehyde_outfolder)
        assert output.xtbopt_log_file is None
        assert output.xtbopt_geometry is None
        assert output.g98_file is not None
        orientations, symbols, energies = output._choose_orientations()
        assert len(orientations) == 1
        assert symbols == list(output.g98_file.symbols)

    def test_tier4_input_geometry_used_for_completed_sp_job(
        self, xtb_p_benzyne_sp_outfolder
    ):
        output = XTBOutput(folder=xtb_p_benzyne_sp_outfolder)
        assert output.xtbopt_log_file is None
        assert output.xtbopt_geometry is None
        assert output.g98_file is None
        assert output.normal_termination is True
        orientations, symbols, energies = output._choose_orientations()
        assert len(orientations) == 1

    def test_no_tier_available_returns_empty(self, xtb_p_benzyne_sp_outfolder):
        """Forcing normal_termination False removes the only tier
        available to this fixture (tier 4), leaving nothing."""
        output = XTBOutput(folder=xtb_p_benzyne_sp_outfolder)
        with patch.object(
            XTBOutput,
            "normal_termination",
            new_callable=PropertyMock,
            return_value=False,
        ):
            orientations, symbols, energies = output._choose_orientations()
            assert output.all_structures == []
        assert orientations == []
        assert symbols is None
        assert energies is None


class TestXTBOutputIsOptimizedStructureList:
    """Covers _compute_is_optimized_structure_list's hess-job branch
    (opt-job branch already covered by the ohess fixtures elsewhere)."""

    def test_hess_job_marks_last_structure_optimized(
        self, xtb_acetaldehyde_outfolder
    ):
        output = XTBOutput(folder=xtb_acetaldehyde_outfolder)
        assert output.jobtype == "hess"
        flags = [m.is_optimized_structure for m in output.all_structures]
        assert flags == [True]

    def test_neither_opt_nor_hess_job_marks_nothing_optimized(
        self, xtb_p_benzyne_sp_outfolder
    ):
        output = XTBOutput(folder=xtb_p_benzyne_sp_outfolder)
        assert output.jobtype == "sp"
        flags = [m.is_optimized_structure for m in output.all_structures]
        assert flags == [False]


class TestXTBOutputNumAtomsAndVibrationalFrequenciesNaturalPaths:
    """Covers the "primary source available" branches for num_atoms
    (engrad_file present) and vibrational_frequencies (main_out
    present), which the fallback-focused tests above always bypassed
    via __dict__ overrides."""

    def test_num_atoms_from_engrad_file(self, xtb_co2_outfolder):
        output = XTBOutput(folder=xtb_co2_outfolder)
        assert output.engrad_file is not None
        assert output.num_atoms == 3

    def test_vibrational_frequencies_from_main_out(self, xtb_co2_outfolder):
        output = XTBOutput(folder=xtb_co2_outfolder)
        assert output.main_out is not None
        assert (
            output.vibrational_frequencies
            == output.main_out.vibrational_frequencies
        )


class TestXTBOutputOptimizedStructureNaturalAndEdgeCases:
    def test_returns_none_when_not_converged(self, xtb_p_benzyne_sp_outfolder):
        output = XTBOutput(folder=xtb_p_benzyne_sp_outfolder)
        assert output.geometry_optimization_converged is False
        assert output.optimized_structure is None

    def test_returns_xtbopt_geometry_when_available(self, xtb_co2_outfolder):
        output = XTBOutput(folder=xtb_co2_outfolder)
        assert output.xtbopt_geometry is not None
        assert output.optimized_structure is output.xtbopt_geometry

    def test_returns_none_when_converged_but_no_structures_available(
        self, xtb_co2_outfolder
    ):
        output = XTBOutput(folder=xtb_co2_outfolder)
        output.__dict__["xtbopt_geometry"] = None
        output.__dict__["all_structures"] = []
        with patch.object(
            XTBOutput,
            "geometry_optimization_converged",
            new_callable=PropertyMock,
            return_value=True,
        ):
            assert output.optimized_structure is None


class TestXTBOutputGetMolecule:
    def test_returns_molecule_at_given_index(self, xtb_co2_outfolder):
        output = XTBOutput(folder=xtb_co2_outfolder)
        assert output.get_molecule("-1") == output.all_structures[-1]
        assert output.get_molecule("1") == output.all_structures[0]

    def test_raises_when_no_structures_available(
        self, xtb_p_benzyne_sp_outfolder
    ):
        output = XTBOutput(folder=xtb_p_benzyne_sp_outfolder)
        with patch.object(
            XTBOutput,
            "normal_termination",
            new_callable=PropertyMock,
            return_value=False,
        ):
            with pytest.raises(ValueError, match="No molecule could be found"):
                output.get_molecule()
