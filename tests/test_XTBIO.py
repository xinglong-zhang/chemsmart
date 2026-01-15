import os.path

from chemsmart.io.xtb.file import XTBChargesFile, XTBMainOut
from chemsmart.io.xtb.input import XTBInput
from chemsmart.io.xtb.output import XTBOutput


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

    def test_main_out_opt(self, xtb_co2_outfolder):
        """Test parsing main output from CO2 opt calculation."""
        xtb_main_out_file = os.path.join(xtb_co2_outfolder, "co2.out")
        assert os.path.exists(xtb_main_out_file)
        co2_main_out = XTBMainOut(xtb_main_out_file)
        assert co2_main_out.xtb_version == "6.7.1"
        assert co2_main_out.normal_termination
        assert (
            co2_main_out.route_string
            == "/Users/xinglongzhang/miniconda3/envs/chemsmart/bin/xtb /Users/xinglongzhang/chemsmart_test/xtb/co2.xyz "
            "--gfn2 --opt vtight --chrg 0 --uhf 0"
        )
        assert co2_main_out.num_basis_functions == 12
        assert co2_main_out.num_atomic_orbital == 12
        assert co2_main_out.num_shells == 6
        assert co2_main_out.num_electrons == 16
        assert co2_main_out.max_iter == 250
        assert co2_main_out.hamiltonian == "GFN2-xTB"
        assert co2_main_out.restart is True
        assert co2_main_out.solvent_on is False
        assert co2_main_out.pc_potential is False
        assert co2_main_out.electronic_temperature == 300.0
        assert co2_main_out.accuracy == 1.0
        assert co2_main_out.integral_cutoff == 25.0
        assert co2_main_out.integral_neglect == 1e-8
        assert co2_main_out.scf_convergence == 1.0e-6
        assert co2_main_out.wf_convergence == 1.0e-4
        assert co2_main_out.broyden_damping == 0.4
        assert co2_main_out.net_charge == 0
        assert co2_main_out.unpaired_electrons == 0
        assert co2_main_out.optimization_level == "verytight"
        # assert xtb_output.max_optcycles is None
        # assert xtb_output.anc_microcycles is None
        # assert xtb_output.degrees_of_freedom is None
        # assert xtb_output.rf_solver is None
        # assert not xtb_output.write_all_intermediate_geometries
        # assert not xtb_output.is_linear
        # assert xtb_output.energy_convergence is None
        # assert xtb_output.gradient_convergence is None
        # assert xtb_output.max_rf_displacement is None
        # assert xtb_output.low_frequency_cutoff is None
        # assert xtb_output.max_frequency_cutoff is None
        # assert xtb_output.s6_in_model_hessian is None
        # assert xtb_output.homo_energy == -12.1467
        # assert xtb_output.lumo_energy == 2.2442
        # assert xtb_output.c6_coefficient == 44.535326
        # assert xtb_output.c8_coefficient == 795.739567
        # assert xtb_output.alpha_coefficient == 9.429122
        # assert xtb_output.total_energy_without_gsasa_hb is None
        # assert xtb_output.scc_energy == -5.104925504312
        # assert xtb_output.isotropic_es == 0.031459394051
        # assert xtb_output.anisotropic_es == 0.000394673573
        # assert xtb_output.anisotropic_xc == -0.000882256681
        # assert xtb_output.dispersion_energy == -0.000141082937
        # assert xtb_output.solvation_energy_gsolv is None
        # assert xtb_output.electronic_solvation_energy_gelec is None
        # assert xtb_output.surface_area_solvation_energy_gsasa is None
        # assert xtb_output.hydrogen_bonding_solvation_energy_ghb is None
        # assert xtb_output.empirical_shift_correction_gshift is None
        # assert xtb_output.repulsion_energy == 0.034381060848
        # assert xtb_output.additional_restraining_energy == 0.0
        # assert not xtb_output.numfreq
        # assert xtb_output.hessian_step_length is None
        # assert xtb_output.scc_accuracy is None
        # assert xtb_output.hessian_scale_factor is None
        # assert xtb_output.rms_gradient is None
        # assert xtb_output.total_charge == 0
        # assert np.allclose(
        #     xtb_output.qonly_molecular_dipole, [-0.0, 0.0, 0.607]
        # )
        # assert np.allclose(
        #     xtb_output.full_molecular_dipole, [-0.0, -0.0, 0.872]
        # )
        # assert xtb_output.total_molecular_dipole_moment == 2.217
        # assert np.allclose(
        #     xtb_output.qonly_molecular_quadrupole,
        #     [[1.311, 0.0, 0.0], [0.0, -0.492, 0.0], [0.0, 0.0, -0.819]],
        # )
        # assert np.allclose(
        #     xtb_output.q_dip_molecular_quadrupole,
        #     [[1.747, 0.0, 0.0], [0.0, -0.572, 0.0], [0.0, 0.0, -1.176]],
        # )
        # assert np.allclose(
        #     xtb_output.full_molecular_quadrupole,
        #     [[1.951, 0.0, 0.0], [0.0, -0.831, 0.0], [0.0, 0.0, -1.121]],
        # )
        # assert xtb_output.total_energy == -5.070544443464
        # assert xtb_output.gradient_norm == 0.000075164743
        # assert xtb_output.fmo_gap == 14.390891673350
        # assert xtb_output.accuracy == 1.0
        # assert xtb_output.integral_cutoff == 25.0
        # assert xtb_output.integral_neglect == 1e-8
        # assert xtb_output.scf_convergence == 1.0e-6
        # assert xtb_output.wf_convergence == 1.0e-4
        # assert xtb_output.broyden_damping == 0.4
        # assert xtb_output.net_charge == 0
        # assert xtb_output.unpaired_electrons == 0
        # assert xtb_output.optimization_level == "normal"
        # assert xtb_output.max_optcycles == 200
        # assert xtb_output.anc_microcycles == 20
        # assert xtb_output.degrees_of_freedom == 3
        # assert xtb_output.rf_solver == "davidson"
        # assert xtb_output.write_all_intermediate_geometries
        # assert not xtb_output.is_linear
        # assert xtb_output.energy_convergence == 5.0e-6
        # assert xtb_output.gradient_convergence == 1.0e-3
        # assert xtb_output.max_rf_displacement == 1.0
        # assert xtb_output.low_frequency_cutoff == 0.01
        # assert xtb_output.max_frequency_cutoff == 5.0
        # assert xtb_output.s6_in_model_hessian == 20.0
        # assert xtb_output.geometry_optimization_converged
        # assert xtb_output.route_string == "xtb coord --opt"
        # assert xtb_output.homo_energy == -12.1467
        # assert xtb_output.lumo_energy == 2.2442
        # assert xtb_output.optimization_level == "normal"
        # assert xtb_output.degrees_of_freedom == 3
        # assert xtb_output.optimized_structure_block == [
        #     "$coord",
        #     "0.00000000011942       -0.00000000000000       -0.71677520925432      o",
        #     "1.45926122846511       -0.00000000000000        0.35838760458144      h",
        #     "-1.45926122858453        0.00000000000000        0.35838760467288      h",
        #     "$end",
        #     "",
        # ]
        # assert xtb_output.molecular_mass == 18.0152864
        # assert xtb_output.center_of_mass == [0.0, 0.0, -0.3156364]
        # assert xtb_output.moments_of_inertia == [
        #     0.5795334e00,
        #     0.1202080e01,
        #     0.1781614e01,
        # ]
        # assert xtb_output.rotational_constants == [
        #     0.2908828e02,
        #     0.1402372e02,
        #     0.9462003e01,
        # ]
        # assert xtb_output.total_energy == -5.070544443465
        # assert xtb_output.gradient_norm == 0.000074994303
        # assert xtb_output.fmo_gap == 14.390898452735
        #
        # all_summary_blocks = xtb_output.get_all_summary_blocks()
        # assert len(all_summary_blocks) == 2
        # assert len(all_summary_blocks[0]) == 11
        #
        # assert not xtb_output.solvent_on
        # assert xtb_output.total_energy_without_gsasa_hb is None


class TestXTBChargesFile:
    """Tests for XTBChargesFile class."""

    def test_partial_charges_opt(self, xtb_co2_outfolder):
        """Test parsing charges from CO2 opt calculation."""
        charges_file = os.path.join(xtb_co2_outfolder, "charges")
        assert os.path.exists(charges_file)
        co2_charges = XTBChargesFile(charges_file)
        assert co2_charges.partial_charges == [
            -0.23213972,
            -0.23213972,
            0.46427944,
        ]
        assert co2_charges.total_charge == 0

    def test_charges_ohess(self, xtb_water_outfolder):
        """Test parsing charges from water ohess calculation."""
        charges_file = os.path.join(xtb_water_outfolder, "charges")
        assert os.path.exists(charges_file)
        water_charges = XTBChargesFile(charges_file)
        assert water_charges.partial_charges == [
            -0.56472698,
            0.28236349,
            0.28236349,
        ]
        assert water_charges.total_charge == 0


class TestXTBOutput:

    def test_opt_output(self, xtb_co2_outfolder):
        assert os.path.exists(xtb_co2_outfolder)
        xtb_co2_output = XTBOutput(folder=xtb_co2_outfolder)
        assert xtb_co2_output.normal_termination
        assert xtb_co2_output.partial_charges == {
            "O1": -0.23213972,
            "O2": -0.23213972,
            "C1": 0.46427944,
        }
