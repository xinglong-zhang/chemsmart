import os.path

import numpy as np

from chemsmart.io.xtb.file import (
    XTBChargesFile,
    XTBEnergyFile,
    XTBEngradFile,
    XTBG98File,
    XTBMainOut,
)
from chemsmart.io.xtb.folder import XTBFolder
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

    def test_main_out_co2(self, xtb_co2_outfolder):
        """Test parsing main output from CO2 ohess calculation."""
        xtb_main_out_file = os.path.join(xtb_co2_outfolder, "co2_ohess.out")
        assert os.path.exists(xtb_main_out_file)
        co2_main_out = XTBMainOut(xtb_main_out_file)
        assert co2_main_out.xtb_version == "6.7.1"
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
        assert not co2_main_out.restart
        assert not co2_main_out.solvent_on
        assert not co2_main_out.pc_potential
        assert co2_main_out.electronic_temperature == 300.0
        assert co2_main_out.accuracy == 1.0
        assert co2_main_out.integral_cutoff == 25.0
        assert co2_main_out.integral_neglect == 1e-8
        assert co2_main_out.scf_convergence == 1.0e-6
        assert co2_main_out.wf_convergence == 1.0e-4
        assert co2_main_out.broyden_damping == 0.4
        assert co2_main_out.net_charge == 0
        assert co2_main_out.unpaired_electrons == 0
        # Geometry Optimization Setup
        assert co2_main_out.optimization_level == "verytight"
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
            co2_main_out.qonly_molecular_dipole, [0.0, -0.0, -0.0]
        )
        assert np.allclose(
            co2_main_out.full_molecular_dipole, [0.0, 0.0, -0.0]
        )
        assert co2_main_out.total_molecular_dipole_moment == 0.0
        assert np.allclose(
            co2_main_out.qonly_molecular_quadrupole,
            [[-2.169, 0.0, 0.0], [0.0, 1.084, -0.0], [0.0, -0.0, 1.084]],
        )
        assert np.allclose(
            co2_main_out.q_dip_molecular_quadrupole,
            [[-3.107, 0.0, 0.0], [0.0, 1.553, -0.0], [0.0, -0.0, 1.553]],
        )
        assert np.allclose(
            co2_main_out.full_molecular_quadrupole,
            [[-4.360, 0.0, 0.0], [0.0, 2.180, -0.0], [0.0, -0.0, 2.180]],
        )
        assert co2_main_out.molecular_mass == 44.0095457
        assert co2_main_out.center_of_mass == [-0.0, 0.0, 0.0]
        assert co2_main_out.moments_of_inertia == [
            -0.3040259e-14,
            0.4185248e02,
            0.4185248e02,
        ]
        assert co2_main_out.rotational_constants == [
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
        assert co2_main_out.total_enthalpy == -10.292932643737
        assert co2_main_out.total_free_energy == -10.317252637172
        assert co2_main_out.gradient_norm == 0.000000274582
        assert co2_main_out.fmo_gap == 8.448655866329


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


class TestXTBG98File:
    """Tests for XTBG98File class."""

    def test_g98_co2(self, xtb_co2_outfolder):
        """Test parsing G98 output from co2 ohess calculation."""
        g98_file = os.path.join(xtb_co2_outfolder, "g98.out")
        assert os.path.exists(g98_file)
        co2_g98 = XTBG98File(g98_file)
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


class TestXTBFolder:
    """Tests for XTBFolder class."""

    def test_folder_co2(self, xtb_co2_outfolder):
        """Test XTBFolder with CO2 ohess calculation output."""
        assert os.path.exists(xtb_co2_outfolder)
        co2_folder = XTBFolder(xtb_co2_outfolder)

        assert co2_folder._xtb_out() is not None
        assert os.path.basename(co2_folder._xtb_out()) == "co2_ohess.out"

        assert co2_folder._xtbopt_log() is not None
        assert os.path.basename(co2_folder._xtbopt_log()) == "xtbopt.log"

        assert co2_folder._charges() is not None
        assert os.path.basename(co2_folder._charges()) == "charges"

        assert co2_folder._energy() is not None
        assert os.path.basename(co2_folder._energy()) == "energy"

        assert co2_folder._engrad() is not None
        assert os.path.basename(co2_folder._engrad()) == "co2.engrad"

        assert co2_folder._g98_out() is not None
        assert os.path.basename(co2_folder._g98_out()) == "g98.out"

        assert co2_folder._gradient() is not None
        assert os.path.basename(co2_folder._gradient()) == "gradient"

        assert co2_folder._hessian() is not None
        assert os.path.basename(co2_folder._hessian()) == "hessian"

        assert co2_folder._vibspectrum() is not None
        assert os.path.basename(co2_folder._vibspectrum()) == "vibspectrum"

        assert co2_folder._wbo() is not None
        assert os.path.basename(co2_folder._wbo()) == "wbo"

        assert co2_folder._xtbopt_geometry() is not None
        assert os.path.basename(co2_folder._xtbopt_geometry()) == "xtbopt.xyz"

        assert co2_folder._xtbtopo_mol() is not None
        assert os.path.basename(co2_folder._xtbtopo_mol()) == "xtbtopo.mol"

    def test_folder_cyclopentadienyl_anion(
        self, xtb_cyclopentadienyl_anion_outfolder
    ):
        """Test XTBFolder with cyclopentadienyl anion opt calculation output."""
        assert os.path.exists(xtb_cyclopentadienyl_anion_outfolder)
        cyclopentadienyl_anion_folder = XTBFolder(
            xtb_cyclopentadienyl_anion_outfolder
        )

        assert cyclopentadienyl_anion_folder._xtb_out() is not None
        assert (
            os.path.basename(cyclopentadienyl_anion_folder._xtb_out())
            == "cyclopentadienyl_anion_opt.out"
        )

        assert cyclopentadienyl_anion_folder._xtbopt_log() is not None
        assert cyclopentadienyl_anion_folder._charges() is not None
        assert (
            cyclopentadienyl_anion_folder._energy() is None
        )  # --grad calculation is not enabled
        assert (
            cyclopentadienyl_anion_folder._engrad() is None
        )  # --grad calculation is not enabled
        assert (
            cyclopentadienyl_anion_folder._g98_out() is None
        )  # --hess calculation is not enabled
        assert (
            cyclopentadienyl_anion_folder._gradient() is None
        )  # --grad calculation is not enabled
        assert (
            cyclopentadienyl_anion_folder._hessian() is None
        )  # --hess calculation is not enabled
        assert (
            cyclopentadienyl_anion_folder._vibspectrum() is None
        )  # --hess calculation is not enabled
        assert cyclopentadienyl_anion_folder._wbo() is not None
        assert (
            os.path.basename(cyclopentadienyl_anion_folder._xtbopt_geometry())
            == "xtbopt.coord"
        )
        assert cyclopentadienyl_anion_folder._xtbtopo_mol() is not None

    def test_folder_p_benzyne_sp(self, xtb_p_benzyne_sp_outfolder):
        """Test XTBFolder with p-benzyne sp calculation output."""
        assert os.path.exists(xtb_p_benzyne_sp_outfolder)
        p_benzyne_sp_folder = XTBFolder(xtb_p_benzyne_sp_outfolder)

        assert p_benzyne_sp_folder._xtb_out() is not None
        assert (
            os.path.basename(p_benzyne_sp_folder._xtb_out())
            == "p_benzyne_sp_alpb_toluene.out"
        )

        assert (
            p_benzyne_sp_folder._xtbopt_log() is None
        )  # no optimization performed
        assert p_benzyne_sp_folder._charges() is not None
        assert (
            p_benzyne_sp_folder._energy() is None
        )  # --grad calculation is not enabled
        assert (
            p_benzyne_sp_folder._engrad() is None
        )  # --grad calculation is not enabled
        assert (
            p_benzyne_sp_folder._g98_out() is None
        )  # --hess calculation is not enabled
        assert (
            p_benzyne_sp_folder._gradient() is None
        )  # --grad calculation is not enabled
        assert (
            p_benzyne_sp_folder._hessian() is None
        )  # --hess calculation is not enabled
        assert (
            p_benzyne_sp_folder._vibspectrum() is None
        )  # --hess calculation is not enabled
        assert p_benzyne_sp_folder._wbo() is not None
        assert (
            p_benzyne_sp_folder._xtbopt_geometry() is None
        )  # no optimization performed
        assert p_benzyne_sp_folder._xtbtopo_mol() is not None


class TestXTBOutput:

    def test_ohess_output(self, xtb_co2_outfolder):
        assert os.path.exists(xtb_co2_outfolder)
        xtb_co2_output = XTBOutput(folder=xtb_co2_outfolder)
        assert xtb_co2_output.normal_termination
        assert xtb_co2_output.partial_charges == {
            "O1": -0.23213972,
            "O2": -0.23213972,
            "C1": 0.46427944,
        }
