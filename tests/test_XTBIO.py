import os.path

import numpy as np

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
