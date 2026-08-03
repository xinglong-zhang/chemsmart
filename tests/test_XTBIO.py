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
from chemsmart.utils.constants import kcal_per_mol_to_hartree


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

        assert co2_folder._input_geometry() is not None
        assert os.path.basename(co2_folder._input_geometry()) == "co2.xyz"

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
        assert cyclopentadienyl_anion_folder.is_xtb_calculation_directory

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
            os.path.basename(cyclopentadienyl_anion_folder._input_geometry())
            == "cyclopentadienyl_anion.coord"
        )
        assert (
            os.path.basename(cyclopentadienyl_anion_folder._xtbopt_geometry())
            == "xtbopt.coord"
        )
        assert cyclopentadienyl_anion_folder._xtbtopo_mol() is not None

    def test_folder_p_benzyne_sp(self, xtb_p_benzyne_sp_outfolder):
        """Test XTBFolder with p-benzyne sp calculation output."""
        assert os.path.exists(xtb_p_benzyne_sp_outfolder)
        p_benzyne_sp_folder = XTBFolder(xtb_p_benzyne_sp_outfolder)
        assert p_benzyne_sp_folder.is_xtb_calculation_directory

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
        assert p_benzyne_sp_folder._input_geometry() is not None
        assert p_benzyne_sp_folder._xtbtopo_mol() is not None

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
