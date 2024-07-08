import math
import numpy as np
import pytest

# from pyatoms.io.ase.atoms import AtomsWrapper
# from pyatoms.io.orca import ORCARefs
# from pyatoms.io.orca.inputs import ORCAInput
# from pyatoms.io.orca.outputs import ORCAEngradFile, ORCAOutput
# from pyatoms.io.orca.route import ORCARoute

from chemsmart.io.orca.route import ORCARoute



class TestORCARoute:
    def test_read_route(self):
        s1 = "!HF DEF2-SVP"
        r1 = ORCARoute(route_string=s1)
        assert r1.route_keywords == ["hf", "def2-svp"]
        assert r1.functional is None
        assert r1.ab_initio == "hf"
        assert r1.basis == "def2-svp"

        s2 = "! Opt Freq M062X def2-TZVP"
        r2 = ORCARoute(route_string=s2)
        assert r2.route_keywords == ["opt", "freq", "m062x", "def2-tzvp"]
        assert r2.functional == "m062x"
        assert r2.basis == "def2-tzvp"

        s3 = "! CPCM(toluene) DLPNO-CCSD(T) Extrapolate(2/3,cc) AutoAux DEFGRID2 TightSCF KDIIS"
        r3 = ORCARoute(route_string=s3)
        assert r3.route_keywords == [
            "cpcm(toluene)",
            "dlpno-ccsd(t)",
            "extrapolate(2/3,cc)",
            "autoaux",
            "defgrid2",
            "tightscf",
            "kdiis",
        ]
        assert r3.functional is None
        assert r3.ab_initio == "dlpno-ccsd(t)"
        assert r3.basis is None
        assert r3.extrapolation_basis == "extrapolate(2/3,cc)"
        assert r3.auxiliary_basis == "autoaux"
        assert r3.defgrid == "defgrid2"
        assert r3.scf_tol == "tight"
        assert r3.scf_algorithm == "kdiis"

        s4 = "! Engrad BP86 D4 STO-3G\n! MiniPrint NoPrintMos NoPop"  # from orca testsuite ORCA_Test_1637.inp
        r4 = ORCARoute(route_string=s4)
        assert r4.route_keywords == [
            "engrad",
            "bp86",
            "d4",
            "sto-3g",
            "miniprint",
            "noprintmos",
            "nopop",
        ]
        assert r4.functional == "bp86"
        assert r4.ab_initio is None
        assert r4.basis == "sto-3g"
        assert r4.extrapolation_basis is None
        assert r4.auxiliary_basis is None

        s5 = "! RHF def2-QZVP SP VeryTightSCF PModel UseSym"
        r5 = ORCARoute(route_string=s5)
        assert r5.route_keywords == [
            "rhf",
            "def2-qzvp",
            "sp",
            "verytightscf",
            "pmodel",
            "usesym",
        ]
        assert r5.functional is None
        assert r5.ab_initio == "rhf"
        assert r5.basis == "def2-qzvp"
        assert r5.extrapolation_basis is None
        assert r5.auxiliary_basis is None
        assert r5.scf_tol == "verytight"

        s6 = "! B3LYP D3ZERO def2-TZVP "
        r6 = ORCARoute(route_string=s6)
        assert r6.route_keywords == ["b3lyp", "d3zero", "def2-tzvp"]
        assert r6.functional == "b3lyp"
        assert r6.dispersion == "d3zero"
        assert r6.ab_initio is None
        assert r6.basis == "def2-tzvp"
        assert r6.extrapolation_basis is None
        assert r6.auxiliary_basis is None


class TestORCABasis:
    def test_orca_all_auxiliary_basis_sets(self):
        orca_ref = ORCARefs()
        orca_bases = orca_ref.orca_all_auxiliary_basis_sets
        assert "def2/j" in orca_bases
        assert "def2/jk" in orca_bases


class TestORCAInput:
    def test_read_water_sp_input(self, water_sp_input_path):
        orca_inp = ORCAInput(inpfile=water_sp_input_path)
        assert orca_inp.route_string == "! hf def2-svp"
        assert orca_inp.route_keywords == ["hf", "def2-svp"]
        assert orca_inp.functional is None
        assert orca_inp.ab_initio == "hf"
        assert orca_inp.basis == "def2-svp"
        assert orca_inp.coordinate_type == "xyz"
        assert orca_inp.charge == 0
        assert orca_inp.multiplicity == 1
        assert orca_inp.natoms == 3
        assert orca_inp.list_of_symbols == ["O", "H", "H"]
        assert isinstance(orca_inp.atoms, AtomsWrapper)
        assert all(orca_inp.atoms.symbols == ["O", "H", "H"])
        assert orca_inp.empirical_formula == "H2O"
        assert orca_inp.scf_maxiter == 500

    def test_read_water_opt_input(self, water_opt_input_path):
        orca_out = ORCAInput(inpfile=water_opt_input_path)
        route_keywords = orca_out.route_keywords
        assert route_keywords == ["opt", "freq", "m062x", "def2-svp"]
        assert orca_out.functional == "m062x"
        assert orca_out.basis == "def2-svp"
        assert orca_out.coordinate_type == "xyz"
        assert orca_out.charge == 0
        assert orca_out.multiplicity == 1
        assert orca_out.natoms == 3
        assert orca_out.list_of_symbols == ["O", "H", "H"]
        assert isinstance(orca_out.atoms, AtomsWrapper)
        assert all(orca_out.atoms.symbols == ["O", "H", "H"])
        assert orca_out.empirical_formula == "H2O"

    def test_read_solvent(self, orca_epr_solv):
        orca_inp = ORCAInput(inpfile=orca_epr_solv)
        assert orca_inp.functional == "b3lyp"
        assert orca_inp.basis == "6-311++g(2d,2p)"
        assert orca_inp.aux_basis == "def2/jk"
        assert orca_inp.scf_tol == "extreme"
        assert orca_inp.solvent_model == "smd"
        assert orca_inp.solvent_id == "water"

    def test_orca_faulty_solvent(self, orca_faulty_solv):
        orca_inp = ORCAInput(inpfile=orca_faulty_solv)
        assert orca_inp.functional == "b3lyp"
        assert orca_inp.basis == "6-311++g(2d,2p)"
        assert orca_inp.aux_basis == "def2/jk"
        assert orca_inp.scf_tol == "extreme"
        assert orca_inp.solvent_model == "smd"
        with pytest.raises(
            Exception,
            match="Your input file specifies solvent but solvent is not in quotes, "
            "thus, your input file is not valid to run for ORCA!",
        ):
            orca_inp.solvent_id  # noqa: B018


class TestORCAOutput:
    def test_read_water_output(self, water_output_gas_path):
        orca_out = ORCAOutput(outputfile=water_output_gas_path)
        assert isinstance(orca_out.atoms, AtomsWrapper)
        assert orca_out.route_string == "! opt freq m062x def2-svp"
        assert orca_out.functional == "m062x"
        assert orca_out.basis == "def2-svp"
        assert orca_out.ab_initio is None
        assert orca_out.aux_basis is None
        assert orca_out.extrapolation_basis is None
        assert orca_out.natoms == 3
        assert orca_out.num_basis_functions == 24
        assert orca_out.num_shells == 12
        assert orca_out.max_ang_mom == 2
        assert orca_out.contraction_scheme.lower() == "segmented contraction"
        assert orca_out.coulomb_range_seperation.lower() == "not used"
        assert orca_out.exchange_range_seperation.lower() == "not used"
        assert orca_out.finite_nucleus_model.lower() == "not used"
        assert orca_out.aux_j_fitting_basis.lower() == "available"
        assert orca_out.aux_j_num_basis_functions == 71
        assert orca_out.aux_j_num_shells == 25
        assert orca_out.aux_j_max_ang_mom == 4
        assert orca_out.aux_jk_fitting_basis.lower() == "not available"
        assert orca_out.aux_k_fitting_basis.lower() == "not available"
        assert orca_out.aux_external_fitting_basis.lower() == "not available"
        assert orca_out.integral_threshold == 2.5e-11
        assert orca_out.primitive_cutoff == 2.5e-12
        assert orca_out.primitive_pair_threshold == 2.5e-12
        assert orca_out.ri_approx is True
        assert orca_out.rij_cosx is True
        assert orca_out.charge == 0
        assert orca_out.multiplicity == 1
        assert orca_out.num_electrons == 10
        assert orca_out.basis_dim == 24
        assert orca_out.diis_acceleration is True
        assert orca_out.scf_maxiter == 125
        assert orca_out.converged is True
        assert orca_out.normal_termination is True

    def test_water_optimized_output(self, water_output_gas_path):  # noqa: PLR0915
        orca_out = ORCAOutput(outputfile=water_output_gas_path)
        optimized_geometry = orca_out.get_optimized_parameters()
        assert optimized_geometry == {
            "B(H1,O0)": 0.9627,
            "B(H2,O0)": 0.9627,
            "A(H1,O0,H2)": 103.35,
        }
        atoms = orca_out.final_structure
        assert isinstance(atoms, AtomsWrapper)
        assert all(atoms.symbols == ["O", "H", "H"])
        assert orca_out.empirical_formula == "H2O"
        assert np.allclose(
            orca_out.optimized_geometry,
            np.array(
                [
                    [0.0, 0.0, 0.087341],
                    [-0.755205, 0.0, -0.50967],
                    [0.755205, 0.0, -0.50967],
                ]
            ),
            rtol=1e-03,
        )
        assert math.isclose(orca_out.final_energy, -2076.86288, rel_tol=1e-4)
        assert math.isclose(orca_out.final_nuclear_repulsion, 248.85900, rel_tol=1e-4)
        assert math.isclose(orca_out.final_electronic_energy, -2325.72188, rel_tol=1e-4)
        assert math.isclose(orca_out.one_electron_energy, -3346.95900, rel_tol=1e-4)
        assert math.isclose(orca_out.two_electron_energy, 1021.23712, rel_tol=1e-4)
        assert math.isclose(orca_out.max_cosx_asymmetry_energy, 0.00023, rel_tol=1e-2)
        assert math.isclose(orca_out.potential_energy, -4140.24210, rel_tol=1e-4)
        assert math.isclose(orca_out.kinetic_energy, 2063.37922, rel_tol=1e-4)
        assert math.isclose(orca_out.virial_ratio, 2.00653475, rel_tol=1e-4)
        assert math.isclose(orca_out.xc_energy, -121.97818321, rel_tol=1e-4)
        assert orca_out.dfet_embed_energy == 0.0
        assert orca_out.orbital_occupancy == [
            2,
            2,
            2,
            2,
            2,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
        orbital_energies = np.array(
            [
                -533.6944922825586,
                -29.80792369019839,
                -15.936838030119114,
                -12.224062098182399,
                -10.252651603489042,
                2.419962981919028,
                4.519647950331253,
                16.859657764363483,
                18.723882609506855,
                27.008008125221124,
                27.26657071522466,
                29.185381600732917,
                31.307570385387297,
                37.97574774204451,
                39.35882085950502,
                43.269233088136716,
                49.80502705178539,
                62.011537805982066,
                63.17398100555701,
                82.40358722432674,
                83.45748420505049,
                88.40867752634217,
                96.38705590868665,
                105.15657812830358,
            ]
        )
        assert np.allclose(orca_out.orbital_energies, orbital_energies, rtol=1e-6)

        assert orca_out.mulliken_atomic_charges == {
            "O0": -0.32926,
            "H1": 0.16463,
            "H2": 0.16463,
        }
        assert orca_out.loewdin_atomic_charges == {
            "O0": -0.155184,
            "H1": 0.077592,
            "H2": 0.077592,
        }
        assert orca_out.mayer_mulliken_gross_atomic_population == {
            "O0": 8.3293,
            "H1": 0.8354,
            "H2": 0.8354,
        }
        assert orca_out.mayer_total_nuclear_charge == {"O0": 8.0, "H1": 1.0, "H2": 1.0}
        assert orca_out.mayer_mulliken_gross_atomic_charge == {
            "O0": -0.3293,
            "H1": 0.1646,
            "H2": 0.1646,
        }
        assert orca_out.mayer_total_valence == {
            "O0": 2.0072,
            "H1": 1.0104,
            "H2": 1.0104,
        }
        assert orca_out.mayer_bonded_valence == {
            "O0": 2.0072,
            "H1": 1.0104,
            "H2": 1.0104,
        }
        assert orca_out.mayer_free_valence == {"H1": 0.0, "H2": 0.0, "O0": 0.0}
        assert orca_out.mayer_bond_orders_larger_than_zero_point_one == {
            "B(O0,H1)": 1.0036,
            "B(O0,H2)": 1.0036,
        }
        assert np.allclose(
            orca_out.dipole_moment_electric_contribution,
            np.array([0.0, 0.0, 0.18466]),
            rtol=1e-4,
        )
        assert np.allclose(
            orca_out.dipole_moment_nuclear_contribution,
            np.array([0.0, 0.0, -0.99386]),
            rtol=1e-4,
        )
        assert np.allclose(
            orca_out.total_dipole_moment, np.array([0.0, 0.0, -0.8092]), rtol=1e-4
        )
        assert orca_out.dipole_moment_in_au == 0.80920
        assert orca_out.dipole_moment_in_debye == 2.05682
        assert np.allclose(
            orca_out.dipole_moment_along_axis_in_au,
            np.array([0.0, -0.809197, 0.0]),
            rtol=1e-4,
        )
        assert np.allclose(
            orca_out.dipole_moment_along_axis_in_debye,
            np.array([0.0, -2.056815, 0.0]),
            rtol=1e-4,
        )
        assert orca_out.rotational_constants_in_wavenumbers == [
            26.416987,
            14.661432,
            9.428573,
        ]
        assert orca_out.rotational_constants_in_MHz == [
            791961.33697,
            439538.666271,
            282661.493198,
        ]
        assert orca_out.vibrational_frequencies == [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1625.35,
            3875.61,
            3971.9,
        ]
        assert orca_out.vib_freq_scale_factor == 1.0
        assert orca_out.molar_absorption_coefficients == [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.012719,
            0.002968,
            0.009899,
        ]
        assert orca_out.integrated_absorption_coefficients == [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            64.27,
            15.0,
            50.03,
        ]
        assert orca_out.transition_dipole_deriv_norm == [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.002442,
            0.000239,
            0.000778,
        ]
        assert orca_out.num_translation_and_rotation_modes == 6
        assert orca_out.num_vibration_modes == 3
        assert orca_out.temperature_in_K == 298.15
        assert orca_out.pressure_in_atm == 1.0
        assert orca_out.total_mass_in_amu == 18.02
        assert math.isclose(orca_out.internal_energy, -2076.198522987, rel_tol=1e-4)
        assert math.isclose(orca_out.electronic_energy, -2076.8629218525, rel_tol=1e-4)
        assert math.isclose(orca_out.zero_point_energy, 0.587242346752, rel_tol=1e-4)
        assert math.isclose(
            orca_out.thermal_vibration_correction, 7.91851273564e-5, rel_tol=1e-6
        )
        assert math.isclose(
            orca_out.thermal_rotation_correction, 0.038538666777, rel_tol=1e-4
        )
        assert math.isclose(
            orca_out.thermal_translation_correction, 0.038538666777, rel_tol=1e-4
        )
        assert math.isclose(orca_out.zero_point_energy, 0.587242346752, rel_tol=1e-4)
        assert math.isclose(orca_out.zero_point_energy, 0.587242346752, rel_tol=1e-4)
        assert math.isclose(orca_out.enthalpy, -2076.1728297262, rel_tol=1e-4)
        assert math.isclose(
            orca_out.thermal_enthalpy_correction, 0.02569326085952, rel_tol=1e-4
        )
        assert orca_out.electronic_entropy_no_temperature_in_SI == 0.0
        assert math.isclose(
            orca_out.vibrational_entropy_no_temperature_in_SI, 0.028883578, rel_tol=1e-4
        )
        assert math.isclose(
            orca_out.rotational_entropy_no_temperature_in_SI, 43.887276051, rel_tol=1e-4
        )
        assert math.isclose(
            orca_out.translational_entropy_no_temperature_in_SI,
            144.8035920,
            rel_tol=1e-4,
        )
        assert math.isclose(orca_out.entropy_TS, 0.5831641766362, rel_tol=1e-4)

        assert orca_out.mulliken_atomic_charges == {
            "O0": -0.32926,
            "H1": 0.16463,
            "H2": 0.16463,
        }
        assert orca_out.loewdin_atomic_charges == {
            "O0": -0.155184,
            "H1": 0.077592,
            "H2": 0.077592,
        }
        assert orca_out.mayer_mulliken_gross_atomic_population == {
            "O0": 8.3293,
            "H1": 0.8354,
            "H2": 0.8354,
        }
        assert orca_out.mayer_total_nuclear_charge == {"O0": 8.0, "H1": 1.0, "H2": 1.0}
        assert orca_out.mayer_mulliken_gross_atomic_charge == {
            "O0": -0.3293,
            "H1": 0.1646,
            "H2": 0.1646,
        }
        assert orca_out.mayer_total_valence == {
            "O0": 2.0072,
            "H1": 1.0104,
            "H2": 1.0104,
        }
        assert orca_out.mayer_bonded_valence == {
            "O0": 2.0072,
            "H1": 1.0104,
            "H2": 1.0104,
        }
        assert orca_out.mayer_free_valence == {"H1": 0.0, "H2": 0.0, "O0": 0.0}
        assert orca_out.mayer_bond_orders_larger_than_zero_point_one == {
            "B(O0,H1)": 1.0036,
            "B(O0,H2)": 1.0036,
        }
        assert np.allclose(
            orca_out.dipole_moment_electric_contribution,
            np.array([0.0, 0.0, 0.18466]),
            rtol=1e-4,
        )
        assert np.allclose(
            orca_out.dipole_moment_nuclear_contribution,
            np.array([0.0, 0.0, -0.99386]),
            rtol=1e-4,
        )
        assert np.allclose(
            orca_out.total_dipole_moment, np.array([0.0, 0.0, -0.8092]), rtol=1e-4
        )
        assert orca_out.dipole_moment_in_au == 0.80920
        assert orca_out.dipole_moment_in_debye == 2.05682
        assert np.allclose(
            orca_out.dipole_moment_along_axis_in_au,
            np.array([0.0, -0.809197, 0.0]),
            rtol=1e-4,
        )
        assert np.allclose(
            orca_out.dipole_moment_along_axis_in_debye,
            np.array([0.0, -2.056815, 0.0]),
            rtol=1e-4,
        )
        assert orca_out.rotational_constants_in_wavenumbers == [
            26.416987,
            14.661432,
            9.428573,
        ]
        assert orca_out.rotational_constants_in_MHz == [
            791961.33697,
            439538.666271,
            282661.493198,
        ]
        assert orca_out.vibrational_frequencies == [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1625.35,
            3875.61,
            3971.9,
        ]
        assert orca_out.vib_freq_scale_factor == 1.0
        assert orca_out.molar_absorption_coefficients == [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.012719,
            0.002968,
            0.009899,
        ]
        assert orca_out.integrated_absorption_coefficients == [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            64.27,
            15.0,
            50.03,
        ]
        assert orca_out.transition_dipole_deriv_norm == [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.002442,
            0.000239,
            0.000778,
        ]
        assert orca_out.num_translation_and_rotation_modes == 6
        assert orca_out.num_vibration_modes == 3
        assert orca_out.temperature_in_K == 298.15
        assert orca_out.pressure_in_atm == 1.0
        assert orca_out.total_mass_in_amu == 18.02
        assert math.isclose(orca_out.internal_energy, -2076.198522987, rel_tol=1e-4)
        assert math.isclose(orca_out.electronic_energy, -2076.8629218525, rel_tol=1e-4)
        assert math.isclose(orca_out.zero_point_energy, 0.587242346752, rel_tol=1e-4)
        assert math.isclose(
            orca_out.thermal_vibration_correction, 7.91851273564e-5, rel_tol=1e-6
        )
        assert math.isclose(
            orca_out.thermal_rotation_correction, 0.038538666777, rel_tol=1e-4
        )
        assert math.isclose(
            orca_out.thermal_translation_correction, 0.038538666777, rel_tol=1e-4
        )
        assert math.isclose(orca_out.zero_point_energy, 0.587242346752, rel_tol=1e-4)
        assert math.isclose(orca_out.zero_point_energy, 0.587242346752, rel_tol=1e-4)
        assert math.isclose(orca_out.enthalpy, -2076.1728297262, rel_tol=1e-4)
        assert math.isclose(
            orca_out.thermal_enthalpy_correction, 0.02569326085952, rel_tol=1e-4
        )
        assert orca_out.electronic_entropy_no_temperature_in_SI == 0.0
        assert math.isclose(
            orca_out.vibrational_entropy_no_temperature_in_SI, 0.028883578, rel_tol=1e-4
        )
        assert math.isclose(
            orca_out.rotational_entropy_no_temperature_in_SI, 43.887276051, rel_tol=1e-4
        )
        assert math.isclose(
            orca_out.translational_entropy_no_temperature_in_SI,
            144.8035920,
            rel_tol=1e-4,
        )
        assert math.isclose(orca_out.entropy_TS, 0.5831641766362, rel_tol=1e-4)

        entropy_TS_in_J_per_mol = 56266.786841951627
        TS = orca_out.entropy_in_J_per_mol_per_K * 298.15  # 56266.794 J/mol
        # converted to 13.448 kcal/mol, as expected from output
        assert math.isclose(entropy_TS_in_J_per_mol, TS, rel_tol=1e-4)
        assert orca_out.rotational_entropy_symmetry_correction_J_per_mol_per_K == {
            1: 49.65043,
            2: 43.887276,
            3: 40.516087,
            4: 38.124122,
            5: 36.268792,
            6: 34.752933,
            7: 33.471224,
            8: 32.361055,
            9: 31.381743,
            10: 30.505726,
            11: 29.713277,
            12: 28.989778,
        }

        assert math.isclose(
            orca_out.gibbs_free_energy, -2076.7559939028233202, rel_tol=1e-6
        )
        assert isinstance(orca_out.atoms, AtomsWrapper)
        assert orca_out.total_run_time_hours == 0.0028

    def test_read_sp_output(self, water_sp_gas_path):
        orca_out = ORCAOutput(outputfile=water_sp_gas_path)
        assert (
            orca_out.route_string
            == "!  DLPNO-CCSD(T) Extrapolate(2/3,cc) AutoAux DEFGRID3 TightSCF KDIIS".lower()
        )
        assert orca_out.functional is None
        assert orca_out.basis is None
        assert orca_out.ab_initio == "DLPNO-CCSD(T)".lower()
        assert orca_out.aux_basis == "AutoAux".lower()
        assert orca_out.extrapolation_basis == "Extrapolate(2/3,cc)".lower()
        assert orca_out.natoms == 3
        assert orca_out.num_basis_functions == 30
        assert orca_out.num_shells == 18
        assert orca_out.max_ang_mom == 2
        assert orca_out.contraction_scheme == "PARTIAL GENERAL contraction"
        assert orca_out.coulomb_range_seperation.lower() == "not used"
        assert orca_out.exchange_range_seperation.lower() == "not used"
        assert orca_out.finite_nucleus_model.lower() == "not used"
        assert orca_out.aux_j_fitting_basis.lower() == "available"
        assert orca_out.aux_j_num_basis_functions == 113
        assert orca_out.aux_j_num_shells == 37
        assert orca_out.aux_j_max_ang_mom == 4
        assert orca_out.aux_jk_fitting_basis.lower() == "available"
        assert orca_out.aux_k_fitting_basis.lower() == "available"
        assert orca_out.aux_external_fitting_basis.lower() == "not available"
        assert orca_out.integral_threshold == 2.5e-11
        assert orca_out.primitive_cutoff == 2.5e-12
        assert orca_out.primitive_pair_threshold == 2.5e-12
        assert orca_out.ri_approx is None
        assert orca_out.rij_cosx is None
        assert orca_out.charge == 0
        assert orca_out.multiplicity == 1
        assert orca_out.num_electrons == 10
        assert orca_out.basis_dim == 24
        assert orca_out.diis_acceleration is False
        assert orca_out.scf_maxiter == 500
        assert orca_out.converged is None
        assert isinstance(orca_out.atoms, AtomsWrapper)
        assert orca_out.normal_termination is True

    def test_gtoint_errfile(self, gtoint_errfile):
        orca_out = ORCAOutput(outputfile=gtoint_errfile)
        assert orca_out.route_string == "! m062x def2-svp opt freq defgrid3"
        assert orca_out.functional == "m062x"
        assert orca_out.basis == "def2-svp"
        assert orca_out.defgrid == "defgrid3"
        assert orca_out.ab_initio is None
        assert orca_out.aux_basis is None
        assert orca_out.extrapolation_basis is None
        assert orca_out.natoms == 27
        assert orca_out.normal_termination is False


class TestORCAEngrad:
    def test_read_water_output(self, water_engrad_path):
        orca_engrad = ORCAEngradFile(engrad_file=water_engrad_path)
        assert orca_engrad.natoms == 3
        assert math.isclose(
            orca_engrad.energy, -2076.86288, rel_tol=1e-4
        )  # energy in eV
        gradient = np.array(
            [
                [3.15988602e-07, -2.28910474e-06, -2.40955630e-04],
                [7.40004117e-04, 7.35335559e-08, 4.97007134e-04],
                [-7.38385711e-04, 6.79285506e-08, 4.96422516e-04],
            ]
        )
        assert np.allclose(orca_engrad.gradient, gradient, rtol=1e-6)
        assert isinstance(orca_engrad.atoms, AtomsWrapper)
        assert orca_engrad.atoms.chemical_symbols == ["O", "H", "H"]
        coordinates = np.array(
            [
                [-0.0, 0.0, 0.08734059],
                [-0.75520525, 0.0, -0.50967031],
                [0.75520525, 0.0, -0.50967031],
            ]
        )
        assert np.allclose(orca_engrad.atoms.positions, coordinates, rtol=1e-6)
