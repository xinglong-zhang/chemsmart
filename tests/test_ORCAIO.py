import math
import os
from shutil import copy

import numpy as np
import pytest
from ase import units

from chemsmart.io.molecules.structure import CoordinateBlock, Molecule
from chemsmart.io.orca import ORCARefs
from chemsmart.io.orca.input import ORCAInput
from chemsmart.io.orca.output import ORCAEngradFile, ORCAOutput
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
        assert os.path.exists(water_sp_input_path)
        orca_inp = ORCAInput(filename=water_sp_input_path)
        assert orca_inp.route_string == "! hf def2-svp"
        assert orca_inp.route_keywords == ["hf", "def2-svp"]
        assert orca_inp.functional is None
        assert orca_inp.ab_initio == "hf"
        assert orca_inp.basis == "def2-svp"
        assert orca_inp.coordinate_type == "xyz"
        assert orca_inp.charge == 0
        assert orca_inp.multiplicity == 1
        assert orca_inp.molecule.num_atoms == 3
        assert orca_inp.molecule.chemical_symbols == ["O", "H", "H"]
        assert isinstance(orca_inp.molecule, Molecule)
        assert all(orca_inp.molecule.symbols == ["O", "H", "H"])
        assert orca_inp.molecule.empirical_formula == "H2O"
        assert orca_inp.scf_maxiter == 500

    def test_read_water_opt_input(self, water_opt_input_path):
        orca_inp = ORCAInput(filename=water_opt_input_path)
        route_keywords = orca_inp.route_keywords
        assert route_keywords == ["opt", "freq", "m062x", "def2-svp"]
        assert orca_inp.functional == "m062x"
        assert orca_inp.basis == "def2-svp"
        assert orca_inp.coordinate_type == "xyz"
        assert orca_inp.charge == 0
        assert orca_inp.multiplicity == 1
        assert orca_inp.molecule.num_atoms == 3
        assert orca_inp.molecule.chemical_symbols == ["O", "H", "H"]
        assert isinstance(orca_inp.molecule, Molecule)
        assert all(orca_inp.molecule.symbols == ["O", "H", "H"])
        assert orca_inp.molecule.empirical_formula == "H2O"

    def test_read_solvent(self, orca_epr_solv):
        orca_inp = ORCAInput(filename=orca_epr_solv)
        assert orca_inp.functional == "b3lyp"
        assert orca_inp.basis == "6-311++g(2d,2p)"
        assert orca_inp.aux_basis == "def2/jk"
        assert orca_inp.scf_tol == "extreme"
        assert orca_inp.solvent_model == "smd"
        assert orca_inp.solvent_id == "water"

    def test_orca_faulty_solvent(self, orca_faulty_solv):
        orca_inp = ORCAInput(filename=orca_faulty_solv)
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
            orca_inp.solvent_id

    def test_orca_input_with_xyz_files_specified(
        self,
        tmpdir,
        orca_input_nebts_file,
        orca_input_nebts_reactant_xyz_file,
        orca_input_nebts_product_xyz_file,
        orca_input_nebts_ts_xyz_file,
    ):
        # copy all files to tmpdir
        orca_input_nebts_file_tmp = os.path.join(
            tmpdir, "orca_input_nebts.inp"
        )
        orca_input_nebts_reactant_xyz_file_tmp = os.path.join(
            tmpdir, "reactant.xyz"
        )
        orca_input_nebts_product_xyz_file_tmp = os.path.join(
            tmpdir, "product.xyz"
        )
        orca_input_nebts_ts_xyz_file_tmp = os.path.join(tmpdir, "ts.xyz")
        copy(orca_input_nebts_file, orca_input_nebts_file_tmp)
        copy(
            orca_input_nebts_reactant_xyz_file,
            orca_input_nebts_reactant_xyz_file_tmp,
        )
        copy(
            orca_input_nebts_product_xyz_file,
            orca_input_nebts_product_xyz_file_tmp,
        )
        copy(orca_input_nebts_ts_xyz_file, orca_input_nebts_ts_xyz_file_tmp)
        assert os.path.exists(orca_input_nebts_file_tmp)
        assert os.path.exists(orca_input_nebts_reactant_xyz_file_tmp)
        assert os.path.exists(orca_input_nebts_product_xyz_file_tmp)
        assert os.path.exists(orca_input_nebts_ts_xyz_file_tmp)

        orca_inp = ORCAInput(filename=orca_input_nebts_file)
        assert orca_inp.route_string == "!  GFN2-xTB NEB-TS Freq".lower()
        assert orca_inp.functional is None
        assert orca_inp.basis is None
        assert orca_inp.coordinate_type == "xyzfile"  # xyzfile is specified
        assert orca_inp.charge == 0
        assert orca_inp.multiplicity == 1
        assert orca_inp.molecule.num_atoms == 40
        assert isinstance(orca_inp.molecule, Molecule)
        assert orca_inp.molecule.empirical_formula == "C23H15NO"


class TestORCAOutput:
    def test_read_water_output(self, water_output_gas_path):
        orca_out = ORCAOutput(filename=water_output_gas_path)
        assert isinstance(orca_out.molecule, Molecule)
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
        assert orca_out.has_forces
        orca_out_first_forces = np.array(
            [
                [-0.000000001, 0.000000112, -0.002606324],
                [-0.013427516, 0.000000001, 0.001404818],
                [0.013427527, 0.000000001, 0.001404820],
            ]
        )

        assert np.allclose(
            orca_out.forces[0], orca_out_first_forces, rtol=1e-6
        )
        orca_out_last_forces = np.array(
            [
                [0.000000006, -0.000000045, -0.000004686],
                [0.000014391, 0.000000001, 0.000009665],
                [-0.000014359, 0.000000001, 0.000009654],
            ]
        )
        assert np.allclose(
            orca_out.forces[-1], orca_out_last_forces, rtol=1e-6
        )

        orca_out_first_forces_in_eV_per_angstrom = np.array(
            [
                np.array([-0.000000001, 0.000000112, -0.002606324])
                * units.Hartree
                / units.Bohr,
                np.array([-0.013427516, 0.000000001, 0.001404818])
                * units.Hartree
                / units.Bohr,
                np.array([0.013427527, 0.000000001, 0.001404820])
                * units.Hartree
                / units.Bohr,
            ]
        )
        assert np.allclose(
            orca_out.forces_in_eV_per_angstrom[0],
            orca_out_first_forces_in_eV_per_angstrom,
            rtol=1e-6,
        )

        assert isinstance(orca_out.input_coordinates_block, CoordinateBlock)
        molecule = orca_out.input_coordinates_block.molecule
        assert isinstance(molecule, Molecule)
        assert all(molecule.symbols == ["O", "H", "H"])
        expected_coordinate_block = [
            "O  0.0000  0.0000  0.0626",
            "H  -0.7920  0.0000  -0.4973",
            "H  0.7920  0.0000  -0.4973",
        ]
        orca_coordinate_block = (
            orca_out.input_coordinates_block.coordinate_block
        )
        for i, line in enumerate(orca_coordinate_block):
            assert len(line.split()) == len(
                expected_coordinate_block[i].split()
            )
            assert line.split()[0] == expected_coordinate_block[i].split()[0]
            assert math.isclose(
                float(line.split()[1]),
                float(expected_coordinate_block[i].split()[1]),
                rel_tol=1e-4,
            )
            assert math.isclose(
                float(line.split()[2]),
                float(expected_coordinate_block[i].split()[2]),
                rel_tol=1e-4,
            )
            assert math.isclose(
                float(line.split()[3]),
                float(expected_coordinate_block[i].split()[3]),
                rel_tol=1e-4,
            )

        assert orca_out.molecule.empirical_formula == "H2O"
        assert len(orca_out.energies) == 6
        assert orca_out.energies[0] == -76.322282695198
        assert orca_out.energies[-1] == -76.323311011349
        assert orca_out.total_core_hours == orca_out.total_service_unit == 0.0
        assert orca_out.total_elapsed_walltime == 0.0

    def test_water_optimized_output(self, water_output_gas_path):
        orca_out = ORCAOutput(filename=water_output_gas_path)
        assert orca_out.forces is not None
        optimized_geometry = orca_out.get_optimized_parameters()
        assert optimized_geometry == {
            "B(H2,O1)": 0.9627,
            "B(H3,O1)": 0.9627,
            "A(H2,O1,H3)": 103.35,
        }
        molecule = orca_out.final_structure
        assert isinstance(molecule, Molecule)
        assert molecule.symbols == ["O", "H", "H"]
        assert orca_out.molecule.empirical_formula == "H2O"
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
        assert math.isclose(orca_out.final_energy, -76.32331101, rel_tol=1e-4)
        assert math.isclose(
            orca_out.final_nuclear_repulsion, 9.14540060, rel_tol=1e-4
        )
        assert math.isclose(
            orca_out.final_electronic_energy, -85.46871161, rel_tol=1e-4
        )
        assert math.isclose(
            orca_out.one_electron_energy, -122.99848749, rel_tol=1e-4
        )
        assert math.isclose(
            orca_out.two_electron_energy, 37.52977588, rel_tol=1e-4
        )
        assert math.isclose(
            orca_out.max_cosx_asymmetry_energy, 0.00000841, rel_tol=1e-2
        )
        assert math.isclose(
            orca_out.potential_energy, -152.15110671, rel_tol=1e-4
        )
        assert math.isclose(orca_out.kinetic_energy, 75.82779570, rel_tol=1e-4)
        assert math.isclose(orca_out.virial_ratio, 2.00653475, rel_tol=1e-4)
        assert math.isclose(orca_out.xc_energy, -4.482615927956, rel_tol=1e-4)
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
        assert np.allclose(
            orca_out.orbital_energies, orbital_energies, rtol=1e-6
        )
        assert orca_out.homo_energy == -10.252651603489042
        assert orca_out.lumo_energy == 2.419962981919028
        assert np.isclose(orca_out.fmo_gap, 12.6726145854, rtol=1e-6)

        # test HOMO LUMO

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
        assert orca_out.mayer_total_nuclear_charge == {
            "O0": 8.0,
            "H1": 1.0,
            "H2": 1.0,
        }
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
            orca_out.total_dipole_moment,
            np.array([0.0, 0.0, -0.8092]),
            rtol=1e-4,
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
            1625.35,
            3875.61,
            3971.9,
        ]
        assert orca_out.vib_freq_scale_factor == 1.0
        assert orca_out.molar_absorption_coefficients == [
            0.012719,
            0.002968,
            0.009899,
        ]
        assert orca_out.integrated_absorption_coefficients == [
            64.27,
            15.0,
            50.03,
        ]
        assert orca_out.transition_dipole_deriv_norm == [
            0.002442,
            0.000239,
            0.000778,
        ]
        assert orca_out.num_translation_and_rotation_modes == 6
        assert orca_out.num_vibration_modes == 3
        assert orca_out.temperature_in_K == 298.15
        assert orca_out.pressure_in_atm == 1.0
        assert orca_out.total_mass_in_amu == 18.02
        assert math.isclose(
            orca_out.internal_energy, -76.29889480, rel_tol=1e-4
        )  # in Hartrees, default unit in ORCA output file
        assert math.isclose(
            orca_out.electronic_energy, -76.32331101, rel_tol=1e-4
        )
        assert math.isclose(
            orca_out.zero_point_energy, 0.02158076, rel_tol=1e-8
        )
        assert math.isclose(
            orca_out.thermal_vibration_correction,
            0.00000291,
            rel_tol=1e-8,
        )
        assert math.isclose(
            orca_out.thermal_rotation_correction, 0.00141627, rel_tol=1e-8
        )
        assert math.isclose(
            orca_out.thermal_translation_correction,
            0.00141627,
            rel_tol=1e-8,
        )
        assert math.isclose(orca_out.enthalpy, -76.29795059, rel_tol=1e-4)
        assert math.isclose(
            orca_out.thermal_enthalpy_correction,
            0.00094421,
            rel_tol=1e-8,
        )
        assert orca_out.electronic_entropy_no_temperature_in_SI == 0.0
        assert math.isclose(
            orca_out.vibrational_entropy_no_temperature_in_SI,
            0.028883578,
            rel_tol=1e-4,
        )
        assert math.isclose(
            orca_out.rotational_entropy_no_temperature_in_SI,
            43.887276051,
            rel_tol=1e-4,
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
        assert orca_out.mayer_total_nuclear_charge == {
            "O0": 8.0,
            "H1": 1.0,
            "H2": 1.0,
        }
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
            orca_out.total_dipole_moment,
            np.array([0.0, 0.0, -0.8092]),
            rtol=1e-4,
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
            1625.35,
            3875.61,
            3971.9,
        ]
        assert orca_out.vib_freq_scale_factor == 1.0
        assert orca_out.molar_absorption_coefficients == [
            0.012719,
            0.002968,
            0.009899,
        ]
        assert orca_out.integrated_absorption_coefficients == [
            64.27,
            15.0,
            50.03,
        ]
        assert orca_out.transition_dipole_deriv_norm == [
            0.002442,
            0.000239,
            0.000778,
        ]
        assert orca_out.num_translation_and_rotation_modes == 6
        assert orca_out.num_vibration_modes == 3
        assert orca_out.temperature_in_K == 298.15
        assert orca_out.pressure_in_atm == 1.0
        assert orca_out.total_mass_in_amu == 18.02

        entropy_TS_in_J_per_mol = 56266.786841951627
        TS = orca_out.entropy_in_J_per_mol_per_K * 298.15  # 56266.794 J/mol
        # converted to 13.448 kcal/mol, as expected from output
        assert math.isclose(entropy_TS_in_J_per_mol, TS, rel_tol=1e-4)
        assert (
            orca_out.rotational_entropy_symmetry_correction_J_per_mol_per_K
            == {
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
        )

        assert math.isclose(
            orca_out.gibbs_free_energy, -76.31938148, rel_tol=1e-8
        )
        assert isinstance(orca_out.molecule, Molecule)
        assert orca_out.total_elapsed_walltime == 0.0

    def test_read_sp_output(self, water_sp_gas_path):
        orca_out = ORCAOutput(filename=water_sp_gas_path)
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
        assert isinstance(orca_out.molecule, Molecule)
        assert orca_out.normal_termination is True
        assert orca_out.total_elapsed_walltime == 0.0

    def test_read_sp_full_print_output(self, dlpno_ccsdt_sp_full_print):
        orca_out = ORCAOutput(filename=dlpno_ccsdt_sp_full_print)
        assert (
            orca_out.route_string
            == "! CPCM DLPNO-CCSD(T) Extrapolate(2/3,cc) AutoAux DEFGRID3 TightSCF KDIIS".lower()
        )
        assert orca_out.functional is None
        assert orca_out.basis is None
        assert orca_out.ab_initio == "DLPNO-CCSD(T)".lower()
        assert orca_out.aux_basis == "AutoAux".lower()
        assert orca_out.extrapolation_basis == "Extrapolate(2/3,cc)".lower()
        assert orca_out.natoms == 78
        assert orca_out.num_basis_functions == 1200
        assert orca_out.num_shells == 720
        assert orca_out.max_ang_mom == 2
        assert orca_out.contraction_scheme == "PARTIAL GENERAL contraction"
        assert orca_out.coulomb_range_seperation.lower() == "not used"
        assert orca_out.exchange_range_seperation.lower() == "not used"
        assert orca_out.finite_nucleus_model.lower() == "not used"
        assert orca_out.aux_j_fitting_basis.lower() == "available"
        assert orca_out.aux_j_num_basis_functions == 4518
        assert orca_out.aux_j_num_shells == 1494
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
        assert orca_out.num_electrons == 396
        assert orca_out.basis_dim == 876
        assert orca_out.diis_acceleration is False
        assert orca_out.scf_maxiter == 500
        assert orca_out.converged is None
        assert isinstance(orca_out.molecule, Molecule)
        assert orca_out.normal_termination is True
        assert np.isclose(orca_out.homo_energy, -9.40597, rtol=1e-5)
        assert np.isclose(orca_out.lumo_energy, 1.56267, rtol=1e-5)
        assert np.isclose(orca_out.fmo_gap, 10.968637, rtol=1e-5)

    def test_read_hirshfeld_full_print_output(self, hirshfeld_full_print):
        orca_out = ORCAOutput(filename=hirshfeld_full_print)
        assert orca_out.route_string == "! Hirshfeld".lower()
        assert orca_out.functional is None
        assert orca_out.basis is None
        assert orca_out.ab_initio is None
        assert orca_out.aux_basis is None
        assert orca_out.extrapolation_basis is None
        assert orca_out.natoms == 78
        assert orca_out.num_basis_functions == 876
        assert orca_out.num_shells == 396
        assert orca_out.max_ang_mom == 2
        assert orca_out.contraction_scheme == "SEGMENTED contraction"
        assert orca_out.coulomb_range_seperation.lower() == "not used"
        assert orca_out.exchange_range_seperation.lower() == "not used"
        assert orca_out.finite_nucleus_model.lower() == "not used"
        assert orca_out.aux_j_fitting_basis.lower() == "not available"
        assert orca_out.aux_j_num_basis_functions is None
        assert orca_out.aux_j_num_shells is None
        assert orca_out.aux_j_max_ang_mom is None
        assert orca_out.aux_jk_fitting_basis.lower() == "not available"
        assert orca_out.aux_k_fitting_basis.lower() == "not available"
        assert orca_out.aux_external_fitting_basis.lower() == "not available"
        assert orca_out.integral_threshold == 1e-10
        assert orca_out.primitive_cutoff == 1e-11
        assert orca_out.primitive_pair_threshold == 1e-11
        assert orca_out.ri_approx is None
        assert orca_out.rij_cosx is None
        assert orca_out.charge == 0
        assert orca_out.multiplicity == 1
        assert orca_out.num_electrons == 396
        assert orca_out.basis_dim == 876
        assert orca_out.diis_acceleration is True
        assert orca_out.scf_maxiter == 125
        assert orca_out.converged is None
        assert isinstance(orca_out.molecule, Molecule)
        assert orca_out.normal_termination is True
        assert orca_out.dfet_embed_energy is None
        assert np.isclose(orca_out.homo_energy, -9.382, rtol=1e-4)
        assert np.isclose(orca_out.lumo_energy, 1.3054, rtol=1e-4)
        assert np.isclose(orca_out.fmo_gap, 10.6874, rtol=1e-4)
        assert np.isclose(orca_out.final_energy, -3017.439958087227, rtol=1e-4)
        assert np.isclose(orca_out.final_energy, -3017.44162461, rtol=1e-4)
        assert orca_out.mulliken_atomic_charges == {
            "O0": -0.544377,
            "O1": -0.409137,
            "O2": -0.606993,
            "C3": 0.677261,
            "C4": 0.867621,
            "O5": -0.528054,
            "O6": -0.567098,
            "O7": -0.379967,
            "C8": 0.817666,
            "C9": 0.693779,
            "C10": -0.203223,
            "C11": 0.011123,
            "C12": 0.046179,
            "C13": -0.002551,
            "H14": 0.063852,
            "C15": -0.00488,
            "H16": 0.084111,
            "C17": -0.208516,
            "H18": 0.042411,
            "H19": 0.040349,
            "C20": -0.197355,
            "C21": 0.029296,
            "C22": 0.059636,
            "C23": -0.003394,
            "H24": 0.047614,
            "C25": -0.006339,
            "H26": 0.069309,
            "C27": -0.208924,
            "H28": 0.039999,
            "H29": 0.042263,
            "C30": -0.265534,
            "C31": 0.061251,
            "C32": 0.057398,
            "C33": -0.002789,
            "H34": 0.078851,
            "C35": -0.002288,
            "H36": 0.069767,
            "C37": -0.19645,
            "H38": 0.051254,
            "H39": 0.05065,
            "C40": -0.294612,
            "C41": 0.040797,
            "C42": 0.052521,
            "C43": 0.000933,
            "H44": 0.073947,
            "C45": 0.000287,
            "H46": 0.076551,
            "C47": -0.20139,
            "H48": 0.048462,
            "H49": 0.049734,
            "C50": -0.029678,
            "H51": 0.079891,
            "H52": 0.081113,
            "C53": 0.898756,
            "C54": -0.027638,
            "H55": 0.074581,
            "H56": 0.075537,
            "C57": 0.898419,
            "C58": -0.0271,
            "H59": 0.0765,
            "H60": 0.075212,
            "C61": 0.897581,
            "C62": -0.03138,
            "H63": 0.081663,
            "H64": 0.08177,
            "C65": 0.899537,
            "F66": -0.300158,
            "F67": -0.298385,
            "F68": -0.300066,
            "F69": -0.301256,
            "F70": -0.299519,
            "F71": -0.301222,
            "F72": -0.301847,
            "F73": -0.30405,
            "F74": -0.30104,
            "F75": -0.303988,
            "F76": -0.300916,
            "F77": -0.30332,
        }
        assert orca_out.loewdin_atomic_charges == {
            "O0": -0.153581,
            "O1": -0.246133,
            "O2": -0.333881,
            "C3": 0.232162,
            "C4": 0.303606,
            "O5": -0.169518,
            "O6": -0.204746,
            "O7": -0.271108,
            "C8": 0.294594,
            "C9": 0.315946,
            "C10": -0.075941,
            "C11": 0.016906,
            "C12": 0.029775,
            "C13": -0.025404,
            "H14": 0.040526,
            "C15": -0.024275,
            "H16": 0.041787,
            "C17": -0.03112,
            "H18": 0.028985,
            "H19": 0.028792,
            "C20": -0.069475,
            "C21": 0.010548,
            "C22": 0.032577,
            "C23": -0.025207,
            "H24": 0.031748,
            "C25": -0.019802,
            "H26": 0.040965,
            "C27": -0.032088,
            "H28": 0.028071,
            "H29": 0.029927,
            "C30": -0.09081,
            "C31": 0.042271,
            "C32": 0.040749,
            "C33": -0.0233,
            "H34": 0.045632,
            "C35": -0.023676,
            "H36": 0.041708,
            "C37": -0.009317,
            "H38": 0.033282,
            "H39": 0.033126,
            "C40": -0.097899,
            "C41": 0.033678,
            "C42": 0.026727,
            "C43": -0.022871,
            "H44": 0.043848,
            "C45": -0.022604,
            "H46": 0.042688,
            "C47": -0.016927,
            "H48": 0.032007,
            "H49": 0.032482,
            "C50": -0.02776,
            "H51": 0.045321,
            "H52": 0.045782,
            "C53": 0.294505,
            "C54": -0.028567,
            "H55": 0.043014,
            "H56": 0.043381,
            "C57": 0.294417,
            "C58": -0.028161,
            "H59": 0.044058,
            "H60": 0.043151,
            "C61": 0.294271,
            "C62": -0.027717,
            "H63": 0.046041,
            "H64": 0.046117,
            "C65": 0.294476,
            "F66": -0.114027,
            "F67": -0.113991,
            "F68": -0.113914,
            "F69": -0.115055,
            "F70": -0.115092,
            "F71": -0.114978,
            "F72": -0.11563,
            "F73": -0.119358,
            "F74": -0.114803,
            "F75": -0.117527,
            "F76": -0.114716,
            "F77": -0.118666,
        }
        assert orca_out.mayer_total_nuclear_charge == {
            "O0": 8.0,
            "O1": 8.0,
            "O2": 8.0,
            "C3": 6.0,
            "C4": 6.0,
            "O5": 8.0,
            "O6": 8.0,
            "O7": 8.0,
            "C8": 6.0,
            "C9": 6.0,
            "C10": 6.0,
            "C11": 6.0,
            "C12": 6.0,
            "C13": 6.0,
            "H14": 1.0,
            "C15": 6.0,
            "H16": 1.0,
            "C17": 6.0,
            "H18": 1.0,
            "H19": 1.0,
            "C20": 6.0,
            "C21": 6.0,
            "C22": 6.0,
            "C23": 6.0,
            "H24": 1.0,
            "C25": 6.0,
            "H26": 1.0,
            "C27": 6.0,
            "H28": 1.0,
            "H29": 1.0,
            "C30": 6.0,
            "C31": 6.0,
            "C32": 6.0,
            "C33": 6.0,
            "H34": 1.0,
            "C35": 6.0,
            "H36": 1.0,
            "C37": 6.0,
            "H38": 1.0,
            "H39": 1.0,
            "C40": 6.0,
            "C41": 6.0,
            "C42": 6.0,
            "C43": 6.0,
            "H44": 1.0,
            "C45": 6.0,
            "H46": 1.0,
            "C47": 6.0,
            "H48": 1.0,
            "H49": 1.0,
            "C50": 6.0,
            "H51": 1.0,
            "H52": 1.0,
            "C53": 6.0,
            "C54": 6.0,
            "H55": 1.0,
            "H56": 1.0,
            "C57": 6.0,
            "C58": 6.0,
            "H59": 1.0,
            "H60": 1.0,
            "C61": 6.0,
            "C62": 6.0,
            "H63": 1.0,
            "H64": 1.0,
            "C65": 6.0,
            "F66": 9.0,
            "F67": 9.0,
            "F68": 9.0,
            "F69": 9.0,
            "F70": 9.0,
            "F71": 9.0,
            "F72": 9.0,
            "F73": 9.0,
            "F74": 9.0,
            "F75": 9.0,
            "F76": 9.0,
            "F77": 9.0,
        }
        assert orca_out.mayer_mulliken_gross_atomic_charge == {
            "O0": -0.5444,
            "O1": -0.4091,
            "O2": -0.607,
            "C3": 0.6773,
            "C4": 0.8676,
            "O5": -0.5281,
            "O6": -0.5671,
            "O7": -0.38,
            "C8": 0.8177,
            "C9": 0.6938,
            "C10": -0.2032,
            "C11": 0.0111,
            "C12": 0.0462,
            "C13": -0.0026,
            "H14": 0.0639,
            "C15": -0.0049,
            "H16": 0.0841,
            "C17": -0.2085,
            "H18": 0.0424,
            "H19": 0.0403,
            "C20": -0.1974,
            "C21": 0.0293,
            "C22": 0.0596,
            "C23": -0.0034,
            "H24": 0.0476,
            "C25": -0.0063,
            "H26": 0.0693,
            "C27": -0.2089,
            "H28": 0.04,
            "H29": 0.0423,
            "C30": -0.2655,
            "C31": 0.0613,
            "C32": 0.0574,
            "C33": -0.0028,
            "H34": 0.0789,
            "C35": -0.0023,
            "H36": 0.0698,
            "C37": -0.1965,
            "H38": 0.0513,
            "H39": 0.0507,
            "C40": -0.2946,
            "C41": 0.0408,
            "C42": 0.0525,
            "C43": 0.0009,
            "H44": 0.0739,
            "C45": 0.0003,
            "H46": 0.0766,
            "C47": -0.2014,
            "H48": 0.0485,
            "H49": 0.0497,
            "C50": -0.0297,
            "H51": 0.0799,
            "H52": 0.0811,
            "C53": 0.8988,
            "C54": -0.0276,
            "H55": 0.0746,
            "H56": 0.0755,
            "C57": 0.8984,
            "C58": -0.0271,
            "H59": 0.0765,
            "H60": 0.0752,
            "C61": 0.8976,
            "C62": -0.0314,
            "H63": 0.0817,
            "H64": 0.0818,
            "C65": 0.8995,
            "F66": -0.3002,
            "F67": -0.2984,
            "F68": -0.3001,
            "F69": -0.3013,
            "F70": -0.2995,
            "F71": -0.3012,
            "F72": -0.3018,
            "F73": -0.304,
            "F74": -0.301,
            "F75": -0.304,
            "F76": -0.3009,
            "F77": -0.3033,
        }
        assert orca_out.mayer_total_valence == {
            "O0": 1.9052,
            "O1": 2.0594,
            "O2": 1.7735,
            "C3": 3.9521,
            "C4": 3.7926,
            "O5": 1.9317,
            "O6": 1.9052,
            "O7": 2.1188,
            "C8": 3.787,
            "C9": 3.976,
            "C10": 3.7922,
            "C11": 3.9566,
            "C12": 3.9143,
            "C13": 3.89,
            "H14": 1.0025,
            "C15": 3.8948,
            "H16": 1.0121,
            "C17": 3.8309,
            "H18": 0.9961,
            "H19": 0.9958,
            "C20": 3.8656,
            "C21": 3.9155,
            "C22": 3.9648,
            "C23": 3.9025,
            "H24": 0.999,
            "C25": 3.9074,
            "H26": 1.0074,
            "C27": 3.8315,
            "H28": 0.9961,
            "H29": 0.9954,
            "C30": 3.8275,
            "C31": 3.9196,
            "C32": 3.934,
            "C33": 3.8981,
            "H34": 0.9991,
            "C35": 3.9001,
            "H36": 1.0001,
            "C37": 3.8314,
            "H38": 0.9962,
            "H39": 0.9962,
            "C40": 3.7756,
            "C41": 3.9398,
            "C42": 3.9018,
            "C43": 3.8946,
            "H44": 1.0026,
            "C45": 3.8959,
            "H46": 1.0001,
            "C47": 3.8302,
            "H48": 0.996,
            "H49": 0.9962,
            "C50": 3.8917,
            "H51": 0.9901,
            "H52": 0.9896,
            "C53": 3.9948,
            "C54": 3.8931,
            "H55": 0.9907,
            "H56": 0.9899,
            "C57": 3.9901,
            "C58": 3.8928,
            "H59": 0.9908,
            "H60": 0.9896,
            "C61": 3.989,
            "C62": 3.8918,
            "H63": 0.9897,
            "H64": 0.9898,
            "C65": 3.9981,
            "F66": 0.9947,
            "F67": 0.9985,
            "F68": 0.9948,
            "F69": 0.9935,
            "F70": 0.9972,
            "F71": 0.9934,
            "F72": 0.9931,
            "F73": 0.9916,
            "F74": 0.9938,
            "F75": 0.9905,
            "F76": 0.9939,
            "F77": 0.9925,
        }
        assert orca_out.mayer_bonded_valence == {
            "O0": 1.9052,
            "O1": 2.0594,
            "O2": 1.7735,
            "C3": 3.9521,
            "C4": 3.7926,
            "O5": 1.9317,
            "O6": 1.9052,
            "O7": 2.1188,
            "C8": 3.787,
            "C9": 3.976,
            "C10": 3.7922,
            "C11": 3.9566,
            "C12": 3.9143,
            "C13": 3.89,
            "H14": 1.0025,
            "C15": 3.8948,
            "H16": 1.0121,
            "C17": 3.8309,
            "H18": 0.9961,
            "H19": 0.9958,
            "C20": 3.8656,
            "C21": 3.9155,
            "C22": 3.9648,
            "C23": 3.9025,
            "H24": 0.999,
            "C25": 3.9074,
            "H26": 1.0074,
            "C27": 3.8315,
            "H28": 0.9961,
            "H29": 0.9954,
            "C30": 3.8275,
            "C31": 3.9196,
            "C32": 3.934,
            "C33": 3.8981,
            "H34": 0.9991,
            "C35": 3.9001,
            "H36": 1.0001,
            "C37": 3.8314,
            "H38": 0.9962,
            "H39": 0.9962,
            "C40": 3.7756,
            "C41": 3.9398,
            "C42": 3.9018,
            "C43": 3.8946,
            "H44": 1.0026,
            "C45": 3.8959,
            "H46": 1.0001,
            "C47": 3.8302,
            "H48": 0.996,
            "H49": 0.9962,
            "C50": 3.8917,
            "H51": 0.9901,
            "H52": 0.9896,
            "C53": 3.9948,
            "C54": 3.8931,
            "H55": 0.9907,
            "H56": 0.9899,
            "C57": 3.9901,
            "C58": 3.8928,
            "H59": 0.9908,
            "H60": 0.9896,
            "C61": 3.989,
            "C62": 3.8918,
            "H63": 0.9897,
            "H64": 0.9898,
            "C65": 3.9981,
            "F66": 0.9947,
            "F67": 0.9985,
            "F68": 0.9948,
            "F69": 0.9935,
            "F70": 0.9972,
            "F71": 0.9934,
            "F72": 0.9931,
            "F73": 0.9916,
            "F74": 0.9938,
            "F75": 0.9905,
            "F76": 0.9939,
            "F77": 0.9925,
        }
        assert orca_out.mayer_free_valence == {
            "O0": -0.0,
            "O1": -0.0,
            "O2": -0.0,
            "C3": 0.0,
            "C4": -0.0,
            "O5": 0.0,
            "O6": 0.0,
            "O7": 0.0,
            "C8": 0.0,
            "C9": 0.0,
            "C10": 0.0,
            "C11": -0.0,
            "C12": -0.0,
            "C13": 0.0,
            "H14": 0.0,
            "C15": -0.0,
            "H16": -0.0,
            "C17": -0.0,
            "H18": 0.0,
            "H19": -0.0,
            "C20": -0.0,
            "C21": -0.0,
            "C22": 0.0,
            "C23": 0.0,
            "H24": -0.0,
            "C25": -0.0,
            "H26": 0.0,
            "C27": 0.0,
            "H28": -0.0,
            "H29": 0.0,
            "C30": 0.0,
            "C31": 0.0,
            "C32": -0.0,
            "C33": -0.0,
            "H34": -0.0,
            "C35": -0.0,
            "H36": 0.0,
            "C37": -0.0,
            "H38": -0.0,
            "H39": -0.0,
            "C40": -0.0,
            "C41": -0.0,
            "C42": 0.0,
            "C43": -0.0,
            "H44": -0.0,
            "C45": -0.0,
            "H46": 0.0,
            "C47": -0.0,
            "H48": -0.0,
            "H49": -0.0,
            "C50": 0.0,
            "H51": 0.0,
            "H52": -0.0,
            "C53": -0.0,
            "C54": 0.0,
            "H55": 0.0,
            "H56": -0.0,
            "C57": -0.0,
            "C58": -0.0,
            "H59": -0.0,
            "H60": -0.0,
            "C61": 0.0,
            "C62": 0.0,
            "H63": -0.0,
            "H64": 0.0,
            "C65": 0.0,
            "F66": 0.0,
            "F67": 0.0,
            "F68": 0.0,
            "F69": -0.0,
            "F70": -0.0,
            "F71": 0.0,
            "F72": -0.0,
            "F73": -0.0,
            "F74": 0.0,
            "F75": 0.0,
            "F76": -0.0,
            "F77": -0.0,
        }
        assert orca_out.mayer_bond_orders_larger_than_zero_point_one == {
            "B(O0,C3)": 1.0348,
            "B(O0,C4)": 0.8825,
            "B(O1,C3)": 1.949,
            "B(O2,C4)": 1.3485,
            "B(O2,C9)": 0.4023,
            "B(C3,C40)": 0.9896,
            "B(C4,O6)": 0.523,
            "B(C4,C20)": 0.9928,
            "B(O5,C8)": 1.3616,
            "B(O5,C9)": 0.5414,
            "B(O6,C8)": 1.3566,
            "B(O7,C9)": 2.0474,
            "B(C8,C30)": 1.0372,
            "B(C9,C10)": 0.9694,
            "B(C10,C11)": 1.4585,
            "B(C10,C12)": 1.3794,
            "B(C11,C13)": 1.4076,
            "B(C11,H14)": 0.9993,
            "B(C12,C15)": 1.4513,
            "B(C12,H16)": 0.992,
            "B(C13,C17)": 1.4249,
            "B(C13,H18)": 1.0007,
            "B(C15,C17)": 1.3845,
            "B(C15,H19)": 1.0017,
            "B(C17,C54)": 1.0073,
            "B(C20,C21)": 1.3659,
            "B(C20,C22)": 1.4737,
            "B(C21,C23)": 1.4569,
            "B(C21,H24)": 1.0024,
            "B(C22,C25)": 1.4132,
            "B(C22,H26)": 0.9894,
            "B(C23,C27)": 1.3856,
            "B(C23,H28)": 0.9999,
            "B(C25,C27)": 1.4263,
            "B(C25,H29)": 1.0007,
            "B(C27,C58)": 1.0067,
            "B(C30,C31)": 1.3915,
            "B(C30,C32)": 1.3908,
            "B(C31,C33)": 1.4404,
            "B(C31,H34)": 0.995,
            "B(C32,C35)": 1.4461,
            "B(C32,H36)": 0.9956,
            "B(C33,C37)": 1.4007,
            "B(C33,H38)": 0.9997,
            "B(C35,C37)": 1.3964,
            "B(C35,H39)": 0.9999,
            "B(C37,C62)": 1.0052,
            "B(C40,C41)": 1.4046,
            "B(C40,C42)": 1.3898,
            "B(C41,C43)": 1.4423,
            "B(C41,H44)": 0.9961,
            "B(C42,C45)": 1.4366,
            "B(C42,H46)": 0.996,
            "B(C43,C47)": 1.3967,
            "B(C43,H48)": 1.0006,
            "B(C45,C47)": 1.4039,
            "B(C45,H49)": 0.9998,
            "B(C47,C50)": 1.0057,
            "B(C50,H51)": 0.962,
            "B(C50,H52)": 0.9619,
            "B(C50,C53)": 1.0319,
            "B(C53,F69)": 0.9887,
            "B(C53,F70)": 0.9915,
            "B(C53,F71)": 0.9887,
            "B(C54,H55)": 0.9626,
            "B(C54,H56)": 0.9623,
            "B(C54,C57)": 1.0315,
            "B(C57,F72)": 0.9887,
            "B(C57,F73)": 0.9862,
            "B(C57,F74)": 0.9898,
            "B(C58,H59)": 0.9621,
            "B(C58,H60)": 0.9622,
            "B(C58,C61)": 1.0317,
            "B(C61,F75)": 0.9857,
            "B(C61,F76)": 0.9904,
            "B(C61,F77)": 0.9871,
            "B(C62,H63)": 0.962,
            "B(C62,H64)": 0.9619,
            "B(C62,C65)": 1.0319,
            "B(C65,F66)": 0.9897,
            "B(C65,F67)": 0.9927,
            "B(C65,F68)": 0.9897,
        }
        assert all(
            np.isclose(
                orca_out.dipole_moment_electric_contribution,
                np.array([-9.52242, 11.27995, 14.97503]),
                rtol=1e-4,
            )
        )
        assert all(
            np.isclose(
                orca_out.dipole_moment_nuclear_contribution,
                np.array([8.32449, -12.31917, -16.05228]),
                rtol=1e-4,
            )
        )
        assert all(
            np.isclose(
                orca_out.dipole_moment_nuclear_contribution,
                np.array([8.32449, -12.31917, -16.05228]),
                rtol=1e-4,
            )
        )
        assert all(
            np.isclose(
                orca_out.total_dipole_moment,
                np.array([-1.19793, -1.03921, -1.07725]),
                rtol=1e-4,
            )
        )
        assert all(
            np.isclose(
                orca_out.dipole_moment_nuclear_contribution,
                np.array([8.32449, -12.31917, -16.05228]),
                rtol=1e-4,
            )
        )
        assert orca_out.dipole_moment_in_au == 1.91715
        assert orca_out.dipole_moment_in_debye == 4.87300
        assert all(
            np.isclose(
                orca_out.dipole_moment_along_axis_in_au,
                np.array([1.230072, -1.002711, -1.075613]),
                rtol=1e-4,
            )
        )

        assert all(
            np.isclose(
                orca_out.dipole_moment_along_axis_in_debye,
                np.array([3.126594, -2.54869, -2.733992]),
                rtol=1e-4,
            )
        )

        assert all(
            np.isclose(
                orca_out.rotational_constants_in_wavenumbers,
                np.array([0.001029, 0.000946, 0.000779]),
                rtol=1e-4,
            )
        )

        assert all(
            np.isclose(
                orca_out.rotational_constants_in_MHz,
                np.array([30.855324, 28.357580, 23.346614]),
                rtol=1e-4,
            )
        )

        assert orca_out.total_integrated_alpha_density == 198.000148686
        assert orca_out.total_integrated_beta_density == 198.000148686
        assert orca_out.hirshfeld_charges == {
            "O1": -0.177356,
            "O2": -0.318637,
            "O3": -0.327626,
            "C4": 0.307497,
            "C5": 0.323403,
            "O6": -0.188261,
            "O7": -0.219327,
            "O8": -0.356855,
            "C9": 0.341168,
            "C10": 0.337474,
            "C11": -0.006514,
            "C12": -0.009178,
            "C13": -0.005508,
            "C14": -0.032259,
            "H15": 0.040593,
            "C16": -0.030158,
            "H17": 0.034323,
            "C18": 0.009491,
            "H19": 0.031173,
            "H20": 0.031144,
            "C21": -0.014393,
            "C22": -0.01634,
            "C23": -0.002196,
            "C24": -0.032917,
            "H25": 0.032388,
            "C26": -0.026389,
            "H27": 0.038283,
            "C28": 0.009334,
            "H29": 0.030111,
            "H30": 0.032852,
            "C31": -0.014945,
            "C32": 0.006681,
            "C33": 0.005171,
            "C34": -0.026901,
            "H35": 0.045827,
            "C36": -0.027402,
            "H37": 0.042418,
            "C38": 0.025729,
            "H39": 0.036044,
            "H40": 0.035829,
            "C41": -0.019938,
            "C42": -0.000832,
            "C43": -0.004276,
            "C44": -0.02817,
            "H45": 0.043913,
            "C46": -0.028343,
            "H47": 0.040816,
            "C48": 0.020241,
            "H49": 0.034702,
            "H50": 0.035052,
            "C51": -0.025018,
            "H52": 0.041582,
            "H53": 0.041853,
            "C54": 0.381108,
            "C55": -0.026646,
            "H56": 0.039017,
            "H57": 0.03902,
            "C58": 0.380289,
            "C59": -0.026345,
            "H60": 0.040108,
            "H61": 0.038841,
            "C62": 0.380091,
            "C63": -0.024559,
            "H64": 0.042274,
            "H65": 0.042368,
            "C66": 0.381591,
            "F67": -0.148676,
            "F68": -0.148678,
            "F69": -0.148572,
            "F70": -0.149508,
            "F71": -0.149659,
            "F72": -0.149747,
            "F73": -0.150144,
            "F74": -0.153264,
            "F75": -0.149999,
            "F76": -0.151498,
            "F77": -0.150262,
            "F78": -0.152801,
        }
        assert orca_out.hirshfeld_spin_densities == {
            "O1": 0.0,
            "O2": 0.0,
            "O3": 0.0,
            "C4": 0.0,
            "C5": 0.0,
            "O6": 0.0,
            "O7": 0.0,
            "O8": 0.0,
            "C9": 0.0,
            "C10": 0.0,
            "C11": 0.0,
            "C12": 0.0,
            "C13": 0.0,
            "C14": 0.0,
            "H15": 0.0,
            "C16": 0.0,
            "H17": 0.0,
            "C18": 0.0,
            "H19": 0.0,
            "H20": 0.0,
            "C21": 0.0,
            "C22": 0.0,
            "C23": 0.0,
            "C24": 0.0,
            "H25": 0.0,
            "C26": 0.0,
            "H27": 0.0,
            "C28": 0.0,
            "H29": 0.0,
            "H30": 0.0,
            "C31": 0.0,
            "C32": 0.0,
            "C33": 0.0,
            "C34": 0.0,
            "H35": 0.0,
            "C36": 0.0,
            "H37": 0.0,
            "C38": 0.0,
            "H39": 0.0,
            "H40": 0.0,
            "C41": 0.0,
            "C42": 0.0,
            "C43": 0.0,
            "C44": 0.0,
            "H45": 0.0,
            "C46": 0.0,
            "H47": 0.0,
            "C48": 0.0,
            "H49": 0.0,
            "H50": 0.0,
            "C51": 0.0,
            "H52": 0.0,
            "H53": 0.0,
            "C54": 0.0,
            "C55": 0.0,
            "H56": 0.0,
            "H57": 0.0,
            "C58": 0.0,
            "C59": 0.0,
            "H60": 0.0,
            "H61": 0.0,
            "C62": 0.0,
            "C63": 0.0,
            "H64": 0.0,
            "H65": 0.0,
            "C66": 0.0,
            "F67": 0.0,
            "F68": 0.0,
            "F69": 0.0,
            "F70": 0.0,
            "F71": 0.0,
            "F72": 0.0,
            "F73": 0.0,
            "F74": 0.0,
            "F75": 0.0,
            "F76": 0.0,
            "F77": 0.0,
            "F78": 0.0,
        }

    def test_gtoint_errfile(self, gtoint_errfile):
        orca_out = ORCAOutput(filename=gtoint_errfile)
        assert orca_out.route_string == "! m062x def2-svp opt freq defgrid3"
        assert orca_out.functional == "m062x"
        assert orca_out.basis == "def2-svp"
        assert orca_out.defgrid == "defgrid3"
        assert orca_out.ab_initio is None
        assert orca_out.aux_basis is None
        assert orca_out.extrapolation_basis is None
        assert orca_out.natoms == 27
        assert orca_out.normal_termination is False

    def test_get_constrained_atoms(
        self,
        orca_fixed_atoms,
        orca_fixed_bonds_and_angles,
        orca_fixed_dihedral,
    ):
        fixed_atoms = ORCAOutput(filename=orca_fixed_atoms)
        assert fixed_atoms.frozen_atoms == [3, 12]  # atoms 3 and 12 are frozen

    def test_get_constrained_bond_lengths_and_angles(
        self, orca_fixed_bonds_and_angles
    ):
        fixed_bond = ORCAOutput(filename=orca_fixed_bonds_and_angles)
        assert fixed_bond.constrained_bond_lengths == {
            "B(H10,H9)": 2.4714,
        }
        assert fixed_bond.constrained_bond_angles == {
            "A(C2,C6,H9)": 69.0631,
        }

    def test_get_constrained_dihedral_angles(self, orca_fixed_dihedral):
        fixed_dihedral = ORCAOutput(filename=orca_fixed_dihedral)
        assert fixed_dihedral.constrained_dihedral_angles == {
            "D(O18,H14,C13,C4)": -125.9028,
        }


class TestORCAEngrad:
    def test_read_water_output(self, water_engrad_path):
        orca_engrad = ORCAEngradFile(filename=water_engrad_path)
        assert orca_engrad.natoms == 3
        assert math.isclose(
            orca_engrad.energy, -76.323311011349, rel_tol=1e-4
        )  # energy in Hartree
        gradient = np.array(
            [
                [0.000000006145, -0.000000044516, -0.000004685841],
                [0.000014390789, 0.000000001430, 0.000009665250],
                [-0.000014359316, 0.000000001321, 0.000009653881],
            ]
        )
        assert np.allclose(orca_engrad.gradient, gradient, rtol=1e-6)
        assert isinstance(orca_engrad.molecule, Molecule)
        assert orca_engrad.molecule.chemical_symbols == ["O", "H", "H"]
        coordinates = np.array(
            [
                [-0.0, 0.0, 0.08734059],
                [-0.75520525, 0.0, -0.50967031],
                [0.75520525, 0.0, -0.50967031],
            ]
        )
        assert np.allclose(
            orca_engrad.molecule.positions, coordinates, rtol=1e-6
        )
