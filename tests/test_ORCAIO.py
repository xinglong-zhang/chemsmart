import math
import os
from shutil import copy

import numpy as np
import pytest
from ase import units

from chemsmart.io.molecules.structure import CoordinateBlock, Molecule
from chemsmart.io.orca import ORCARefs
from chemsmart.io.orca.input import ORCAInput, ORCAQMMMInput
from chemsmart.io.orca.output import (
    ORCAEngradFile,
    ORCANEBFile,
    ORCAOutput,
    ORCAQMMMOutput,
)
from chemsmart.io.orca.route import ORCARoute
from chemsmart.jobs.orca.settings import ORCANEBJobSettings


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

        # two-layer ONIOM
        s7 = "!QM/XTB BP86 def2-TZVP def2/J"
        r7 = ORCARoute(route_string=s7)
        assert r7.route_keywords == ["qm/xtb", "bp86", "def2-tzvp", "def2/j"]
        assert r7.qm_functional == "bp86"
        assert r7.qm_basis == "def2-tzvp"
        assert r7.auxiliary_basis == "def2/j"
        assert r7.qm2_method == "xtb"
        assert r7.qmmm_jobtype == "qm/xtb"

        # three-layer ONIOM
        s8 = "!QM/HF-3c/MM Opt B3LYP def2-TZVP def2/J NumFreq CPCM(water)"
        r8 = ORCARoute(route_string=s8)
        assert r8.route_keywords == [
            "qm/hf-3c/mm",
            "opt",
            "b3lyp",
            "def2-tzvp",
            "def2/j",
            "numfreq",
            "cpcm(water)",
        ]
        assert r8.qm_functional == "b3lyp"
        assert r8.qm_basis == "def2-tzvp"
        assert r8.auxiliary_basis == "def2/j"
        assert r8.qm2_method == "hf-3c"
        assert r8.qmmm_jobtype == "qm/hf-3c/mm"

        # MOL-CRYSTAL-QMMM route
        s9 = "! MOL-CRYSTAL-QMMM PBE def2-SVP Opt NumFreq"
        r9 = ORCARoute(route_string=s9)
        assert r9.route_keywords == [
            "mol-crystal-qmmm",
            "pbe",
            "def2-svp",
            "opt",
            "numfreq",
        ]
        assert r9.qm_functional == "pbe"
        assert r9.qm_basis == "def2-svp"
        assert r9.qmmm_jobtype == "mol-crystal-qmmm"

        # IONIC-CRYSTAL-QMMM route
        s10 = "! IONIC-CRYSTAL-QMMM"
        r10 = ORCARoute(route_string=s10)
        assert r10.qmmm_jobtype == "ionic-crystal-qmmm"


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

    def test_orca_neb_input_with_xyz_files_specified(
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
        assert orca_inp.route_string == "! gfn2-xtb neb-ts freq".lower()
        assert orca_inp.functional is None
        assert orca_inp.basis is None
        assert orca_inp.coordinate_type == "xyzfile"  # xyzfile is specified
        assert orca_inp.charge == 0
        assert orca_inp.multiplicity == 2
        assert orca_inp.molecule.num_atoms == 40
        assert isinstance(orca_inp.molecule, Molecule)
        assert orca_inp.molecule.empirical_formula == "C23H15NO"

    def test_orca_qmmm_input(self, orca_inputs_directory):
        orca_inp1 = os.path.join(orca_inputs_directory, "dna_qmmm1.inp")
        orca_inp1 = ORCAQMMMInput(filename=orca_inp1)
        # charge and multiplicity of QM region (instead of real system in regular input)
        assert orca_inp1.qm_charge == 2
        assert orca_inp1.qm_multiplicity == 1
        assert orca_inp1.qm_atoms == [
            "54",
            "124:133",
            "209",
            "210",
            "259:263",
            "271",
            "272",
            "326:340",
            "424:476",
            "488:516",
        ]
        assert orca_inp1.qm_active_atoms == ["0:5", "16", "21:30"]
        # assert orca_inp.qm_force_field
        assert orca_inp1.qm_h_bond_length == [
            ("c", "hla", "1.09"),
            ("o", "hla", "0.98"),
            ("n", "hla", "0.99"),
        ]
        assert orca_inp1.qm_boundary_interaction == (
            "Will neglect bends at QM2-QM1-MM1 and torsions at QM3-QM2-QM1-MM1 boundary.\n"
            "Will include bonds at QM1-MM1 boundary.\n"
        )
        assert orca_inp1.qm_embedding_type == "electrostatic"
        assert orca_inp1.qm2_functional.strip('"') == "b3lyp"
        assert orca_inp1.qm2_basis.strip('"') == "def2-svp def2/j"

        orca_inp2 = os.path.join(orca_inputs_directory, "dna_qmmm2.inp")
        orca_inp2 = ORCAQMMMInput(filename=orca_inp2)
        assert orca_inp2.qm2_level_of_theory.strip('"') == "myqm2method.txt"
        assert orca_inp2.qm_qm2_boundary_treatment == "pbeh3c"
        assert orca_inp2.qm2_atoms == ["5:22"]
        assert orca_inp2.qm2_charge == 0
        assert orca_inp2.qm2_multiplicity == 3

        # todo:tests for crystal QMMM
        # orca_inp3 = os.path.join(orca_inputs_directory, "ionic_crystal_qmmm.inp")


class TestORCAOutput:
    def test_read_water_output(self, water_output_gas_path):
        orca_out = ORCAOutput(filename=water_output_gas_path)
        assert isinstance(orca_out.molecule, Molecule)
        assert orca_out.route_string == "! opt freq m062x def2-svp"
        assert orca_out.functional == "m062x"
        assert orca_out.basis == "def2-svp"
        assert orca_out.spin == "restricted"
        assert orca_out.ab_initio is None
        assert orca_out.aux_basis is None
        assert orca_out.extrapolation_basis is None
        assert orca_out.num_atoms == 3
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
        assert orca_out.alpha_occ_eigenvalues == [
            -19.612911 * units.Hartree,
            -1.095421 * units.Hartree,
            -0.585668 * units.Hartree,
            -0.449226 * units.Hartree,
            -0.376778 * units.Hartree,
        ]
        assert orca_out.beta_occ_eigenvalues == [
            -19.612911 * units.Hartree,
            -1.095421 * units.Hartree,
            -0.585668 * units.Hartree,
            -0.449226 * units.Hartree,
            -0.376778 * units.Hartree,
        ]
        assert orca_out.alpha_virtual_eigenvalues == [
            0.088932 * units.Hartree,
            0.166094 * units.Hartree,
            0.619581 * units.Hartree,
            0.688090 * units.Hartree,
            0.992526 * units.Hartree,
            1.002028 * units.Hartree,
            1.072543 * units.Hartree,
            1.150532 * units.Hartree,
            1.395583 * units.Hartree,
            1.446410 * units.Hartree,
            1.590115 * units.Hartree,
            1.830301 * units.Hartree,
            2.278882 * units.Hartree,
            2.321601 * units.Hartree,
            3.028276 * units.Hartree,
            3.067006 * units.Hartree,
            3.248959 * units.Hartree,
            3.542159 * units.Hartree,
            3.864433 * units.Hartree,
        ]
        assert orca_out.beta_virtual_eigenvalues == [
            0.088932 * units.Hartree,
            0.166094 * units.Hartree,
            0.619581 * units.Hartree,
            0.688090 * units.Hartree,
            0.992526 * units.Hartree,
            1.002028 * units.Hartree,
            1.072543 * units.Hartree,
            1.150532 * units.Hartree,
            1.395583 * units.Hartree,
            1.446410 * units.Hartree,
            1.590115 * units.Hartree,
            1.830301 * units.Hartree,
            2.278882 * units.Hartree,
            2.321601 * units.Hartree,
            3.028276 * units.Hartree,
            3.067006 * units.Hartree,
            3.248959 * units.Hartree,
            3.542159 * units.Hartree,
            3.864433 * units.Hartree,
        ]
        assert orca_out.num_unpaired_electrons == 0
        # For closed-shell systems, somo_energies should be None
        assert orca_out.somo_energies is None
        assert orca_out.lowest_somo_energy is None
        assert orca_out.highest_somo_energy is None
        assert orca_out.alpha_homo_energy == -0.376778 * units.Hartree
        assert orca_out.beta_homo_energy == -0.376778 * units.Hartree
        assert orca_out.alpha_lumo_energy == 0.088932 * units.Hartree
        assert orca_out.beta_lumo_energy == 0.088932 * units.Hartree
        assert orca_out.homo_energy == -0.376778 * units.Hartree
        assert orca_out.lumo_energy == 0.088932 * units.Hartree
        assert np.isclose(
            orca_out.fmo_gap,
            (0.088932 - (-0.376778)) * units.Hartree,
            rtol=1e-6,
        )

        assert orca_out.mulliken_atomic_charges == {
            "O1": -0.32926,
            "H2": 0.16463,
            "H3": 0.16463,
        }
        assert orca_out.loewdin_atomic_charges == {
            "O1": -0.155184,
            "H2": 0.077592,
            "H3": 0.077592,
        }
        assert orca_out.mayer_mulliken_gross_atomic_population == {
            "O1": 8.3293,
            "H2": 0.8354,
            "H3": 0.8354,
        }
        assert orca_out.mayer_total_nuclear_charge == {
            "O1": 8.0,
            "H2": 1.0,
            "H3": 1.0,
        }
        assert orca_out.mayer_mulliken_gross_atomic_charge == {
            "O1": -0.3293,
            "H2": 0.1646,
            "H3": 0.1646,
        }
        assert orca_out.mayer_total_valence == {
            "O1": 2.0072,
            "H2": 1.0104,
            "H3": 1.0104,
        }
        assert orca_out.mayer_bonded_valence == {
            "O1": 2.0072,
            "H2": 1.0104,
            "H3": 1.0104,
        }
        assert orca_out.mayer_free_valence == {"H2": 0.0, "H3": 0.0, "O1": 0.0}
        assert orca_out.mayer_bond_orders_larger_than_zero_point_one == {
            "B(O1,H2)": 1.0036,
            "B(O1,H3)": 1.0036,
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
        mode1 = np.array(
            [
                [-0.0, 0.0, -0.069893],
                [-0.43577, -0.0, 0.554673],
                [0.43577, -0.0, 0.554673],
            ]
        )
        mode2 = np.array(
            [
                [0.0, -0.0, -0.050918],
                [0.579153, -0.0, 0.404085],
                [-0.579154, 0.0, 0.404085],
            ]
        )
        mode3 = np.array(
            [
                [
                    [-0.069728, -0.0, -0.0],
                    [0.553362, 0.0, 0.437448],
                    [0.553361, 0.0, -0.437447],
                ],
            ]
        )
        assert len(orca_out.vibrational_modes) == 3
        assert np.allclose(orca_out.vibrational_modes[0], mode1, rtol=1e-4)
        assert np.allclose(orca_out.vibrational_modes[1], mode2, rtol=1e-4)
        assert np.allclose(orca_out.vibrational_modes[2], mode3, rtol=1e-4)

        # test for molecule obtained from output file
        mol = orca_out.molecule
        assert len(mol.vibrational_modes) == 3
        assert np.allclose(mol.vibrational_modes[0], mode1, rtol=1e-4)
        assert np.allclose(mol.vibrational_modes[1], mode2, rtol=1e-4)
        assert np.allclose(mol.vibrational_modes[2], mode3, rtol=1e-4)

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
        transition_dipole1 = np.array([0.000000, -0.000000, 0.049416])
        transition_dipole2 = np.array([-0.000000, -0.000000, 0.015460])
        transition_dipole3 = np.array([0.027888, -0.000000, 0.000000])
        assert np.allclose(
            orca_out.transition_dipoles[0], transition_dipole1, rtol=1e-4
        )
        assert np.allclose(
            orca_out.transition_dipoles[1], transition_dipole2, rtol=1e-4
        )
        assert np.allclose(
            orca_out.transition_dipoles[2], transition_dipole3, rtol=1e-4
        )

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
        assert math.isclose(orca_out.entropy_TS, 0.02143089, rel_tol=1e-4)

        assert orca_out.mulliken_atomic_charges == {
            "O1": -0.32926,
            "H2": 0.16463,
            "H3": 0.16463,
        }
        assert orca_out.loewdin_atomic_charges == {
            "O1": -0.155184,
            "H2": 0.077592,
            "H3": 0.077592,
        }
        assert orca_out.mayer_mulliken_gross_atomic_population == {
            "O1": 8.3293,
            "H2": 0.8354,
            "H3": 0.8354,
        }
        assert orca_out.mayer_total_nuclear_charge == {
            "O1": 8.0,
            "H2": 1.0,
            "H3": 1.0,
        }
        assert orca_out.mayer_mulliken_gross_atomic_charge == {
            "O1": -0.3293,
            "H2": 0.1646,
            "H3": 0.1646,
        }
        assert orca_out.mayer_total_valence == {
            "O1": 2.0072,
            "H2": 1.0104,
            "H3": 1.0104,
        }
        assert orca_out.mayer_bonded_valence == {
            "O1": 2.0072,
            "H2": 1.0104,
            "H3": 1.0104,
        }
        assert orca_out.mayer_free_valence == {"H2": 0.0, "H3": 0.0, "O1": 0.0}
        assert orca_out.mayer_bond_orders_larger_than_zero_point_one == {
            "B(O1,H2)": 1.0036,
            "B(O1,H3)": 1.0036,
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
        assert orca_out.rotational_symmetry_number == 2
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
                "sn=1": 49.65043,
                "sn=2": 43.887276,
                "sn=3": 40.516087,
                "sn=4": 38.124122,
                "sn=5": 36.268792,
                "sn=6": 34.752933,
                "sn=7": 33.471224,
                "sn=8": 32.361055,
                "sn=9": 31.381743,
                "sn=10": 30.505726,
                "sn=11": 29.713277,
                "sn=12": 28.989778,
            }
        )

        assert math.isclose(
            orca_out.gibbs_free_energy, -76.31938148, rel_tol=1e-8
        )
        assert isinstance(orca_out.molecule, Molecule)
        assert orca_out.total_elapsed_walltime == 0.0

    def test_he_freq_output(self, orca_he_output_freq):
        orca_out = ORCAOutput(filename=orca_he_output_freq)
        assert isinstance(orca_out.molecule, Molecule)
        assert orca_out.spin == "restricted"
        assert orca_out.rotational_symmetry_number == 1
        assert orca_out.rotational_constants_in_wavenumbers == [0, 0, 0]
        assert orca_out.rotational_constants_in_MHz == [0, 0, 0]
        assert orca_out.num_vibration_modes == 0
        assert orca_out.temperature_in_K == 298.15
        assert orca_out.pressure_in_atm == 1.00
        assert orca_out.total_mass_in_amu == 4.00

    def test_read_sp_output(self, water_sp_gas_path):
        orca_out = ORCAOutput(filename=water_sp_gas_path)
        assert orca_out.spin == "restricted"
        assert (
            orca_out.route_string
            == "!  DLPNO-CCSD(T) Extrapolate(2/3,cc) AutoAux DEFGRID3 TightSCF KDIIS".lower()
        )
        assert orca_out.functional is None
        assert orca_out.basis is None
        assert orca_out.ab_initio == "DLPNO-CCSD(T)".lower()
        assert orca_out.aux_basis == "AutoAux".lower()
        assert orca_out.extrapolation_basis == "Extrapolate(2/3,cc)".lower()
        assert orca_out.num_atoms == 3
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
        assert orca_out.spin == "restricted"
        assert orca_out.functional is None
        assert orca_out.basis is None
        assert orca_out.ab_initio == "DLPNO-CCSD(T)".lower()
        assert orca_out.aux_basis == "AutoAux".lower()
        assert orca_out.extrapolation_basis == "Extrapolate(2/3,cc)".lower()
        assert orca_out.num_atoms == 78
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
        assert orca_out.spin == "restricted"
        assert orca_out.functional is None
        assert orca_out.basis is None
        assert orca_out.ab_initio is None
        assert orca_out.aux_basis is None
        assert orca_out.extrapolation_basis is None
        assert orca_out.num_atoms == 78
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
            "O1": -0.544377,
            "O2": -0.409137,
            "O3": -0.606993,
            "C4": 0.677261,
            "C5": 0.867621,
            "O6": -0.528054,
            "O7": -0.567098,
            "O8": -0.379967,
            "C9": 0.817666,
            "C10": 0.693779,
            "C11": -0.203223,
            "C12": 0.011123,
            "C13": 0.046179,
            "C14": -0.002551,
            "H15": 0.063852,
            "C16": -0.00488,
            "H17": 0.084111,
            "C18": -0.208516,
            "H19": 0.042411,
            "H20": 0.040349,
            "C21": -0.197355,
            "C22": 0.029296,
            "C23": 0.059636,
            "C24": -0.003394,
            "H25": 0.047614,
            "C26": -0.006339,
            "H27": 0.069309,
            "C28": -0.208924,
            "H29": 0.039999,
            "H30": 0.042263,
            "C31": -0.265534,
            "C32": 0.061251,
            "C33": 0.057398,
            "C34": -0.002789,
            "H35": 0.078851,
            "C36": -0.002288,
            "H37": 0.069767,
            "C38": -0.19645,
            "H39": 0.051254,
            "H40": 0.05065,
            "C41": -0.294612,
            "C42": 0.040797,
            "C43": 0.052521,
            "C44": 0.000933,
            "H45": 0.073947,
            "C46": 0.000287,
            "H47": 0.076551,
            "C48": -0.20139,
            "H49": 0.048462,
            "H50": 0.049734,
            "C51": -0.029678,
            "H52": 0.079891,
            "H53": 0.081113,
            "C54": 0.898756,
            "C55": -0.027638,
            "H56": 0.074581,
            "H57": 0.075537,
            "C58": 0.898419,
            "C59": -0.0271,
            "H60": 0.0765,
            "H61": 0.075212,
            "C62": 0.897581,
            "C63": -0.03138,
            "H64": 0.081663,
            "H65": 0.08177,
            "C66": 0.899537,
            "F67": -0.300158,
            "F68": -0.298385,
            "F69": -0.300066,
            "F70": -0.301256,
            "F71": -0.299519,
            "F72": -0.301222,
            "F73": -0.301847,
            "F74": -0.30405,
            "F75": -0.30104,
            "F76": -0.303988,
            "F77": -0.300916,
            "F78": -0.30332,
        }
        assert orca_out.loewdin_atomic_charges == {
            "O1": -0.153581,
            "O2": -0.246133,
            "O3": -0.333881,
            "C4": 0.232162,
            "C5": 0.303606,
            "O6": -0.169518,
            "O7": -0.204746,
            "O8": -0.271108,
            "C9": 0.294594,
            "C10": 0.315946,
            "C11": -0.075941,
            "C12": 0.016906,
            "C13": 0.029775,
            "C14": -0.025404,
            "H15": 0.040526,
            "C16": -0.024275,
            "H17": 0.041787,
            "C18": -0.03112,
            "H19": 0.028985,
            "H20": 0.028792,
            "C21": -0.069475,
            "C22": 0.010548,
            "C23": 0.032577,
            "C24": -0.025207,
            "H25": 0.031748,
            "C26": -0.019802,
            "H27": 0.040965,
            "C28": -0.032088,
            "H29": 0.028071,
            "H30": 0.029927,
            "C31": -0.09081,
            "C32": 0.042271,
            "C33": 0.040749,
            "C34": -0.0233,
            "H35": 0.045632,
            "C36": -0.023676,
            "H37": 0.041708,
            "C38": -0.009317,
            "H39": 0.033282,
            "H40": 0.033126,
            "C41": -0.097899,
            "C42": 0.033678,
            "C43": 0.026727,
            "C44": -0.022871,
            "H45": 0.043848,
            "C46": -0.022604,
            "H47": 0.042688,
            "C48": -0.016927,
            "H49": 0.032007,
            "H50": 0.032482,
            "C51": -0.02776,
            "H52": 0.045321,
            "H53": 0.045782,
            "C54": 0.294505,
            "C55": -0.028567,
            "H56": 0.043014,
            "H57": 0.043381,
            "C58": 0.294417,
            "C59": -0.028161,
            "H60": 0.044058,
            "H61": 0.043151,
            "C62": 0.294271,
            "C63": -0.027717,
            "H64": 0.046041,
            "H65": 0.046117,
            "C66": 0.294476,
            "F67": -0.114027,
            "F68": -0.113991,
            "F69": -0.113914,
            "F70": -0.115055,
            "F71": -0.115092,
            "F72": -0.114978,
            "F73": -0.11563,
            "F74": -0.119358,
            "F75": -0.114803,
            "F76": -0.117527,
            "F77": -0.114716,
            "F78": -0.118666,
        }
        assert orca_out.mayer_total_nuclear_charge == {
            "O1": 8.0,
            "O2": 8.0,
            "O3": 8.0,
            "C4": 6.0,
            "C5": 6.0,
            "O6": 8.0,
            "O7": 8.0,
            "O8": 8.0,
            "C9": 6.0,
            "C10": 6.0,
            "C11": 6.0,
            "C12": 6.0,
            "C13": 6.0,
            "C14": 6.0,
            "H15": 1.0,
            "C16": 6.0,
            "H17": 1.0,
            "C18": 6.0,
            "H19": 1.0,
            "H20": 1.0,
            "C21": 6.0,
            "C22": 6.0,
            "C23": 6.0,
            "C24": 6.0,
            "H25": 1.0,
            "C26": 6.0,
            "H27": 1.0,
            "C28": 6.0,
            "H29": 1.0,
            "H30": 1.0,
            "C31": 6.0,
            "C32": 6.0,
            "C33": 6.0,
            "C34": 6.0,
            "H35": 1.0,
            "C36": 6.0,
            "H37": 1.0,
            "C38": 6.0,
            "H39": 1.0,
            "H40": 1.0,
            "C41": 6.0,
            "C42": 6.0,
            "C43": 6.0,
            "C44": 6.0,
            "H45": 1.0,
            "C46": 6.0,
            "H47": 1.0,
            "C48": 6.0,
            "H49": 1.0,
            "H50": 1.0,
            "C51": 6.0,
            "H52": 1.0,
            "H53": 1.0,
            "C54": 6.0,
            "C55": 6.0,
            "H56": 1.0,
            "H57": 1.0,
            "C58": 6.0,
            "C59": 6.0,
            "H60": 1.0,
            "H61": 1.0,
            "C62": 6.0,
            "C63": 6.0,
            "H64": 1.0,
            "H65": 1.0,
            "C66": 6.0,
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
            "F78": 9.0,
        }
        assert orca_out.mayer_mulliken_gross_atomic_charge == {
            "O1": -0.5444,
            "O2": -0.4091,
            "O3": -0.607,
            "C4": 0.6773,
            "C5": 0.8676,
            "O6": -0.5281,
            "O7": -0.5671,
            "O8": -0.38,
            "C9": 0.8177,
            "C10": 0.6938,
            "C11": -0.2032,
            "C12": 0.0111,
            "C13": 0.0462,
            "C14": -0.0026,
            "H15": 0.0639,
            "C16": -0.0049,
            "H17": 0.0841,
            "C18": -0.2085,
            "H19": 0.0424,
            "H20": 0.0403,
            "C21": -0.1974,
            "C22": 0.0293,
            "C23": 0.0596,
            "C24": -0.0034,
            "H25": 0.0476,
            "C26": -0.0063,
            "H27": 0.0693,
            "C28": -0.2089,
            "H29": 0.04,
            "H30": 0.0423,
            "C31": -0.2655,
            "C32": 0.0613,
            "C33": 0.0574,
            "C34": -0.0028,
            "H35": 0.0789,
            "C36": -0.0023,
            "H37": 0.0698,
            "C38": -0.1965,
            "H39": 0.0513,
            "H40": 0.0507,
            "C41": -0.2946,
            "C42": 0.0408,
            "C43": 0.0525,
            "C44": 0.0009,
            "H45": 0.0739,
            "C46": 0.0003,
            "H47": 0.0766,
            "C48": -0.2014,
            "H49": 0.0485,
            "H50": 0.0497,
            "C51": -0.0297,
            "H52": 0.0799,
            "H53": 0.0811,
            "C54": 0.8988,
            "C55": -0.0276,
            "H56": 0.0746,
            "H57": 0.0755,
            "C58": 0.8984,
            "C59": -0.0271,
            "H60": 0.0765,
            "H61": 0.0752,
            "C62": 0.8976,
            "C63": -0.0314,
            "H64": 0.0817,
            "H65": 0.0818,
            "C66": 0.8995,
            "F67": -0.3002,
            "F68": -0.2984,
            "F69": -0.3001,
            "F70": -0.3013,
            "F71": -0.2995,
            "F72": -0.3012,
            "F73": -0.3018,
            "F74": -0.304,
            "F75": -0.301,
            "F76": -0.304,
            "F77": -0.3009,
            "F78": -0.3033,
        }
        assert orca_out.mayer_total_valence == {
            "O1": 1.9052,
            "O2": 2.0594,
            "O3": 1.7735,
            "C4": 3.9521,
            "C5": 3.7926,
            "O6": 1.9317,
            "O7": 1.9052,
            "O8": 2.1188,
            "C9": 3.787,
            "C10": 3.976,
            "C11": 3.7922,
            "C12": 3.9566,
            "C13": 3.9143,
            "C14": 3.89,
            "H15": 1.0025,
            "C16": 3.8948,
            "H17": 1.0121,
            "C18": 3.8309,
            "H19": 0.9961,
            "H20": 0.9958,
            "C21": 3.8656,
            "C22": 3.9155,
            "C23": 3.9648,
            "C24": 3.9025,
            "H25": 0.999,
            "C26": 3.9074,
            "H27": 1.0074,
            "C28": 3.8315,
            "H29": 0.9961,
            "H30": 0.9954,
            "C31": 3.8275,
            "C32": 3.9196,
            "C33": 3.934,
            "C34": 3.8981,
            "H35": 0.9991,
            "C36": 3.9001,
            "H37": 1.0001,
            "C38": 3.8314,
            "H39": 0.9962,
            "H40": 0.9962,
            "C41": 3.7756,
            "C42": 3.9398,
            "C43": 3.9018,
            "C44": 3.8946,
            "H45": 1.0026,
            "C46": 3.8959,
            "H47": 1.0001,
            "C48": 3.8302,
            "H49": 0.996,
            "H50": 0.9962,
            "C51": 3.8917,
            "H52": 0.9901,
            "H53": 0.9896,
            "C54": 3.9948,
            "C55": 3.8931,
            "H56": 0.9907,
            "H57": 0.9899,
            "C58": 3.9901,
            "C59": 3.8928,
            "H60": 0.9908,
            "H61": 0.9896,
            "C62": 3.989,
            "C63": 3.8918,
            "H64": 0.9897,
            "H65": 0.9898,
            "C66": 3.9981,
            "F67": 0.9947,
            "F68": 0.9985,
            "F69": 0.9948,
            "F70": 0.9935,
            "F71": 0.9972,
            "F72": 0.9934,
            "F73": 0.9931,
            "F74": 0.9916,
            "F75": 0.9938,
            "F76": 0.9905,
            "F77": 0.9939,
            "F78": 0.9925,
        }
        assert orca_out.mayer_bonded_valence == {
            "O1": 1.9052,
            "O2": 2.0594,
            "O3": 1.7735,
            "C4": 3.9521,
            "C5": 3.7926,
            "O6": 1.9317,
            "O7": 1.9052,
            "O8": 2.1188,
            "C9": 3.787,
            "C10": 3.976,
            "C11": 3.7922,
            "C12": 3.9566,
            "C13": 3.9143,
            "C14": 3.89,
            "H15": 1.0025,
            "C16": 3.8948,
            "H17": 1.0121,
            "C18": 3.8309,
            "H19": 0.9961,
            "H20": 0.9958,
            "C21": 3.8656,
            "C22": 3.9155,
            "C23": 3.9648,
            "C24": 3.9025,
            "H25": 0.999,
            "C26": 3.9074,
            "H27": 1.0074,
            "C28": 3.8315,
            "H29": 0.9961,
            "H30": 0.9954,
            "C31": 3.8275,
            "C32": 3.9196,
            "C33": 3.934,
            "C34": 3.8981,
            "H35": 0.9991,
            "C36": 3.9001,
            "H37": 1.0001,
            "C38": 3.8314,
            "H39": 0.9962,
            "H40": 0.9962,
            "C41": 3.7756,
            "C42": 3.9398,
            "C43": 3.9018,
            "C44": 3.8946,
            "H45": 1.0026,
            "C46": 3.8959,
            "H47": 1.0001,
            "C48": 3.8302,
            "H49": 0.996,
            "H50": 0.9962,
            "C51": 3.8917,
            "H52": 0.9901,
            "H53": 0.9896,
            "C54": 3.9948,
            "C55": 3.8931,
            "H56": 0.9907,
            "H57": 0.9899,
            "C58": 3.9901,
            "C59": 3.8928,
            "H60": 0.9908,
            "H61": 0.9896,
            "C62": 3.989,
            "C63": 3.8918,
            "H64": 0.9897,
            "H65": 0.9898,
            "C66": 3.9981,
            "F67": 0.9947,
            "F68": 0.9985,
            "F69": 0.9948,
            "F70": 0.9935,
            "F71": 0.9972,
            "F72": 0.9934,
            "F73": 0.9931,
            "F74": 0.9916,
            "F75": 0.9938,
            "F76": 0.9905,
            "F77": 0.9939,
            "F78": 0.9925,
        }
        assert orca_out.mayer_free_valence == {
            "O1": -0.0,
            "O2": -0.0,
            "O3": -0.0,
            "C4": 0.0,
            "C5": -0.0,
            "O6": 0.0,
            "O7": 0.0,
            "O8": 0.0,
            "C9": 0.0,
            "C10": 0.0,
            "C11": 0.0,
            "C12": -0.0,
            "C13": -0.0,
            "C14": 0.0,
            "H15": 0.0,
            "C16": -0.0,
            "H17": -0.0,
            "C18": -0.0,
            "H19": 0.0,
            "H20": -0.0,
            "C21": -0.0,
            "C22": -0.0,
            "C23": 0.0,
            "C24": 0.0,
            "H25": -0.0,
            "C26": -0.0,
            "H27": 0.0,
            "C28": 0.0,
            "H29": -0.0,
            "H30": 0.0,
            "C31": 0.0,
            "C32": 0.0,
            "C33": -0.0,
            "C34": -0.0,
            "H35": -0.0,
            "C36": -0.0,
            "H37": 0.0,
            "C38": -0.0,
            "H39": -0.0,
            "H40": -0.0,
            "C41": -0.0,
            "C42": -0.0,
            "C43": 0.0,
            "C44": -0.0,
            "H45": -0.0,
            "C46": -0.0,
            "H47": 0.0,
            "C48": -0.0,
            "H49": -0.0,
            "H50": -0.0,
            "C51": 0.0,
            "H52": 0.0,
            "H53": -0.0,
            "C54": -0.0,
            "C55": 0.0,
            "H56": 0.0,
            "H57": -0.0,
            "C58": -0.0,
            "C59": -0.0,
            "H60": -0.0,
            "H61": -0.0,
            "C62": 0.0,
            "C63": 0.0,
            "H64": -0.0,
            "H65": 0.0,
            "C66": 0.0,
            "F67": 0.0,
            "F68": 0.0,
            "F69": 0.0,
            "F70": -0.0,
            "F71": -0.0,
            "F72": 0.0,
            "F73": -0.0,
            "F74": -0.0,
            "F75": 0.0,
            "F76": 0.0,
            "F77": -0.0,
            "F78": -0.0,
        }
        assert orca_out.mayer_bond_orders_larger_than_zero_point_one == {
            "B(O1,C4)": 1.0348,
            "B(O1,C5)": 0.8825,
            "B(O2,C4)": 1.949,
            "B(O3,C5)": 1.3485,
            "B(O3,C10)": 0.4023,
            "B(C4,C41)": 0.9896,
            "B(C5,O7)": 0.523,
            "B(C5,C21)": 0.9928,
            "B(O6,C9)": 1.3616,
            "B(O6,C10)": 0.5414,
            "B(O7,C9)": 1.3566,
            "B(O8,C10)": 2.0474,
            "B(C9,C31)": 1.0372,
            "B(C10,C11)": 0.9694,
            "B(C11,C12)": 1.4585,
            "B(C11,C13)": 1.3794,
            "B(C12,C14)": 1.4076,
            "B(C12,H15)": 0.9993,
            "B(C13,C16)": 1.4513,
            "B(C13,H17)": 0.992,
            "B(C14,C18)": 1.4249,
            "B(C14,H19)": 1.0007,
            "B(C16,C18)": 1.3845,
            "B(C16,H20)": 1.0017,
            "B(C18,C55)": 1.0073,
            "B(C21,C22)": 1.3659,
            "B(C21,C23)": 1.4737,
            "B(C22,C24)": 1.4569,
            "B(C22,H25)": 1.0024,
            "B(C23,C26)": 1.4132,
            "B(C23,H27)": 0.9894,
            "B(C24,C28)": 1.3856,
            "B(C24,H29)": 0.9999,
            "B(C26,C28)": 1.4263,
            "B(C26,H30)": 1.0007,
            "B(C28,C59)": 1.0067,
            "B(C31,C32)": 1.3915,
            "B(C31,C33)": 1.3908,
            "B(C32,C34)": 1.4404,
            "B(C32,H35)": 0.995,
            "B(C33,C36)": 1.4461,
            "B(C33,H37)": 0.9956,
            "B(C34,C38)": 1.4007,
            "B(C34,H39)": 0.9997,
            "B(C36,C38)": 1.3964,
            "B(C36,H40)": 0.9999,
            "B(C38,C63)": 1.0052,
            "B(C41,C42)": 1.4046,
            "B(C41,C43)": 1.3898,
            "B(C42,C44)": 1.4423,
            "B(C42,H45)": 0.9961,
            "B(C43,C46)": 1.4366,
            "B(C43,H47)": 0.996,
            "B(C44,C48)": 1.3967,
            "B(C44,H49)": 1.0006,
            "B(C46,C48)": 1.4039,
            "B(C46,H50)": 0.9998,
            "B(C48,C51)": 1.0057,
            "B(C51,H52)": 0.962,
            "B(C51,H53)": 0.9619,
            "B(C51,C54)": 1.0319,
            "B(C54,F70)": 0.9887,
            "B(C54,F71)": 0.9915,
            "B(C54,F72)": 0.9887,
            "B(C55,H56)": 0.9626,
            "B(C55,H57)": 0.9623,
            "B(C55,C58)": 1.0315,
            "B(C58,F73)": 0.9887,
            "B(C58,F74)": 0.9862,
            "B(C58,F75)": 0.9898,
            "B(C59,H60)": 0.9621,
            "B(C59,H61)": 0.9622,
            "B(C59,C62)": 1.0317,
            "B(C62,F76)": 0.9857,
            "B(C62,F77)": 0.9904,
            "B(C62,F78)": 0.9871,
            "B(C63,H64)": 0.962,
            "B(C63,H65)": 0.9619,
            "B(C63,C66)": 1.0319,
            "B(C66,F67)": 0.9897,
            "B(C66,F68)": 0.9927,
            "B(C66,F69)": 0.9897,
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
        assert orca_out.num_atoms == 27
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

    def test_sn2_ts_orca_output(self, orca_sn2_ts_output):
        orca_out = ORCAOutput(filename=orca_sn2_ts_output)
        assert orca_out.spin == "restricted"
        assert orca_out.forces is not None
        optimized_geometry = orca_out.get_optimized_parameters()
        assert isinstance(optimized_geometry, dict)
        assert optimized_geometry["B(Cl2,C1)"] == 2.0226
        assert optimized_geometry["B(F6,C1)"] == 2.1106
        molecule = orca_out.molecule
        assert isinstance(molecule, Molecule)
        assert molecule.symbols == ["C", "Cl", "H", "H", "H", "F"]
        assert molecule.empirical_formula == "CH3ClF"
        assert np.allclose(
            molecule.positions,
            np.array(
                [
                    [
                        [-0.000000e00, -1.080000e-03, -4.434500e-02],
                        [0.000000e00, 1.833000e-03, 1.978187e00],
                        [-0.000000e00, 1.057065e00, -2.682420e-01],
                        [9.166220e-01, -5.311370e-01, -2.657160e-01],
                        [-9.166220e-01, -5.311370e-01, -2.657160e-01],
                        [1.000000e-06, 4.457000e-03, -2.155338e00],
                    ]
                ]
            ),
            rtol=1e-4,
        )
        assert math.isclose(molecule.energy, -599.599020937503, rel_tol=1e-4)
        assert math.isclose(
            orca_out.final_nuclear_repulsion, 86.95944182533995, rel_tol=1e-4
        )
        assert math.isclose(
            orca_out.final_electronic_energy, -686.55844243196202, rel_tol=1e-4
        )
        assert math.isclose(
            orca_out.one_electron_energy, -1010.54352403257508, rel_tol=1e-4
        )
        assert math.isclose(
            orca_out.two_electron_energy, 323.98508160061300, rel_tol=1e-4
        )
        assert orca_out.max_cosx_asymmetry_energy is None
        assert math.isclose(
            orca_out.potential_energy, -1197.54777156261616, rel_tol=1e-4
        )
        assert math.isclose(
            orca_out.kinetic_energy, 597.94875062511289, rel_tol=1e-4
        )
        assert math.isclose(
            orca_out.virial_ratio, 2.00275988587762, rel_tol=1e-4
        )
        assert math.isclose(orca_out.xc_energy, -21.652294220893, rel_tol=1e-4)
        assert orca_out.dfet_embed_energy is None
        assert orca_out.orbital_occupancy == [
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
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
        ]
        orbital_energies = np.array(
            [
                -2786.2790,
                -674.8933,
                -283.4531,
                -260.2004,
                -198.3255,
                -198.2058,
                -198.2058,
                -22.5055,
                -18.7812,
                -15.0033,
                -8.3983,
                -8.3973,
                -6.1101,
                -4.2133,
                -4.2132,
                -2.2556,
                -2.2539,
                -1.8648,
                7.2164,
                8.2028,
                9.7068,
                9.7163,
                18.0067,
                18.0069,
                20.3332,
                22.5511,
                24.8903,
                24.8961,
                25.9780,
            ]
        )
        assert np.allclose(
            orca_out.orbital_energies, orbital_energies, rtol=1e-4
        )

        # test HOMO LUMO

        assert np.isclose(orca_out.homo_energy, -1.8648, rtol=1e-4)
        assert np.isclose(orca_out.lumo_energy, 7.2164, rtol=1e-4)
        assert np.isclose(orca_out.fmo_gap, 9.08117, rtol=1e-4)

        assert orca_out.mulliken_atomic_charges == {
            "C1": 0.074477,
            "Cl2": -0.503787,
            "H3": 0.059299,
            "H4": 0.058549,
            "H5": 0.058549,
            "F6": -0.747088,
        }
        assert orca_out.loewdin_atomic_charges == {
            "C1": 0.007991,
            "Cl2": -0.333494,
            "H3": 0.015310,
            "H4": 0.015231,
            "H5": 0.015231,
            "F6": -0.720269,
        }
        assert orca_out.mayer_mulliken_gross_atomic_population == {
            "C1": 5.9255,
            "Cl2": 17.5038,
            "H3": 0.9407,
            "H4": 0.9415,
            "H5": 0.9415,
            "F6": 9.7471,
        }
        assert orca_out.mayer_total_nuclear_charge == {
            "C1": 6.0000,
            "Cl2": 17.0000,
            "H3": 1.0000,
            "H4": 1.0000,
            "H5": 1.0000,
            "F6": 9.0000,
        }
        assert orca_out.mayer_mulliken_gross_atomic_charge == {
            "C1": 0.0745,
            "Cl2": -0.5038,
            "H3": 0.0593,
            "H4": 0.0585,
            "H5": 0.0585,
            "F6": -0.7471,
        }
        assert orca_out.mayer_total_valence == {
            "C1": 3.8739,
            "Cl2": 0.7224,
            "H3": 1.0198,
            "H4": 1.0190,
            "H5": 1.0190,
            "F6": 0.4499,
        }
        assert orca_out.mayer_bonded_valence == {
            "C1": 3.8739,
            "Cl2": 0.7224,
            "H3": 1.0198,
            "H4": 1.0190,
            "H5": 1.0190,
            "F6": 0.4499,
        }
        assert orca_out.mayer_free_valence == {
            "C1": 0.0000,
            "Cl2": 0.0000,
            "H3": 0.0000,
            "H4": 0.0000,
            "H5": 0.0000,
            "F6": 0.0000,
        }
        assert orca_out.mayer_bond_orders_larger_than_zero_point_one == {
            "B(C1,H3)": 0.9819,
            "B(C1,H4)": 0.9823,
            "B(C1,H5)": 0.9823,
            "B(C1,F6)": 0.2871,
        }
        assert np.allclose(
            orca_out.dipole_moment_electric_contribution,
            np.array([0.000000936, 0.007520938, 3.440735059]),
            rtol=1e-9,
        )
        assert np.allclose(
            orca_out.dipole_moment_nuclear_contribution,
            np.array([-0.000001417, -0.01222165, -1.626726359]),
            rtol=1e-9,
        )
        assert np.allclose(
            orca_out.total_dipole_moment,
            np.array([-0.000000481, -0.004700712, 1.814008700]),
            rtol=1e-9,
        )
        assert orca_out.dipole_moment_in_au == 1.814014790
        assert orca_out.dipole_moment_in_debye == 4.610859166
        assert np.allclose(
            orca_out.dipole_moment_along_axis_in_au,
            np.array([1.814011, -0.003782, 0.000001]),
            rtol=1e-6,
        )
        assert np.allclose(
            orca_out.dipole_moment_along_axis_in_debye,
            np.array([4.610849, -0.009613, 0.000003]),
            rtol=1e-6,
        )
        assert orca_out.rotational_constants_in_wavenumbers == [
            4.973563,
            0.077423,
            0.077423,
        ]
        assert orca_out.rotational_constants_in_MHz == [
            149103.667270,
            2321.080687,
            2321.072706,
        ]
        assert molecule.vibrational_frequencies == [
            -407.58,
            212.75,
            214.67,
            308.27,
            882.84,
            884.88,
            1062.79,
            1361.79,
            1363.05,
            3209.04,
            3397.13,
            3398.56,
        ]
        mode1 = np.array(
            [
                [-0.00000e00, -4.79000e-03, 8.59726e-01],
                [-0.00000e00, 7.86000e-04, -1.18060e-01],
                [1.00000e-06, -1.97735e-01, -1.49549e-01],
                [
                    -1.62414e-01,
                    8.02020e-02,
                    -1.18726e-01,
                ],
                [1.62414e-01, 8.02010e-02, -1.18726e-01],
                [0.00000e00, 3.54200e-03, -3.02689e-01],
            ],
            dtype=float,
        )
        mode2 = np.array(
            [
                [4.28110e-01, 3.00000e-06, -1.00000e-06],
                [-9.02300e-02, -1.00000e-06, 1.00000e-06],
                [4.82113e-01, 4.00000e-06, 2.00000e-06],
                [4.79999e-01, -1.44490e-02, 2.03671e-01],
                [4.79998e-01, 1.44560e-02, -2.03671e-01],
                [-1.78795e-01, -1.00000e-06, -2.00000e-06],
            ]
        )
        last_mode = np.array(
            [
                [
                    [-1.01365e-01, -2.20000e-05, -0.00000e00],
                    [4.64000e-04, 0.00000e00, -0.00000e00],
                    [-9.63800e-03, 1.71000e-04, -3.50000e-05],
                    [5.94810e-01, -3.49328e-01, -1.37333e-01],
                    [5.94964e-01, 3.49413e-01, 1.37368e-01],
                    [6.04000e-04, 0.00000e00, 0.00000e00],
                ],
            ]
        )
        assert len(orca_out.vibrational_modes) == 12
        assert np.allclose(orca_out.vibrational_modes[0], mode1, rtol=1e-4)
        assert np.allclose(orca_out.vibrational_modes[1], mode2, rtol=1e-4)
        assert np.allclose(
            orca_out.vibrational_modes[-1], last_mode, rtol=1e-4
        )

        # test for molecule obtained from output file
        mol = orca_out.molecule
        assert len(mol.vibrational_modes) == 12
        assert np.allclose(mol.vibrational_modes[0], mode1, rtol=1e-4)
        assert np.allclose(mol.vibrational_modes[1], mode2, rtol=1e-4)
        assert np.allclose(mol.vibrational_modes[-1], last_mode, rtol=1e-4)

        assert orca_out.vib_freq_scale_factor == 1.0
        assert orca_out.molar_absorption_coefficients == [
            0.003487,
            0.003475,
            0.000119,
            0.000054,
            0.000062,
            0.001142,
            0.001459,
            0.001641,
            0.000142,
            0.000125,
            0.000119,
        ]
        assert orca_out.integrated_absorption_coefficients == [
            17.62,
            17.56,
            0.6,
            0.27,
            0.31,
            5.77,
            7.38,
            8.29,
            0.72,
            0.63,
            0.6,
        ]
        assert orca_out.transition_dipole_deriv_norm == [
            0.005115,
            0.005052,
            0.00012,
            1.9e-05,
            2.2e-05,
            0.000335,
            0.000334,
            0.000376,
            1.4e-05,
            1.2e-05,
            1.1e-05,
        ]
        transition_dipole1 = np.array([0.071520, 0.000001, -0.000000])
        transition_dipole2 = np.array([0.000001, -0.070961, -0.004067])
        transition_dipole3 = np.array([0.000000, 0.000216, 0.010974])
        assert np.allclose(
            orca_out.transition_dipoles[0], transition_dipole1, rtol=1e-4
        )
        assert np.allclose(
            orca_out.transition_dipoles[1], transition_dipole2, rtol=1e-4
        )
        assert np.allclose(
            orca_out.transition_dipoles[2], transition_dipole3, rtol=1e-4
        )

        assert orca_out.num_translation_and_rotation_modes == 6
        assert orca_out.num_vibration_modes == 12
        assert orca_out.temperature_in_K == 298.15
        assert orca_out.pressure_in_atm == 1.0
        assert orca_out.total_mass_in_amu == 69.49
        assert math.isclose(
            orca_out.internal_energy, -599.55741380, rel_tol=1e-4
        )  # in Hartrees, default unit in ORCA output file
        assert math.isclose(
            orca_out.electronic_energy, -599.59902094, rel_tol=1e-4
        )
        assert math.isclose(
            orca_out.zero_point_energy, 0.03712451, rel_tol=1e-8
        )
        assert math.isclose(
            orca_out.thermal_vibration_correction,
            0.00165009,
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
        assert math.isclose(orca_out.enthalpy, -599.55646959, rel_tol=1e-4)
        assert math.isclose(
            orca_out.thermal_enthalpy_correction,
            0.00094421,
            rel_tol=1e-8,
        )
        assert orca_out.electronic_entropy_no_temperature_in_SI == 0.0
        assert math.isclose(
            orca_out.vibrational_entropy_no_temperature_in_SI,
            0.00276525 * units.Hartree / (units.J / units.mol) / 298.15,
            rel_tol=1e-4,
        )
        assert math.isclose(
            orca_out.rotational_entropy_no_temperature_in_SI,
            0.01116917 * units.Hartree / (units.J / units.mol),
            rel_tol=1,
        )
        assert math.isclose(
            orca_out.translational_entropy_no_temperature_in_SI,
            0.01835566 * units.Hartree / (units.J / units.mol),
            rel_tol=1,
        )
        assert math.isclose(orca_out.entropy_TS, 0.03229008, rel_tol=1e-4)

        entropy_TS_in_J_per_mol = (
            0.03229008 * units.Hartree / (units.J / units.mol)
        )
        TS = orca_out.entropy_in_J_per_mol_per_K * 298.15  # 56266.794 J/mol
        # converted to 13.448 kcal/mol, as expected from output
        assert math.isclose(entropy_TS_in_J_per_mol, TS, rel_tol=1e-4)
        assert (
            orca_out.rotational_entropy_symmetry_correction_J_per_mol_per_K
            == {
                "sn=1": 98.355364,
                "sn=2": 92.592298,
                "sn=3": 89.221109,
                "sn=4": 86.829143,
                "sn=5": 84.973814,
                "sn=6": 83.457954,
                "sn=7": 82.176245,
                "sn=8": 81.065989,
                "sn=9": 80.086765,
                "sn=10": 79.210747,
                "sn=11": 78.418298,
                "sn=12": 77.6948,
                "Dinfh": 68.171144,
                "Cinfv": 73.934299,
            }
        )

        assert math.isclose(
            orca_out.gibbs_free_energy, -599.58875967, rel_tol=1e-8
        )
        assert isinstance(orca_out.molecule, Molecule)
        assert orca_out.total_elapsed_walltime == 0.0

    def test_fe2_singlet_orbital_properties(self, fe2_singlet_output):
        """Test HOMO/LUMO properties for Fe2 singlet state."""
        orca_out = ORCAOutput(filename=fe2_singlet_output)
        assert orca_out.multiplicity == 1  # Singlet state
        assert orca_out.spin == "restricted"
        assert orca_out.num_unpaired_electrons == 0

        # SOMO properties (should be None for closed-shell)
        assert orca_out.somo_energies is None
        assert orca_out.lowest_somo_energy is None
        assert orca_out.highest_somo_energy is None

        # HOMO/LUMO properties
        assert orca_out.homo_energy == -0.738698 * units.Hartree
        assert orca_out.lumo_energy == -0.403510 * units.Hartree
        assert orca_out.fmo_gap == (-0.403510 - (-0.738698)) * units.Hartree

    def test_fe2_triplet_orbital_properties(self, fe2_triplet_output):
        """Test HOMO/LUMO/SOMO properties for Fe2 triplet state."""
        orca_out = ORCAOutput(filename=fe2_triplet_output)
        assert orca_out.multiplicity == 3  # Triplet state
        assert orca_out.spin == "unrestricted"
        assert orca_out.num_unpaired_electrons == 2

        # SOMO properties
        assert orca_out.somo_energies == [
            -0.767986 * units.Hartree,
            -0.761816 * units.Hartree,
        ]
        assert len(orca_out.somo_energies) == 2  # 2 unpaired electrons
        assert orca_out.lowest_somo_energy == -0.767986 * units.Hartree
        assert orca_out.highest_somo_energy == -0.761816 * units.Hartree

        # HOMO/LUMO properties for open-shell
        assert orca_out.alpha_homo_energy == -0.761816 * units.Hartree
        assert orca_out.beta_homo_energy == -0.721020 * units.Hartree
        assert orca_out.alpha_lumo_energy == -0.346513 * units.Hartree
        assert orca_out.beta_lumo_energy == -0.381082 * units.Hartree
        assert (
            orca_out.fmo_gap
            == (min(-0.346513, -0.381082) - (-0.761816)) * units.Hartree
        )
        assert np.isclose(
            orca_out.alpha_fmo_gap,
            (-0.346513 - (-0.761816)) * units.Hartree,
            rtol=1e-6,
        )
        assert np.isclose(
            orca_out.beta_fmo_gap,
            (-0.381082 - (-0.721020)) * units.Hartree,
            rtol=1e-6,
        )

    def test_fe2_quintet_orbital_properties(self, fe2_quintet_output):
        """Test HOMO/LUMO/SOMO properties for Fe2 quintet state."""
        orca_out = ORCAOutput(filename=fe2_quintet_output)
        assert orca_out.multiplicity == 5  # Quintet state
        assert orca_out.spin == "unrestricted"
        assert orca_out.num_unpaired_electrons == 4

        # SOMO properties
        assert orca_out.somo_energies == [
            -0.819921 * units.Hartree,
            -0.806437 * units.Hartree,
            -0.802013 * units.Hartree,
            -0.705655 * units.Hartree,
        ]
        assert len(orca_out.somo_energies) == 4  # 4 unpaired electrons
        assert orca_out.lowest_somo_energy == -0.819921 * units.Hartree
        assert orca_out.highest_somo_energy == -0.705655 * units.Hartree

        # HOMO/LUMO properties
        assert orca_out.alpha_homo_energy == -0.705655 * units.Hartree
        assert orca_out.beta_homo_energy == -0.688110 * units.Hartree
        assert orca_out.alpha_lumo_energy == -0.336304 * units.Hartree
        assert orca_out.beta_lumo_energy == -0.328840 * units.Hartree
        assert np.isclose(
            orca_out.fmo_gap,
            (min(-0.336304, -0.328840) - (-0.705655)) * units.Hartree,
            # 0.369351
        )
        assert np.isclose(
            orca_out.alpha_fmo_gap,
            (-0.336304 - (-0.705655)) * units.Hartree,
            rtol=1e-6,
        )
        assert np.isclose(
            orca_out.beta_fmo_gap,
            (-0.328840 - (-0.688110)) * units.Hartree,
            rtol=1e-6,
        )

    def test_fe3_doublet_orbital_properties(self, fe3_doublet_output):
        """Test HOMO/LUMO/SOMO properties for Fe3 doublet state."""
        orca_out = ORCAOutput(filename=fe3_doublet_output)
        assert orca_out.multiplicity == 2  # Doublet state
        assert orca_out.spin == "unrestricted"
        assert orca_out.num_unpaired_electrons == 1

        # SOMO properties
        assert orca_out.somo_energies == [-1.060149 * units.Hartree]
        assert len(orca_out.somo_energies) == 1  # 1 unpaired electron
        assert orca_out.lowest_somo_energy == -1.060149 * units.Hartree
        assert orca_out.highest_somo_energy == -1.060149 * units.Hartree

        # HOMO/LUMO properties
        assert orca_out.alpha_homo_energy == -1.060149 * units.Hartree
        assert orca_out.beta_homo_energy == -1.061789 * units.Hartree
        assert orca_out.alpha_lumo_energy == -0.764634 * units.Hartree
        assert orca_out.beta_lumo_energy == -0.744986 * units.Hartree
        assert np.isclose(
            orca_out.fmo_gap,
            (min(-0.764634, -0.744986) - (-1.060149)) * units.Hartree,
        )
        assert np.isclose(
            orca_out.alpha_fmo_gap,
            (-0.764634 - (-1.060149)) * units.Hartree,
            rtol=1e-6,
        )
        assert np.isclose(
            orca_out.beta_fmo_gap,
            (-0.744986 - (-1.061789)) * units.Hartree,
            rtol=1e-6,
        )

    def test_fe3_quartet_orbital_properties(self, fe3_quartet_output):
        """Test HOMO/LUMO/SOMO properties for Fe3 quartet state."""
        orca_out = ORCAOutput(filename=fe3_quartet_output)
        assert orca_out.multiplicity == 4  # Quartet state
        assert orca_out.spin == "unrestricted"
        assert orca_out.num_unpaired_electrons == 3

        # SOMO properties
        assert orca_out.somo_energies == [
            -1.040564 * units.Hartree,
            -1.037095 * units.Hartree,
            -1.034892 * units.Hartree,
        ]
        assert len(orca_out.somo_energies) == 3  # 3 unpaired electrons
        assert orca_out.lowest_somo_energy == -1.040564 * units.Hartree
        assert orca_out.highest_somo_energy == -1.034892 * units.Hartree

        # HOMO/LUMO properties
        assert orca_out.alpha_homo_energy == -1.034892 * units.Hartree
        assert orca_out.beta_homo_energy == -1.024711 * units.Hartree
        assert orca_out.alpha_lumo_energy == -0.749293 * units.Hartree
        assert orca_out.beta_lumo_energy == -0.721360 * units.Hartree
        assert np.isclose(
            orca_out.fmo_gap,
            (min(-0.749293, -0.721360) - (-1.034892)) * units.Hartree,
        )
        assert np.isclose(
            orca_out.alpha_fmo_gap,
            (-0.749293 - (-1.034892)) * units.Hartree,
            rtol=1e-6,
        )
        assert np.isclose(
            orca_out.beta_fmo_gap,
            (-0.721360 - (-1.024711)) * units.Hartree,
            rtol=1e-6,
        )

    def test_fe3_sextet_orbital_properties(self, fe3_sextet_output):
        """Test HOMO/LUMO/SOMO properties for Fe3 sextet state."""
        orca_out = ORCAOutput(filename=fe3_sextet_output)
        assert orca_out.multiplicity == 6  # Sextet state
        assert orca_out.spin == "unrestricted"
        assert orca_out.num_unpaired_electrons == 5

        # SOMO properties
        assert orca_out.somo_energies == [
            -1.114592 * units.Hartree,
            -1.093714 * units.Hartree,
            -1.080409 * units.Hartree,
            -1.066772 * units.Hartree,
            -1.030076 * units.Hartree,
        ]
        assert len(orca_out.somo_energies) == 5  # 5 unpaired electrons
        assert orca_out.lowest_somo_energy == -1.114592 * units.Hartree
        assert orca_out.highest_somo_energy == -1.030076 * units.Hartree

        # HOMO/LUMO properties
        assert orca_out.alpha_homo_energy == -1.030076 * units.Hartree
        assert orca_out.beta_homo_energy == -1.045382 * units.Hartree
        assert orca_out.alpha_lumo_energy == -0.553561 * units.Hartree
        assert orca_out.beta_lumo_energy == -0.775129 * units.Hartree
        assert np.isclose(
            orca_out.fmo_gap,
            (min(-0.553561, -0.775129) - (-1.030076)) * units.Hartree,
        )
        assert np.isclose(
            orca_out.alpha_fmo_gap,
            (-0.553561 - (-1.030076)) * units.Hartree,
            rtol=1e-6,
        )
        assert np.isclose(
            orca_out.beta_fmo_gap,
            (-0.775129 - (-1.045382)) * units.Hartree,
            rtol=1e-6,
        )


class TestORCAEngrad:
    def test_read_water_output(self, water_engrad_path):
        orca_engrad = ORCAEngradFile(filename=water_engrad_path)
        assert orca_engrad.num_atoms == 3
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


class TestORCAQMMM:
    def test_read_qmmm_output(self, orca_two_layer_qmmmm_output_file):
        orca_qmmm1 = ORCAQMMMOutput(filename=orca_two_layer_qmmmm_output_file)
        assert orca_qmmm1.multiscale_model == "QM1/QM2"
        assert orca_qmmm1.qm2_method == "XTB2"
        assert orca_qmmm1.total_charge == 0
        assert orca_qmmm1.scaling_factor_qm2 == 1.0
        assert orca_qmmm1.point_charges_in_qm_from_mm == 24
        assert orca_qmmm1.point_charges_in_qm_from_charge_shift == 0
        assert orca_qmmm1.total_system_size == 36
        assert orca_qmmm1.qm_system_size == 12
        assert orca_qmmm1.qm2_system_size == 24
        assert orca_qmmm1.number_of_link_atoms == 0
        assert orca_qmmm1.qm_plus_link_atoms_size == 12
        assert orca_qmmm1.qm_region == ["1-12"]
        assert orca_qmmm1.qm2_energy_of_large_region == -994.9374837306615
        assert orca_qmmm1.qm2_energy_of_small_region == -396.0605045891306
        assert orca_qmmm1.qm_qm2_energy == -5889.533884098047
        assert orca_qmmm1.qm_energy == -5290.656904956516


class TestORCAQMMMJobSettings:
    def test_partition_string_single_and_list_input(self):
        """Partition string should accept '1-15,37,39' or list and compress ranges."""
        from chemsmart.jobs.orca.settings import ORCAQMMMJobSettings

        s = ORCAQMMMJobSettings()
        s.high_level_atoms = "1-15,37,39"
        out = s._get_partition_string()
        assert out.strip() == "QMAtoms {0:14 36 38} end"

        # list input should produce the same output
        s.high_level_atoms = [1, *range(2, 16), 37, 39]
        out2 = s._get_partition_string()
        assert out2.strip() == "QMAtoms {0:14 36 38} end"

    def test_partition_string_qm_and_qm2(self):
        """When both high_level_atoms and medium_level_atoms provided, both lines should be returned."""
        from chemsmart.jobs.orca.settings import ORCAQMMMJobSettings

        s = ORCAQMMMJobSettings()
        s.high_level_atoms = "1-3,5"
        s.intermediate_level_atoms = "7-9,12"
        out = s._get_partition_string()
        # order: QMAtoms then QM2Atoms
        lines = [ln for ln in out.splitlines() if ln.strip()]
        assert lines[0].strip() == "QMAtoms {0:2 4} end"
        assert lines[1].strip() == "QM2Atoms {6:8 11} end"

    def test_charge_and_multiplicity_population(self):
        """ORCAQMMMJobSettings should populate .charge and .multiplicity from intermediate or high fields."""
        from chemsmart.jobs.orca.settings import ORCAQMMMJobSettings

        # when both high and intermediate specified, high-region takes precedence
        s1 = ORCAQMMMJobSettings(
            charge_intermediate=0,
            mult_intermediate=1,
            charge_high=2,
            mult_high=3,
        )
        assert s1.charge == 2
        assert s1.multiplicity == 3

        # intermediate missing -> fall back to high
        s2 = ORCAQMMMJobSettings(charge_high=-1, mult_high=2)
        assert s2.charge == -1
        assert s2.multiplicity == 2

    def test_partition_string_empty_and_none(self):
        """Empty string or None should return empty partition block."""
        from chemsmart.jobs.orca.settings import ORCAQMMMJobSettings

        s = ORCAQMMMJobSettings()
        s.high_level_atoms = ""
        assert s._get_partition_string() == ""

        s.high_level_atoms = None
        assert s._get_partition_string() == ""


class TestORCANEB:
    def test_read_neb_output(self, orca_neb_output_file):
        import pathlib

        src = pathlib.Path(orca_neb_output_file)

        # Read with tolerant UTF-8 decoding and inject into the parser to avoid locale decoding errors
        data = src.read_text(encoding="utf-8", errors="replace")
        lines = [ln.strip() for ln in data.splitlines()]

        orca_neb = ORCANEBFile(filename=str(src))
        orca_neb.__dict__["contents"] = lines
        orca_neb.__dict__["content_lines_string"] = data
        assert orca_neb.nimages == 10
        assert orca_neb.num_atoms == 148
        assert orca_neb.ci_converged is True
        assert orca_neb.ts_converged is True
        assert orca_neb.ci == "Climbing Image:  image 4."
        assert orca_neb.ci_energy == -219.0833212
        assert (
            orca_neb.reactant.empirical_formula
            == orca_neb.product.empirical_formula
            == "C72H68NO6P"
        )
        assert orca_neb.ci_max_abs_force == 0.001963
        assert orca_neb.ts_delta_energy == 4.02
        assert orca_neb.ts_rms_force == 0.00034
        assert orca_neb.ts_max_abs_force == 0.00543
        assert orca_neb.ts_energy == -219.09056
        assert orca_neb.preopt_ends


class TestORCANEBJobSettings:
    """Test suite for ORCANEBJobSettings class."""

    def test_init_default(self):
        """Test default initialization."""
        settings = ORCANEBJobSettings()
        assert settings.joboption is None
        assert settings.nimages is None
        assert settings.preopt_ends is False

    def test_init_with_parameters(self):
        """Test initialization with parameters."""
        settings = ORCANEBJobSettings(
            joboption="NEB-TS", nimages=8, semiempirical="XTB2"
        )
        assert settings.joboption == "NEB-TS"
        assert settings.nimages == 8
        assert settings.semiempirical == "XTB2"

    def test_route_string_generation(self):
        """Test route string generation."""
        settings = ORCANEBJobSettings(joboption="NEB-CI", semiempirical="XTB2")
        assert settings.route_string == "!  XTB2 NEB-CI"

    def test_neb_block_basic(self):
        """Test basic NEB block generation."""
        settings = ORCANEBJobSettings(
            nimages=5, starting_xyz="start.xyz", ending_xyzfile="end.xyz"
        )
        neb_block = settings.neb_block
        assert "%neb" in neb_block
        assert "NImages 5" in neb_block
        assert 'NEB_END_XYZFile "end.xyz"' in neb_block

    def test_inheritance(self):
        """Test inheritance from ORCAJobSettings."""
        settings = ORCANEBJobSettings(functional="B3LYP", basis="def2-SVP")
        assert isinstance(settings, ORCANEBJobSettings)
        assert settings.functional == "B3LYP"
        assert settings.basis == "def2-SVP"

    def test_validation_errors(self):
        """Test validation raises appropriate errors."""
        settings = ORCANEBJobSettings(starting_xyz="start.xyz")
        with pytest.raises(
            AssertionError, match="The number of images is missing"
        ):
            _ = settings.neb_block

    def test_equality_identical_settings(self):
        """Test that identical NEB settings are equal."""
        settings1 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            intermediate_xyzfile="ts_guess.xyz",
            preopt_ends=True,
            semiempirical="XTB2",
            functional="B3LYP",
            basis="def2-SVP",
            charge=0,
            multiplicity=1,
        )
        settings2 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            intermediate_xyzfile="ts_guess.xyz",
            preopt_ends=True,
            semiempirical="XTB2",
            functional="B3LYP",
            basis="def2-SVP",
            charge=0,
            multiplicity=1,
        )
        assert settings1 == settings2
        assert not (settings1 != settings2)

    def test_equality_different_joboption(self):
        """Test that settings with different joboption are not equal."""
        settings1 = ORCANEBJobSettings(
            joboption="NEB-TS", nimages=8, ending_xyzfile="product.xyz"
        )
        settings2 = ORCANEBJobSettings(
            joboption="NEB-CI", nimages=8, ending_xyzfile="product.xyz"
        )
        assert settings1 != settings2
        assert not (settings1 == settings2)

    def test_equality_different_nimages(self):
        """Test that settings with different nimages are not equal."""
        settings1 = ORCANEBJobSettings(
            joboption="NEB-TS", nimages=8, ending_xyzfile="product.xyz"
        )
        settings2 = ORCANEBJobSettings(
            joboption="NEB-TS", nimages=12, ending_xyzfile="product.xyz"
        )
        assert settings1 != settings2

    def test_equality_different_intermediate_xyzfile(self):
        """Test that settings with different intermediate file are not equal."""
        settings1 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            intermediate_xyzfile="ts1.xyz",
        )
        settings2 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            intermediate_xyzfile="ts2.xyz",
        )
        assert settings1 != settings2

    def test_equality_different_restarting_xyzfile(self):
        """Test that settings with different restart file are not equal."""
        settings1 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            restarting_xyzfile="restart1.allxyz",
        )
        settings2 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            restarting_xyzfile="restart2.allxyz",
        )
        assert settings1 != settings2

    def test_equality_different_preopt_ends(self):
        """Test that settings with different preopt_ends are not equal."""
        settings1 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            preopt_ends=True,
        )
        settings2 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            preopt_ends=False,
        )
        assert settings1 != settings2

    def test_equality_different_semiempirical(self):
        """Test that settings with different semiempirical method are not equal."""
        settings1 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            semiempirical="XTB2",
        )
        settings2 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            semiempirical="XTB1",
        )
        assert settings1 != settings2

    def test_equality_different_parent_attributes(self):
        """Test that settings with different parent class attributes are not equal."""
        settings1 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            functional="B3LYP",
            basis="def2-SVP",
        )
        settings2 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            functional="PBE0",
            basis="def2-SVP",
        )
        assert settings1 != settings2

    def test_equality_includes_all_neb_attributes(self):
        """Test that all 7 NEB-specific attributes are included in equality check."""
        # Create two settings that differ only in each NEB attribute
        base_kwargs = {
            "joboption": "NEB-TS",
            "nimages": 8,
            "ending_xyzfile": "product.xyz",
            "intermediate_xyzfile": "ts.xyz",
            "restarting_xyzfile": "restart.allxyz",
            "preopt_ends": True,
            "semiempirical": "XTB2",
        }

        # Test each attribute individually
        for attr in [
            "joboption",
            "nimages",
            "ending_xyzfile",
            "intermediate_xyzfile",
            "restarting_xyzfile",
            "preopt_ends",
            "semiempirical",
        ]:
            settings1 = ORCANEBJobSettings(**base_kwargs)
            modified_kwargs = base_kwargs.copy()

            # Change the attribute to a different value
            if attr == "joboption":
                modified_kwargs[attr] = "NEB-CI"
            elif attr == "nimages":
                modified_kwargs[attr] = 12
            elif attr == "preopt_ends":
                modified_kwargs[attr] = False
            elif attr == "semiempirical":
                modified_kwargs[attr] = "XTB1"
            else:  # file paths
                modified_kwargs[attr] = "different_file.xyz"

            settings2 = ORCANEBJobSettings(**modified_kwargs)
            assert (
                settings1 != settings2
            ), f"Equality failed for attribute: {attr}"
