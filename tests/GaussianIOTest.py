import os.path
import numpy as np
from chemsmart.io.gaussian.inputs import Gaussian16Input
from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.io.gaussian.output import Gaussian16WBIOutput
from chemsmart.io.gaussian.cube import GaussianCubeFile
from ase.symbols import Symbols
from ase import units


class TestGaussian16Input:
    def test_read_gaussian_input(self, gaussian_opt_inputfile):
        assert os.path.exists(gaussian_opt_inputfile)
        g16_input = Gaussian16Input(filename=gaussian_opt_inputfile)
        assert g16_input.molecule.chemical_symbols == [
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "H",
            "H",
            "H",
            "H",
            "C",
            "O",
            "H",
            "Cl",
        ]  # list of chemical symbols
        assert isinstance(g16_input.molecule.symbols, Symbols)
        assert g16_input.molecule.symbols.formula == "C6H4COHCl"
        assert g16_input.molecule.natoms == 14
        assert g16_input.molecule.empirical_formula == "C7H5ClO"
        assert all(
            np.isclose(
                g16_input.molecule.positions[0],
                [-0.5448210000, -1.1694570000, 0.0001270000],
                atol=10e-5,
            )
        )
        assert g16_input.additional_opt_options_in_route is None
        assert g16_input.additional_route_parameters is None
        assert g16_input.job_type == "opt"
        assert g16_input.functional == "m062x"
        assert g16_input.basis == "def2svp"
        assert g16_input.molecule.frozen_atoms is None

    def test_read_frozen_coords(self, gaussian_frozen_opt_inputfile):
        assert os.path.exists(gaussian_frozen_opt_inputfile)
        g16_frozen = Gaussian16Input(filename=gaussian_frozen_opt_inputfile)
        assert g16_frozen.molecule.symbols.formula == "C6H4COHCl"
        assert g16_frozen.molecule.empirical_formula == "C7H5ClO"
        assert g16_frozen.additional_opt_options_in_route is None
        assert g16_frozen.additional_route_parameters is None
        assert g16_frozen.job_type == "opt"

    def test_read_modred_inputfile(self, gaussian_modred_inputfile):
        assert os.path.exists(gaussian_modred_inputfile)
        g16_modred = Gaussian16Input(filename=gaussian_modred_inputfile)
        assert g16_modred.molecule.chemical_symbols == [
            "O",
            "N",
            "C",
            "C",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "C",
            "O",
            "O",
        ]  # list of chemical symbols
        assert g16_modred.molecule.symbols.formula == "ONC2H7CO2"
        assert g16_modred.molecule.empirical_formula == "C3H7NO3"
        assert g16_modred.additional_opt_options_in_route is None
        assert g16_modred.additional_route_parameters is None
        assert g16_modred.job_type == "modred"
        assert g16_modred.modredundant == [[2, 12], [9, 2]]
        assert g16_modred.functional == "m062x"
        assert g16_modred.basis == "def2svp"

    def test_read_scan_inputfile(self, gaussian_scan_inputfile):
        assert os.path.exists(gaussian_scan_inputfile)
        g16_scan = Gaussian16Input(filename=gaussian_scan_inputfile)
        assert g16_scan.molecule.chemical_symbols == [
            "O",
            "N",
            "C",
            "C",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "C",
            "O",
            "O",
        ]  # list of chemical symbols
        assert g16_scan.molecule.symbols.formula == "ONC2H7CO2"
        assert g16_scan.molecule.empirical_formula == "C3H7NO3"
        assert g16_scan.additional_opt_options_in_route is None
        assert g16_scan.additional_route_parameters is None
        assert g16_scan.job_type == "modred"
        assert g16_scan.modredundant == {
            "coords": [[2, 12], [9, 2]],
            "num_steps": 10,
            "step_size": 0.05,
        }
        assert g16_scan.functional == "m062x"
        assert g16_scan.basis == "def2svp"

    def test_pbc_1d_input(self, gaussian_pbc_1d_inputfile):
        assert os.path.exists(gaussian_pbc_1d_inputfile)
        g16_pbc_1d = Gaussian16Input(filename=gaussian_pbc_1d_inputfile)
        assert g16_pbc_1d.molecule.symbols.formula == "CH2CHC2H2Cl"
        assert g16_pbc_1d.molecule.empirical_formula == "C4H5Cl"
        assert all(
            np.isclose(
                g16_pbc_1d.molecule.positions[-1],
                [0.62098257, 0.98609446, -1.78763987],
                atol=1e-5,
            )
        )
        assert g16_pbc_1d.additional_opt_options_in_route is None
        assert g16_pbc_1d.additional_route_parameters is None
        assert g16_pbc_1d.job_type == "sp"
        assert g16_pbc_1d.modredundant is None
        assert g16_pbc_1d.functional == "pbepbe"
        assert g16_pbc_1d.basis == "6-31g(d,p)/auto"


class TestGaussian16Output:
    def test_normal_termination_with_forces_and_frequencies(
        self, td_outputfile
    ):
        assert os.path.exists(td_outputfile)
        g16_output = Gaussian16Output(filename=td_outputfile)
        assert g16_output.num_atoms == 49
        assert g16_output.tddft_transitions[0] == (0.7744, 1601.13, 0.0084)
        assert g16_output.tddft_transitions[1] == (1.0201, 1215.37, 0.0632)
        assert g16_output.excitation_energies_eV == [
            0.7744,
            1.0201,
            1.502,
            2.052,
            2.1157,
            2.4471,
            2.6665,
            2.8332,
            3.0814,
            3.2134,
            3.2777,
            3.3555,
            3.3963,
            3.5764,
            3.604,
            3.6596,
            3.6907,
            3.697,
            3.8718,
            3.9218,
            3.9461,
            3.9949,
            4.0171,
            4.0813,
            4.0981,
            4.1212,
            4.2337,
            4.3012,
            4.3178,
            4.3324,
            4.3623,
            4.4078,
            4.4256,
            4.4396,
            4.4734,
            4.486,
            4.4893,
            4.5261,
            4.5624,
            4.6544,
            4.6823,
            4.7346,
            4.7521,
            4.7704,
            4.798,
            4.8059,
            4.8211,
            4.8303,
            4.8511,
            4.8561,
        ]
        assert len(g16_output.excitation_energies_eV) == 50
        assert len(g16_output.transitions) == 50
        assert len(g16_output.contribution_coefficients) == 50
        assert g16_output.transitions[0] == [
            "104A ->108A",
            "105A ->107A",
            "106A ->107A",
            "106A ->108A",
            "105B ->106B",
            "106A <-107A",
        ]
        assert g16_output.contribution_coefficients[0] == [
            0.15573,
            -0.1244,
            0.93545,
            -0.10308,
            0.26021,
            0.12114,
        ]
        assert g16_output.contribution_coefficients[-1] == [
            -0.17274,
            0.14866,
            0.13926,
            -0.31107,
            0.79088,
            0.17825,
        ]

    def test_singlet_opt_output(self, gaussian_singlet_opt_outfile):
        assert os.path.exists(gaussian_singlet_opt_outfile)
        g16_output = Gaussian16Output(filename=gaussian_singlet_opt_outfile)
        assert g16_output.normal_termination
        assert g16_output.tddft_transitions == []  # no tddft calcs
        assert len(g16_output.alpha_occ_eigenvalues) == 116
        assert g16_output.alpha_occ_eigenvalues[0] == -25.29096 * units.Hartree
        assert g16_output.alpha_occ_eigenvalues[-1] == -0.29814 * units.Hartree
        assert len(g16_output.alpha_virtual_eigenvalues) == 378
        assert (
            g16_output.alpha_virtual_eigenvalues[0] == -0.02917 * units.Hartree
        )
        assert (
            g16_output.alpha_virtual_eigenvalues[-1]
            == 56.20437 * units.Hartree
        )
        assert g16_output.beta_occ_eigenvalues is None
        assert g16_output.beta_virtual_eigenvalues is None
        assert g16_output.homo_energy == -0.29814 * units.Hartree
        assert g16_output.lumo_energy == -0.02917 * units.Hartree
        assert np.isclose(g16_output.fmo_gap, 0.26897 * units.Hartree)

    def test_triplet_opt_output(self, gaussian_triplet_opt_outfile):
        assert os.path.exists(gaussian_triplet_opt_outfile)
        g16_output = Gaussian16Output(filename=gaussian_triplet_opt_outfile)
        assert g16_output.normal_termination
        assert g16_output.tddft_transitions == []  # no tddft calcs
        assert len(g16_output.alpha_occ_eigenvalues) == 215
        assert (
            g16_output.alpha_occ_eigenvalues[0] == -482.71377 * units.Hartree
        )
        assert g16_output.alpha_occ_eigenvalues[-1] == -0.15673 * units.Hartree
        assert len(g16_output.alpha_virtual_eigenvalues) == 750
        assert (
            g16_output.alpha_virtual_eigenvalues[0] == -0.07423 * units.Hartree
        )
        assert (
            g16_output.alpha_virtual_eigenvalues[-1] == 4.23682 * units.Hartree
        )
        assert len(g16_output.beta_occ_eigenvalues) == 213
        assert g16_output.beta_occ_eigenvalues[0] == -482.71362 * units.Hartree
        assert g16_output.beta_occ_eigenvalues[-1] == -0.18923 * units.Hartree
        assert len(g16_output.beta_virtual_eigenvalues) == 752
        assert (
            g16_output.beta_virtual_eigenvalues[0] == -0.05025 * units.Hartree
        )
        assert (
            g16_output.beta_virtual_eigenvalues[-1] == 4.26643 * units.Hartree
        )
        assert g16_output.homo_energy is None
        assert g16_output.lumo_energy is None
        assert g16_output.somo_energy == -0.15673 * units.Hartree

    def test_quintet_opt_output(self, gaussian_quintet_opt_outfile):
        assert os.path.exists(gaussian_quintet_opt_outfile)
        g16_output = Gaussian16Output(filename=gaussian_quintet_opt_outfile)
        assert g16_output.tddft_transitions == []  # no tddft calcs
        assert len(g16_output.alpha_occ_eigenvalues) == 216
        assert (
            g16_output.alpha_occ_eigenvalues[0] == -482.71572 * units.Hartree
        )
        assert g16_output.alpha_occ_eigenvalues[-1] == -0.18764 * units.Hartree
        assert len(g16_output.alpha_virtual_eigenvalues) == 749
        assert (
            g16_output.alpha_virtual_eigenvalues[0] == -0.03881 * units.Hartree
        )
        assert (
            g16_output.alpha_virtual_eigenvalues[-1] == 4.23318 * units.Hartree
        )
        assert len(g16_output.beta_occ_eigenvalues) == 212
        assert g16_output.beta_occ_eigenvalues[0] == -482.71538 * units.Hartree
        assert g16_output.beta_occ_eigenvalues[-1] == -0.19564 * units.Hartree
        assert len(g16_output.beta_virtual_eigenvalues) == 753
        assert (
            g16_output.beta_virtual_eigenvalues[0] == -0.06116 * units.Hartree
        )
        assert (
            g16_output.beta_virtual_eigenvalues[-1] == 4.23626 * units.Hartree
        )
        assert g16_output.somo_energy == -0.18764 * units.Hartree
        assert g16_output.fmo_gap is None


class TestGaussianWBIOutput:
    def test_normal_termination_with_forces_and_frequencies(
        self, wbi_outputfile
    ):
        assert os.path.exists(wbi_outputfile)
        g16_output = Gaussian16WBIOutput(filename=wbi_outputfile)
        assert g16_output.nbo_version == "3.1"
        assert len(g16_output.natural_atomic_orbitals) == 128
        assert len(g16_output.natural_atomic_orbitals["Ni1"]) == 31
        assert len(g16_output.natural_atomic_orbitals["P2"]) == 18
        assert len(g16_output.natural_atomic_orbitals["H128"]) == 5
        assert (
            g16_output.natural_atomic_orbitals["Ni1"]["NAO_Ni10"]["nao_type"]
            == "3py"
        )
        assert (
            g16_output.natural_atomic_orbitals["Ni1"]["NAO_Ni10"][
                "electron_type"
            ]
            == "Cor"
        )
        assert (
            g16_output.natural_atomic_orbitals["Ni1"]["NAO_Ni10"]["occupancy"]
            == 1.99858
        )
        assert (
            g16_output.natural_atomic_orbitals["Ni1"]["NAO_Ni10"]["energy"]
            == -2.68937
        )
        assert g16_output.get_num_naos("Ni1") == 31
        assert np.isclose(
            g16_output.get_total_electron_occ("Ni1"), 27.47171, rtol=1e-4
        )
        assert np.isclose(
            g16_output.get_total_electron_occ("H17"), 0.78631, rtol=1e-4
        )
        # import pprint
        # pprint.pprint(g16_output.natural_atomic_orbitals['Ni1'])
        assert len(g16_output.natural_population_analysis) == 128
        assert (
            g16_output.natural_population_analysis["Ni1"]["natural_charge"]
            == 0.52827
        )
        assert (
            g16_output.natural_population_analysis["C100"]["natural_charge"]
            == -0.42062
        )
        assert g16_output.natural_charges["Ni1"] == 0.52827
        assert g16_output.natural_charges["C100"] == -0.42062
        assert g16_output.total_electrons["Ni1"] == 27.47173
        assert g16_output.total_electrons["C100"] == 6.42062
        assert (
            g16_output.electronic_configuration["Ni1"]
            == "[core]4S(0.27)3d(8.70)4p(0.51)"
        )
        assert (
            g16_output.electronic_configuration["C100"]
            == "[core]2S(0.95)2p(3.44)3S(0.01)3p(0.02)"
        )
        assert g16_output.electronic_configuration["H128"] == "1S(0.80)"
        assert (
            g16_output.get_electronic_configuration("Ni1")
            == "[core]4S(0.27)3d(8.70)4p(0.51)"
        )


class TestGaussianCubeFile:
    def test_read_file_content(self, spin_cube_file):
        spin_cube = GaussianCubeFile(filename=spin_cube_file)
        assert spin_cube.cube_job_title == "Gaussian job density"
        assert (
            spin_cube.cube_job_description
            == "Electron density from Total SCF Density"
        )
        assert spin_cube.num_atoms == 2
        assert spin_cube.coordinate_origin == (-5.483229, -5.483229, -6.522947)
        assert type(spin_cube.coordinate_origin) is tuple
        assert spin_cube.grid_points == (9, 9, 11)
        assert type(spin_cube.grid_points) is tuple
        assert spin_cube.grid_increment_vector == (
            (1.2911, 0.0, 0.0),
            (0.0, 1.2911, 0.0),
            (0.0, 0.0, 1.2911),
        )
