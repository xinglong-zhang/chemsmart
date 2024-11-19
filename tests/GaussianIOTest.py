import os.path
import numpy as np
from chemsmart.io.gaussian.output import Gaussian16TDDFTOutput
from chemsmart.io.gaussian.output import Gaussian16WBIOutput
from chemsmart.io.gaussian.cube import GaussianCubeFile


class TestGaussian16Output:
    def test_normal_termination_with_forces_and_frequencies(self, td_outputfile):
        assert os.path.exists(td_outputfile)
        g16_output = Gaussian16TDDFTOutput(filename=td_outputfile)
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


class TestGaussianWBIOutput:
    def test_normal_termination_with_forces_and_frequencies(self, wbi_outputfile):
        assert os.path.exists(wbi_outputfile)
        g16_output = Gaussian16WBIOutput(filename=wbi_outputfile)
        assert g16_output.nbo_version == "3.1"
        assert len(g16_output.natural_atomic_orbitals) == 128
        assert len(g16_output.natural_atomic_orbitals["Ni1"]) == 31
        assert len(g16_output.natural_atomic_orbitals["P2"]) == 18
        assert len(g16_output.natural_atomic_orbitals["H128"]) == 5
        assert (
            g16_output.natural_atomic_orbitals["Ni1"]["NAO_Ni10"]["nao_type"] == "3py"
        )
        assert (
            g16_output.natural_atomic_orbitals["Ni1"]["NAO_Ni10"]["electron_type"]
            == "Cor"
        )
        assert (
            g16_output.natural_atomic_orbitals["Ni1"]["NAO_Ni10"]["occupancy"]
            == 1.99858
        )
        assert (
            g16_output.natural_atomic_orbitals["Ni1"]["NAO_Ni10"]["energy"] == -2.68937
        )
        assert g16_output.get_num_naos("Ni1") == 31
        assert np.isclose(g16_output.get_total_electron_occ("Ni1"), 27.47171, rtol=1e-4)
        assert np.isclose(g16_output.get_total_electron_occ("H17"), 0.78631, rtol=1e-4)
        # import pprint
        # pprint.pprint(g16_output.natural_atomic_orbitals['Ni1'])
        assert len(g16_output.natural_population_analysis) == 128
        assert (
            g16_output.natural_population_analysis["Ni1"]["natural_charge"] == 0.52827
        )
        assert (
            g16_output.natural_population_analysis["C100"]["natural_charge"] == -0.42062
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
            spin_cube.cube_job_description == "Electron density from Total SCF Density"
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
