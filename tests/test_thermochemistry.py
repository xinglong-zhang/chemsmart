import os.path

import numpy as np

from chemsmart.analysis.thermochemistry import Thermochemistry
from chemsmart.io.gaussian.output import (
    Gaussian16Output,
)
from chemsmart.io.molecules.structure import Molecule


class TestThermochemistry:
    def test_thermochemistry_from_gaussian_output(
        self, gaussian_singlet_opt_outfile
    ):
        """Values from Goodvices, as a reference:
                 Structure                                           E        ZPE             H        T.S     T.qh-S          G(T)       qh-G(T)
           ********************************************************************************************************************************
        o  nhc_neutral_singlet                      -1864.040180   0.284336  -1863.732135   0.079236   0.072784  -1863.811371  -1863.804919
        """
        assert os.path.exists(gaussian_singlet_opt_outfile)
        g16_output = Gaussian16Output(filename=gaussian_singlet_opt_outfile)
        assert g16_output.normal_termination
        mol = g16_output.molecule
        assert mol.empirical_formula == "C19H12F3I2N3O"
        assert np.isclose(mol.mass, 609.1230, rtol=1e-3)

        expected_E = (
            -1864.040180
        )  # in Hartree, from Gaussian output last SCF Done value:
        # SCF Done:  E(RM062X) =  -1864.04018044     A.U. after    1 cycles
        assert np.isclose(g16_output.energies[-1], expected_E, rtol=10e-6)

        expected_ZPE = 0.284336  # in Hartree, from Gaussian output:
        # Zero-point correction=                           0.284336 (Hartree/Particle)
        assert np.isclose(
            g16_output.zero_point_energy, expected_ZPE, rtol=10e-6
        )

        thermochem1 = Thermochemistry(
            filename=gaussian_singlet_opt_outfile,
            temperature=298.15,  # in Kelvin
            pressure=1,  # in atm
        )

        # q_t = (2 * pi * m * k_B * T / h^2)^(3/2) * (k_B * T / P)
        # here T = 298.15 K, P = 1 atm = 101325 Pa
        # k_B = 1.380649 * 10^-23 J/K
        # h = 6.62607015 * 10^-34 J s
        # A_N = 6.0221408e+23 mol^-1
        # using these constants, we got 8.732609925501497e+48
        expected_translational_partition_function = (
            (
                2
                * np.pi
                * (mol.mass / (6.0221408 * 1e23 * 1000))
                * 1.380649
                * 1e-23
                * 298.15
                / (6.62607015 * 1e-34 * 6.62607015 * 1e-34)
            )
        ) ** (3 / 2) * (1.380649 * 1e-23 * 298.15 / 101325)
        assert np.isclose(
            thermochem1.translational_partition_function,
            expected_translational_partition_function,
        )
        # print(expected_translational_partition_function)
        # print(thermochem1.translational_partition_function)
        # # assert np.isclose(thermochem1.translational_internal_energy()

        assert g16_output.freq

        thermochem2 = Thermochemistry(
            filename=gaussian_singlet_opt_outfile,
            temperature=598.15,  # in Kelvin
            pressure=1.2,  # in atm
        )

        expected_translational_partition_function2 = (
            (
                2
                * np.pi
                * (mol.mass / (6.0221408 * 1e23 * 1000))
                * 1.380649
                * 1e-23
                * 598.15
                / (6.62607015 * 1e-34 * 6.62607015 * 1e-34)
            )
        ) ** (3 / 2) * (1.380649 * 1e-23 * 298.15 / (1.2 * 101325))

        assert np.isclose(
            thermochem2.translational_partition_function,
            expected_translational_partition_function2,
            rtol=1e5,
        )

        # tests on thermochem2

    def test_thermochemistry_co2(self, tmpdir):
        mol = Molecule.from_pubchem("280")
        tmp_path = tmpdir.join("co2.com")
        mol.write_com(tmp_path)
        assert os.path.exists(tmp_path)
        assert np.isclose(mol.mass, 44.01, rtol=1e-2)
        thermochem1 = Thermochemistry(
            filename=tmp_path,
            temperature=298.15,  # in Kelvin
            pressure=1,  # in atm
        )

        assert np.isclose(
            thermochem1.translational_partition_function, 1.15e7, rtol=1e5
        )
