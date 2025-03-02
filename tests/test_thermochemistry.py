import math
import os.path

import numpy as np

from chemsmart.analysis.thermochemistry import Thermochemistry
from chemsmart.io.gaussian.output import (
    Gaussian16Output,
)
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian import GaussianOptJob
from chemsmart.settings.gaussian import GaussianProjectSettings
from chemsmart.utils.constants import cal_to_joules


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

        mol_as_ase_atoms = mol.to_ase()
        (
            mol_as_ase_atoms_moi_along_principal_axes,
            mol_as_ase_atoms_moi_principal_axes,
        ) = mol_as_ase_atoms.get_moments_of_inertia(vectors=True)

        assert np.allclose(mol.masses, mol_as_ase_atoms.get_masses())
        assert np.allclose(mol.positions, mol_as_ase_atoms.get_positions())
        assert np.allclose(
            mol.center_of_mass, mol_as_ase_atoms.get_center_of_mass()
        )

        # same moments of inertia along principal axes as from ase
        assert np.allclose(
            mol.moments_of_inertia,
            mol_as_ase_atoms_moi_along_principal_axes,
        )

        assert np.allclose(
            mol.moments_of_inertia,
            np.array([3111.56720644, 6850.18931668, 9483.18883579]),
            atol=1e-4,
        )

        # same moments of inertia principal axes as from ase
        assert np.allclose(
            mol.moments_of_inertia_principal_axes,
            mol_as_ase_atoms_moi_principal_axes,
        )

        # test moments of inertia from molecular structure and from Gaussian output
        # can't get moments of inertia from Gaussian output due to printing *****
        assert np.allclose(
            g16_output.moments_of_inertia, np.array([np.inf, np.inf, np.inf])
        )

        # test moments of inertia principal axes from molecular structure and from Gaussian output
        assert np.allclose(
            mol.moments_of_inertia_principal_axes,
            g16_output.moments_of_inertia_principal_axes,
            atol=1e-4,
        )

        assert np.allclose(
            g16_output.moments_of_inertia_principal_axes[0],
            [0.99974, -0.02291, 0.00147],
        )

        # test rotational temperatures from molecule (direct calc) same as from Gaussian output
        assert np.allclose(
            mol.rotational_temperatures,
            g16_output.rotational_temperatures,
            atol=1e-4,
        )  # [0.0078, 0.00354, 0.00256]

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

    def test_thermochemistry_co2(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        jobrunner_scratch,
    ):
        mol = Molecule.from_pubchem("280")
        tmp_path = tmpdir.join("co2.com")

        os.chdir(tmpdir)
        mol.write_com(tmp_path)

        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1

        # create gaussian job

        job = GaussianOptJob.from_filename(
            filename=tmp_path,
            settings=settings,
        )
        assert isinstance(job, GaussianOptJob)

        # set scratch
        jobrunner_scratch.scratch_dir = tmpdir

        # run the job
        job.run(jobrunner=jobrunner_scratch)
        assert job.is_complete()

        # check that the job finished
        tmp_log_path = tmpdir.join("CO2_fake.log")
        assert os.path.exists(tmp_log_path)
        assert np.isclose(mol.mass, 44.01, rtol=1e-2)

        # get thermochemistry for CO2 molecule
        thermochem1 = Thermochemistry(
            filename=tmp_log_path,
            temperature=298.15,  # in Kelvin
            pressure=1,  # in atm
        )

        assert np.isclose(
            thermochem1.translational_partition_function, 1.15e7, rtol=1e5
        )

    def test_thermochemistry_co2_gaussian_output(
        self, gaussian_co2_opt_outfile
    ):
        """Values from Goodvices, as a reference:
                goodvibes -f 100 -c 1.0 -t 298.15 co2.log
        Structure                                           E        ZPE             H        T.S     T.qh-S          G(T)       qh-G(T)
           ********************************************************************************************************************************
        o  co2                                       -188.444680   0.011776   -188.429325   0.021262   0.021262   -188.450587   -188.450588
           ********************************************************************************************************************************
        """
        assert os.path.exists(gaussian_co2_opt_outfile)
        g16_output = Gaussian16Output(filename=gaussian_co2_opt_outfile)
        assert g16_output.normal_termination
        assert g16_output.num_atoms == 3
        mol = g16_output.molecule
        mol_as_ase_atoms = mol.to_ase()
        (
            mol_as_ase_atoms_moi_along_principal_axes,
            mol_as_ase_atoms_moi_principal_axes,
        ) = mol_as_ase_atoms.get_moments_of_inertia(vectors=True)

        assert np.allclose(mol.masses, mol_as_ase_atoms.get_masses())
        assert np.allclose(mol.positions, mol_as_ase_atoms.get_positions())
        assert np.allclose(
            mol.center_of_mass, mol_as_ase_atoms.get_center_of_mass()
        )
        assert np.allclose(
            mol.moments_of_inertia,
            mol_as_ase_atoms_moi_along_principal_axes,
        )
        assert np.allclose(
            mol.moments_of_inertia_principal_axes,
            mol_as_ase_atoms_moi_principal_axes,
        )
        # test moments of inertia from molecular structure and from gaussian output
        assert np.allclose(
            mol.moments_of_inertia, g16_output.moments_of_inertia, rtol=1e-2
        )
        assert np.isclose(g16_output.moments_of_inertia[-1], 43.27307045, atol=1e-4)
        assert np.allclose(
            mol.moments_of_inertia_principal_axes[0],
            g16_output.moments_of_inertia_principal_axes[0],
        )

        # components X and Y are swapped but they are physically the same (same moment of inertia values)
        assert np.allclose(
            mol.moments_of_inertia_principal_axes[1],
            g16_output.moments_of_inertia_principal_axes[2],
        )
        assert np.allclose(
            mol.moments_of_inertia_principal_axes[2],
            g16_output.moments_of_inertia_principal_axes[1],
        )

        assert mol.empirical_formula == "CO2"
        assert np.isclose(mol.mass, 44.01, rtol=1e-2)
        assert g16_output.multiplicity == 1
        assert np.isclose(g16_output.energies[-1], -188.444680, rtol=1e-6)
        assert np.isclose(g16_output.zero_point_energy, 0.011776, rtol=1e-6)
        assert g16_output.rotational_temperatures == [0.56050]
        assert g16_output.rotational_symmetry_number == 2
        assert g16_output.rotational_constants_in_Hz == [11.678834 * 1e9]
        assert g16_output.vibrational_frequencies == [
            653.7393,
            653.7393,
            1389.0254,
            2472.7660,
        ]
        assert mol.is_linear
        assert (
            g16_output.num_vib_frequencies == 4
        )  # vDOF = 3 * num_atoms - 5 since CO2 is a linear molecule
        thermochem1 = Thermochemistry(
            filename=gaussian_co2_opt_outfile,
            temperature=298.15,  # in Kelvin
            pressure=1,  # in atm
        )

        """Values from Gaussian output
        Temperature   298.150 Kelvin.  Pressure   1.00000 Atm.
        
                            E (Thermal)             CV                S
                             KCal/Mol        Cal/Mol-Kelvin    Cal/Mol-Kelvin
        Total                    9.043              6.920             51.103
        Electronic               0.000              0.000              0.000
        Translational            0.889              2.981             37.270
        Rotational               0.592              1.987             13.083
        Vibrational              7.561              1.952              0.751
                              Q            Log10(Q)             Ln(Q)
        Total Bot       0.127619D+05          4.105917          9.454223
        Total V=0       0.333203D+10          9.522709         21.926848
        Vib (Bot)       0.418410D-05         -5.378398        -12.384219
        Vib (V=0)       0.109243D+01          0.038394          0.088406
        Electronic      0.100000D+01          0.000000          0.000000
        Translational   0.114679D+08          7.059482         16.255059
        Rotational      0.265970D+03          2.424833          5.583383
        """

        # q_t = (2 * pi * m * k_B * T / h^2)^(3/2) * (k_B * T / P)
        # here T = 298.15 K, P = 1 atm = 101325 Pa
        # k_B = 1.38064852 * 10^-23 J/K
        # h = 6.62607015 * 10^-34 J s
        # N_A = 6.022140857 * 10^23 mol^-1
        # using these constants, we got 11475738.782739257
        expected_translational_partition_function = (
            2
            * np.pi
            * (mol.mass / (6.022140857 * 1e23 * 1000))
            * 1.38064852
            * 1e-23
            * 298.15
            / (6.62607015 * 1e-34) ** 2
        ) ** (3 / 2) * (1.38064852 * 1e-23 * 298.15 / 101325)
        assert np.isclose(
            thermochem1.translational_partition_function,
            expected_translational_partition_function,
        )
        assert np.isclose(
            thermochem1.translational_partition_function, 0.114679e08, rtol=1e4
        )

        # S_t = R * [ln(q_t) + 1 + d / 2]
        # R = 8.314459861448581 J mol^-1 K^-1
        # d = 2 for linear molecules
        # using these constants, we got 151.78666481171763 J mol^-1 K^-1
        expected_translational_entropy = 8.314459861448581 * (
            np.log(expected_translational_partition_function) + 1 + 2 / 2
        )
        assert np.isclose(
            thermochem1.translational_entropy,
            expected_translational_entropy,
        )
        assert np.isclose(
            thermochem1.translational_entropy, 37.270 * cal_to_joules, rtol=1e2
        )

        # E_t = 3/2 * R * T
        # we got 3718.4343115363413 J mol^-1
        expected_translational_internal_energy = (
            3 / 2 * 8.314459861448581 * 298.15
        )
        assert np.isclose(
            thermochem1.translational_internal_energy,
            expected_translational_internal_energy,
        )
        assert np.isclose(
            thermochem1.translational_internal_energy,
            0.889 * cal_to_joules * 1000,
            rtol=1e1,
        )

        # C_t = 3/2 * R
        # we got 12.47168979217287 J mol^-1 K^-1
        expected_translational_heat_capacity = 3 / 2 * 8.314459861448581
        assert np.isclose(
            thermochem1.translational_heat_capacity,
            expected_translational_heat_capacity,
        )
        assert np.isclose(
            thermochem1.translational_heat_capacity,
            2.981 * cal_to_joules,
            rtol=1e-1,
        )

        # q_e = ω_0
        # ω_0 = multiplicity
        # we got 1
        expected_electronic_partition_function = g16_output.multiplicity
        assert np.isclose(
            thermochem1.electronic_partition_function,
            expected_electronic_partition_function,
        )
        assert np.isclose(
            thermochem1.electronic_partition_function, 0.100000e01
        )

        # S_e = R * ln(q_e)
        # we got 0 J mol^-1 K^-1
        expected_electronic_entropy = 8.314459861448581 * np.log(
            expected_electronic_partition_function
        )
        assert np.isclose(
            thermochem1.electronic_entropy,
            expected_electronic_entropy,
        )
        assert np.isclose(
            thermochem1.electronic_entropy, 0.000 * cal_to_joules
        )

        # E_e = 0 J mol^-1
        expected_electronic_internal_energy = 0
        assert np.isclose(
            thermochem1.electronic_internal_energy,
            expected_electronic_internal_energy,
        )
        assert np.isclose(
            thermochem1.electronic_internal_energy,
            0.000 * cal_to_joules * 1000,
        )

        # C_e = 0 J mol^-1 K^-1
        expected_electronic_heat_capacity = 0
        assert np.isclose(
            thermochem1.electronic_heat_capacity,
            expected_electronic_heat_capacity,
        )
        assert np.isclose(
            thermochem1.electronic_heat_capacity, 0.000 * cal_to_joules
        )

        # q_r = 1 / σ_r * (T / Θ_r)
        # Θ_r = h^2 / (8 * pi^2 * I * k_B)
        # 1 Bohr = 0.5291772105638411 Angstrom, 1 Angstrom = 1 * 10^-10 m
        # using these constants, we got 265.9698605567915
        expected_rotational_partition_function = (
            1
            / g16_output.rotational_symmetry_number
            * (
                298.15
                / (
                    (6.62607015 * 1e-34) ** 2
                    / (
                        8
                        * np.pi**2
                        * (
                            mol.moments_of_inertia
                            / (6.022140857 * 1e23 * 1000)
                            * (0.5291772105638411 * 1e-10) ** 2
                        )
                        * 1.38064852
                        * 1e-23
                    )
                )
            )
        )
        assert np.isclose(
            thermochem1.rotational_partition_function,
            expected_rotational_partition_function,
        )
        assert np.isclose(
            thermochem1.rotational_partition_function, 0.265970e03, rtol=1e-1
        )

        # S_r = R * (ln(q_r) + 1)
        # we got 54.7372736743199 J mol^-1 K^-1
        expected_rotational_entropy = 8.314459861448581 * (
            np.log(expected_rotational_partition_function) + 1
        )
        assert np.isclose(
            thermochem1.rotational_entropy,
            expected_rotational_entropy,
        )
        assert np.isclose(
            thermochem1.rotational_entropy, 13.083 * cal_to_joules, rtol=1e-1
        )

        # E_r = R * T
        # we got 2478.956207690894 J mol^-1
        expected_rotational_internal_energy = 8.314459861448581 * 298.15
        assert np.isclose(
            thermochem1.rotational_internal_energy,
            expected_rotational_internal_energy,
        )
        assert np.isclose(
            thermochem1.rotational_internal_energy,
            0.592 * cal_to_joules * 1000,
            rtol=1e1,
        )

        # C_r = R
        expected_rotational_heat_capacity = 8.314459861448581
        assert np.isclose(
            thermochem1.rotational_heat_capacity,
            expected_rotational_heat_capacity,
        )
        assert np.isclose(
            thermochem1.rotational_heat_capacity,
            1.987 * cal_to_joules,
            rtol=1e-2,
        )

        # q_v = q_1 * q_2 * ... * q_vDOF
        # For the zero reference point is the bottom of the well (BOT)
        # q_v,K = exp(-Θ_v,K / (2 * T)) / (1 - exp(-Θ_v,K / T))
        # Θ_v,K = h * v_K / k_B
        # c = 2.99792458 * 10^10 cm s^-1
        expected_theta_1 = (
            6.62607015
            * 1e-34
            * g16_output.vibrational_frequencies[0]
            * 2.99792458
            * 1e10
            / (1.38064852 * 1e-23)
        )
        expected_theta_2 = (
            6.62607015
            * 1e-34
            * g16_output.vibrational_frequencies[1]
            * 2.99792458
            * 1e10
            / (1.38064852 * 1e-23)
        )
        expected_theta_3 = (
            6.62607015
            * 1e-34
            * g16_output.vibrational_frequencies[2]
            * 2.99792458
            * 1e10
            / (1.38064852 * 1e-23)
        )
        expected_theta_4 = (
            6.62607015
            * 1e-34
            * g16_output.vibrational_frequencies[3]
            * 2.99792458
            * 1e10
            / (1.38064852 * 1e-23)
        )
        expected_q_1_BOT = math.exp(-expected_theta_1 / (2 * 298.15)) / (
            1 - math.exp(-expected_theta_1 / 298.15)
        )
        expected_q_2_BOT = math.exp(-expected_theta_2 / (2 * 298.15)) / (
            1 - math.exp(-expected_theta_2 / 298.15)
        )
        expected_q_3_BOT = math.exp(-expected_theta_3 / (2 * 298.15)) / (
            1 - math.exp(-expected_theta_3 / 298.15)
        )
        expected_q_4_BOT = math.exp(-expected_theta_4 / (2 * 298.15)) / (
            1 - math.exp(-expected_theta_4 / 298.15)
        )
        # we got [0.21571794219328444, 0.21571794219328444, 0.03507487817334666, 0.002563489617325113]
        expected_vibrational_partition_function_by_mode_BOT = [
            expected_q_1_BOT,
            expected_q_2_BOT,
            expected_q_3_BOT,
            expected_q_4_BOT,
        ]
        assert np.allclose(
            thermochem1.vibrational_partition_function_by_mode_BOT,
            expected_vibrational_partition_function_by_mode_BOT,
        )
        # we got 4.184082811907725e-06
        expected_vibrational_partition_function_BOT = (
            expected_q_1_BOT
            * expected_q_2_BOT
            * expected_q_3_BOT
            * expected_q_4_BOT
        )
        assert np.isclose(
            thermochem1.vibrational_partition_function_BOT,
            expected_vibrational_partition_function_BOT,
        )
        assert np.isclose(
            thermochem1.vibrational_partition_function_BOT,
            0.418410e-05,
            rtol=1e-8,
        )

        # For the zero reference point is the first vibrational energy level (V=0)
        # q_v,K = 1 / (1 - exp(-Θ_v,K / T))
        expected_q_1_V0 = 1 / (1 - math.exp(-expected_theta_1 / 298.15))
        expected_q_2_V0 = 1 / (1 - math.exp(-expected_theta_2 / 298.15))
        expected_q_3_V0 = 1 / (1 - math.exp(-expected_theta_3 / 298.15))
        expected_q_4_V0 = 1 / (1 - math.exp(-expected_theta_4 / 298.15))
        # we got [1.044549566691688, 1.044549566691688, 1.001228737283563, 1.0000065714358344]
        expected_vibrational_partition_function_by_mode_V0 = [
            expected_q_1_V0,
            expected_q_2_V0,
            expected_q_3_V0,
            expected_q_4_V0,
        ]
        assert np.allclose(
            thermochem1.vibrational_partition_function_by_mode_V0,
            expected_vibrational_partition_function_by_mode_V0,
        )
        # we got 1.0924316314141918
        expected_vibrational_partition_function_V0 = (
            expected_q_1_V0
            * expected_q_2_V0
            * expected_q_3_V0
            * expected_q_4_V0
        )
        assert np.isclose(
            thermochem1.vibrational_partition_function_V0,
            expected_vibrational_partition_function_V0,
        )
        assert np.isclose(
            thermochem1.vibrational_partition_function_V0,
            0.109243e01,
            rtol=1e-4,
        )

        # S_v = R * Σ((Θ_v,K / T) / (exp(Θ_v,K / T) - 1) - ln(1 - exp(-Θ_v,K / T)))
        # we got 3.1412459957707966 J mol^-1 K^-1
        expected_vibrational_entropy = 8.314459861448581 * (
            (
                (expected_theta_1 / 298.15)
                / (math.exp(expected_theta_1 / 298.15) - 1)
                - np.log(1 - math.exp(-expected_theta_1 / 298.15))
            )
            + (
                (expected_theta_2 / 298.15)
                / (math.exp(expected_theta_2 / 298.15) - 1)
                - np.log(1 - math.exp(-expected_theta_2 / 298.15))
            )
            + (
                (expected_theta_3 / 298.15)
                / (math.exp(expected_theta_3 / 298.15) - 1)
                - np.log(1 - math.exp(-expected_theta_3 / 298.15))
            )
            + (
                (expected_theta_4 / 298.15)
                / (math.exp(expected_theta_4 / 298.15) - 1)
                - np.log(1 - math.exp(-expected_theta_4 / 298.15))
            )
        )
        assert np.isclose(
            thermochem1.vibrational_entropy,
            expected_vibrational_entropy,
        )
        assert np.isclose(
            thermochem1.vibrational_entropy, 0.751 * cal_to_joules, rtol=1e-4
        )

        # E_v = R * Σ(Θ_v,K * (1/2 + 1 / (exp(Θ_v,K / T) - 1)))
        # we got 31636.50907329069 J mol^-1
        expected_vibrational_internal_energy = 8.314459861448581 * (
            (
                expected_theta_1
                * (1 / 2 + 1 / (math.exp(expected_theta_1 / 298.15) - 1))
            )
            + (
                expected_theta_2
                * (1 / 2 + 1 / (math.exp(expected_theta_2 / 298.15) - 1))
            )
            + (
                expected_theta_3
                * (1 / 2 + 1 / (math.exp(expected_theta_3 / 298.15) - 1))
            )
            + (
                expected_theta_4
                * (1 / 2 + 1 / (math.exp(expected_theta_4 / 298.15) - 1))
            )
        )
        assert np.isclose(
            thermochem1.vibrational_internal_energy,
            expected_vibrational_internal_energy,
        )
        assert np.isclose(
            thermochem1.vibrational_internal_energy,
            7.561 * cal_to_joules * 1000,
            rtol=1e1,
        )

        # C_v = R * Σ(exp(-Θ_v,K / T) * ((Θ_v,K / T) / (exp(-Θ_v,K / T) - 1))^2)
        # we got 8.168650902687183 J mol^-1 K^-1
        expected_vibrational_heat_capacity = 8.314459861448581 * (
            (
                math.exp(-expected_theta_1 / 298.15)
                * (
                    (expected_theta_1 / 298.15)
                    / (math.exp(-expected_theta_1 / 298.15) - 1)
                )
                ** 2
            )
            + (
                math.exp(-expected_theta_2 / 298.15)
                * (
                    (expected_theta_2 / 298.15)
                    / (math.exp(-expected_theta_2 / 298.15) - 1)
                )
                ** 2
            )
            + (
                math.exp(-expected_theta_3 / 298.15)
                * (
                    (expected_theta_3 / 298.15)
                    / (math.exp(-expected_theta_3 / 298.15) - 1)
                )
                ** 2
            )
            + (
                math.exp(-expected_theta_4 / 298.15)
                * (
                    (expected_theta_4 / 298.15)
                    / (math.exp(-expected_theta_4 / 298.15) - 1)
                )
                ** 2
            )
        )
        assert np.isclose(
            thermochem1.vibrational_heat_capacity,
            expected_vibrational_heat_capacity,
        )
        assert np.isclose(
            thermochem1.vibrational_heat_capacity,
            1.952 * cal_to_joules,
            rtol=1e-2,
        )

        # q_tot = q_t * q_r * q_v * q_e
        # we got 3334320528.7441
        expected_total_partition_function = (
            expected_translational_partition_function
            * expected_rotational_partition_function
            * expected_vibrational_partition_function_V0
            * expected_electronic_partition_function
        )
        assert np.isclose(
            thermochem1.total_partition_function,
            expected_total_partition_function,
        )
        assert np.isclose(
            thermochem1.total_partition_function, 0.333203e10, rtol=1e7
        )

        # S_tot = S_t + S_r + S_v + S_e
        # we got 209.66518448180832 J mol^-1 K^-1
        expected_total_entropy = (
            expected_translational_entropy
            + expected_rotational_entropy
            + expected_vibrational_entropy
            + expected_electronic_entropy
        )
        assert np.isclose(
            thermochem1.total_entropy,
            expected_total_entropy,
        )
        assert np.isclose(
            thermochem1.total_entropy, 51.103 * cal_to_joules, rtol=1e1
        )

        # E_tot = E_t + E_r + E_v + E_e
        # we got 37833.899592517926 J mol^-1
        expected_total_internal_energy = (
            expected_translational_internal_energy
            + expected_rotational_internal_energy
            + expected_vibrational_internal_energy
            + expected_electronic_internal_energy
        )
        assert np.isclose(
            thermochem1.total_internal_energy,
            expected_total_internal_energy,
        )
        assert np.isclose(
            thermochem1.total_internal_energy,
            9.043 * cal_to_joules * 1000,
            rtol=1e1,
        )

        # C_tot = C_t + C_r + C_v + C_e
        # we got 28.954800556308633 J mol^-1 K^-1
        expected_total_heat_capacity = (
            expected_translational_heat_capacity
            + expected_rotational_heat_capacity
            + expected_vibrational_heat_capacity
            + expected_electronic_heat_capacity
        )
        assert np.isclose(
            thermochem1.total_heat_capacity,
            expected_total_heat_capacity,
        )
        assert np.isclose(
            thermochem1.total_heat_capacity, 6.920 * cal_to_joules, rtol=1e-2
        )