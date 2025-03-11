import os.path

import numpy as np

from chemsmart.analysis.thermochemistry import (
    GaussianThermochemistry,
    Thermochemistry,
    qRRHOThermochemistry,
)
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
        # h = 6.62606957 * 10^-34 J s
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
                / (6.62606957 * 1e-34 * 6.62606957 * 1e-34)
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
                / (6.62606957 * 1e-34 * 6.62606957 * 1e-34)
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
        """Values from Gaussian output
        Temperature   298.150 Kelvin.  Pressure   1.00000 Atm.
        Zero-point vibrational energy      30919.1 (Joules/Mol)
                                           7.38984 (Kcal/Mol)
        Vibrational temperatures:    940.59   940.59  1998.50  3557.76
               (Kelvin)

        Zero-point correction=                           0.011776 (Hartree/Particle)
        Thermal correction to Energy=                    0.014410
        Thermal correction to Enthalpy=                  0.015354
        Thermal correction to Gibbs Free Energy=        -0.008927
        Sum of electronic and zero-point Energies=           -188.432903
        Sum of electronic and thermal Energies=              -188.430269
        Sum of electronic and thermal Enthalpies=            -188.429325
        Sum of electronic and thermal Free Energies=         -188.453606

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
        assert np.isclose(g16_output.moments_of_inertia[-1], 43.27307045)
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
        assert np.isclose(g16_output.mass, 43.98983)
        assert g16_output.multiplicity == 1
        assert np.isclose(g16_output.energies[-1], -188.444679593)
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
        # q_t = (2 * pi * m * k_B * T / h^2)^(3/2) * (k_B * T / P)
        # here T = 298.15 K, P = 1 atm = 101325 Pa
        # k_B = 1.3806488 * 10^-23 J/K
        # h = 6.62606957 * 10^-34 J s
        # N_A = 6.02214129 * 10^23 mol^-1
        # using these constants, we got 11467858.194120048
        expected_translational_partition_function = (
            2
            * np.pi
            * (g16_output.mass / (6.02214129 * 1e23 * 1000))
            * 1.3806488
            * 1e-23
            * 298.15
            / (6.62606957 * 1e-34) ** 2
        ) ** (3 / 2) * (1.3806488 * 1e-23 * 298.15 / 101325)
        assert np.isclose(
            thermochem1.translational_partition_function,
            expected_translational_partition_function,
        )
        assert np.isclose(
            thermochem1.translational_partition_function,
            0.114679e08,
        )

        # S_t = R * [ln(q_t) + 1 + 3 / 2]
        # R = 8.314462145468951 J mol^-1 K^-1
        # using these constants, we got 155.9382259343905 J mol^-1 K^-1
        expected_translational_entropy = 8.314462145468951 * (
            np.log(expected_translational_partition_function) + 1 + 3 / 2
        )
        assert np.isclose(
            thermochem1.translational_entropy,
            expected_translational_entropy,
        )
        assert np.isclose(
            thermochem1.translational_entropy,
            37.270 * cal_to_joules,
            atol=1e-2,
        )

        # E_t = 3/2 * R * T
        # we got 3718.4353330073513 J mol^-1
        expected_translational_internal_energy = (
            3 / 2 * 8.314462145468951 * 298.15
        )
        assert np.isclose(
            thermochem1.translational_internal_energy,
            expected_translational_internal_energy,
        )
        assert np.isclose(
            thermochem1.translational_internal_energy,
            0.889 * cal_to_joules * 1000,
            atol=1e1,
        )

        # C_t = 3/2 * R
        # we got 12.471693218203427 J mol^-1 K^-1
        expected_translational_heat_capacity = 3 / 2 * 8.314462145468951
        assert np.isclose(
            thermochem1.translational_heat_capacity,
            expected_translational_heat_capacity,
        )
        assert np.isclose(
            thermochem1.translational_heat_capacity,
            2.981 * cal_to_joules,
            atol=1e-2,
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
        expected_electronic_entropy = 8.314462145468951 * np.log(
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
        # 1 Angstrom = 1 * 10^-10 m
        # using these constants, we got 265.9699419350582
        expected_rotational_partition_function = (
            1
            / g16_output.rotational_symmetry_number
            * (
                298.15
                / (
                    (6.62606957 * 1e-34) ** 2
                    / (
                        8
                        * np.pi**2
                        * (
                            g16_output.moments_of_inertia[-1]
                            / (6.02214129 * 1e23 * 1000)
                            * 1e-10**2
                        )
                        * 1.3806488
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
            thermochem1.rotational_partition_function,
            0.265970e03,
        )

        # S_r = R * (ln(q_r) + 1)
        # we got 54.73729125485975 J mol^-1 K^-1
        expected_rotational_entropy = 8.314462145468951 * (
            np.log(expected_rotational_partition_function) + 1
        )
        assert np.isclose(
            thermochem1.rotational_entropy,
            expected_rotational_entropy,
        )
        assert np.isclose(
            thermochem1.rotational_entropy, 13.083 * cal_to_joules, atol=1e-2
        )

        # E_r = R * T
        # we got 2478.956888671568 J mol^-1
        expected_rotational_internal_energy = 8.314462145468951 * 298.15
        assert np.isclose(
            thermochem1.rotational_internal_energy,
            expected_rotational_internal_energy,
        )
        assert np.isclose(
            thermochem1.rotational_internal_energy,
            0.592 * cal_to_joules * 1000,
            atol=1e1,
        )

        # C_r = R
        expected_rotational_heat_capacity = 8.314462145468951
        assert np.isclose(
            thermochem1.rotational_heat_capacity,
            expected_rotational_heat_capacity,
        )
        assert np.isclose(
            thermochem1.rotational_heat_capacity,
            1.987 * cal_to_joules,
            atol=1e-2,
        )

        # q_v = q_1 * q_2 * ... * q_vDOF
        # For the zero reference point is the bottom of the well (BOT)
        # q_v,K = exp(-Θ_v,K / (2 * T)) / (1 - exp(-Θ_v,K / T))
        # Θ_v,K = h * v_K / k_B
        # c = 2.99792458 * 10^10 cm s^-1
        vibrational_frequencies = np.array(g16_output.vibrational_frequencies)
        expected_theta = (
            6.62606957
            * 1e-34
            * vibrational_frequencies
            * 2.99792458
            * 1e10
            / (1.3806488 * 1e-23)
        )
        # we got [0.21571805, 0.21571805, 0.03507491, 0.00256349]
        expected_vibrational_partition_function_by_mode_BOT = np.exp(
            -expected_theta / (2 * 298.15)
        ) / (1 - np.exp(-expected_theta / 298.15))
        assert np.allclose(
            thermochem1.vibrational_partition_function_by_mode_BOT,
            expected_vibrational_partition_function_by_mode_BOT,
        )
        # we got 4.184098315130738e-06
        expected_vibrational_partition_function_BOT = np.prod(
            expected_vibrational_partition_function_by_mode_BOT
        )
        assert np.isclose(
            thermochem1.vibrational_partition_function_BOT,
            expected_vibrational_partition_function_BOT,
        )
        assert np.isclose(
            thermochem1.vibrational_partition_function_BOT,
            0.418410e-05,
        )

        # For the zero reference point is the first vibrational energy level (V=0)
        # q_v,K = 1 / (1 - exp(-Θ_v,K / T))
        # we got [1.04454961, 1.04454961, 1.00122874, 1.00000657]
        expected_vibrational_partition_function_by_mode_V0 = 1 / (
            1 - np.exp(-expected_theta / 298.15)
        )
        assert np.allclose(
            thermochem1.vibrational_partition_function_by_mode_V0,
            expected_vibrational_partition_function_by_mode_V0,
        )
        # we got 1.0924317232036713
        expected_vibrational_partition_function_V0 = np.prod(
            expected_vibrational_partition_function_by_mode_V0
        )
        assert np.isclose(
            thermochem1.vibrational_partition_function_V0,
            expected_vibrational_partition_function_V0,
        )
        assert np.isclose(
            thermochem1.vibrational_partition_function_V0,
            0.109243e01,
        )

        # S_v = R * Σ((Θ_v,K / T) / (exp(Θ_v,K / T) - 1) - ln(1 - exp(-Θ_v,K / T)))
        # we got 3.1412492303422708 J mol^-1 K^-1
        expected_vibrational_entropy = 8.314462145468951 * np.sum(
            (expected_theta / 298.15) / (np.exp(expected_theta / 298.15) - 1)
            - np.log(1 - np.exp(-expected_theta / 298.15))
        )
        assert np.isclose(
            thermochem1.vibrational_entropy,
            expected_vibrational_entropy,
        )
        assert np.isclose(
            thermochem1.vibrational_entropy, 0.751 * cal_to_joules, atol=1e-2
        )

        # E_v = R * Σ(Θ_v,K * (1/2 + 1 / (exp(Θ_v,K / T) - 1)))
        # we got 31636.50928586775 J mol^-1
        expected_vibrational_internal_energy = 8.314462145468951 * np.sum(
            expected_theta
            * (1 / 2 + 1 / (np.exp(expected_theta / 298.15) - 1))
        )
        assert np.isclose(
            thermochem1.vibrational_internal_energy,
            expected_vibrational_internal_energy,
        )
        assert np.isclose(
            thermochem1.vibrational_internal_energy,
            7.561 * cal_to_joules * 1000,
            atol=1e1,
        )

        # C_v = R * Σ(exp(-Θ_v,K / T) * ((Θ_v,K / T) / (exp(-Θ_v,K / T) - 1))^2)
        # we got 8.16865700927472 J mol^-1 K^-1
        expected_vibrational_heat_capacity = 8.314462145468951 * np.sum(
            np.exp(-expected_theta / 298.15)
            * (
                (expected_theta / 298.15)
                / (np.exp(-expected_theta / 298.15) - 1)
            )
            ** 2
        )
        assert np.isclose(
            thermochem1.vibrational_heat_capacity,
            expected_vibrational_heat_capacity,
        )
        assert np.isclose(
            thermochem1.vibrational_heat_capacity,
            1.952 * cal_to_joules,
            atol=1e-2,
        )

        # q_tot = q_t * q_r * q_v * q_e
        # we got 3332032092.51935
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
        assert np.isclose(thermochem1.total_partition_function, 0.333203e10)

        # S_tot = S_t + S_r + S_v + S_e
        # we got 213.8167664195925 J mol^-1 K^-1
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
            thermochem1.total_entropy, 51.103 * cal_to_joules, atol=1e-2
        )

        # E_tot = E_t + E_r + E_v + E_e
        # we got 37833.90150754667 J mol^-1
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
            atol=1e1,
        )

        # C_tot = C_t + C_r + C_v + C_e
        # we got 28.954812372947096 J mol^-1 K^-1
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
            thermochem1.total_heat_capacity, 6.920 * cal_to_joules, atol=1e-2
        )

        gaussian_thermochem1 = GaussianThermochemistry(
            filename=gaussian_co2_opt_outfile,
            temperature=298.15,
            pressure=1,
        )
        assert np.isclose(
            gaussian_thermochem1.zero_point_vibrational_energy,
            30919.1,
            atol=1e-1,
        )
        assert np.allclose(
            gaussian_thermochem1.vibrational_temperatures,
            [940.59, 940.59, 1998.50, 3557.76],
            atol=1e-2,
        )
        assert np.isclose(
            gaussian_thermochem1.zero_point_correction, 0.011776, atol=1e-6
        )
        assert np.isclose(
            gaussian_thermochem1.thermal_correction_energy, 0.014410, atol=1e-6
        )
        assert np.isclose(
            gaussian_thermochem1.thermal_correction_enthalpy,
            0.015354,
            atol=1e-6,
        )
        assert np.isclose(
            gaussian_thermochem1.thermal_correction_free_energy,
            -0.008927,
            atol=1e-6,
        )
        assert np.isclose(
            gaussian_thermochem1.sum_of_electronic_and_zero_point_energies,
            -188.432903,
            atol=1e-6,
        )
        assert np.isclose(
            gaussian_thermochem1.sum_of_electronic_and_thermal_energies,
            -188.430269,
            atol=1e-6,
        )
        assert np.isclose(
            gaussian_thermochem1.sum_of_electronic_and_thermal_enthalpies,
            -188.429325,
            atol=1e-6,
        )
        assert np.isclose(
            gaussian_thermochem1.sum_of_electronic_and_thermal_free_energies,
            -188.453606,
            atol=1e-6,
        )

    def test_thermochemistry_co2_qrrho(self, gaussian_co2_opt_outfile):
        """Values from Goodvibes, as a reference:
                goodvibes -f 100 -c 1.0 -t 298.15 co2.log
        Structure                                           E        ZPE             H        T.S     T.qh-S          G(T)       qh-G(T)
           ********************************************************************************************************************************
        o  co2                                       -188.444680   0.011776   -188.429325   0.021262   0.021262   -188.450587   -188.450588
           ********************************************************************************************************************************
        """
        assert os.path.exists(gaussian_co2_opt_outfile)
        g16_output = Gaussian16Output(filename=gaussian_co2_opt_outfile)
        vibrational_frequencies = np.array(g16_output.vibrational_frequencies)
        assert g16_output.normal_termination
        qrrho_thermochem1 = qRRHOThermochemistry(
            filename=gaussian_co2_opt_outfile,
            temperature=298.15,  # in Kelvin
            concentration=1.0,  # in mol/L
        )

        # when arguments are not specified, the quasi-rrho calculation use
        # default cutoff frequency of 100 cm^-1 for entropy
        # and default alpha of 4
        # we got [0.9994528, 0.9994528, 0.99997314, 0.99999733]
        expected_entropy_damping_function = 1 / (
            1 + (100 / vibrational_frequencies) ** 4
        )
        assert np.allclose(
            qrrho_thermochem1.entropy_dumping_function,
            expected_entropy_damping_function,
        )

        # S_R,K = R * (1/2 + ln((8 * pi^3 * u'_K * k_B * T / h^2)^(1/2)))
        # u'_K = u_K * B_av / (u_K + B_av)
        # u_K = h / (8 * pi^2 * v_K)
        # B_av = 1 * 10^44 kg m^2
        expected_mu = (
            6.62606957
            * 1e-34
            / (8 * np.pi**2 * vibrational_frequencies * 2.99792458 * 1e10)
        )
        expected_mu_prime = expected_mu * 1.00e-44 / (expected_mu + 1.00e-44)
        # we got [4.13969469, 4.13969469, 1.00669595, -1.39088804] in J mol^-1 K^-1
        expected_freerot_entropy = 8.314462145468951 * (
            1 / 2
            + np.log(
                (
                    8
                    * np.pi**3
                    * expected_mu_prime
                    * 1.3806488
                    * 1e-23
                    * 298.15
                    / (6.62606957 * 1e-34) ** 2
                )
                ** (1 / 2)
            )
        )
        assert np.allclose(
            qrrho_thermochem1.freerot_entropy,
            expected_freerot_entropy,
        )

        # S^rrho_v,K = R * [(Θ_v,K / T) / (exp(Θ_v,K / T) - 1) - ln(1 - exp(-Θ_v,K / T))]
        # Θ_v,K = h * v_K / k_B
        expected_theta = (
            6.62606957
            * 1e-34
            * vibrational_frequencies
            * 2.99792458
            * 1e10
            / (1.3806488 * 1e-23)
        )
        # we got [1.53092635e+00, 1.53092635e+00, 7.86899024e-02, 7.06622985e-04] in J mol^-1 K^-1
        expected_rrho_entropy = 8.314462145468951 * (
            (expected_theta / 298.15) / (np.exp(expected_theta / 298.15) - 1)
            - np.log(1 - np.exp(-expected_theta / 298.15))
        )
        assert np.allclose(
            qrrho_thermochem1.rrho_entropy,
            expected_rrho_entropy,
        )

        # S^qrrho_v = Σ(w(v_K) * S^rrho_v,K + (1 - w(v_K)) * S_R,K)
        # we got 3.144125458725274 J mol^-1 K^-1
        expected_qrrho_vibrational_entropy = np.sum(
            expected_entropy_damping_function * expected_rrho_entropy
            + (1 - expected_entropy_damping_function)
            * expected_freerot_entropy
        )
        assert np.isclose(
            qrrho_thermochem1.qrrho_vibrational_entropy,
            expected_qrrho_vibrational_entropy,
        )

        # q_t,c = (2 * pi * m * k_B * T / h^2)^(3/2) * (1 / c)
        # we got 468737.7730646621
        expected_translational_partition_function_concentration = (
            2
            * np.pi
            * (g16_output.mass / (6.02214129 * 1e23 * 1000))
            * 1.3806488
            * 1e-23
            * 298.15
            / (6.62606957 * 1e-34) ** 2
        ) ** (3 / 2) * (1 / (1.0 * 6.02214129 * 1e23 * 1000))
        assert np.isclose(
            qrrho_thermochem1.translational_partition_function_concentration,
            expected_translational_partition_function_concentration,
        )

        # S_t,c = R * [ln(q_t,c) + 1 + 3/2]
        # we got 129.3547289549365 J mol^-1 K^-1
        expected_translational_entropy_concentration = 8.314462145468951 * (
            np.log(expected_translational_partition_function_concentration)
            + 1
            + 3 / 2
        )
        # S^qrrho_tot = S_t,c + S_r + S^qrrho_v + S_e
        # we got 129.3547289549365 + 54.73729125485975 + 3.144125458725274 + 0 = 187.23614566852152 J mol^-1 K^-1
        expected_rotational_entropy = 8.314462145468951 * (
            np.log(
                1
                / g16_output.rotational_symmetry_number
                * (
                    298.15
                    / (
                        (6.62606957 * 1e-34) ** 2
                        / (
                            8
                            * np.pi**2
                            * (
                                g16_output.moments_of_inertia[-1]
                                / (6.02214129 * 1e23 * 1000)
                                * 1e-10**2
                            )
                            * 1.3806488
                            * 1e-23
                        )
                    )
                )
            )
            + 1
        )
        expected_electronic_entropy = 8.314462145468951 * np.log(
            g16_output.multiplicity
        )
        expected_qrrho_total_entropy = (
            expected_translational_entropy_concentration
            + expected_rotational_entropy
            + expected_qrrho_vibrational_entropy
            + expected_electronic_entropy
        )
        assert np.isclose(
            qrrho_thermochem1.qrrho_total_entropy,
            expected_qrrho_total_entropy,
        )

        # E0 in Hartree
        assert np.isclose(qrrho_thermochem1.energies, -188.444680, atol=1e-6)

        # ZPE in Hartree
        assert np.isclose(
            qrrho_thermochem1.zero_point_energy_hartree, 0.011776, atol=1e-6
        )

        # E_tot = E_t + E_r + E_v + E_e
        # we got 3718.4353330073513 + 2478.956888671568 + 31636.50928586775 + 0 = 37833.90150754667 J mol^-1
        expected_translational_internal_energy = (
            3 / 2 * 8.314462145468951 * 298.15
        )
        expected_rotational_internal_energy = 8.314462145468951 * 298.15
        expected_vibrational_internal_energy = 8.314462145468951 * np.sum(
            expected_theta
            * (1 / 2 + 1 / (np.exp(expected_theta / 298.15) - 1))
        )
        expected_electronic_internal_energy = 0
        expected_total_internal_energy = (
            expected_translational_internal_energy
            + expected_rotational_internal_energy
            + expected_vibrational_internal_energy
            + expected_electronic_internal_energy
        )
        # H = E0 + E_tot + k_B * T in Hartree
        # we got -188.4293252361468 Hartree
        expected_enthalpy = g16_output.energies[-1] + (
            expected_total_internal_energy + 8.314462145468951 * 298.15
        ) / (4.35974434e-18 * 6.02214129 * 1e23)
        assert np.isclose(
            qrrho_thermochem1.enthalpy,
            expected_enthalpy,
        )
        assert np.isclose(qrrho_thermochem1.enthalpy, -188.429325, atol=1e-6)

        # T * S_tot in Hartree
        assert np.isclose(
            qrrho_thermochem1.entropy_times_temperature, 0.021262, atol=1e-6
        )

        # T * S^qrrho_tot in Hartree
        # 1 Hartree = 4.35974434 × 10^-18 Joules
        # 1 mol = 6.02214129 * 10^23 Particle
        # we got 0.021262412674747003 Hartree
        expected_qrrho_entropy_times_temperture = (
            298.15 * expected_qrrho_total_entropy
        ) / (4.35974434e-18 * 6.02214129 * 1e23)
        assert np.isclose(
            qrrho_thermochem1.qrrho_entropy_times_temperture,
            expected_qrrho_entropy_times_temperture,
        )
        assert np.isclose(
            qrrho_thermochem1.qrrho_entropy_times_temperture,
            0.021262,
            atol=1e-6,
        )

        # G = H - T * S_tot
        assert np.isclose(
            qrrho_thermochem1.gibbs_free_energy, -188.450587, atol=1e-6
        )

        # G^qrrho_qs = H - T * S^qrrho_tot
        # we got -188.45058764882157 Hartree
        expected_qrrho_gibbs_free_energy_qs = (
            expected_enthalpy - expected_qrrho_entropy_times_temperture
        )
        assert np.isclose(
            qrrho_thermochem1.qrrho_gibbs_free_energy,
            expected_qrrho_gibbs_free_energy_qs,
        )
        assert np.isclose(
            qrrho_thermochem1.qrrho_gibbs_free_energy_qs,
            -188.450588,
            atol=1e-6,
        )

        """Values from Goodvibes, as a reference:
                goodvibes -f 100 -c 1.0 -t 298.15 -q co2.log
           Structure                                           E        ZPE             H          qh-H        T.S     T.qh-S          G(T)       qh-G(T)
           **********************************************************************************************************************************************
        o  co2                                       -188.444680   0.011776   -188.429325   -188.429327   0.021262   0.021262   -188.450587   -188.450589
           **********************************************************************************************************************************************
        """
        # quasi-rrho correction is turned for both entropy and enthalpy,
        # with default cutoff frequencies of 100 cm^-1
        expected_enthalpy_damping_function = 1 / (
            1 + (100 / vibrational_frequencies) ** 4
        )
        assert np.allclose(
            qrrho_thermochem1.enthalpy_dumping_function,
            expected_enthalpy_damping_function,
        )

        # E^rrho_v,K = R * Θ_v,K * (1/2 + 1 / (exp(Θ_v,K / T) - 1))
        # we got [4258.62774713, 4258.62774713, 8328.63418484, 14790.61960677] in J mol^-1
        expected_rrho_internal_energy = (
            8.314462145468951
            * expected_theta
            * (1 / 2 + 1 / (np.exp(expected_theta / 298.15) - 1))
        )
        assert np.allclose(
            qrrho_thermochem1.rrho_internal_energy,
            expected_rrho_internal_energy,
        )

        # E^qrrho_v = Σ(w(v_K) * E^rrho_v,K + (1 - w(v_K)) * 1/2 * R * T)
        # we got 31632.978467909088 J mol^-1
        expected_qrrho_vibrational_internal_energy = np.sum(
            expected_enthalpy_damping_function * expected_rrho_internal_energy
            + (1 - expected_enthalpy_damping_function)
            * 1
            / 2
            * 8.314462145468951
            * 298.15
        )
        assert np.isclose(
            qrrho_thermochem1.qrrho_vibrational_internal_energy,
            expected_qrrho_vibrational_internal_energy,
        )

        # E^qrrho_tot = E_t + E_r + E^qrrho_v + E_e
        # we got 3718.4353330073513 + 2478.956888671568 + 31632.978467909088 + 0 = 37830.37068958801 J mol^-1
        expected_qrrho_total_internal_energy = (
            expected_translational_internal_energy
            + expected_rotational_internal_energy
            + expected_qrrho_vibrational_internal_energy
            + expected_electronic_internal_energy
        )
        assert np.isclose(
            qrrho_thermochem1.qrrho_total_internal_energy,
            expected_qrrho_total_internal_energy,
        )
        # H^qrrho = E0 + E^qrrho_tot + R * T
        # we got -188.42932658096436 Hartree
        expected_qrrho_enthalpy = g16_output.energies[-1] + (
            expected_qrrho_total_internal_energy + 8.314462145468951 * 298.15
        ) / (4.35974434e-18 * 6.02214129 * 1e23)
        assert np.isclose(
            qrrho_thermochem1.qrrho_enthalpy,
            expected_qrrho_enthalpy,
        )
        assert np.isclose(
            qrrho_thermochem1.qrrho_enthalpy, -188.429327, atol=1e-6
        )

        # G^qrrho = H^qrrho - T * S^qrrho_tot
        # we got -188.45058899363912 Hartree
        expected_qrrho_gibbs_free_energy = (
            expected_qrrho_enthalpy - expected_qrrho_entropy_times_temperture
        )
        assert np.isclose(
            qrrho_thermochem1.qrrho_gibbs_free_energy,
            expected_qrrho_gibbs_free_energy,
        )
        assert np.isclose(
            qrrho_thermochem1.qrrho_gibbs_free_energy, -188.450589, atol=1e-6
        )

        """Values from Goodvibes, as a reference:
                goodvibes -f 100 -c 0.5 -t 598.15 co2.log
        Structure                                           E        ZPE             H        T.S     T.qh-S          G(T)       qh-G(T)
           ********************************************************************************************************************************
        o  co2                                       -188.444680   0.011776   -188.424452   0.049327   0.049327   -188.473778   -188.473779
           ********************************************************************************************************************************
        """
        qrrho_thermochem2 = qRRHOThermochemistry(
            filename=gaussian_co2_opt_outfile,
            temperature=598.15,  # in Kelvin
            concentration=0.5,  # in mol/L
        )
        assert np.isclose(qrrho_thermochem2.energies, -188.444680, atol=1e-6)
        assert np.isclose(
            qrrho_thermochem2.zero_point_energy_hartree, 0.011776, atol=1e-6
        )
        assert np.isclose(qrrho_thermochem2.enthalpy, -188.424452, atol=1e-6)
        assert np.isclose(
            qrrho_thermochem2.entropy_times_temperature, 0.049327, atol=1e-6
        )
        assert np.isclose(
            qrrho_thermochem2.qrrho_entropy_times_temperture,
            0.049327,
            atol=1e-6,
        )
        assert np.isclose(
            qrrho_thermochem2.gibbs_free_energy, -188.473778, atol=1e-6
        )
        assert np.isclose(
            qrrho_thermochem2.qrrho_gibbs_free_energy_qs,
            -188.473779,
            atol=1e-6,
        )

        """Values from Goodvibes, as a reference:
                goodvibes -f 1000 -c 1.0 -t 298.15 -q co2.log
        Structure                                           E        ZPE             H          qh-H        T.S     T.qh-S          G(T)       qh-G(T)
           **********************************************************************************************************************************************
        o  co2                                       -188.444680   0.011776   -188.429325   -188.431976   0.021262   0.021781   -188.450587   -188.453757
           **********************************************************************************************************************************************
        """
        qrrho_thermochem3 = qRRHOThermochemistry(
            filename=gaussian_co2_opt_outfile,
            temperature=298.15,  # in Kelvin
            concentration=1.0,  # in mol/L
            s_freq_cutoff=1000,  # in cm^-1
            h_freq_cutoff=1000,  # in cm^-1
        )
        # the cutoff frequency for both entropy and enthalpy is specified as 1000 cm^-1
        # we got [0.15444091, 0.15444091, 0.78825007, 0.97395018]
        expected_damping_function = 1 / (
            1 + (1000 / vibrational_frequencies) ** 4
        )
        assert np.allclose(
            qrrho_thermochem3.entropy_dumping_function,
            expected_damping_function,
        )
        assert np.allclose(
            qrrho_thermochem3.enthalpy_dumping_function,
            expected_damping_function,
        )
        assert np.isclose(qrrho_thermochem3.energies, -188.444680, atol=1e-6)
        assert np.isclose(
            qrrho_thermochem3.zero_point_energy_hartree, 0.011776, atol=1e-6
        )
        assert np.isclose(qrrho_thermochem3.enthalpy, -188.429325, atol=1e-6)
        assert np.isclose(
            qrrho_thermochem3.qrrho_enthalpy, -188.431976, atol=1e-6
        )
        assert np.isclose(
            qrrho_thermochem3.entropy_times_temperature, 0.021262, atol=1e-6
        )
        assert np.isclose(
            qrrho_thermochem3.qrrho_entropy_times_temperture,
            0.021781,
            atol=1e-6,
        )
        assert np.isclose(
            qrrho_thermochem3.gibbs_free_energy, -188.450587, atol=1e-6
        )
        assert np.isclose(
            qrrho_thermochem3.qrrho_gibbs_free_energy, -188.453757, atol=1e-6
        )

    def test_thermochemistry_he_gaussian_output(self, gaussian_he_opt_outfile):
        """Values from Gaussian output
        Temperature   298.150 Kelvin.  Pressure   1.00000 Atm.
        Zero-point vibrational energy          0.0 (Joules/Mol)
                                           0.00000 (Kcal/Mol)
        Vibrational temperatures:
                 (Kelvin)

        Zero-point correction=                           0.000000 (Hartree/Particle)
        Thermal correction to Energy=                    0.001416
        Thermal correction to Enthalpy=                  0.002360
        Thermal correction to Gibbs Free Energy=        -0.011953
        Sum of electronic and zero-point Energies=             -2.915130
        Sum of electronic and thermal Energies=                -2.913713
        Sum of electronic and thermal Enthalpies=              -2.912769
        Sum of electronic and thermal Free Energies=           -2.927083

                            E (Thermal)             CV                S
                             KCal/Mol        Cal/Mol-Kelvin    Cal/Mol-Kelvin
        Total                    0.889              2.981             30.125
        Electronic               0.000              0.000              0.000
        Translational            0.889              2.981             30.125
        Rotational               0.000              0.000              0.000
        Vibrational              0.000              0.000              0.000
                              Q            Log10(Q)             Ln(Q)
        Total Bot       0.314751D+06          5.497968         12.659538
        Total V=0       0.314751D+06          5.497968         12.659538
        Vib (Bot)       0.100000D+01          0.000000          0.000000
        Vib (V=0)       0.100000D+01          0.000000          0.000000
        Electronic      0.100000D+01          0.000000          0.000000
        Translational   0.314751D+06          5.497968         12.659538
        Rotational      0.100000D+01          0.000000          0.000000
        """
        assert os.path.exists(gaussian_he_opt_outfile)
        g16_output = Gaussian16Output(filename=gaussian_he_opt_outfile)
        assert g16_output.normal_termination
        assert g16_output.num_atoms == 1
        mol = g16_output.molecule
        assert mol.empirical_formula == "He"
        assert np.isclose(g16_output.mass, 4.00260)
        assert g16_output.multiplicity == 1
        assert np.isclose(g16_output.energies[-1], -2.91512971456)
        assert mol.is_monoatomic

        thermochem2 = GaussianThermochemistry(
            filename=gaussian_he_opt_outfile,
            temperature=298.15,
            pressure=1,
        )
        assert np.isclose(
            thermochem2.translational_partition_function, 0.314751e06
        )
        assert np.isclose(
            thermochem2.translational_entropy,
            30.125 * cal_to_joules,
            atol=1e-2,
        )
        assert np.isclose(
            thermochem2.translational_internal_energy,
            0.889 * cal_to_joules * 1000,
            atol=1e1,
        )
        assert np.isclose(
            thermochem2.translational_heat_capacity,
            2.981 * cal_to_joules,
            atol=1e-2,
        )
        assert np.isclose(
            thermochem2.electronic_partition_function, 0.100000e01
        )
        assert np.isclose(
            thermochem2.electronic_entropy, 0.000 * cal_to_joules, atol=1e-2
        )
        assert np.isclose(
            thermochem2.electronic_internal_energy,
            0.000 * cal_to_joules * 1000,
            atol=1e1,
        )
        assert np.isclose(
            thermochem2.electronic_heat_capacity,
            0.000 * cal_to_joules,
            atol=1e-2,
        )
        assert np.isclose(
            thermochem2.rotational_partition_function, 0.100000e01
        )
        assert np.isclose(
            thermochem2.rotational_entropy, 0.000 * cal_to_joules, atol=1e-2
        )
        assert np.isclose(
            thermochem2.rotational_internal_energy,
            0.000 * cal_to_joules * 1000,
            atol=1e1,
        )
        assert np.isclose(
            thermochem2.rotational_heat_capacity,
            0.000 * cal_to_joules,
            atol=1e-2,
        )
        assert np.isclose(
            thermochem2.vibrational_partition_function_BOT, 0.100000e01
        )
        assert np.isclose(
            thermochem2.vibrational_partition_function_V0, 0.100000e01
        )
        assert np.isclose(
            thermochem2.vibrational_entropy, 0.000 * cal_to_joules, atol=1e-2
        )
        assert np.isclose(
            thermochem2.vibrational_internal_energy,
            0.000 * cal_to_joules * 1000,
            atol=1e1,
        )
        assert np.isclose(
            thermochem2.vibrational_heat_capacity,
            0.000 * cal_to_joules,
            atol=1e-2,
        )
        assert np.isclose(thermochem2.total_partition_function, 0.314751e06)
        assert np.isclose(
            thermochem2.total_entropy, 30.125 * cal_to_joules, atol=1e-2
        )
        assert np.isclose(
            thermochem2.total_internal_energy,
            0.889 * cal_to_joules * 1000,
            atol=1e1,
        )
        assert np.isclose(
            thermochem2.total_heat_capacity, 2.981 * cal_to_joules, atol=1e-2
        )

        gaussian_thermochem2 = GaussianThermochemistry(
            filename=gaussian_he_opt_outfile,
            temperature=298.15,
            pressure=1,
        )
        assert np.isclose(
            gaussian_thermochem2.zero_point_vibrational_energy,
            0.00000,
            atol=1e-1,
        )
        assert np.allclose(
            gaussian_thermochem2.vibrational_temperatures,
            [],
        )
        assert np.isclose(
            gaussian_thermochem2.zero_point_correction, 0.000000, atol=1e-6
        )
        assert np.isclose(
            gaussian_thermochem2.thermal_correction_energy, 0.001416, atol=1e-6
        )
        assert np.isclose(
            gaussian_thermochem2.thermal_correction_enthalpy,
            0.002360,
            atol=1e-6,
        )
        assert np.isclose(
            gaussian_thermochem2.thermal_correction_free_energy,
            -0.011953,
            atol=1e-6,
        )
        assert np.isclose(
            gaussian_thermochem2.sum_of_electronic_and_zero_point_energies,
            -2.915130,
            atol=1e-6,
        )
        assert np.isclose(
            gaussian_thermochem2.sum_of_electronic_and_thermal_energies,
            -2.913713,
            atol=1e-6,
        )
        assert np.isclose(
            gaussian_thermochem2.sum_of_electronic_and_thermal_enthalpies,
            -2.912769,
            atol=1e-6,
        )
        assert np.isclose(
            gaussian_thermochem2.sum_of_electronic_and_thermal_free_energies,
            -2.927083,
            atol=1e-6,
        )

    def test_thermochemistry_he_qrrho(self, gaussian_he_opt_outfile):
        """Values from Goodvibes, as a reference:
                goodvibes -f 1000 -c 0.5 -t 598.15 -q he.log
        Structure                                           E        ZPE             H          qh-H        T.S     T.qh-S          G(T)       qh-G(T)
           **********************************************************************************************************************************************
        o  He                                          -2.915130   0.000000     -2.910394     -2.910394   0.025951   0.025951     -2.936345     -2.936345
           **********************************************************************************************************************************************
        """
        assert os.path.exists(gaussian_he_opt_outfile)
        g16_output = Gaussian16Output(filename=gaussian_he_opt_outfile)
        assert g16_output.normal_termination
        qrrho_thermochem_he = qRRHOThermochemistry(
            filename=gaussian_he_opt_outfile,
            temperature=598.15,  # in Kelvin
            concentration=0.5,  # in mol/L
            s_freq_cutoff=1000,  # in cm^-1
            h_freq_cutoff=1000,  # in cm^-1
        )
        assert np.isclose(qrrho_thermochem_he.energies, -2.915130, atol=1e-6)
        assert np.isclose(
            qrrho_thermochem_he.zero_point_energy_hartree, 0.000000, atol=1e-6
        )
        assert np.isclose(qrrho_thermochem_he.enthalpy, -2.910394, atol=1e-6)
        assert np.isclose(
            qrrho_thermochem_he.qrrho_enthalpy, -2.910394, atol=1e-6
        )
        assert np.isclose(
            qrrho_thermochem_he.entropy_times_temperature, 0.025951, atol=1e-6
        )
        assert np.isclose(
            qrrho_thermochem_he.qrrho_entropy_times_temperture,
            0.025951,
            atol=1e-6,
        )
        assert np.isclose(
            qrrho_thermochem_he.gibbs_free_energy, -2.936345, atol=1e-6
        )
        assert np.isclose(
            qrrho_thermochem_he.qrrho_gibbs_free_energy, -2.936345, atol=1e-6
        )

    def test_thermochemistry_water_gaussian_output(
        self, gaussian_mp2_outputfile
    ):
        """Values from Gaussian output
        Temperature   298.150 Kelvin.  Pressure   1.00000 Atm.
        Zero-point vibrational energy      56211.1 (Joules/Mol)
                                           13.43478 (Kcal/Mol)
        Vibrational temperatures:   2342.81  5498.69  5679.79
                 (Kelvin)

        Zero-point correction=                           0.021410 (Hartree/Particle)
        Thermal correction to Energy=                    0.024245
        Thermal correction to Enthalpy=                  0.025189
        Thermal correction to Gibbs Free Energy=         0.003766
        Sum of electronic and zero-point Energies=            -76.307583
        Sum of electronic and thermal Energies=               -76.304747
        Sum of electronic and thermal Enthalpies=             -76.303803
        Sum of electronic and thermal Free Energies=          -76.325227

                            E (Thermal)             CV                S
                             KCal/Mol        Cal/Mol-Kelvin    Cal/Mol-Kelvin
        Total                   15.214              6.009             45.090
        Electronic               0.000              0.000              0.000
        Translational            0.889              2.981             34.608
        Rotational               0.889              2.981             10.475
        Vibrational             13.437              0.047              0.007
                              Q            Log10(Q)             Ln(Q)
        Total Bot       0.185336D-01         -1.732040         -3.988170
        Total V=0       0.130534D+09          8.115722         18.687141
        Vib (Bot)       0.142038D-09         -9.847595        -22.674925
        Vib (V=0)       0.100039D+01          0.000168          0.000387
        Electronic      0.100000D+01          0.000000          0.000000
        Translational   0.300431D+07          6.477745         14.915559
        Rotational      0.434320D+02          1.637809          3.771196
        """
        assert os.path.exists(gaussian_mp2_outputfile)
        g16_output = Gaussian16Output(filename=gaussian_mp2_outputfile)
        assert g16_output.normal_termination
        assert g16_output.num_atoms == 3
        mol = g16_output.molecule
        assert np.allclose(
            g16_output.moments_of_inertia, [0.62549131, 1.15863761, 1.78412891]
        )
        assert mol.empirical_formula == "H2O"
        assert np.isclose(g16_output.mass, 18.01056)
        assert g16_output.multiplicity == 1
        assert np.isclose(g16_output.energies[-1], -76.328992324258)
        assert g16_output.rotational_temperatures == [
            38.77653,
            20.93353,
            13.59452,
        ]
        assert g16_output.rotational_symmetry_number == 2
        assert g16_output.rotational_constants_in_Hz == [
            807.97175 * 1e9,
            436.18400 * 1e9,
            283.26385 * 1e9,
        ]
        assert g16_output.vibrational_frequencies == [
            1628.3334,
            3821.7812,
            3947.6507,
        ]
        assert not mol.is_monoatomic
        assert not mol.is_linear
        assert (
            g16_output.num_vib_frequencies == 3
        )  # vDOF = 3 * num_atoms - 6 since H2O is a non-linear molecule

        thermochem3 = GaussianThermochemistry(
            filename=gaussian_mp2_outputfile,
            temperature=298.15,
            pressure=1,
        )
        assert np.isclose(
            thermochem3.translational_partition_function, 0.300431e07
        )
        assert np.isclose(
            thermochem3.translational_entropy,
            34.608 * cal_to_joules,
            atol=1e-2,
        )
        assert np.isclose(
            thermochem3.translational_internal_energy,
            0.889 * cal_to_joules * 1000,
            atol=1e1,
        )
        assert np.isclose(
            thermochem3.translational_heat_capacity,
            2.981 * cal_to_joules,
            atol=1e-2,
        )
        assert np.isclose(
            thermochem3.electronic_partition_function, 0.100000e01
        )
        assert np.isclose(
            thermochem3.electronic_entropy, 0.000 * cal_to_joules, atol=1e-2
        )
        assert np.isclose(
            thermochem3.electronic_internal_energy,
            0.000 * cal_to_joules * 1000,
            atol=1e1,
        )
        assert np.isclose(
            thermochem3.electronic_heat_capacity,
            0.000 * cal_to_joules,
            atol=1e-2,
        )
        assert np.isclose(
            thermochem3.rotational_partition_function, 0.434320e02
        )
        assert np.isclose(
            thermochem3.rotational_entropy, 10.475 * cal_to_joules, atol=1e-2
        )
        assert np.isclose(
            thermochem3.rotational_internal_energy,
            0.889 * cal_to_joules * 1000,
            atol=1e1,
        )
        assert np.isclose(
            thermochem3.rotational_heat_capacity,
            2.981 * cal_to_joules,
            atol=1e-2,
        )
        assert np.isclose(
            thermochem3.vibrational_partition_function_BOT, 0.142038e-09
        )
        assert np.isclose(
            thermochem3.vibrational_partition_function_V0, 0.100039e01
        )
        assert np.isclose(
            thermochem3.vibrational_entropy, 0.007 * cal_to_joules, atol=1e-2
        )
        assert np.isclose(
            thermochem3.vibrational_internal_energy,
            13.437 * cal_to_joules * 1000,
            atol=1e1,
        )
        assert np.isclose(
            thermochem3.vibrational_heat_capacity,
            0.047 * cal_to_joules,
            atol=1e-2,
        )
        assert np.isclose(thermochem3.total_partition_function, 0.130534e09)
        assert np.isclose(
            thermochem3.total_entropy, 45.090 * cal_to_joules, atol=1e-2
        )
        assert np.isclose(
            thermochem3.total_internal_energy,
            15.214 * cal_to_joules * 1000,
            atol=1e1,
        )
        assert np.isclose(
            thermochem3.total_heat_capacity, 6.009 * cal_to_joules, atol=1e-2
        )

        gaussian_thermochem3 = GaussianThermochemistry(
            filename=gaussian_mp2_outputfile,
            temperature=298.15,
            pressure=1,
        )
        assert np.isclose(
            gaussian_thermochem3.zero_point_vibrational_energy,
            56211.1,
            atol=1e-1,
        )
        assert np.allclose(
            gaussian_thermochem3.vibrational_temperatures,
            [2342.81, 5498.69, 5679.79],
        )
        assert np.isclose(
            gaussian_thermochem3.zero_point_correction, 0.021410, atol=1e-6
        )
        assert np.isclose(
            gaussian_thermochem3.thermal_correction_energy, 0.024245, atol=1e-6
        )
        assert np.isclose(
            gaussian_thermochem3.thermal_correction_enthalpy,
            0.025189,
            atol=1e-6,
        )
        assert np.isclose(
            gaussian_thermochem3.thermal_correction_free_energy,
            0.003766,
            atol=1e-6,
        )
        assert np.isclose(
            gaussian_thermochem3.sum_of_electronic_and_zero_point_energies,
            -76.307583,
            atol=1e-6,
        )
        assert np.isclose(
            gaussian_thermochem3.sum_of_electronic_and_thermal_energies,
            -76.304747,
            atol=1e-6,
        )
        assert np.isclose(
            gaussian_thermochem3.sum_of_electronic_and_thermal_enthalpies,
            -76.303803,
            atol=1e-6,
        )
        assert np.isclose(
            gaussian_thermochem3.sum_of_electronic_and_thermal_free_energies,
            -76.325227,
            atol=1e-6,
        )

    def test_thermochemistry_water_qrrho(self, gaussian_mp2_outputfile):
        """Values from Goodvibes, as a reference:
                goodvibes -f 500 -c 2.0 -t 1298.15 -q water_mp2.log
        Structure                                           E        ZPE             H          qh-H        T.S     T.qh-S          G(T)       qh-G(T)
           **********************************************************************************************************************************************
        o  water_mp2                                  -76.328992   0.021410    -76.289193    -76.289224   0.098212   0.098221    -76.387404    -76.387445
           **********************************************************************************************************************************************
        """
        assert os.path.exists(gaussian_mp2_outputfile)
        g16_output = Gaussian16Output(filename=gaussian_mp2_outputfile)
        assert g16_output.normal_termination
        qrrho_thermochem_water = qRRHOThermochemistry(
            filename=gaussian_mp2_outputfile,
            temperature=1298.15,  # in Kelvin
            concentration=2.0,  # in mol/L
            s_freq_cutoff=500,  # in cm^-1
            h_freq_cutoff=500,  # in cm^-1
        )
        assert np.isclose(
            qrrho_thermochem_water.energies, -76.328992, atol=1e-6
        )
        assert np.isclose(
            qrrho_thermochem_water.zero_point_energy_hartree,
            0.021410,
            atol=1e-6,
        )
        assert np.isclose(
            qrrho_thermochem_water.enthalpy, -76.289193, atol=1e-6
        )
        assert np.isclose(
            qrrho_thermochem_water.qrrho_enthalpy, -76.289224, atol=1e-6
        )
        assert np.isclose(
            qrrho_thermochem_water.entropy_times_temperature,
            0.098212,
            atol=1e-6,
        )
        assert np.isclose(
            qrrho_thermochem_water.qrrho_entropy_times_temperture,
            0.098221,
            atol=1e-6,
        )
        assert np.isclose(
            qrrho_thermochem_water.gibbs_free_energy, -76.387404, atol=1e-6
        )
        assert np.isclose(
            qrrho_thermochem_water.qrrho_gibbs_free_energy,
            -76.387445,
            atol=1e-6,
        )