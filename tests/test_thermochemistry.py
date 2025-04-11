import os.path

import numpy as np
import pytest

from chemsmart.analysis.thermochemistry import (
    Thermochemistry,
)
from chemsmart.io.gaussian.output import (
    Gaussian16Output,
)
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.orca.output import ORCAOutput
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
        # CO2_fake.log file has no thermochemistry data
        # since it was run with fake gaussian
        with pytest.raises(TypeError):
            thermochem1 = Thermochemistry(
                filename=tmp_log_path,
                temperature=298.15,  # in Kelvin
                pressure=1,  # in atm
            )
            assert np.isclose(
                thermochem1.translational_partition_function, 1.15e7, rtol=1e5
            )


class TestThermochemistryCO2:
    """CO2 is used as an example to test the thermochemical properties of linear molecules."""

    def test_thermochemistry_co2_gaussian_output(
        self, gaussian_co2_opt_outfile
    ):
        """Values from Gaussian output
        Temperature   298.150 Kelvin.  Pressure   1.00000 Atm.
        Atom     1 has atomic number  8 and mass  15.99491
        Atom     2 has atomic number  8 and mass  15.99491
        Atom     3 has atomic number  6 and mass  12.00000
        Molecular mass:    43.98983 amu.
        Principal axes and moments of inertia in atomic units:
                                  1         2         3
            Eigenvalues --     0.00000 154.53094 154.53094
                  X           -0.00000   1.00000  -0.00000
                  Y            0.00000   0.00000   1.00000
                  Z            1.00000   0.00000  -0.00000
        This molecule is a prolate symmetric top.
        Rotational symmetry number  2.
        Rotational temperature (Kelvin)      0.56050
        Rotational constant (GHZ):          11.678834
        """
        assert os.path.exists(gaussian_co2_opt_outfile)
        g16_output = Gaussian16Output(filename=gaussian_co2_opt_outfile)
        assert g16_output.normal_termination
        assert g16_output.job_type == "opt"
        assert g16_output.num_atoms == 3
        mol = g16_output.molecule
        assert mol.empirical_formula == "CO2"
        assert g16_output.multiplicity == 1
        assert np.isclose(g16_output.energies[-1], -188.444679593)
        assert np.isclose(mol.mass, 44.009)  # weighted_atomic_mass=True
        assert np.isclose(g16_output.mass, 43.98983)
        assert np.isclose(
            g16_output.moments_of_inertia[-1], 43.27307045
        )  # in amu Å^2
        assert np.isclose(mol.moments_of_inertia[-1], 43.28411748057631)
        assert g16_output.rotational_symmetry_number == 2
        assert g16_output.rotational_temperatures == [0.56050]
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

        thermochem0 = Thermochemistry(
            filename=gaussian_co2_opt_outfile,
            temperature=298.15,  # in Kelvin
            pressure=1,  # in atm
            weighted_atomic_mass=True,  # use natural abundance weighted masses
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
            * (mol.mass / (6.02214129 * 1e23 * 1000))
            * 1.3806488
            * 1e-23
            * 298.15
            / (6.62606957 * 1e-34) ** 2
        ) ** (3 / 2) * (1.3806488 * 1e-23 * 298.15 / 101325)
        assert np.isclose(
            thermochem0.translational_partition_function,
            expected_translational_partition_function,
        )

        # S_t = R * [ln(q_t) + 1 + 3 / 2]
        # R = 8.314462145468951 J mol^-1 K^-1
        # using these constants, we got 155.9382259343905 J mol^-1 K^-1
        expected_translational_entropy = 8.314462145468951 * (
            np.log(expected_translational_partition_function) + 1 + 3 / 2
        )
        assert np.isclose(
            thermochem0.translational_entropy,
            expected_translational_entropy,
        )

        # E_t = 3/2 * R * T
        # we got 3718.4353330073513 J mol^-1
        expected_translational_internal_energy = (
            3 / 2 * 8.314462145468951 * 298.15
        )
        assert np.isclose(
            thermochem0.translational_internal_energy,
            expected_translational_internal_energy,
        )

        # C_t = 3/2 * R
        # we got 12.471693218203427 J mol^-1 K^-1
        expected_translational_heat_capacity = 3 / 2 * 8.314462145468951
        assert np.isclose(
            thermochem0.translational_heat_capacity,
            expected_translational_heat_capacity,
        )

        # q_e = ω_0
        # ω_0 = multiplicity
        # we got 1
        expected_electronic_partition_function = g16_output.multiplicity
        assert np.isclose(
            thermochem0.electronic_partition_function,
            expected_electronic_partition_function,
        )

        # S_e = R * ln(q_e)
        # we got 0 J mol^-1 K^-1
        expected_electronic_entropy = 8.314462145468951 * np.log(
            expected_electronic_partition_function
        )
        assert np.isclose(
            thermochem0.electronic_entropy,
            expected_electronic_entropy,
        )

        # E_e = 0 J mol^-1
        expected_electronic_internal_energy = 0
        assert np.isclose(
            thermochem0.electronic_internal_energy,
            expected_electronic_internal_energy,
        )

        # C_e = 0 J mol^-1 K^-1
        expected_electronic_heat_capacity = 0
        assert np.isclose(
            thermochem0.electronic_heat_capacity,
            expected_electronic_heat_capacity,
        )

        # q_r = 1 / σ_r * (T / Θ_r)
        # Θ_r = h^2 / (8 * pi^2 * I * k_B)
        # 1 Angstrom = 1 * 10^-10 m
        # using these constants, we got 265.9699419335578
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
                            mol.moments_of_inertia[-1]
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
            thermochem0.rotational_partition_function,
            expected_rotational_partition_function,
        )

        # S_r = R * (ln(q_r) + 1)
        # we got 54.737291254812845 J mol^-1 K^-1
        expected_rotational_entropy = 8.314462145468951 * (
            np.log(expected_rotational_partition_function) + 1
        )
        assert np.isclose(
            thermochem0.rotational_entropy,
            expected_rotational_entropy,
        )

        # E_r = R * T
        # we got 2478.956888671568 J mol^-1
        expected_rotational_internal_energy = 8.314462145468951 * 298.15
        assert np.isclose(
            thermochem0.rotational_internal_energy,
            expected_rotational_internal_energy,
        )

        # C_r = R
        expected_rotational_heat_capacity = 8.314462145468951
        assert np.isclose(
            thermochem0.rotational_heat_capacity,
            expected_rotational_heat_capacity,
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
        expected_vibrational_partition_function_by_mode_bot = np.exp(
            -expected_theta / (2 * 298.15)
        ) / (1 - np.exp(-expected_theta / 298.15))
        assert np.allclose(
            thermochem0.vibrational_partition_function_by_mode_bot,
            expected_vibrational_partition_function_by_mode_bot,
        )
        # we got 4.184098315130738e-06
        expected_vibrational_partition_function_bot = np.prod(
            expected_vibrational_partition_function_by_mode_bot
        )
        assert np.isclose(
            thermochem0.vibrational_partition_function_bot,
            expected_vibrational_partition_function_bot,
        )

        # For the zero reference point is the first vibrational energy level (V=0)
        # q_v,K = 1 / (1 - exp(-Θ_v,K / T))
        # we got [1.04454961, 1.04454961, 1.00122874, 1.00000657]
        expected_vibrational_partition_function_by_mode_v0 = 1 / (
            1 - np.exp(-expected_theta / 298.15)
        )
        assert np.allclose(
            thermochem0.vibrational_partition_function_by_mode_v0,
            expected_vibrational_partition_function_by_mode_v0,
        )
        # we got 1.0924317232036713
        expected_vibrational_partition_function_v0 = np.prod(
            expected_vibrational_partition_function_by_mode_v0
        )
        assert np.isclose(
            thermochem0.vibrational_partition_function_v0,
            expected_vibrational_partition_function_v0,
        )

        # S_v = R * Σ((Θ_v,K / T) / (exp(Θ_v,K / T) - 1) - ln(1 - exp(-Θ_v,K / T)))
        # we got 3.1412492303422708 J mol^-1 K^-1
        expected_vibrational_entropy = 8.314462145468951 * np.sum(
            (expected_theta / 298.15) / (np.exp(expected_theta / 298.15) - 1)
            - np.log(1 - np.exp(-expected_theta / 298.15))
        )
        assert np.isclose(
            thermochem0.vibrational_entropy,
            expected_vibrational_entropy,
        )

        # E_v = R * Σ(Θ_v,K * (1/2 + 1 / (exp(Θ_v,K / T) - 1)))
        # we got 31636.50928586775 J mol^-1
        expected_vibrational_internal_energy = 8.314462145468951 * np.sum(
            expected_theta
            * (1 / 2 + 1 / (np.exp(expected_theta / 298.15) - 1))
        )
        assert np.isclose(
            thermochem0.vibrational_internal_energy,
            expected_vibrational_internal_energy,
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
            thermochem0.vibrational_heat_capacity,
            expected_vibrational_heat_capacity,
        )

        # q_tot = q_t * q_r * q_v * q_e
        # we got 3332032092.51935
        expected_total_partition_function = (
            expected_translational_partition_function
            * expected_rotational_partition_function
            * expected_vibrational_partition_function_v0
            * expected_electronic_partition_function
        )
        assert np.isclose(
            thermochem0.total_partition_function,
            expected_total_partition_function,
        )

        # S_tot = S_t + S_r + S_v + S_e
        # we got 213.81676641954562 J mol^-1 K^-1
        expected_total_entropy = (
            expected_translational_entropy
            + expected_rotational_entropy
            + expected_vibrational_entropy
            + expected_electronic_entropy
        )
        assert np.isclose(
            thermochem0.total_entropy,
            expected_total_entropy,
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
            thermochem0.total_internal_energy,
            expected_total_internal_energy,
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
            thermochem0.total_heat_capacity,
            expected_total_heat_capacity,
        )

        thermochem1 = Thermochemistry(
            filename=gaussian_co2_opt_outfile,
            temperature=298.15,  # in Kelvin
            pressure=1,  # in atm
            weighted_atomic_mass=False,  # use single isotope masses
        )
        """Values from Gaussian output
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
        assert np.isclose(
            thermochem1.translational_partition_function,
            0.114679e08,
        )
        assert np.isclose(
            thermochem1.translational_entropy,
            37.270 * cal_to_joules,
            atol=1e-2,
        )
        assert np.isclose(
            thermochem1.translational_internal_energy,
            0.889 * cal_to_joules * 1000,
            atol=1e1,
        )
        assert np.isclose(
            thermochem1.translational_heat_capacity,
            2.981 * cal_to_joules,
            atol=1e-2,
        )
        assert np.isclose(
            thermochem1.electronic_partition_function, 0.100000e01
        )
        assert np.isclose(
            thermochem1.electronic_entropy, 0.000 * cal_to_joules
        )
        assert np.isclose(
            thermochem1.electronic_internal_energy,
            0.000 * cal_to_joules * 1000,
        )
        assert np.isclose(
            thermochem1.electronic_heat_capacity, 0.000 * cal_to_joules
        )
        assert np.isclose(
            thermochem1.rotational_partition_function, 0.265970e03, atol=1e0
        )
        assert np.isclose(
            thermochem1.rotational_entropy, 13.083 * cal_to_joules, atol=1e-2
        )
        assert np.isclose(
            thermochem1.rotational_internal_energy,
            0.592 * cal_to_joules * 1000,
            atol=1e1,
        )
        assert np.isclose(
            thermochem1.rotational_heat_capacity,
            1.987 * cal_to_joules,
            atol=1e-2,
        )
        assert np.isclose(
            thermochem1.vibrational_partition_function_bot,
            0.418410e-05,
        )
        assert np.isclose(
            thermochem1.vibrational_partition_function_v0,
            0.109243e01,
        )
        assert np.isclose(
            thermochem1.vibrational_entropy, 0.751 * cal_to_joules, atol=1e-2
        )
        assert np.isclose(
            thermochem1.vibrational_internal_energy,
            7.561 * cal_to_joules * 1000,
            atol=1e1,
        )
        assert np.isclose(
            thermochem1.vibrational_heat_capacity,
            1.952 * cal_to_joules,
            atol=1e-2,
        )
        assert np.isclose(
            thermochem1.total_partition_function, 0.333203e10, atol=1e6
        )
        assert np.isclose(
            thermochem1.total_entropy, 51.103 * cal_to_joules, atol=1e-2
        )
        assert np.isclose(
            thermochem1.total_internal_energy,
            9.043 * cal_to_joules * 1000,
            atol=1e1,
        )
        assert np.isclose(
            thermochem1.total_heat_capacity, 6.920 * cal_to_joules, atol=1e-2
        )

    def test_thermochemistry_co2_orca_output(self, orca_co2_output):
        """Values from ORCA output
        --------------------------
        THERMOCHEMISTRY AT 298.15K
        --------------------------

        Temperature         ...   298.15 K
        Pressure            ...     1.00 atm
        Total Mass          ...    44.01 AMU
        ...
        The molecule is recognized as being linear
        ...
        freq.     679.10  E(vib)   ...       0.08
        freq.     679.10  E(vib)   ...       0.08
        freq.    1440.27  E(vib)   ...       0.00
        freq.    2532.19  E(vib)   ...       0.00
        ...
        Point Group:  D(inf)h, Symmetry Number:   1
        Rotational constants in cm-1:     0.000000     0.394105     0.394105
        """
        assert os.path.exists(orca_co2_output)
        orca_out = ORCAOutput(filename=orca_co2_output)
        assert orca_out.normal_termination
        assert orca_out.job_type == "opt"
        assert orca_out.natoms == 3
        mol = orca_out.molecule
        assert mol.empirical_formula == "CO2"
        assert orca_out.multiplicity == 1
        assert np.isclose(orca_out.energies[-1], -188.370538039014)
        assert np.isclose(mol.mass, 44.009)  # weighted_atomic_mass=True
        assert np.isclose(orca_out.mass, 44.01)
        assert orca_out.rotational_symmetry_number == 1
        assert orca_out.rotational_constants_in_wavenumbers == [
            0.000000,
            0.394105,
            0.394105,
        ]
        assert orca_out.vibrational_frequencies == [
            679.10,
            679.10,
            1440.27,
            2532.19,
        ]
        assert mol.is_linear
        assert orca_out.num_vibration_modes == 4

    def test_thermochemistry_co2_qrrho(self, gaussian_co2_opt_outfile):
        """Values from Goodvibes, as a reference:
                goodvibes -f 100 -c 1.0 -t 298.15 --qs grimme --bav "conf" co2.log
        Structure                                           E        ZPE             H        T.S     T.qh-S          G(T)       qh-G(T)
           ********************************************************************************************************************************
        o  co2                                       -188.444680   0.011776   -188.429325   0.021262   0.021262   -188.450587   -188.450588
           ********************************************************************************************************************************
        """
        assert os.path.exists(gaussian_co2_opt_outfile)
        g16_output = Gaussian16Output(filename=gaussian_co2_opt_outfile)
        vibrational_frequencies = np.array(g16_output.vibrational_frequencies)
        assert g16_output.normal_termination
        mol = g16_output.molecule
        qrrho_thermochem_co2_1 = Thermochemistry(
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
            qrrho_thermochem_co2_1.entropy_damping_function,
            expected_entropy_damping_function,
        )

        # S_R,K = R * (1/2 + ln((8 * pi^3 * u'_K * k_B * T / h^2)^(1/2)))
        # u'_K = u_K * B_av / (u_K + B_av)
        # u_K = h / (8 * pi^2 * v_K)
        # B_av = h / (8 * pi^2 * I)
        expected_i = (
            mol.moments_of_inertia[-1] / (6.02214129 * 1e23 * 1000) * 1e-10**2
        )
        # we got B_av = 5.675019509661021 * 10^-44 kg m^2
        expected_bav = (
            6.62606957
            * 1e-34
            / (6.62606957 * 1e-34 / (8 * np.pi**2 * expected_i))
        )
        expected_mu = (
            6.62606957
            * 1e-34
            / (8 * np.pi**2 * vibrational_frequencies * 2.99792458 * 1e10)
        )
        expected_mu_prime = (
            expected_mu * expected_bav / (expected_mu + expected_bav)
        )
        # we got S_R,K = [4.13984133, 4.13984133, 1.00676497, -1.39084927] in J mol^-1 K^-1
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
            qrrho_thermochem_co2_1.freerot_entropy,
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
            qrrho_thermochem_co2_1.rrho_entropy,
            expected_rrho_entropy,
        )

        # S^qrrho_v = Σ(w(v_K) * S^rrho_v,K + (1 - w(v_K)) * S_R,K)
        # we got S^qrrho_v = 3.1441256211641195 J mol^-1 K^-1
        expected_qrrho_vibrational_entropy = np.sum(
            expected_entropy_damping_function * expected_rrho_entropy
            + (1 - expected_entropy_damping_function)
            * expected_freerot_entropy
        )
        assert np.isclose(
            qrrho_thermochem_co2_1.qrrho_vibrational_entropy,
            expected_qrrho_vibrational_entropy,
        )

        # q_t,c = (2 * pi * m * k_B * T / h^2)^(3/2) * (1 / c)
        # we got 469044.20805173856
        expected_translational_partition_function_concentration = (
            2
            * np.pi
            * (mol.mass / (6.02214129 * 1e23 * 1000))
            * 1.3806488
            * 1e-23
            * 298.15
            / (6.62606957 * 1e-34) ** 2
        ) ** (3 / 2) * (1 / (1.0 * 6.02214129 * 1e23 * 1000))
        assert np.isclose(
            qrrho_thermochem_co2_1.translational_partition_function_concentration,
            expected_translational_partition_function_concentration,
        )

        # S_t,c = R * [ln(q_t,c) + 1 + 3/2]
        # we got 129.3601627172439 J mol^-1 K^-1
        expected_translational_entropy_concentration = 8.314462145468951 * (
            np.log(expected_translational_partition_function_concentration)
            + 1
            + 3 / 2
        )
        # S^qrrho_tot = S_t,c + S_r + S^qrrho_v + S_e
        # we got 129.3601627172439 + 54.737291254812845 + 3.1441256211641195 + 0 = 187.24157959322085 J mol^-1 K^-1
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
                                mol.moments_of_inertia[-1]
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
            qrrho_thermochem_co2_1.qrrho_total_entropy,
            expected_qrrho_total_entropy,
        )

        # E0 in Hartree
        assert np.isclose(
            qrrho_thermochem_co2_1.energies, -188.444680, atol=1e-6
        )

        # ZPE in Hartree
        assert np.isclose(
            qrrho_thermochem_co2_1.zero_point_energy_hartree,
            0.011776,
            atol=1e-6,
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
            qrrho_thermochem_co2_1.enthalpy,
            expected_enthalpy,
        )
        assert np.isclose(
            qrrho_thermochem_co2_1.enthalpy, -188.429325, atol=1e-6
        )

        # T * S_tot in Hartree
        assert np.isclose(
            qrrho_thermochem_co2_1.entropy_times_temperature,
            0.021262,
            atol=1e-6,
        )

        # T * S^qrrho_tot in Hartree
        # 1 Hartree = 4.35974434 × 10^-18 Joules
        # 1 mol = 6.02214129 * 10^23 Particle
        # we got 0.021263029747636365 Hartree
        expected_qrrho_entropy_times_temperature = (
            298.15 * expected_qrrho_total_entropy
        ) / (4.35974434e-18 * 6.02214129 * 1e23)
        assert np.isclose(
            qrrho_thermochem_co2_1.qrrho_entropy_times_temperature,
            expected_qrrho_entropy_times_temperature,
        )
        assert np.isclose(
            qrrho_thermochem_co2_1.qrrho_entropy_times_temperature,
            0.021262,
            atol=1e-5,
        )

        # G = H - T * S_tot
        assert np.isclose(
            qrrho_thermochem_co2_1.gibbs_free_energy, -188.450587, atol=1e-6
        )

        # G^qrrho_qs = H - T * S^qrrho_tot
        # we got -188.45058826589445 Hartree
        expected_qrrho_gibbs_free_energy_qs = (
            expected_enthalpy - expected_qrrho_entropy_times_temperature
        )
        assert np.isclose(
            qrrho_thermochem_co2_1.qrrho_gibbs_free_energy,
            expected_qrrho_gibbs_free_energy_qs,
        )
        assert np.isclose(
            qrrho_thermochem_co2_1.qrrho_gibbs_free_energy_qs,
            -188.450588,
            atol=1e-6,
        )

        """Values from Goodvibes, as a reference:
                goodvibes -f 100 -c 1.0 -t 298.15 -q --bav "conf" co2.log
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
            qrrho_thermochem_co2_1.enthalpy_damping_function,
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
            qrrho_thermochem_co2_1.rrho_internal_energy,
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
            qrrho_thermochem_co2_1.qrrho_vibrational_internal_energy,
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
            qrrho_thermochem_co2_1.qrrho_total_internal_energy,
            expected_qrrho_total_internal_energy,
        )
        # H^qrrho = E0 + E^qrrho_tot + R * T
        # we got -188.42932658096436 Hartree
        expected_qrrho_enthalpy = g16_output.energies[-1] + (
            expected_qrrho_total_internal_energy + 8.314462145468951 * 298.15
        ) / (4.35974434e-18 * 6.02214129 * 1e23)
        assert np.isclose(
            qrrho_thermochem_co2_1.qrrho_enthalpy,
            expected_qrrho_enthalpy,
        )
        assert np.isclose(
            qrrho_thermochem_co2_1.qrrho_enthalpy, -188.429327, atol=1e-6
        )

        # G^qrrho_qh = H - T * S^qrrho_tot
        # we got -188.45058826589445 Hartree
        expected_qrrho_gibbs_free_energy_qh = (
            expected_enthalpy - expected_qrrho_entropy_times_temperature
        )
        assert np.isclose(
            qrrho_thermochem_co2_1.qrrho_gibbs_free_energy_qh,
            expected_qrrho_gibbs_free_energy_qh,
        )

        # G^qrrho_q = H^qrrho - T * S^qrrho_tot
        # we got -188.450589610712 Hartree
        expected_qrrho_gibbs_free_energy = (
            expected_qrrho_enthalpy - expected_qrrho_entropy_times_temperature
        )
        assert np.isclose(
            qrrho_thermochem_co2_1.qrrho_gibbs_free_energy,
            expected_qrrho_gibbs_free_energy,
        )
        assert np.isclose(
            qrrho_thermochem_co2_1.qrrho_gibbs_free_energy,
            -188.450589,
            atol=1e-6,
        )

        """Values from Goodvibes, as a reference:
                goodvibes -f 100 -c 0.5 -t 598.15 --qs grimme --bav "conf" co2.log
        Structure                                           E        ZPE             H        T.S     T.qh-S          G(T)       qh-G(T)
           ********************************************************************************************************************************
        o  co2                                       -188.444680   0.011776   -188.424452   0.049327   0.049327   -188.473778   -188.473779
           ********************************************************************************************************************************
        """
        qrrho_thermochem_co2_2 = Thermochemistry(
            filename=gaussian_co2_opt_outfile,
            temperature=598.15,  # in Kelvin
            concentration=0.5,  # in mol/L
        )
        assert np.isclose(
            qrrho_thermochem_co2_2.energies, -188.444680, atol=1e-6
        )
        assert np.isclose(
            qrrho_thermochem_co2_2.zero_point_energy_hartree,
            0.011776,
            atol=1e-6,
        )
        assert np.isclose(
            qrrho_thermochem_co2_2.enthalpy, -188.424452, atol=1e-6
        )
        assert np.isclose(
            qrrho_thermochem_co2_2.entropy_times_temperature,
            0.049327,
            atol=1e-6,
        )
        assert np.isclose(
            qrrho_thermochem_co2_2.qrrho_entropy_times_temperature,
            0.049327,
            atol=1e-5,
        )
        assert np.isclose(
            qrrho_thermochem_co2_2.gibbs_free_energy, -188.473778, atol=1e-6
        )
        assert np.isclose(
            qrrho_thermochem_co2_2.qrrho_gibbs_free_energy_qs,
            -188.473779,
            atol=1e-6,
        )

        """Values from Goodvibes, as a reference:
                goodvibes -f 1000 -c 1.0 -t 298.15 -q --bav "conf" co2.log
        Structure                                           E        ZPE             H          qh-H        T.S     T.qh-S          G(T)       qh-G(T)
           **********************************************************************************************************************************************
        o  co2                                       -188.444680   0.011776   -188.429325   -188.431976   0.021262   0.021781   -188.450587   -188.453757
           **********************************************************************************************************************************************
        """
        qrrho_thermochem_co2_3 = Thermochemistry(
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
            qrrho_thermochem_co2_3.entropy_damping_function,
            expected_damping_function,
        )
        assert np.allclose(
            qrrho_thermochem_co2_3.enthalpy_damping_function,
            expected_damping_function,
        )
        assert np.isclose(
            qrrho_thermochem_co2_3.energies, -188.444680, atol=1e-6
        )
        assert np.isclose(
            qrrho_thermochem_co2_3.zero_point_energy_hartree,
            0.011776,
            atol=1e-6,
        )
        assert np.isclose(
            qrrho_thermochem_co2_3.enthalpy, -188.429325, atol=1e-6
        )
        assert np.isclose(
            qrrho_thermochem_co2_3.qrrho_enthalpy, -188.431976, atol=1e-6
        )
        assert np.isclose(
            qrrho_thermochem_co2_3.entropy_times_temperature,
            0.021262,
            atol=1e-6,
        )
        assert np.isclose(
            qrrho_thermochem_co2_3.qrrho_entropy_times_temperature,
            0.021781,
            atol=1e-6,
        )
        assert np.isclose(
            qrrho_thermochem_co2_3.gibbs_free_energy, -188.450587, atol=1e-6
        )
        assert np.isclose(
            qrrho_thermochem_co2_3.qrrho_gibbs_free_energy,
            -188.453757,
            atol=1e-6,
        )


class TestThermochemistryHe:
    """He is used as an example to test the thermochemical properties of monoatomic molecules."""

    def test_thermochemistry_he_gaussian_output(self, gaussian_he_opt_outfile):
        """Values from Gaussian output
        Temperature   298.150 Kelvin.  Pressure   1.00000 Atm.
        Atom     1 has atomic number  2 and mass   4.00260
        Molecular mass:     4.00260 amu
        """
        assert os.path.exists(gaussian_he_opt_outfile)
        g16_output = Gaussian16Output(filename=gaussian_he_opt_outfile)
        assert g16_output.normal_termination
        assert g16_output.job_type == "opt"
        assert g16_output.num_atoms == 1
        mol = g16_output.molecule
        assert mol.empirical_formula == "He"
        assert g16_output.multiplicity == 1
        assert np.isclose(g16_output.energies[-1], -2.91512971456)
        assert np.isclose(mol.mass, 4.002602)  # weighted_atomic_mass=True
        assert np.isclose(g16_output.mass, 4.00260)
        assert mol.is_monoatomic

        thermochem2 = Thermochemistry(
            filename=gaussian_he_opt_outfile,
            temperature=298.15,
            pressure=1,
            weighted_atomic_mass=False,
        )
        """Values from Gaussian output
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
            thermochem2.vibrational_partition_function_bot, 0.100000e01
        )
        assert np.isclose(
            thermochem2.vibrational_partition_function_v0, 0.100000e01
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

    def test_thermochemistry_he_orca_output(self, orca_he_output_freq):
        """Values from ORCA output
        --------------------------
        THERMOCHEMISTRY AT 298.15K
        --------------------------

        Temperature         ...   298.15 K
        Pressure            ...     1.00 atm
        Total Mass          ...     4.00 AMU
        ...
        Point Group:  Kh, Symmetry Number:   1
        Rotational constants in cm-1:     0.000000     0.000000     0.000000
        """
        assert os.path.exists(orca_he_output_freq)
        orca_out = ORCAOutput(filename=orca_he_output_freq)
        assert orca_out.normal_termination
        #        assert orca_out.job_type == "opt"
        assert orca_out.natoms == 1
        mol = orca_out.molecule
        assert mol.empirical_formula == "He"
        assert orca_out.multiplicity == 1
        assert np.isclose(orca_out.energies[-1], -2.899160731389)
        assert np.isclose(mol.mass, 4.002602)  # weighted_atomic_mass=True
        assert np.isclose(orca_out.mass, 4.00)
        assert orca_out.rotational_symmetry_number == 1
        assert orca_out.rotational_constants_in_wavenumbers == [
            0,
            0,
            0,
        ]
        assert orca_out.vibrational_frequencies == []
        assert mol.is_monoatomic
        assert orca_out.num_vibration_modes == 0

    def test_thermochemistry_he_qrrho(self, gaussian_he_opt_outfile):
        """Values from Goodvibes, as a reference:
                goodvibes -f 1000 -c 0.5 -t 598.15 -q --bav "conf" he.log
        Structure                                           E        ZPE             H          qh-H        T.S     T.qh-S          G(T)       qh-G(T)
           **********************************************************************************************************************************************
        o  He                                          -2.915130   0.000000     -2.910394     -2.910394   0.025951   0.025951     -2.936345     -2.936345
           **********************************************************************************************************************************************
        """
        assert os.path.exists(gaussian_he_opt_outfile)
        g16_output = Gaussian16Output(filename=gaussian_he_opt_outfile)
        assert g16_output.normal_termination
        qrrho_thermochem_he = Thermochemistry(
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
            qrrho_thermochem_he.qrrho_entropy_times_temperature,
            0.025951,
            atol=1e-6,
        )
        assert np.isclose(
            qrrho_thermochem_he.gibbs_free_energy, -2.936345, atol=1e-6
        )
        assert np.isclose(
            qrrho_thermochem_he.qrrho_gibbs_free_energy, -2.936345, atol=1e-6
        )


class TestThermochemistryH2O:
    """Water is used as an example to test the thermochemical properties of non-linear polyatomic molecules."""

    def test_thermochemistry_water_gaussian_output(
        self, gaussian_mp2_outputfile
    ):
        """Values from Gaussian output
        Temperature   298.150 Kelvin.  Pressure   1.00000 Atm.
        Atom     1 has atomic number  8 and mass  15.99491
        Atom     2 has atomic number  1 and mass   1.00783
        Atom     3 has atomic number  1 and mass   1.00783
        Molecular mass:    18.01056 amu.
        Principal axes and moments of inertia in atomic units:
                                  1         2         3
            Eigenvalues --     2.23367   4.13757   6.37124
                  X           -0.00000   0.00000   1.00000
                  Y            1.00000   0.00000   0.00000
                  Z           -0.00000   1.00000   0.00000
        This molecule is an asymmetric top.
        Rotational symmetry number  2.
        Rotational temperatures (Kelvin)     38.77653    20.93353    13.59452
        Rotational constants (GHZ):         807.97175   436.18400   283.26385
        """
        assert os.path.exists(gaussian_mp2_outputfile)
        g16_output = Gaussian16Output(filename=gaussian_mp2_outputfile)
        assert g16_output.normal_termination
        assert g16_output.job_type == "opt"
        assert g16_output.num_atoms == 3
        mol = g16_output.molecule
        assert mol.empirical_formula == "H2O"
        assert g16_output.multiplicity == 1
        assert np.isclose(g16_output.energies[-1], -76.328992324258)
        assert np.isclose(mol.mass, 18.015)  # weighted_atomic_mass=True
        assert np.isclose(g16_output.mass, 18.01056)
        assert np.allclose(
            g16_output.moments_of_inertia, [0.62549131, 1.15863761, 1.78412891]
        )  # in amu Å^2
        assert np.allclose(
            mol.moments_of_inertia, [0.62560528, 1.15883759, 1.78444287]
        )  # in amu Å^2
        assert g16_output.rotational_symmetry_number == 2
        assert g16_output.rotational_temperatures == [
            38.77653,
            20.93353,
            13.59452,
        ]
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

        thermochem3 = Thermochemistry(
            filename=gaussian_mp2_outputfile,
            temperature=298.15,
            pressure=1,
            weighted_atomic_mass=False,
        )
        """Values from Gaussian output
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
            thermochem3.rotational_partition_function, 0.434320e02, atol=1e-1
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
            thermochem3.vibrational_partition_function_bot, 0.142038e-09
        )
        assert np.isclose(
            thermochem3.vibrational_partition_function_v0, 0.100039e01
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
        assert np.isclose(
            thermochem3.total_partition_function, 0.130534e09, atol=1e5
        )
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

    def test_thermochemistry_water_orca_output(self, water_output_gas_path):
        """Values from ORCA output
        --------------------------
        THERMOCHEMISTRY AT 298.15K
        --------------------------

        Temperature         ... 298.15 K
        Pressure            ... 1.00 atm
        Total Mass          ... 18.02 AMU
        ...
        freq.    1625.35  E(vib)   ...       0.00
        freq.    3875.61  E(vib)   ...       0.00
        freq.    3971.90  E(vib)   ...       0.00
        ...
        Point Group:  C2v, Symmetry Number:   2
        Rotational constants in cm-1:    26.416987    14.661432     9.428573
        """
        assert os.path.exists(water_output_gas_path)
        orca_out = ORCAOutput(filename=water_output_gas_path)
        assert orca_out.normal_termination
        assert orca_out.job_type == "opt"
        assert orca_out.natoms == 3
        mol = orca_out.molecule
        assert mol.empirical_formula == "H2O"
        assert orca_out.multiplicity == 1
        assert np.isclose(orca_out.energies[-1], -76.323311011349)
        assert np.isclose(mol.mass, 18.015)  # weighted_atomic_mass=True
        assert np.isclose(orca_out.mass, 18.02)
        assert orca_out.rotational_symmetry_number == 2
        assert orca_out.rotational_constants_in_wavenumbers == [
            26.416987,
            14.661432,
            9.428573,
        ]
        assert orca_out.vibrational_frequencies == [
            1625.35,
            3875.61,
            3971.90,
        ]
        assert not mol.is_monoatomic
        assert not mol.is_linear
        assert orca_out.num_vibration_modes == 3

    def test_thermochemistry_water_qrrho(self, gaussian_mp2_outputfile):
        """Values from Goodvibes, as a reference:
                goodvibes -f 500 -c 2.0 -t 1298.15 -q --bav "conf" water_mp2.log
        Structure                                           E        ZPE             H          qh-H        T.S     T.qh-S          G(T)       qh-G(T)
           **********************************************************************************************************************************************
        o  water_mp2                                  -76.328992   0.021410    -76.289193    -76.289224   0.098212   0.098221    -76.387404    -76.387445
           **********************************************************************************************************************************************
        """
        assert os.path.exists(gaussian_mp2_outputfile)
        g16_output = Gaussian16Output(filename=gaussian_mp2_outputfile)
        assert g16_output.normal_termination
        mol = g16_output.molecule
        qrrho_thermochem_water = Thermochemistry(
            filename=gaussian_mp2_outputfile,
            temperature=1298.15,  # in Kelvin
            concentration=2.0,  # in mol/L
            s_freq_cutoff=500,  # in cm^-1
            h_freq_cutoff=500,  # in cm^-1
        )
        vibrational_frequencies = np.array(g16_output.vibrational_frequencies)
        moments_of_inertia = np.array(mol.moments_of_inertia)
        expected_mu = (
            6.62606957
            * 1e-34
            / (8 * np.pi**2 * vibrational_frequencies * 2.99792458 * 1e10)
        )
        expected_i = moments_of_inertia / (6.02214129 * 1e23 * 1000) * 1e-10**2
        expected_b = 6.62606957 * 1e-34 / (8 * np.pi**2 * expected_i)
        # we got 509048848271.2939 Hz
        expected_average_rotational_constant = sum(expected_b) / len(
            expected_b
        )
        # B_av = (B_x + B_y + B_z) / 3
        # B_i = h / (8 * pi^2 * I_i) for i = x, y, z
        # we got B_av = 1.3016569220226748 * 10^-45 kg m^2
        expected_bav = (
            6.62606957 * 1e-34 / expected_average_rotational_constant
        )
        expected_mu_prime = (
            expected_mu * expected_bav / (expected_mu + expected_bav)
        )
        # we got S_R,K = [6.46111674, 2.9146504, 2.77994675] in J mol^-1 K^-1
        expected_freerot_entropy = 8.314462145468951 * (
            1 / 2
            + np.log(
                (
                    8
                    * np.pi**3
                    * expected_mu_prime
                    * 1.3806488
                    * 1e-23
                    * 1298.15
                    / (6.62606957 * 1e-34) ** 2
                )
                ** (1 / 2)
            )
        )
        assert np.allclose(
            qrrho_thermochem_water.freerot_entropy,
            expected_freerot_entropy,
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
            atol=1e-5,
        )
        assert np.isclose(
            qrrho_thermochem_water.qrrho_entropy_times_temperature,
            0.098221,
            atol=1e-5,
        )
        assert np.isclose(
            qrrho_thermochem_water.gibbs_free_energy, -76.387404, atol=1e-6
        )
        assert np.isclose(
            qrrho_thermochem_water.qrrho_gibbs_free_energy,
            -76.387445,
            atol=1e-6,
        )
