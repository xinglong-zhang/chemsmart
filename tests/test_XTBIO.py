import os.path

from chemsmart.io.xtb.output import XTBOutput


class TestXTBOutput:
    def test_sp_output(self, xtb_sp_outfile):
        assert os.path.exists(xtb_sp_outfile)
        xtb_output = XTBOutput(filename=xtb_sp_outfile)
        assert xtb_output.normal_termination
        assert xtb_output.route_string == "xtb xtbopt.coord"
        assert xtb_output.num_basis_functions == 6
        assert xtb_output.num_atomic_orbital == 6
        assert xtb_output.num_shells == 4
        assert xtb_output.num_electrons == 8
        assert xtb_output.hamiltonian == "GFN2-xTB"
        assert xtb_output.restart == False
        assert xtb_output.solvation == False
        assert xtb_output.pc_potential == False
        assert xtb_output.electronic_temperature == 300.0
        assert xtb_output.broyden_damping == 0.4
        assert xtb_output.net_charge == 0
        assert xtb_output.unpaired_electrons == 0
        assert xtb_output.homo_energy == -12.1467
        assert xtb_output.lumo_energy == 2.2442
        assert xtb_output.scc_energy == -5.104925504312
        assert xtb_output.isotropic_es == 0.031459394051
        assert xtb_output.anisotropic_es == 0.000394673573
        assert xtb_output.anisotropic_xc == -0.000882256681
        assert xtb_output.dispersion == -0.000141082937
        assert xtb_output.repulsion_energy == 0.034381060848
        assert xtb_output.total_charge == 0
        assert xtb_output.molecular_dipole == {
            "q_only": [-0.0, 0.0, 0.607],
            "full": [-0.0, -0.0, 0.872],
        }
        assert xtb_output.total_dipole == 2.217
        assert xtb_output.molecular_quadrupole == {
            "q_only": [1.311, 0.0, -0.492, 0.0, 0.0, -0.819],
            "q+dip": [1.747, 0.0, -0.572, -0.0, 0.0, -1.176],
            "full": [1.951, 0.0, -0.831, -0.0, 0.0, -1.121],
        }
        assert xtb_output.total_energy == -5.070544443464
        assert xtb_output.gradient_norm == 0.000075164743
        assert xtb_output.fmo_gap == 14.390891673350

    def test_opt_output(self, xtb_opt_outfile):
        assert os.path.exists(xtb_opt_outfile)
        xtb_output = XTBOutput(filename=xtb_opt_outfile)
        assert xtb_output.normal_termination
        assert xtb_output.geometry_optimization_converged
        assert xtb_output.route_string == "xtb coord --opt"
        assert xtb_output.homo_energy == -12.1467
        assert xtb_output.lumo_energy == 2.2442
        assert xtb_output.optimization_level == "normal"
        assert xtb_output.degrees_of_freedom == 3
        assert xtb_output.optimized_structure_block == [
            "$coord",
            "0.00000000011942       -0.00000000000000       -0.71677520925432      o",
            "1.45926122846511       -0.00000000000000        0.35838760458144      h",
            "-1.45926122858453        0.00000000000000        0.35838760467288      h",
            "$end",
            "",
        ]
        assert xtb_output.molecular_mass == 18.0152864
        assert xtb_output.center_of_mass == [0.0, 0.0, -0.3156364]
        assert xtb_output.moments_of_inertia == [
            0.5795334E+00,
            0.1202080E+01,
            0.1781614E+01,
        ]
        assert xtb_output.rotational_constants == [
            0.2908828E+02,
            0.1402372E+02,
            0.9462003E+01,
        ]
        assert xtb_output.total_energy == -5.070544443465
        assert xtb_output.gradient_norm == 0.000074994303
        assert xtb_output.fmo_gap == 14.390898452735

    def test_opt_gbsa_output(self, xtb_opt_gbsa_outfile):
        assert os.path.exists(xtb_opt_gbsa_outfile)
        xtb_output = XTBOutput(filename=xtb_opt_gbsa_outfile)
        assert xtb_output.normal_termination
        assert xtb_output.geometry_optimization_converged
        assert (
            xtb_output.route_string
            == "xtb pyridine_opt.sdf -gbsa acetonitrile -opt"
        )
        assert xtb_output.solvation
        assert xtb_output.solvent_model == "GBSA"
        assert xtb_output.solvent == "acetonitrile"
        assert xtb_output.dielectric_constant == 37.5
        assert xtb_output.free_energy_shift == 0.0020473
        assert xtb_output.temperature == 298.15
        assert xtb_output.density == 0.786
        assert xtb_output.solvent_mass == 41.05
        assert xtb_output.h_bond_correction == True
        assert xtb_output.ion_screening == False
        assert xtb_output.surface_tension == 1.0000E-05
        assert xtb_output.gsolv == -0.008095712200
        assert xtb_output.gelec == -0.000777170108
        assert xtb_output.gsasa == -0.009170771809
        assert xtb_output.ghb == -0.000195045636
        assert xtb_output.gshift == 0.002047275352
        assert xtb_output.degrees_of_freedom == 27
        assert xtb_output.optimized_structure_block == [
            "xtb: 6.7.1 (edcfbbe)",
            "02012514123D",
            "",
            "11 11  0     0  0            999 V2000",
            "-1.3731   -0.0003   -0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0",
            "1.4025    0.0003    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "0.6951    1.1927    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "0.6956   -1.1925    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "-0.6910    1.1338    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "-0.6905   -1.1342   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "2.4828    0.0005   -0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",
            "1.2041    2.1438   -0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",
            "1.2049   -2.1434   -0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",
            "-1.2892    2.0363   -0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",
            "-1.2884   -2.0368   -0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",
            "1  5  2  0  0  0  0",
            "1  6  1  0  0  0  0",
            "2  3  2  0  0  0  0",
            "2  4  1  0  0  0  0",
            "2  7  1  0  0  0  0",
            "3  5  1  0  0  0  0",
            "3  8  1  0  0  0  0",
            "4  6  2  0  0  0  0",
            "4  9  1  0  0  0  0",
            "5 10  1  0  0  0  0",
            "6 11  1  0  0  0  0",
            "M  END",
            "$$$$",
            "",
        ]
        assert xtb_output.total_energy == -16.158054009147
        assert xtb_output.gradient_norm == 0.000467095333
        assert xtb_output.fmo_gap == 3.234266958383

    def test_hess_output(self, xtb_hess_outfile):
        assert os.path.exists(xtb_hess_outfile)
        xtb_output = XTBOutput(filename=xtb_hess_outfile)
        assert xtb_output.normal_termination
        assert xtb_output.vibrational_frequencies == [
            -0.00,
            -0.00,
            -0.00,
            0.00,
            0.00,
            0.00,
            323.92,
            373.48,
            580.44,
            615.66,
            685.50,
            719.86,
            890.60,
            910.17,
            924.21,
            936.09,
            984.08,
            1077.66,
            1094.10,
            1096.43,
            1158.59,
            1203.16,
            1272.11,
            1279.08,
            1403.91,
            1444.90,
            1567.99,
            1578.36,
            3047.39,
            3048.65,
            3072.56,
            3092.02,
            3096.48,
        ]
        assert xtb_output.reduced_masses == [
            10.59,
            10.48,
            11.55,
            11.20,
            10.92,
            11.30,
            8.31,
            10.29,
            11.82,
            11.22,
            9.05,
            4.27,
            3.36,
            3.25,
            4.70,
            3.52,
            11.72,
            11.29,
            4.71,
            4.26,
            1.43,
            2.61,
            8.73,
            6.67,
            8.84,
            8.78,
            11.85,
            11.46,
            1.79,
            1.80,
            1.79,
            1.85,
            1.91,
        ]
        assert xtb_output.ir_intensities == [
            0.64,
            1.05,
            0.66,
            0.47,
            0.96,
            0.18,
            0.00,
            3.57,
            5.57,
            0.34,
            16.72,
            8.46,
            0.00,
            0.82,
            0.00,
            0.02,
            16.88,
            2.78,
            3.32,
            2.95,
            0.02,
            2.28,
            1.18,
            1.04,
            6.58,
            12.10,
            14.76,
            61.31,
            73.32,
            51.38,
            27.05,
            95.18,
            29.94,
        ]
        assert xtb_output.raman_intensities == [
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
        ]
        assert xtb_output.num_frequencies == 27
        assert xtb_output.num_imaginary_frequencies == 0
        assert xtb_output.symmetry == "c2v"
        assert xtb_output.rotational_number == 2
        assert xtb_output.scaling_factor == 1.0
        assert xtb_output.partition_function == {
            "vibrational": 1.99,
            "rotational": 0.417e05,
            "internal": 0.828e05,
            "translational": 0.681e27,
        }
        assert xtb_output.zero_point_energy == 0.085379802829
        assert xtb_output.grrho_without_zpve == -0.026863175517
        assert xtb_output.grrho_contribution == 0.058516627312
        assert xtb_output.total_enthalpy == -16.067220717226
        assert xtb_output.total_free_energy == -16.099537329248
