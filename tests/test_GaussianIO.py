import os.path

import numpy as np
from ase import units
from ase.symbols import Symbols

from chemsmart.io.gaussian.cube import GaussianCubeFile
from chemsmart.io.gaussian.input import Gaussian16Input, Gaussian16QMMMInput
from chemsmart.io.gaussian.output import (
    Gaussian16Output,
    Gaussian16OutputWithPBC,
    Gaussian16WBIOutput,
)
from chemsmart.io.gaussian.route import GaussianRoute
from chemsmart.io.molecules.structure import Molecule


class TestRouteString:
    def test_read_route_string_standard(self):
        s1a = "# opt freq mn15 def2svp"
        r1a = GaussianRoute(s1a)
        assert r1a.functional == "mn15"
        assert r1a.basis == "def2svp"
        assert r1a.jobtype == "opt"
        assert r1a.solv is False
        assert r1a.dieze_tag is None
        assert r1a.additional_opt_options_in_route is None
        assert r1a.additional_route_parameters is None

    def test_read_route_string_standard2(self):
        s1b = "# opt=(ts,calcfc,noeigentest) freq b3lyp/6-311+G(d,p) empiricaldispersion=gd3bj"
        r1b = GaussianRoute(s1b)
        assert r1b.functional == "b3lyp empiricaldispersion=gd3bj"
        assert r1b.basis == "6-311+g(d,p)"
        assert r1b.jobtype == "ts"
        assert r1b.solv is False
        assert r1b.dieze_tag is None
        assert (
            r1b.additional_opt_options_in_route is None
        )  # noeigentest prevents Gaussian from stopping
        # if no negative Hessian eigenvalue was found
        #                                                   # (not additional opt options for geometry opt)
        assert r1b.additional_route_parameters is None

    def test_read_route_string_standard3a(self):
        s1c = "# opt freq mn15 gen"
        r1c = GaussianRoute(s1c)
        assert r1c.functional == "mn15"
        assert r1c.basis == "gen"
        assert r1c.jobtype == "opt"
        assert r1c.solv is False
        assert r1c.dieze_tag is None
        assert r1c.additional_opt_options_in_route is None
        assert r1c.additional_route_parameters is None

    def test_read_route_string_standard3b(self):
        s1c = "# opt freq mn15 genecp"
        r1c = GaussianRoute(s1c)
        assert r1c.functional == "mn15"
        assert r1c.basis == "genecp"
        assert r1c.jobtype == "opt"
        assert r1c.solv is False
        assert r1c.dieze_tag is None
        assert r1c.additional_opt_options_in_route is None
        assert r1c.additional_route_parameters is None

    def test_read_route_string_standard4(self):
        s1d = "#t mn15 def2qzvp scrf=(smd,solvent=generic,read)"
        r1d = GaussianRoute(s1d)
        assert r1d.functional == "mn15"
        assert r1d.basis == "def2qzvp"
        assert r1d.jobtype == "sp"
        assert r1d.solv is True
        assert r1d.solvent_model == "smd"
        assert r1d.solvent_id == "generic,read"
        assert r1d.dieze_tag == "#t"
        assert r1d.additional_opt_options_in_route is None
        assert r1d.additional_route_parameters is None

    def test_read_route_string_standard5(self):
        s1e = "#p opt=modred freq tpsstpss/def2tzvp/fit empiricaldispersion=gd3bj scrf=(cpcm,solvent=toluene)"
        r1e = GaussianRoute(s1e)
        assert r1e.functional == "tpsstpss empiricaldispersion=gd3bj"
        assert (
            r1e.basis == "def2tzvp/fit"
        )  # density fitting basis set (for pure functionals)
        assert r1e.jobtype == "modred"
        assert r1e.solv is True
        assert r1e.solvent_model == "cpcm"
        assert r1e.solvent_id == "toluene"
        assert r1e.dieze_tag == "#p"
        assert r1e.additional_opt_options_in_route is None
        assert r1e.additional_route_parameters is None

    def test_read_route_string_standard6(self):
        s1f = "# mpw1pw91/6-311+G(2d,p) nmr=(GIAO,Mixed)"  # NMR route
        r1f = GaussianRoute(s1f)
        assert r1f.functional == "mpw1pw91"
        assert r1f.basis == "6-311+g(2d,p)"
        assert r1f.jobtype == "sp"
        assert r1f.dieze_tag is None
        assert r1f.additional_opt_options_in_route is None
        assert r1f.additional_route_parameters is None
        # assert r1f.additional_route_parameters == 'nmr=(GIAO,Mixed)'
        # TODO: nmr route to be specified

    def test_read_route_string_standard7(self):
        s1g = "# TD(nstates=30) wB97XD/def2SVP scrf(solvent=dichloroethane)"  # TD-DFT route
        r1g = GaussianRoute(s1g)
        assert r1g.functional == "wb97xd"
        assert r1g.basis == "def2svp"
        assert r1g.jobtype == "sp"
        assert r1g.dieze_tag is None
        assert r1g.solv is True
        assert r1g.solvent_model == "pcm"  # default solvet model in Gaussian
        assert r1g.solvent_id == "dichloroethane"
        assert r1g.additional_opt_options_in_route is None
        # TODO: TD-DFT route to be specified

    def test_read_route_string_nonstandard(self):
        s1 = "# pbepbe 6-31g(d,p)/auto force scrf=(dipole,solvent=water) pbc=gammaonly"
        r1 = GaussianRoute(s1)
        assert r1.solvent_model == "dipole"
        assert r1.solvent_id == "water"
        # TODO: fix nonstandard functional/basis (very rare cases such as this)

    def test_read_route_semiempirical(self):
        s1 = "# opt freq PM6"
        r1 = GaussianRoute(s1)
        assert r1.functional is None
        assert r1.basis is None
        assert r1.ab_initio is None
        assert r1.semiempirical == "PM6"
        assert r1.solv is False
        assert r1.dieze_tag is None
        assert r1.additional_opt_options_in_route is None
        assert r1.additional_route_parameters is None

    def test_read_route_string_opt_options(self):
        s2a = "# opt=(recalcfc=5) freq mn15 def2svp"
        r2a = GaussianRoute(s2a)
        assert r2a.functional == "mn15"
        assert r2a.basis == "def2svp"
        assert r2a.jobtype == "opt"
        assert r2a.solv is False
        assert r2a.dieze_tag is None
        assert r2a.additional_opt_options_in_route == "recalcfc=5"
        assert r2a.additional_route_parameters is None

        s2b = "# opt=(recalcfc=5,MaxStep=3,MaxCycles=128) freq mn15 def2svp"
        r2b = GaussianRoute(s2b)
        assert r2b.jobtype == "opt"
        assert (
            r2b.additional_opt_options_in_route
            == "recalcfc=5,maxstep=3,maxcycles=128"
        )

        s2c = "# opt=(ts,calcfc,noeigentest,recalcfc=5,MaxStep=3,MaxCycles=128) freq mn15 def2svp"
        r2c = GaussianRoute(s2c)
        assert r2c.jobtype == "ts"
        assert (
            r2c.additional_opt_options_in_route
            == "recalcfc=5,maxstep=3,maxcycles=128"
        )

    def test_read_additional_route_parameters(self):
        s3a = "# opt=(recalcfc=5) freq=numer pbepbe/def2svp nosymm guess=mix"
        r3a = GaussianRoute(s3a)
        assert r3a.jobtype == "opt"
        assert r3a.additional_opt_options_in_route == "recalcfc=5"
        assert r3a.freq is True
        assert r3a.numfreq is True
        assert r3a.solv is False
        assert r3a.functional == "pbepbe"
        assert r3a.basis == "def2svp"
        assert r3a.additional_route_parameters == "nosymm guess=mix"

    def test_solvent_in_route(self):
        s4a = (
            "# opt=(recalcfc=5) freq mn15 def2svp scrf=(dipole,solvent=water)"
        )
        r4a = GaussianRoute(s4a)
        assert r4a.additional_opt_options_in_route == "recalcfc=5"
        assert r4a.solvent_model == "dipole"
        assert r4a.solvent_id == "water"
        assert r4a.additional_solvent_options is None

        s4b = "# opt=(recalcfc=5) freq mn15 def2svp scrf=(smd,solvent=generic,read)"
        r4b = GaussianRoute(s4b)
        assert r4b.solvent_model == "smd"
        assert r4b.solvent_id == "generic,read"
        assert r4b.additional_solvent_options is None

        s4c = "# opt=(recalcfc=5) freq mn15 def2svp scrf=(cpcm,solvent=toluene,iterative)"
        r4c = GaussianRoute(s4c)
        assert r4c.solvent_model == "cpcm"
        assert r4c.solvent_id == "toluene"
        assert r4c.additional_solvent_options == "iterative"

        s4d = "# opt=(recalcfc=5) freq mn15 def2svp scrf=(cpcm,iterative,solvent=toluene)"
        r4d = GaussianRoute(s4d)
        assert r4d.solvent_model == "cpcm"
        assert r4d.solvent_id == "toluene"
        assert r4d.additional_solvent_options == "iterative"

        s4e = "# opt=(recalcfc=5) freq mn15 def2svp\n scrf=(cpcm,iterative,solvent=toluene)"
        r4e = GaussianRoute(s4e)
        assert r4e.solvent_model == "cpcm"
        assert r4e.solvent_id == "toluene"
        assert r4e.additional_solvent_options == "iterative"


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
        assert g16_input.molecule.num_atoms == 14
        assert g16_input.num_atoms == 14
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
        assert g16_input.jobtype == "opt"
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
        assert g16_frozen.jobtype == "opt"

    def test_partition(self, gaussian_qmmm_inputfile_2layer):
        assert os.path.exists(gaussian_qmmm_inputfile_2layer)
        g16_oniom = Gaussian16QMMMInput(
            filename=gaussian_qmmm_inputfile_2layer
        )
        assert g16_oniom.molecule.symbols.formula == "CH3CH3"
        assert g16_oniom.partition == {
            "high level atoms": ["2-5"],
            "low level atoms": ["6-9"],
        }

    def test_oniom_charge_multiplicity(self, gaussian_qmmm_inputfile_3layer):
        g16_oniom = Gaussian16QMMMInput(
            filename=gaussian_qmmm_inputfile_3layer
        )
        assert g16_oniom.oniom_charge == {
            "charge_total": "0",
            "int_charge": "0",
            "model_charge": "0",
        }
        assert g16_oniom.oniom_multiplicity == {
            "real_multiplicity": "1",
            "int_multiplicity": "1",
            "model_multiplicity": "1",
        }
        assert g16_oniom.real_charge == 0
        assert g16_oniom.int_charge == 0
        assert g16_oniom.model_charge == 0
        assert g16_oniom.real_multiplicity == 1
        assert g16_oniom.int_multiplicity == 1
        assert g16_oniom.model_multiplicity == 1

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
        assert g16_modred.jobtype == "modred"
        assert g16_modred.modred == [[2, 12], [9, 2]]
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
        assert g16_scan.jobtype == "modred"
        assert g16_scan.modred == {
            "coords": [[2, 12], [9, 2]],
            "num_steps": 10,
            "step_size": 0.05,
        }
        assert g16_scan.functional == "m062x"
        assert g16_scan.basis == "def2svp"

    def test_read_genecp_inputfile(self, gaussian_opt_genecp_inputfile):
        assert os.path.exists(gaussian_opt_genecp_inputfile)
        g16_genecp = Gaussian16Input(filename=gaussian_opt_genecp_inputfile)
        assert g16_genecp.molecule.symbols.formula == "PdC2O2C2O2H6"
        assert g16_genecp.molecule.empirical_formula == "C4H6O4Pd"
        assert g16_genecp.additional_opt_options_in_route is None
        assert g16_genecp.additional_route_parameters is None
        assert g16_genecp.jobtype == "opt"
        assert g16_genecp.functional == "mn15"
        assert g16_genecp.basis == "genecp"
        assert g16_genecp.genecp_section.genecp_type == "genecp"
        assert g16_genecp.genecp_section.light_elements == ["H", "C", "O"]
        assert g16_genecp.genecp_section.heavy_elements == ["Pd"]
        assert g16_genecp.genecp_section.light_elements_basis == "def2svp"
        assert g16_genecp.genecp_section.heavy_elements_basis == "def2-tzvppd"
        assert g16_genecp.molecule.frozen_atoms is None

    def test_read_gaussian_link_opt_input(self, gaussian_link_opt_input):
        assert os.path.exists(gaussian_link_opt_input)
        g16_link_opt = Gaussian16Input(filename=gaussian_link_opt_input)
        assert g16_link_opt.molecule.empirical_formula == "C7H5ClO"
        assert g16_link_opt.is_link
        assert (
            g16_link_opt.route_string
            == "# opt freq um062x def2svp scrf=(smd,solvent=dichloroethane) geom=check guess=read"
        )
        assert (
            g16_link_opt.additional_route_parameters == "geom=check guess=read"
        )
        assert g16_link_opt.additional_opt_options_in_route is None
        assert g16_link_opt.jobtype == "opt"
        assert g16_link_opt.functional == "um062x"
        assert g16_link_opt.basis == "def2svp"
        assert g16_link_opt.molecule.frozen_atoms is None

    def test_read_gaussian_link_ts_input(self, gaussian_link_ts_input):
        assert os.path.exists(gaussian_link_ts_input)
        g16_link_ts = Gaussian16Input(filename=gaussian_link_ts_input)
        assert g16_link_ts.molecule.empirical_formula == "C7H5ClO"
        assert g16_link_ts.is_link
        assert (
            g16_link_ts.route_string
            == "# opt=(ts,calcfc,noeigentest) freq um062x def2svp scrf=(smd,solvent=dichloroethane) geom=check guess=read"
        )
        assert (
            g16_link_ts.additional_route_parameters == "geom=check guess=read"
        )
        assert g16_link_ts.additional_opt_options_in_route is None
        assert g16_link_ts.jobtype == "ts"
        assert g16_link_ts.functional == "um062x"
        assert g16_link_ts.basis == "def2svp"
        assert g16_link_ts.molecule.frozen_atoms is None

    def test_read_gausssian_link_sp_input(self, gaussian_link_sp_input):
        assert os.path.exists(gaussian_link_sp_input)
        g16_link_sp = Gaussian16Input(filename=gaussian_link_sp_input)
        assert g16_link_sp.molecule.empirical_formula == "C7H5ClO"
        assert g16_link_sp.is_link
        assert (
            g16_link_sp.route_string
            == "# um062x def2tzvp scrf=(smd,solvent=dichloroethane) geom=check guess=read"
        )
        assert (
            g16_link_sp.additional_route_parameters == "geom=check guess=read"
        )
        assert g16_link_sp.additional_opt_options_in_route is None
        assert g16_link_sp.jobtype == "sp"
        assert g16_link_sp.functional == "um062x"
        assert g16_link_sp.basis == "def2tzvp"
        assert g16_link_sp.molecule.frozen_atoms is None

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
        assert g16_pbc_1d.additional_route_parameters == "scf=tight"
        assert g16_pbc_1d.jobtype == "sp"
        assert g16_pbc_1d.modred is None
        assert g16_pbc_1d.functional == "pbepbe"
        assert g16_pbc_1d.basis == "6-31g(d,p)/auto"


class TestGaussian16Output:
    def test_normal_termination_with_forces_and_frequencies(
        self, td_outputfile
    ):
        assert os.path.exists(td_outputfile)
        g16_output = Gaussian16Output(filename=td_outputfile)
        assert (
            g16_output.route_string
            == "# cam-b3lyp gen td(singlets,nstates=50,root=1)"
        )
        assert g16_output.spin == "unrestricted"
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

        assert (
            g16_output.total_core_hours
            == g16_output.total_service_unit
            == 361.7
        )
        assert g16_output.total_elapsed_walltime == 6.4
        mol = g16_output.molecule
        assert not mol.has_vibrations

    def test_singlet_opt_output(self, gaussian_singlet_opt_outfile):
        assert os.path.exists(gaussian_singlet_opt_outfile)
        g16_output = Gaussian16Output(filename=gaussian_singlet_opt_outfile)
        assert g16_output.normal_termination
        assert g16_output.molecule.num_atoms == 40
        assert g16_output.spin == "restricted"
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
        assert g16_output.fmo_gap == g16_output.alpha_fmo_gap
        assert np.allclose(
            g16_output.rotational_temperatures, [0.0078, 0.00354, 0.00256]
        )
        assert np.allclose(
            g16_output.rotational_constants_in_Hz,
            [0.16245 * 1e9, 0.07382 * 1e9, 0.05332 * 1e9],
        )
        assert g16_output.rotational_symmetry_number == 1
        mol = g16_output.molecule
        assert mol.has_vibrations
        assert mol.num_vib_frequencies == mol.num_vib_modes == 114
        vibrational_mode1 = [
            [0.0, 0.08, 0.03],
            [0.01, 0.01, 0.04],
            [0.01, 0.03, 0.04],
            [0.02, 0.09, 0.06],
            [0.0, 0.01, 0.03],
            [0.0, 0.06, 0.02],
            [0.01, 0.05, 0.04],
            [0.01, -0.05, 0.06],
            [0.02, -0.06, 0.07],
            [-0.02, -0.08, 0.05],
            [0.06, -0.06, 0.1],
            [0.06, -0.06, 0.09],
            [-0.05, -0.12, 0.0],
            [-0.05, -0.02, -0.04],
            [-0.11, 0.03, -0.14],
            [-0.16, -0.03, -0.21],
            [-0.16, -0.13, -0.16],
            [-0.11, -0.18, -0.06],
            [0.02, -0.15, 0.11],
            [-0.11, 0.11, -0.17],
            [-0.21, -0.18, -0.21],
            [-0.11, -0.26, -0.03],
            [-0.02, -0.25, 0.15],
            [0.08, -0.14, 0.18],
            [0.01, 0.04, 0.01],
            [0.04, 0.03, 0.03],
            [-0.01, 0.01, -0.01],
            [0.05, -0.0, 0.02],
            [-0.0, -0.02, -0.03],
            [0.03, -0.02, -0.01],
            [0.07, -0.01, 0.03],
            [-0.02, -0.04, -0.05],
            [0.14, -0.05, 0.13],
            [-0.21, 0.0, -0.29],
            [0.07, 0.07, 0.05],
            [-0.06, 0.02, -0.03],
            [0.04, -0.06, -0.03],
            [0.02, -0.07, -0.03],
            [0.07, -0.06, -0.04],
            [0.03, -0.08, -0.03],
        ]

        assert np.allclose(
            mol.vibrational_modes[0], vibrational_mode1, atol=1e-4
        )
        assert mol.vibrational_frequencies[0] == 11.9481

    def test_triplet_opt_output(self, gaussian_triplet_opt_outfile):
        assert os.path.exists(gaussian_triplet_opt_outfile)
        g16_output = Gaussian16Output(filename=gaussian_triplet_opt_outfile)
        assert g16_output.normal_termination
        assert g16_output.spin == "unrestricted"
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
        assert g16_output.num_unpaired_electrons == 2
        assert g16_output.multiplicity == 3
        # somo_energies should return list of 2 SOMOs for triplet
        assert len(g16_output.somo_energies) == 2
        assert g16_output.somo_energies == [
            -0.19177 * units.Hartree,
            -0.15673 * units.Hartree,
        ]
        assert g16_output.lowest_somo_energy == -0.19177 * units.Hartree
        assert g16_output.highest_somo_energy == -0.15673 * units.Hartree
        assert g16_output.alpha_homo_energy == -0.15673 * units.Hartree
        assert g16_output.beta_homo_energy == -0.18923 * units.Hartree
        assert g16_output.alpha_lumo_energy == -0.07423 * units.Hartree
        assert g16_output.beta_lumo_energy == -0.05025 * units.Hartree
        assert np.isclose(
            g16_output.fmo_gap,
            (min(-0.07423, -0.05025) - (-0.15673)) * units.Hartree,
        )
        assert np.isclose(
            g16_output.alpha_fmo_gap,
            (-0.07423 - (-0.15673)) * units.Hartree,
            rtol=1e-6,
        )
        assert np.isclose(
            g16_output.beta_fmo_gap,
            (-0.05025 - (-0.18923)) * units.Hartree,
            rtol=1e-6,
        )

    def test_quintet_opt_output(self, gaussian_quintet_opt_outfile):
        assert os.path.exists(gaussian_quintet_opt_outfile)
        g16_output = Gaussian16Output(filename=gaussian_quintet_opt_outfile)
        assert g16_output.spin == "unrestricted"
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
        assert g16_output.num_unpaired_electrons == 4
        assert g16_output.multiplicity == 5
        # somo_energies should return list of 4 SOMOs for quintet
        assert len(g16_output.somo_energies) == 4
        assert g16_output.somo_energies == [
            -0.22065 * units.Hartree,
            -0.21055 * units.Hartree,
            -0.19474 * units.Hartree,
            -0.18764 * units.Hartree,
        ]
        assert g16_output.lowest_somo_energy == -0.22065 * units.Hartree
        assert g16_output.highest_somo_energy == -0.18764 * units.Hartree
        assert g16_output.alpha_homo_energy == -0.18764 * units.Hartree
        assert g16_output.beta_homo_energy == -0.19564 * units.Hartree
        assert g16_output.alpha_lumo_energy == -0.03881 * units.Hartree
        assert g16_output.beta_lumo_energy == -0.06116 * units.Hartree
        assert np.isclose(
            g16_output.fmo_gap,
            (min(-0.03881, -0.06116) - (-0.18764)) * units.Hartree,
        )
        assert np.isclose(
            g16_output.alpha_fmo_gap,
            (-0.03881 - (-0.18764)) * units.Hartree,
            rtol=1e-6,
        )
        assert np.isclose(
            g16_output.beta_fmo_gap,
            (-0.06116 - (-0.19564)) * units.Hartree,
            rtol=1e-6,
        )

    def test_read_gaussian_link_opt_output_file(
        self, gaussian_link_opt_outputfile
    ):
        assert os.path.exists(gaussian_link_opt_outputfile)
        g16_link_opt = Gaussian16Output(filename=gaussian_link_opt_outputfile)
        assert (
            g16_link_opt.route_string
            == "# opt freq um062x def2svp geom=check guess=read"
        )
        assert g16_link_opt.is_link
        assert g16_link_opt.jobtype == "opt"
        assert g16_link_opt.normal_termination
        assert isinstance(g16_link_opt.molecule, Molecule)
        assert g16_link_opt.tddft_transitions == []
        assert len(g16_link_opt.alpha_occ_eigenvalues) == 8
        assert (
            g16_link_opt.alpha_occ_eigenvalues[0] == -19.77692 * units.Hartree
        )
        assert (
            g16_link_opt.alpha_occ_eigenvalues[-1] == -0.36639 * units.Hartree
        )
        assert len(g16_link_opt.alpha_virtual_eigenvalues) == 20
        assert (
            g16_link_opt.alpha_virtual_eigenvalues[0]
            == -0.06479 * units.Hartree
        )
        assert (
            g16_link_opt.alpha_virtual_eigenvalues[-1]
            == 3.87784 * units.Hartree
        )
        assert len(g16_link_opt.beta_occ_eigenvalues) == 8
        assert (
            g16_link_opt.beta_occ_eigenvalues[0] == -19.77692 * units.Hartree
        )
        assert (
            g16_link_opt.beta_occ_eigenvalues[-1] == -0.36639 * units.Hartree
        )
        assert len(g16_link_opt.beta_virtual_eigenvalues) == 20
        assert (
            g16_link_opt.beta_virtual_eigenvalues[0]
            == -0.06479 * units.Hartree
        )
        assert (
            g16_link_opt.beta_virtual_eigenvalues[-1]
            == 3.87784 * units.Hartree
        )
        assert np.isclose(
            g16_link_opt.fmo_gap, 0.3016 * units.Hartree, atol=1e-5
        )

    def test_read_gaussian_link_ts_output_file(
        self, gaussian_link_ts_outputfile
    ):
        assert os.path.exists(gaussian_link_ts_outputfile)
        g16_link_ts = Gaussian16Output(filename=gaussian_link_ts_outputfile)
        assert not g16_link_ts.normal_termination  # Error termination
        assert (
            g16_link_ts.route_string
            == "# opt=(ts,calcfc,noeigentest,maxstep=10) freq um062x def2svp geom=check guess=read"
        )
        assert g16_link_ts.is_link
        assert g16_link_ts.jobtype == "ts"
        assert len(g16_link_ts.vibrational_frequencies) == 0
        assert (
            g16_link_ts.num_vib_modes == g16_link_ts.num_vib_frequencies == 0
        )
        assert len(g16_link_ts.alpha_occ_eigenvalues) == 8
        assert (
            g16_link_ts.alpha_occ_eigenvalues[0] == -19.78334 * units.Hartree
        )
        assert (
            g16_link_ts.alpha_occ_eigenvalues[-1] == -0.38325 * units.Hartree
        )
        assert len(g16_link_ts.alpha_virtual_eigenvalues) == 20
        assert (
            g16_link_ts.alpha_virtual_eigenvalues[0]
            == -0.08312 * units.Hartree
        )
        assert (
            g16_link_ts.alpha_virtual_eigenvalues[-1]
            == 3.85967 * units.Hartree
        )
        assert len(g16_link_ts.beta_occ_eigenvalues) == 8
        assert g16_link_ts.beta_occ_eigenvalues[0] == -19.78334 * units.Hartree
        assert g16_link_ts.beta_occ_eigenvalues[-1] == -0.38325 * units.Hartree
        assert len(g16_link_ts.beta_virtual_eigenvalues) == 20
        assert (
            g16_link_ts.beta_virtual_eigenvalues[0] == -0.08312 * units.Hartree
        )
        assert (
            g16_link_ts.beta_virtual_eigenvalues[-1] == 3.85967 * units.Hartree
        )
        assert np.isclose(
            g16_link_ts.fmo_gap, 0.30013 * units.Hartree, atol=1e-5
        )

    def test_read_gaussian_link_modred_output_file(
        self, gaussian_link_modred_output
    ):
        assert os.path.exists(gaussian_link_modred_output)
        g16_link_modred = Gaussian16Output(
            filename=gaussian_link_modred_output
        )
        assert g16_link_modred.normal_termination
        assert (
            g16_link_modred.route_string
            == "# opt=modredundant freq umn15 def2svp geom=check guess=read"
        )
        assert g16_link_modred.is_link
        assert g16_link_modred.jobtype == "modred"
        assert isinstance(g16_link_modred.molecule, Molecule)
        assert len(g16_link_modred.vibrational_frequencies) == 126
        assert (
            g16_link_modred.num_vib_modes
            == g16_link_modred.num_vib_frequencies
            == 126
        )
        assert len(g16_link_modred.alpha_occ_eigenvalues) == 97
        assert (
            g16_link_modred.alpha_occ_eigenvalues[0]
            == -254.07064 * units.Hartree
        )
        assert (
            g16_link_modred.alpha_occ_eigenvalues[-1]
            == -0.24253 * units.Hartree
        )
        assert len(g16_link_modred.alpha_virtual_eigenvalues) == 339
        assert (
            g16_link_modred.alpha_virtual_eigenvalues[0]
            == 0.01660 * units.Hartree
        )
        assert (
            g16_link_modred.alpha_virtual_eigenvalues[-1]
            == 4.09404 * units.Hartree
        )
        assert len(g16_link_modred.beta_occ_eigenvalues) == 93
        assert (
            g16_link_modred.beta_occ_eigenvalues[0]
            == -254.07343 * units.Hartree
        )
        assert (
            g16_link_modred.beta_occ_eigenvalues[-1]
            == -0.26404 * units.Hartree
        )
        assert len(g16_link_modred.beta_virtual_eigenvalues) == 343
        assert (
            g16_link_modred.beta_virtual_eigenvalues[0]
            == -0.03779 * units.Hartree
        )
        assert (
            g16_link_modred.beta_virtual_eigenvalues[-1]
            == 4.21075 * units.Hartree
        )
        assert g16_link_modred.multiplicity == 5
        assert g16_link_modred.num_unpaired_electrons == 4
        assert g16_link_modred.somo_energies == [
            -0.30450 * units.Hartree,
            -0.29487 * units.Hartree,
            -0.26983 * units.Hartree,
            -0.24253 * units.Hartree,
        ]

    def test_read_gaussian_link_sp_output_file(
        self, gaussian_link_sp_outputfile
    ):
        assert os.path.exists(gaussian_link_sp_outputfile)
        g16_link_sp = Gaussian16Output(filename=gaussian_link_sp_outputfile)
        assert g16_link_sp.normal_termination
        assert (
            g16_link_sp.route_string
            == "# um062x def2tzvp scrf=(smd,solvent=chloroform) geom=check guess=read"
        )
        assert g16_link_sp.is_link
        assert g16_link_sp.jobtype == "sp"
        assert len(g16_link_sp.vibrational_frequencies) == 0
        assert (
            g16_link_sp.num_vib_modes == g16_link_sp.num_vib_frequencies == 0
        )
        assert len(g16_link_sp.alpha_occ_eigenvalues) == 8
        assert (
            g16_link_sp.alpha_occ_eigenvalues[0] == -19.78515 * units.Hartree
        )
        assert (
            g16_link_sp.alpha_occ_eigenvalues[-1] == -0.38742 * units.Hartree
        )
        assert len(g16_link_sp.alpha_virtual_eigenvalues) == 54
        assert (
            g16_link_sp.alpha_virtual_eigenvalues[0]
            == -0.08907 * units.Hartree
        )
        assert (
            g16_link_sp.alpha_virtual_eigenvalues[-1]
            == 43.63078 * units.Hartree
        )
        assert len(g16_link_sp.beta_occ_eigenvalues) == 8
        assert g16_link_sp.beta_occ_eigenvalues[0] == -19.78515 * units.Hartree
        assert g16_link_sp.beta_occ_eigenvalues[-1] == -0.38742 * units.Hartree
        assert len(g16_link_sp.beta_virtual_eigenvalues) == 54
        assert (
            g16_link_sp.beta_virtual_eigenvalues[0] == -0.08907 * units.Hartree
        )
        assert (
            g16_link_sp.beta_virtual_eigenvalues[-1]
            == 43.63078 * units.Hartree
        )
        assert np.isclose(
            g16_link_sp.fmo_gap, 0.29835 * units.Hartree, atol=1e-5
        )

    def test_read_failed_link_job(self, gaussian_failed_link_output):
        assert os.path.exists(gaussian_failed_link_output)
        g16_failed_link = Gaussian16Output(
            filename=gaussian_failed_link_output
        )
        assert not g16_failed_link.normal_termination
        molecule = Molecule.from_filepath(gaussian_failed_link_output)
        assert isinstance(molecule, Molecule)

    def test_read_genecp_outputfile(self, gaussian_ts_genecp_outfile):
        assert os.path.exists(gaussian_ts_genecp_outfile)
        g16_genecp = Gaussian16Output(filename=gaussian_ts_genecp_outfile)
        assert g16_genecp.normal_termination
        assert g16_genecp.gen_genecp == "genecp"
        assert (
            len(g16_genecp.vibrational_frequencies)
            == g16_genecp.num_atoms * 3 - 6
            == 138
        )
        assert g16_genecp.vibrational_frequencies[0] == -1138.1183
        assert g16_genecp.vibrational_frequencies[1] == 19.1625
        assert g16_genecp.vibrational_frequencies[-1] == 3291.3845
        assert (
            len(g16_genecp.reduced_masses)
            == g16_genecp.num_atoms * 3 - 6
            == 138
        )
        assert g16_genecp.reduced_masses[0] == 1.1629
        assert g16_genecp.reduced_masses[1] == 7.3337
        assert g16_genecp.reduced_masses[-1] == 1.0952
        assert (
            len(g16_genecp.force_constants)
            == g16_genecp.num_atoms * 3 - 6
            == 138
        )
        assert g16_genecp.force_constants[0] == 0.8875
        assert g16_genecp.force_constants[1] == 0.0016
        assert g16_genecp.force_constants[-1] == 6.9902
        assert (
            len(g16_genecp.ir_intensities)
            == g16_genecp.num_atoms * 3 - 6
            == 138
        )
        assert g16_genecp.ir_intensities[0] == 3338.6551
        assert g16_genecp.ir_intensities[1] == 0.1952
        assert g16_genecp.ir_intensities[-1] == 2.0786
        assert (
            len(g16_genecp.vibrational_mode_symmetries)
            == g16_genecp.num_atoms * 3 - 6
            == 138
        )
        # all members are "A"
        assert all(
            sym == "A" for sym in g16_genecp.vibrational_mode_symmetries
        )
        assert (
            g16_genecp.num_vib_modes == g16_genecp.num_vib_frequencies == 138
        )
        vibrational_mode1 = np.array(
            [
                [0.0, -0.0, 0.0],
                [0.0, -0.0, 0.0],
                [-0.0, 0.0, 0.0],
                [0.0, -0.0, 0.0],
                [0.0, -0.0, 0.0],
                [-0.0, 0.0, 0.0],
                [0.0, -0.0, 0.0],
                [0.0, -0.0, 0.0],
                [-0.0, 0.0, 0.0],
                [0.0, -0.0, 0.0],
                [0.0, -0.0, 0.0],
                [-0.0, 0.0, 0.0],
                [0.0, -0.0, 0.0],
                [0.0, -0.0, 0.0],
                [-0.0, -0.01, 0.0],
                [-0.0, 0.0, 0.0],
                [0.0, -0.0, 0.0],
                [0.0, -0.0, 0.0],
                [-0.0, 0.0, 0.0],
                [0.0, -0.0, 0.0],
                [0.0, -0.0, 0.0],
                [-0.0, -0.01, 0.0],
                [0.0, -0.0, -0.01],
                [0.0, -0.0, 0.0],
                [-0.0, 0.06, 0.01],
                [-0.03, 0.06, 0.01],
                [0.0, -0.01, 0.0],
                [-0.01, 0.02, 0.0],
                [-0.0, -0.0, 0.01],
                [0.84, -0.47, 0.23],
                [-0.03, 0.05, -0.01],
                [0.0, 0.0, 0.0],
                [0.0, -0.0, 0.0],
                [-0.0, -0.01, -0.0],
                [0.01, -0.0, 0.0],
                [0.01, 0.01, 0.0],
                [-0.06, 0.01, -0.02],
                [0.01, -0.02, -0.0],
                [0.0, 0.0, 0.0],
                [-0.01, 0.01, -0.0],
                [0.02, 0.05, 0.01],
                [-0.01, -0.01, -0.0],
                [-0.0, -0.01, -0.0],
                [-0.01, -0.01, -0.01],
                [-0.01, -0.01, -0.0],
                [0.0, -0.01, 0.0],
                [-0.01, -0.01, -0.0],
                [-0.0, -0.0, -0.0],
            ]
        )
        assert np.allclose(
            g16_genecp.vibrational_modes[0],
            vibrational_mode1,
            rtol=1e-4,
        )
        assert len(g16_genecp.forces) == 11
        assert g16_genecp.forces[0].shape == (g16_genecp.num_atoms, 3)
        assert np.allclose(
            g16_genecp.forces[0][0], [-0.002864142, 0.002344278, -0.003585424]
        )
        assert np.allclose(
            g16_genecp.forces[0][-1], [0.002024907, 0.001926310, 0.008510237]
        )
        assert np.allclose(
            g16_genecp.forces[-1][0], [0.000000455, 0.000001531, 0.000000084]
        )
        assert np.allclose(
            g16_genecp.forces[-1][-1],
            [-0.000000478, 0.000001912, -0.000001255],
        )
        assert np.allclose(
            g16_genecp.forces_in_eV_per_angstrom[0][0],
            [
                -0.002864142 * units.Hartree / units.Bohr,
                0.002344278 * units.Hartree / units.Bohr,
                -0.003585424 * units.Hartree / units.Bohr,
            ],
        )
        assert len(g16_genecp.input_orientations) == 12
        assert np.allclose(
            g16_genecp.input_orientations[0],
            np.array(
                [
                    [3.72556, -0.854649, -0.217208],
                    [4.885749, -1.558052, 0.105027],
                    [6.080932, -1.298797, -0.536227],
                    [6.145274, -0.307411, -1.501841],
                    [5.004473, 0.391844, -1.84732],
                    [3.786861, 0.111065, -1.237449],
                    [4.838296, -2.330332, 0.858004],
                    [6.963547, -1.863235, -0.278271],
                    [7.079418, -0.088285, -1.995881],
                    [5.03615, 1.151681, -2.614332],
                    [2.505974, -1.037436, 0.571314],
                    [2.56567, -1.227338, 1.957192],
                    [1.219386, -0.921784, 0.042266],
                    [3.510114, -1.303531, 2.478662],
                    [1.051683, -0.825776, -1.02621],
                    [0.318539, -1.129958, 2.129746],
                    [-0.566752, -1.159165, 2.762021],
                    [2.648939, 0.741301, -1.667406],
                    [2.483149, 2.125401, -1.34257],
                    [3.438094, 2.57224, -1.048449],
                    [2.115399, 2.616803, -2.249064],
                    [1.461967, 2.211518, -0.235457],
                    [0.142239, 1.869461, -0.512744],
                    [1.832056, 2.477652, 1.076731],
                    [-0.8136, 1.718584, 0.50609],
                    [-0.148409, 1.696957, -1.548689],
                    [0.894194, 2.401597, 2.095144],
                    [2.856474, 2.738449, 1.304251],
                    [-0.404149, 2.005665, 1.819225],
                    [-2.065762, 2.021708, 0.190945],
                    [-1.120902, 1.929883, 2.629083],
                    [1.483719, -1.309963, 2.709837],
                    [0.144064, -0.909348, 0.827524],
                    [-1.622566, -0.266572, 0.148648],
                    [-3.927221, -3.336699, -1.055001],
                    [-2.236825, -2.169488, -0.253601],
                    [-3.265798, 2.409445, -0.178222],
                    [-3.445303, 0.193621, -0.519209],
                    [-4.2222, -0.961582, -0.931433],
                    [-3.417102, -2.281826, -0.74205],
                    [-3.941817, 1.395953, -0.545514],
                    [-5.352071, 1.644776, -1.021478],
                    [-5.145709, -1.05139, -0.343722],
                    [-5.471368, 1.286147, -2.042],
                    [-6.059441, 1.115128, -0.386089],
                    [-5.554991, 2.710827, -0.981945],
                    [-4.503509, -0.893396, -1.990975],
                    [1.190887, 2.624408, 3.109776],
                ]
            ),
        )
        assert np.allclose(
            g16_genecp.input_orientations[-1],
            np.array(
                [
                    [3.785053, -0.796959, -0.248477],
                    [4.936026, -1.536354, 0.075967],
                    [6.165061, -1.283486, -0.527552],
                    [6.270151, -0.270499, -1.482714],
                    [5.138549, 0.458377, -1.840965],
                    [3.900464, 0.196499, -1.246085],
                    [4.857476, -2.339908, 0.811943],
                    [7.037222, -1.881031, -0.256517],
                    [7.22763, -0.061328, -1.964013],
                    [5.183685, 1.234043, -2.60889],
                    [2.52992, -1.02106, 0.506892],
                    [2.553976, -1.381448, 1.862919],
                    [1.248504, -0.829209, -0.025882],
                    [3.511848, -1.51029, 2.378812],
                    [1.099561, -0.558883, -1.070392],
                    [0.293033, -1.32343, 2.013095],
                    [-0.628222, -1.434285, 2.592938],
                    [2.804106, 0.881335, -1.680539],
                    [2.604367, 2.233236, -1.266643],
                    [3.548501, 2.650646, -0.877792],
                    [2.313292, 2.8026, -2.162191],
                    [1.513044, 2.275613, -0.228027],
                    [0.205942, 1.943357, -0.607426],
                    [1.79619, 2.464054, 1.128646],
                    [-0.808259, 1.728058, 0.340403],
                    [-0.011975, 1.783515, -1.670228],
                    [0.787047, 2.334516, 2.090623],
                    [2.817576, 2.708082, 1.436991],
                    [-0.49976, 1.967624, 1.696475],
                    [-2.044445, 1.966017, -0.060408],
                    [-1.279683, 1.83513, 2.452811],
                    [1.455439, -1.564035, 2.594152],
                    [0.154615, -0.931073, 0.737489],
                    [-1.69539, -0.250282, 0.106815],
                    [-4.144333, -3.345632, -0.467348],
                    [-2.354097, -2.100792, 0.012125],
                    [-3.210154, 2.424655, -0.457068],
                    [-3.524079, 0.186639, -0.477053],
                    [-4.358173, -0.979869, -0.6909],
                    [-3.599183, -2.270983, -0.366868],
                    [-3.954507, 1.421173, -0.656572],
                    [-5.37177, 1.656816, -1.110191],
                    [-5.263398, -0.969009, -0.059797],
                    [-5.558552, 1.150352, -2.068467],
                    [-6.081454, 1.241112, -0.379721],
                    [-5.540001, 2.733122, -1.218723],
                    [-4.703494, -1.060835, -1.736338],
                    [1.015959, 2.500207, 3.145815],
                ]
            ),
        )

        assert np.allclose(
            g16_genecp.input_orientations[-1],
            g16_genecp.input_orientations[-2],
        )  # structures for freq calc and the last opt step

        assert np.allclose(
            g16_genecp.input_orientations[-2],
            g16_genecp.input_orientations[-3],
        )  # structures for the second last and last opt steps

        assert len(g16_genecp.standard_orientations) == 12
        assert np.allclose(
            g16_genecp.standard_orientations[0],
            np.array(
                [
                    [3.670165, -0.853719, -0.227367],
                    [4.837265, -1.551256, 0.082421],
                    [6.028496, -1.272847, -0.558157],
                    [6.081717, -0.267891, -1.510325],
                    [4.933981, 0.425985, -1.843459],
                    [3.720429, 0.126285, -1.234499],
                    [4.798477, -2.334095, 0.824918],
                    [6.916657, -1.832993, -0.310071],
                    [7.012647, -0.033894, -2.003613],
                    [4.957067, 1.196435, -2.60012],
                    [2.454241, -1.057887, 0.561567],
                    [2.519139, -1.266096, 1.944578],
                    [1.165344, -0.946317, 0.037286],
                    [3.465544, -1.341105, 2.462655],
                    [0.994075, -0.837253, -1.02937],
                    [0.271694, -1.19075, 2.123928],
                    [-0.611688, -1.236306, 2.757905],
                    [2.575957, 0.752326, -1.653049],
                    [2.398935, 2.130374, -1.308979],
                    [3.350696, 2.581514, -1.011129],
                    [2.024594, 2.630825, -2.207793],
                    [1.379879, 2.192477, -0.198307],
                    [0.062478, 1.842683, -0.477007],
                    [1.75099, 2.443964, 1.116479],
                    [-0.889398, 1.669597, 0.542007],
                    [-0.229305, 1.681751, -1.514494],
                    [0.816437, 2.345853, 2.136048],
                    [2.773674, 2.710601, 1.34503],
                    [-0.479106, 1.94236, 1.857928],
                    [-2.144958, 1.966006, 0.23408],
                    [-1.193094, 1.849294, 2.668426],
                    [1.439881, -1.368422, 2.698666],
                    [0.091968, -0.95398, 0.825263],
                    [-1.681933, -0.317519, 0.159523],
                    [-3.962799, -3.391046, -1.080221],
                    [-2.280601, -2.220091, -0.267119],
                    [-3.349271, 2.348209, -0.126839],
                    [-3.510318, 0.135742, -0.49755],
                    [-4.278163, -1.020503, -0.923579],
                    [-3.4611, -2.336103, -0.754171],
                    [-4.017363, 1.333927, -0.506257],
                    [-5.430948, 1.576846, -0.975342],
                    [-5.199346, -1.126381, -0.334895],
                    [-5.549725, 1.231103, -2.000362],
                    [-6.132043, 1.03243, -0.345505],
                    [-5.643056, 2.640443, -0.92079],
                    [-4.562767, -0.940374, -1.981404],
                    [1.113773, 2.557428, 3.152894],
                ]
            ),
        )

        last_structure_positions = np.array(
            [
                [3.738125, -0.799262, -0.33422],
                [4.88285, -1.590061, -0.131866],
                [6.120627, -1.239348, -0.664417],
                [6.241126, -0.072719, -1.42203],
                [5.11621, 0.713671, -1.659267],
                [3.869481, 0.355676, -1.136943],
                [4.792375, -2.511643, 0.447817],
                [6.987624, -1.879741, -0.492953],
                [7.205597, 0.214468, -1.845689],
                [5.17365, 1.613668, -2.275798],
                [2.472231, -1.148636, 0.352523],
                [2.476846, -1.74526, 1.622674],
                [1.198691, -0.858815, -0.154214],
                [3.427248, -1.968556, 2.119805],
                [1.064719, -0.405795, -1.135583],
                [0.214191, -1.704367, 1.751149],
                [-0.715188, -1.912566, 2.289781],
                [2.780751, 1.112162, -1.456551],
                [2.578422, 2.369462, -0.810716],
                [3.518074, 2.706377, -0.341283],
                [2.300966, 3.090811, -1.594029],
                [1.47304, 2.230989, 0.2044],
                [0.170519, 1.977876, -0.245295],
                [1.737958, 2.173055, 1.576517],
                [-0.857089, 1.601699, 0.635529],
                [-0.033147, 2.01122, -1.322326],
                [0.715409, 1.878715, 2.48663],
                [2.755566, 2.353372, 1.936817],
                [-0.566695, 1.594058, 2.01651],
                [-2.087111, 1.913135, 0.267447],
                [-1.35723, 1.332419, 2.726778],
                [1.367959, -2.050215, 2.29513],
                [0.094194, -1.090098, 0.564326],
                [-1.745423, -0.299056, 0.041066],
                [-4.193282, -3.230814, -1.108266],
                [-2.406973, -2.099882, -0.390937],
                [-3.246217, 2.440636, -0.056249],
                [-3.564925, 0.243581, -0.479381],
                [-4.398654, -0.862143, -0.908862],
                [-3.647121, -2.193895, -0.810501],
                [-3.99004, 1.492341, -0.441352],
                [-5.400401, 1.811767, -0.864176],
                [-5.312433, -0.959801, -0.297868],
                [-5.575158, 1.48526, -1.899811],
                [-6.120992, 1.275751, -0.228973],
                [-5.564679, 2.890946, -0.781104],
                [-4.72977, -0.753691, -1.956403],
                [0.930185, 1.852433, 3.557361],
            ]
        )

        assert np.allclose(
            g16_genecp.standard_orientations[-1],
            last_structure_positions,
        )

        assert np.allclose(
            g16_genecp.standard_orientations[-1],
            g16_genecp.standard_orientations[-2],
        )  # structures for freq calc and the last opt step

        assert np.allclose(
            g16_genecp.standard_orientations[-2],
            g16_genecp.standard_orientations[-3],
        )  # structures for the second last and last opt steps

        assert (
            g16_genecp.optimized_structure.empirical_formula == "C21H19N3O4Pd"
        )
        assert g16_genecp.additional_opt_options_in_route == "maxstep=10"
        assert g16_genecp.additional_route_parameters is None
        assert g16_genecp.jobtype == "ts"
        assert g16_genecp.functional == "mn15"
        assert g16_genecp.basis == "genecp"
        assert g16_genecp.optimized_structure.frozen_atoms is None
        assert (
            len(g16_genecp.all_structures) == 11
        )  # 11 structures altogether, as shown in GaussView
        assert g16_genecp.optimized_structure.positions.shape == (48, 3)
        assert np.allclose(
            g16_genecp.optimized_structure.positions,
            last_structure_positions,
        )
        assert np.allclose(
            g16_genecp.get_molecule().positions, last_structure_positions
        )

        assert len(g16_genecp.get_molecule(index="1:4")) == 3
        assert np.allclose(
            g16_genecp.get_molecule(index=":4")[-1].positions[0],
            [3.69135800, -0.83587500, -0.25754700],
        )
        assert len(g16_genecp.get_molecule(index="4:")) == 8

        mol = g16_genecp.molecule
        assert np.allclose(mol.positions, last_structure_positions, rtol=1e-4)
        mol2 = mol.vibrationally_displaced(mode_idx=1, amp=0.5)
        assert np.allclose(
            mol2.positions[29], [-1.668781, 1.679069, 0.38199], rtol=1e-4
        )

        mol3 = mol.vibrationally_displaced(mode_idx=1, amp=-0.5)
        assert np.allclose(
            mol3.positions[29], [-2.505441, 2.147201, 0.152904], rtol=1e-4
        )

    def test_read_frozen_opt_outputfile(self, gaussian_frozen_opt_outfile):
        assert os.path.exists(gaussian_frozen_opt_outfile)
        g16_frozen = Gaussian16Output(
            filename=gaussian_frozen_opt_outfile, use_frozen=True
        )
        assert g16_frozen.normal_termination
        assert g16_frozen.num_atoms == 14
        assert g16_frozen.tddft_transitions == []
        assert len(g16_frozen.alpha_occ_eigenvalues) == 36
        assert (
            g16_frozen.alpha_occ_eigenvalues[0] == -102.65018 * units.Hartree
        )
        assert g16_frozen.alpha_occ_eigenvalues[-1] == -0.31442 * units.Hartree
        assert len(g16_frozen.alpha_virtual_eigenvalues) == 119
        assert (
            g16_frozen.alpha_virtual_eigenvalues[0] == -0.03944 * units.Hartree
        )
        assert (
            g16_frozen.alpha_virtual_eigenvalues[-1] == 3.66749 * units.Hartree
        )
        assert g16_frozen.has_frozen_coordinates
        assert g16_frozen.frozen_coordinate_indices == [
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
        ]
        assert g16_frozen.frozen_elements == [
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
        ]
        assert g16_frozen.free_elements == ["C", "O", "H", "Cl"]
        assert g16_frozen.frozen_atoms_masks == [
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            0,
            0,
            0,
            0,
        ]
        assert g16_frozen.optimized_structure.frozen_atoms == [
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            0,
            0,
            0,
            0,
        ]
        assert g16_frozen.optimized_structure.energy == -804.614710796
        assert g16_frozen.free_coordinate_indices == [11, 12, 13, 14]
        assert g16_frozen.num_vib_modes == g16_frozen.num_vib_frequencies == 12
        assert np.allclose(
            g16_frozen.vibrational_modes[0],
            np.array(
                [
                    [0.0, -0.0, 0.37],
                    [0.0, -0.0, 0.91],
                    [-0.0, 0.0, 0.18],
                    [-0.0, 0.0, 0.02],
                ]
            ),
            rtol=1e-4,
        )
        assert np.allclose(
            g16_frozen.vibrational_modes[-1],
            np.array(
                [
                    [-0.03, -0.08, -0.00],
                    [0.0, -0.0, 0.0],
                    [0.37, 0.93, 0.00],
                    [-0.00, -0.00, -0.00],
                ]
            ),
            rtol=1e-4,
        )

        g16_hide_frozen = Gaussian16Output(
            filename=gaussian_frozen_opt_outfile, use_frozen=False
        )
        assert g16_hide_frozen.normal_termination
        assert g16_hide_frozen.num_atoms == 14
        assert g16_hide_frozen.tddft_transitions == []
        assert len(g16_hide_frozen.alpha_occ_eigenvalues) == 36
        assert (
            g16_hide_frozen.alpha_occ_eigenvalues[0]
            == -102.65018 * units.Hartree
        )
        assert (
            g16_hide_frozen.alpha_occ_eigenvalues[-1]
            == -0.31442 * units.Hartree
        )
        assert len(g16_hide_frozen.alpha_virtual_eigenvalues) == 119
        assert (
            g16_hide_frozen.alpha_virtual_eigenvalues[0]
            == -0.03944 * units.Hartree
        )
        assert (
            g16_hide_frozen.alpha_virtual_eigenvalues[-1]
            == 3.66749 * units.Hartree
        )
        assert g16_frozen.modred is None

        # has frozen coordinates
        assert g16_hide_frozen.has_frozen_coordinates
        assert g16_hide_frozen.frozen_coordinate_indices == [
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
        ]
        assert g16_hide_frozen.free_coordinate_indices == [11, 12, 13, 14]
        assert g16_hide_frozen.frozen_elements == [
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
        ]
        assert g16_hide_frozen.free_elements == ["C", "O", "H", "Cl"]
        assert g16_hide_frozen.frozen_atoms_masks == [
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            0,
            0,
            0,
            0,
        ]

        # since use_frozen is False, this is not included in the output structure
        assert g16_hide_frozen.optimized_structure.frozen_atoms is None

        assert g16_hide_frozen.optimized_structure.energy == -804.614710796
        assert (
            g16_hide_frozen.num_vib_modes
            == g16_hide_frozen.num_vib_frequencies
            == 12
        )
        assert np.allclose(
            g16_hide_frozen.vibrational_modes[0],
            np.array(
                [
                    [0.0, -0.0, 0.37],
                    [0.0, -0.0, 0.91],
                    [-0.0, 0.0, 0.18],
                    [-0.0, 0.0, 0.02],
                ]
            ),
            rtol=1e-4,
        )
        assert g16_hide_frozen.modred is None

    def test_read_modred_outputfile(self, gaussian_failed_modred_outfile):
        assert os.path.exists(gaussian_failed_modred_outfile)
        g16_modred = Gaussian16Output(
            filename=gaussian_failed_modred_outfile, use_frozen=True
        )
        assert not g16_modred.normal_termination
        assert g16_modred.num_atoms == 10
        assert g16_modred.tddft_transitions == []
        assert len(g16_modred.alpha_occ_eigenvalues) == 23
        assert g16_modred.alpha_occ_eigenvalues[0] == -19.70039 * units.Hartree
        assert g16_modred.alpha_occ_eigenvalues[-1] == -0.30022 * units.Hartree
        assert len(g16_modred.alpha_virtual_eigenvalues) == 81
        assert (
            g16_modred.alpha_virtual_eigenvalues[0] == 0.03581 * units.Hartree
        )
        assert (
            g16_modred.alpha_virtual_eigenvalues[-1] == 3.95271 * units.Hartree
        )
        assert g16_modred.modred == [[4, 8], [5, 8], [4, 6]]

    def test_read_scan_outputfile(self, gaussian_failed_scan_outfile):
        assert os.path.exists(gaussian_failed_scan_outfile)
        g16_scan = Gaussian16Output(
            filename=gaussian_failed_scan_outfile,
            use_frozen=True,
            include_intermediate=False,
        )
        assert not g16_scan.normal_termination
        assert g16_scan.num_atoms == 110
        assert len(g16_scan.alpha_occ_eigenvalues) == 217
        assert g16_scan.alpha_occ_eigenvalues[0] == -19.75707 * units.Hartree
        assert g16_scan.alpha_occ_eigenvalues[-1] == -0.33917 * units.Hartree
        assert len(g16_scan.alpha_virtual_eigenvalues) == 895
        assert (
            g16_scan.alpha_virtual_eigenvalues[0] == -0.15548 * units.Hartree
        )
        assert (
            g16_scan.alpha_virtual_eigenvalues[-1] == 3.80727 * units.Hartree
        )
        assert g16_scan.modred == {
            "coords": [[1, 19]],
            "num_steps": 10,
            "step_size": -0.1,
        }
        assert g16_scan.num_steps == 11
        assert len(g16_scan.all_structures) == 1

        g16_scan_all_int = Gaussian16Output(
            filename=gaussian_failed_scan_outfile,
            use_frozen=True,
            include_intermediate=True,
        )
        assert not g16_scan_all_int.normal_termination
        assert g16_scan_all_int.num_atoms == 110
        assert len(g16_scan_all_int.alpha_occ_eigenvalues) == 217
        assert (
            g16_scan_all_int.alpha_occ_eigenvalues[0]
            == -19.75707 * units.Hartree
        )
        assert (
            g16_scan_all_int.alpha_occ_eigenvalues[-1]
            == -0.33917 * units.Hartree
        )
        assert len(g16_scan_all_int.alpha_virtual_eigenvalues) == 895

        assert len(g16_scan_all_int.all_structures) == 9
        # 10 orientations with last structure (failed job) removed

    def test_read_hirshfeld_charges_outputfile(
        self, gaussian_hirshfeld_outfile
    ):
        assert os.path.exists(gaussian_hirshfeld_outfile)
        g16_hirshfeld = Gaussian16Output(filename=gaussian_hirshfeld_outfile)
        assert g16_hirshfeld.normal_termination
        assert g16_hirshfeld.num_atoms == 33
        assert len(g16_hirshfeld.mulliken_atomic_charges) == 33
        assert g16_hirshfeld.mulliken_atomic_charges["O1"] == -0.359649
        assert g16_hirshfeld.mulliken_atomic_charges["O2"] == -0.317260
        assert g16_hirshfeld.mulliken_atomic_charges["C3"] == -0.090440
        assert g16_hirshfeld.mulliken_atomic_charges["H33"] == 0.183443
        assert len(g16_hirshfeld.mulliken_atomic_charges_heavy_atoms) == 15
        assert (
            g16_hirshfeld.mulliken_atomic_charges_heavy_atoms["O1"]
            == -0.359649
        )
        assert (
            g16_hirshfeld.mulliken_atomic_charges_heavy_atoms["O2"]
            == -0.317260
        )
        assert (
            g16_hirshfeld.mulliken_atomic_charges_heavy_atoms["C3"] == 0.064107
        )

        assert len(g16_hirshfeld.hirshfeld_charges) == 33
        assert g16_hirshfeld.hirshfeld_charges["O1"] == -0.222183
        assert g16_hirshfeld.hirshfeld_charges["O2"] == -0.175602
        assert g16_hirshfeld.hirshfeld_charges["C3"] == -0.030469
        assert g16_hirshfeld.hirshfeld_charges["H33"] == 0.050255
        assert g16_hirshfeld.hirshfeld_spin_densities["O1"] == 0.000000
        assert g16_hirshfeld.hirshfeld_spin_densities["O2"] == 0.000000
        assert g16_hirshfeld.hirshfeld_spin_densities["C3"] == 0.000000
        assert g16_hirshfeld.hirshfeld_spin_densities["H33"] == 0.000000
        assert np.allclose(
            g16_hirshfeld.hirshfeld_dipoles["O1"],
            np.array([-0.121486, -0.118753, -0.104620]),
        )
        assert np.allclose(
            g16_hirshfeld.hirshfeld_dipoles["O2"],
            np.array([0.024882, -0.086174, 0.133652]),
        )
        assert np.allclose(
            g16_hirshfeld.hirshfeld_dipoles["C3"],
            np.array([-0.008461, -0.029311, -0.015572]),
        )
        assert np.allclose(
            g16_hirshfeld.hirshfeld_dipoles["H33"],
            np.array([-0.143072, 0.058847, -0.063056]),
        )
        assert g16_hirshfeld.hirshfeld_cm5_charges["O1"] == -0.309536
        assert g16_hirshfeld.hirshfeld_cm5_charges["O2"] == -0.278764
        assert g16_hirshfeld.hirshfeld_cm5_charges["C3"] == -0.089643
        assert len(g16_hirshfeld.hirshfeld_charges_heavy_atoms) == 15
        assert g16_hirshfeld.hirshfeld_charges_heavy_atoms["O1"] == -0.222183
        assert g16_hirshfeld.hirshfeld_charges_heavy_atoms["O2"] == -0.175602
        assert g16_hirshfeld.hirshfeld_charges_heavy_atoms["C3"] == 0.011726
        assert (
            g16_hirshfeld.hirshfeld_cm5_charges_heavy_atoms["O1"] == -0.309536
        )
        assert (
            g16_hirshfeld.hirshfeld_cm5_charges_heavy_atoms["O2"] == -0.278764
        )
        assert (
            g16_hirshfeld.hirshfeld_cm5_charges_heavy_atoms["C3"] == 0.012018
        )

    def test_read_hirshfeld_rc_charges_outputfile(
        self, gaussian_rc_hirshfeld_outfile
    ):
        assert os.path.exists(gaussian_rc_hirshfeld_outfile)
        g16_rc_hirshfeld = Gaussian16Output(
            filename=gaussian_rc_hirshfeld_outfile
        )
        assert g16_rc_hirshfeld.normal_termination
        assert g16_rc_hirshfeld.charge == 1
        assert g16_rc_hirshfeld.multiplicity == 2
        assert g16_rc_hirshfeld.num_atoms == 33
        assert len(g16_rc_hirshfeld.mulliken_atomic_charges) == 33
        assert g16_rc_hirshfeld.mulliken_atomic_charges["O1"] == 0.020200
        assert g16_rc_hirshfeld.mulliken_atomic_charges["O2"] == -0.317365
        assert g16_rc_hirshfeld.mulliken_atomic_charges["C3"] == -0.087929
        assert g16_rc_hirshfeld.mulliken_atomic_charges["H33"] == 0.183814
        assert len(g16_rc_hirshfeld.mulliken_atomic_charges_heavy_atoms) == 15
        assert (
            g16_rc_hirshfeld.mulliken_atomic_charges_heavy_atoms["O1"]
            == 0.020200
        )
        assert (
            g16_rc_hirshfeld.mulliken_atomic_charges_heavy_atoms["O2"]
            == -0.317365
        )
        assert (
            g16_rc_hirshfeld.mulliken_atomic_charges_heavy_atoms["C3"]
            == 0.112863
        )

        assert len(g16_rc_hirshfeld.mulliken_spin_densities) == 33
        assert g16_rc_hirshfeld.mulliken_spin_densities["O1"] == 0.684808
        assert g16_rc_hirshfeld.mulliken_spin_densities["O2"] == 0.002091
        assert g16_rc_hirshfeld.mulliken_spin_densities["C3"] == 0.013665
        assert g16_rc_hirshfeld.mulliken_spin_densities["H33"] == -0.000003
        assert (
            g16_rc_hirshfeld.mulliken_spin_densities_heavy_atoms["O1"]
            == 0.684808
        )
        assert (
            g16_rc_hirshfeld.mulliken_spin_densities_heavy_atoms["O2"]
            == 0.002091
        )
        assert (
            g16_rc_hirshfeld.mulliken_spin_densities_heavy_atoms["C3"]
            == 0.018679
        )
        assert g16_rc_hirshfeld.hirshfeld_charges["O1"] == 0.100231
        assert g16_rc_hirshfeld.hirshfeld_charges["O2"] == -0.169411
        assert g16_rc_hirshfeld.hirshfeld_charges["C3"] == -0.000709
        assert g16_rc_hirshfeld.hirshfeld_charges["H33"] == 0.050632
        assert g16_rc_hirshfeld.hirshfeld_charges_heavy_atoms["O1"] == 0.100231
        assert (
            g16_rc_hirshfeld.hirshfeld_charges_heavy_atoms["O2"] == -0.169411
        )
        assert g16_rc_hirshfeld.hirshfeld_charges_heavy_atoms["C3"] == 0.078562
        assert g16_rc_hirshfeld.hirshfeld_spin_densities["O1"] == 0.610176
        assert g16_rc_hirshfeld.hirshfeld_spin_densities["O2"] == 0.003115
        assert g16_rc_hirshfeld.hirshfeld_spin_densities["C3"] == 0.019163
        assert g16_rc_hirshfeld.hirshfeld_spin_densities["H33"] == 0.000005
        assert (
            g16_rc_hirshfeld.hirshfeld_spin_densities_heavy_atoms["O1"]
            == 0.610176
        )
        assert (
            g16_rc_hirshfeld.hirshfeld_spin_densities_heavy_atoms["O2"]
            == 0.003115
        )
        assert (
            g16_rc_hirshfeld.hirshfeld_spin_densities_heavy_atoms["C3"]
            == 0.021583
        )

    def test_read_mp2_outputfile(self, gaussian_mp2_outputfile):
        assert os.path.exists(gaussian_mp2_outputfile)
        g16_mp2 = Gaussian16Output(filename=gaussian_mp2_outputfile)
        assert g16_mp2.normal_termination
        assert g16_mp2.num_atoms == 3
        assert g16_mp2.tddft_transitions == []
        assert len(g16_mp2.alpha_occ_eigenvalues) == 5
        assert g16_mp2.alpha_occ_eigenvalues[0] == -20.56810 * units.Hartree
        assert g16_mp2.alpha_occ_eigenvalues[-1] == -0.51014 * units.Hartree
        assert len(g16_mp2.alpha_virtual_eigenvalues) == 87
        assert g16_mp2.alpha_virtual_eigenvalues[0] == 0.02937 * units.Hartree
        assert (
            g16_mp2.alpha_virtual_eigenvalues[-1] == 15.70360 * units.Hartree
        )
        assert len(g16_mp2.mp2_energies) == 5
        assert g16_mp2.mp2_energies[0] == -76.32896706205
        assert g16_mp2.scf_energies[0] == -76.0599359638

    def test_read_oniom_outputfile(self, gaussian_oniom_outputfile):
        assert os.path.exists(gaussian_oniom_outputfile)
        g16_oniom = Gaussian16Output(filename=gaussian_oniom_outputfile)
        assert g16_oniom.normal_termination is False
        assert g16_oniom.oniom_cutting_bonds == {
            (49, 50): (0.700189, 0.700189),
            (80, 81): (0.700189, 0.700189),
            (176, 177): (0.700189, 0.700189),
            (198, 199): (0.700189, 0.700189),
            (217, 218): (0.700189, 0.700189),
            (439, 438): (0.700189, 0.700189),
        }
        assert g16_oniom.oniom_partition == {
            "high level atoms": [
                "50-60",
                "81-89",
                "177-186",
                "199-207",
                "218-225",
                "291-294",
                "308-312",
                "364-367",
                "375-379",
                "387-390",
                "421-438",
                "440-475",
            ],
            "low level atoms": [
                "1-49",
                "61-80",
                "90-176",
                "187-198",
                "208-217",
                "226-290",
                "295-307",
                "313-363",
                "368-374",
                "380-386",
                "391-420",
                "439",
                "476-483",
            ],
        }
        assert g16_oniom.oniom_get_charge_and_multiplicity == {
            "low-level, real system": (1, 2),
            "high-level, model system": (1, 1),
            "low-level, model system": (1, 1),
        }
        assert g16_oniom.oniom_layer_energies == {
            "method:  high, system:  model": -5303.002072980664,
            "method:  low, system:  model": 6.767438788151,
            "method:  low, system:  real": 9.234384095059,
        }
        assert g16_oniom.num_atoms == 483
        assert len(g16_oniom.oniom_energies) == 2
        assert g16_oniom.oniom_energies[0] == -5278.927903743607
        assert g16_oniom.oniom_energies[1] == -5300.535127673756
        assert g16_oniom.scf_energies[0] == -5303.01662026
        assert (
            g16_oniom.energies_in_eV[0] == -5278.927903743607 * units.Hartree
        )

    def test_normal_termination_semiempirical_pm6_output_file(
        self, gaussian_semiempirical_pm6_output_file
    ):
        g16_pm6 = Gaussian16Output(
            filename=gaussian_semiempirical_pm6_output_file
        )
        assert g16_pm6.normal_termination
        assert g16_pm6.molecule.num_atoms == 27
        assert g16_pm6.molecule.empirical_formula == "C9H16N2"
        assert g16_pm6.ab_initio is None
        assert g16_pm6.functional is None
        assert g16_pm6.basis is None
        assert g16_pm6.jobtype == "opt"
        assert g16_pm6.route_string == "# opt freq pm6"
        assert g16_pm6.freq
        assert (
            g16_pm6.semiempirical == "PM6"
        )  # changed to upper case in route_object.semiempirical


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
        assert g16_output.input_orientations is None
        assert g16_output.standard_orientations is not None
        assert len(g16_output.standard_orientations) == 1
        assert len(g16_output.all_structures) == 1

    def test_molecules(self):
        mol = Molecule.from_pubchem("241")  # benzene molecule
        assert mol.is_aromatic


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


class TestGaussianPBCOutputFile:
    def test_read_2d_pbc_output(self, gaussian_pbc_2d_outputfile):
        assert os.path.exists(gaussian_pbc_2d_outputfile)
        g16_pbc_2d = Gaussian16OutputWithPBC(
            filename=gaussian_pbc_2d_outputfile
        )
        assert g16_pbc_2d.normal_termination is False
        assert g16_pbc_2d.num_atoms == 2
        assert np.array_equal(g16_pbc_2d.pbc, [1, 1, 0])
        assert g16_pbc_2d.dim == 2
        assert np.allclose(
            g16_pbc_2d.input_translation_vectors,
            np.array(
                [
                    [2.475315, 0.0, 0.0],
                    [-1.219952, 2.133447, 0.0],
                ]
            ),
        )
        assert np.allclose(
            g16_pbc_2d.final_translation_vector,
            np.array([[2.47554, -0.0, -0.0], [-1.237852, 2.143856, 0.0]]),
        )
        assert len(g16_pbc_2d.energies) == 5
        assert g16_pbc_2d.energies[0] == -76.1487231466
        assert g16_pbc_2d.energies_in_eV[0] == -76.1487231466 * units.Hartree
        assert np.allclose(
            g16_pbc_2d.forces[-1],
            np.array(
                [
                    [1.5884e-05, 6.7630e-06, 0.0000e00],
                    [-1.5884e-05, -6.7630e-06, -0.0000e00],
                ]
            ),
        )

        assert np.allclose(
            g16_pbc_2d.forces_in_eV_per_angstrom[-1],
            np.array(
                [
                    [1.5884e-05, 6.7630e-06, 0.0000e00],
                    [-1.5884e-05, -6.7630e-06, -0.0000e00],
                ]
            )
            * units.Hartree
            / units.Bohr,
        )

        assert np.allclose(
            g16_pbc_2d.last_structure.positions,
            np.array(
                [[-0.001724, -0.714621, -0.0], [0.001724, 0.714621, 0.0]]
            ),
        )  # last structure that has energy and forces

        assert np.isclose(
            g16_pbc_2d.last_structure.energy,
            -76.1490641879,
            rtol=1e-5,
        )
        assert np.allclose(
            g16_pbc_2d.last_structure.forces,
            np.array(
                [
                    [0.000015884, 0.000006763, 0.000000000],
                    [-0.000015884, -0.000006763, -0.000000000],
                ]
            ),
            rtol=1e-5,
        )

        assert g16_pbc_2d.has_forces

        expected_first_pbc_forces = np.array(
            [
                [-0.005794968, -0.000018277, 0.000000000],
                [-0.015305998, 0.009167561, 0.000000000],
            ]
        )
        assert np.allclose(g16_pbc_2d.pbc_forces[0], expected_first_pbc_forces)

        expected_last_pbc_forces = np.array(
            [
                [0.000006901, 0.000014889, -0.000000000],
                [0.000028860, -0.000004081, -0.000000000],
            ]
        )
        assert np.allclose(g16_pbc_2d.pbc_forces[-1], expected_last_pbc_forces)

        expected_last_translation_vector = np.array(
            [
                [2.475540, -0.000000, -0.000000],
                [-1.237852, 2.143856, 0.000000],
            ]
        )
        assert np.allclose(
            g16_pbc_2d.input_orientations_pbc[-1],
            expected_last_translation_vector,
        )

        # this log file tests/data/GaussianTests/pbc/log/graphite_2d_opt.log
        # has only "Input orientation:" but no "Standard orientation:"
        assert g16_pbc_2d.standard_orientations is None
        assert g16_pbc_2d.standard_orientations_pbc is None
