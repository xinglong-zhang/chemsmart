import os.path

import numpy as np
import pytest
from ase import units
from ase.symbols import Symbols

from chemsmart.io.gaussian.cube import GaussianCubeFile
from chemsmart.io.gaussian.input import Gaussian16Input, Gaussian16QMMMInput
from chemsmart.io.gaussian.output import (
    Gaussian16Output,
    Gaussian16OutputWithPBC,
    Gaussian16pKaOutput,
    Gaussian16WBIOutput,
)
from chemsmart.io.gaussian.route import GaussianRoute
from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.constants import kcal_per_mol_to_hartree


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
        assert r1b.functional == "b3lyp-d3bj"
        assert r1b.basis == "6-311+g(d,p)"
        assert r1b.jobtype == "ts"
        assert r1b.solv is False
        assert r1b.dieze_tag is None
        assert (
            r1b.additional_opt_options_in_route is None
        )  # noeigentest prevents Gaussian from stopping
        # if no negative Hessian eigenvalue was found
        # (not additional opt options for geometry opt)
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
        assert r1e.functional == "tpsstpss-d3bj"
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
        s1g = "# TD(nstates=30) wB97XD/def2SVP scrf(solvent=dichloroethane)"
        # TD-DFT route
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

    def test_read_route_string_oniom_layer_methods_and_bases(self):
        s1qmmm = "# oniom(b3lyp/6-31g(d,p):uff) opt"
        r1qmmm = GaussianRoute(s1qmmm)
        assert r1qmmm.functional == "b3lyp:uff"
        assert r1qmmm.basis == "6-31g(d,p):none"
        s2qmmm = "# oniom(mp2/6-31g:hf/6-31g:am1) geom=connectivity"
        r2qmmm = GaussianRoute(s2qmmm)
        assert r2qmmm.functional == "mp2:hf:am1"
        assert r2qmmm.basis == "6-31g:6-31g:none"
        s3qmmm = "# oniom(mp2/6-31g:hf/6-31g) geom=connectivity"
        r3qmmm = GaussianRoute(s3qmmm)
        assert r3qmmm.functional == "mp2:hf"
        assert r3qmmm.basis == "6-31g:6-31g"

    def test_read_route_string_oniom_semiempirical_high_layer(self):
        s1 = "# oniom(am1:uff) opt"
        r1 = GaussianRoute(s1)
        assert r1.semiempirical == "AM1:UFF"

    def test_read_route_string_oniom_high_layer_not_semiempirical(self):
        s1 = "# oniom(b3lyp:uff) opt"
        r1 = GaussianRoute(s1)
        assert r1.semiempirical is None

    def test_read_route_string_oniom_ab_initio_high_layer(self):
        s1 = "# oniom(mp2/6-31g:uff) opt"
        r1 = GaussianRoute(s1)
        assert r1.ab_initio == "mp2:uff"

    def test_read_route_string_oniom_missing_closing_paren(self):
        # No closing ")" anywhere: the char-scan never finds a match at
        # depth 0, so the oniom parser bails out via its for/else clause.
        s1 = "# oniom(b3lyp:uff opt"
        r1 = GaussianRoute(s1)
        assert r1._get_oniom_layer_methods_and_bases() == (None, None)

    def test_read_route_string_oniom_empty_parens(self):
        s1 = "# oniom() opt"
        r1 = GaussianRoute(s1)
        assert r1.functional is None
        assert r1.basis is None

    def test_read_route_string_oniom_skips_empty_layer(self):
        # A bare "/" between colons produces an empty method/basis split
        # for that middle layer, which must be skipped, not appended.
        s1 = "# oniom(am1:/:uff) opt"
        r1 = GaussianRoute(s1)
        assert r1.functional == "am1:uff"

    def test_jobtype_ircf(self):
        s1 = "# irc=forward b3lyp def2svp"
        r1 = GaussianRoute(s1)
        assert r1.jobtype == "ircf"

    def test_jobtype_ircr(self):
        s1 = "# irc=reverse b3lyp def2svp"
        r1 = GaussianRoute(s1)
        assert r1.jobtype == "ircr"

    def test_jobtype_nci(self):
        s1 = "# b3lyp def2svp output=wfn"
        r1 = GaussianRoute(s1)
        assert r1.jobtype == "nci"

    def test_jobtype_resp(self):
        s1 = "# hf 6-31g* pop=mk iop(6/33=2,6/41=10,6/42=17,6/50=1)"
        r1 = GaussianRoute(s1)
        assert r1.jobtype == "resp"

    def test_jobtype_link(self):
        s1 = "# b3lyp def2svp stable=opt"
        r1 = GaussianRoute(s1)
        assert r1.jobtype == "link"

    def test_three_part_functional_basis_fit(self):
        # A trailing token after the 3-part func/basis/fit spec is needed
        # so the parsing loop continues past it to a further iteration.
        s1 = "# opt tpsstpss/def2tzvp/fit nosymm"
        r1 = GaussianRoute(s1)
        assert r1.functional == "tpsstpss"
        assert r1.basis == "def2tzvp/fit"

    def test_slash_separated_token_with_unsupported_part_count_is_ignored(
        self,
    ):
        # Neither the 2-part nor 3-part func/basis split applies for a
        # token with 4 slash-separated parts; it's silently skipped.
        s1 = "# opt a/b/c/d nosymm"
        r1 = GaussianRoute(s1)
        assert r1.functional is None
        assert r1.basis is None

    def test_functional_leading_hash_is_stripped(self):
        s1 = "#mn15/def2svp opt"
        r1 = GaussianRoute(s1)
        assert r1.functional == "mn15"

    def test_dispersion_with_no_matching_suffix_key_leaves_functional_as_is(
        self,
    ):
        s1 = "# b3lyp def2svp empiricaldispersion=foo"
        r1 = GaussianRoute(s1)
        assert r1.functional == "b3lyp"

    def test_dispersion_merges_into_oniom_functional_first_layer_only(self):
        s1 = "# oniom(b3lyp:uff) empiricaldispersion=gd3bj opt"
        r1 = GaussianRoute(s1)
        assert r1.functional == "b3lyp-d3bj:uff"

    def test_solvent_id_none_when_scrf_has_no_solvent_keyword(self):
        s1 = "# opt b3lyp def2svp scrf=(pcm)"
        r1 = GaussianRoute(s1)
        assert r1.solvent_id is None

    def test_solvent_id_appends_read_when_read_precedes_solvent(self):
        s1 = "# opt b3lyp def2svp scrf=(cpcm,read,solvent=toluene)"
        r1 = GaussianRoute(s1)
        assert r1.solvent_id == "toluene,read"

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
        assert g16_link_opt.functional == "m062x"
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
        assert g16_link_ts.functional == "m062x"
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
        assert g16_link_sp.functional == "m062x"
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


class TestGaussian16InputDirectPropertyCoverage:
    """Direct coverage for Gaussian16Input/Gaussian16QMMMInput
    properties that TestGaussian16Input above doesn't reach directly
    (it mostly asserts on higher-level derived values like molecule
    formulas), plus three real bugs discovered along the way."""

    def test_num_content_groups(self, gaussian_opt_genecp_inputfile):
        g16 = Gaussian16Input(filename=gaussian_opt_genecp_inputfile)
        assert g16.num_content_groups == g16.num_content_blocks

    def test_modredundant_group_present_when_modred_section_exists(
        self, gaussian_modred_inputfile
    ):
        g16 = Gaussian16Input(filename=gaussian_modred_inputfile)
        assert g16.modredundant_group == ["B 2 12 F", "B 9 2 F"]

    def test_modredundant_group_none_without_modred(
        self, gaussian_opt_genecp_inputfile
    ):
        g16 = Gaussian16Input(filename=gaussian_opt_genecp_inputfile)
        assert g16.modredundant_group is None

    def test_is_pbc_true_for_translation_vector_input(
        self, gaussian_pbc_1d_inputfile
    ):
        g16 = Gaussian16Input(filename=gaussian_pbc_1d_inputfile)
        assert g16.is_pbc is True
        assert len(g16.translation_vectors) >= 1

    def test_is_pbc_false_for_non_pbc_input(
        self, gaussian_opt_genecp_inputfile
    ):
        g16 = Gaussian16Input(filename=gaussian_opt_genecp_inputfile)
        assert g16.is_pbc is False

    def test_mem_and_nproc(self, gaussian_opt_genecp_inputfile):
        g16 = Gaussian16Input(filename=gaussian_opt_genecp_inputfile)
        assert g16.mem is not None
        assert g16.nproc is not None

    def test_charge_and_multiplicity_fall_back_to_oniom_parsing(
        self, gaussian_qmmm_inputfile_3layer
    ):
        """When the charge/mult line has more than two numbers (as in a
        combined-layer ONIOM line), the base class's charge/multiplicity
        properties fall back to oniom parsing with use_partition=False
        (which never touches self.partition, unlike oniom_charge)."""
        g16 = Gaussian16Input(filename=gaussian_qmmm_inputfile_3layer)
        # the oniom fallback returns the raw parsed string, unlike the
        # normal charge/mult line path which returns an int
        assert g16.charge == "0"
        assert g16.multiplicity == "1"

    def test_oniom_charge_crashes_on_base_class_non_qmmm_input(
        self, gaussian_qmmm_inputfile_2layer
    ):
        """Documents BUGS_FOUND.md #44: oniom_charge/oniom_multiplicity
        call _get_oniom_charge_and_multiplicity with the default
        use_partition=True, which accesses self.partition -- an
        attribute that only exists on the Gaussian16QMMMInput subclass.
        On the plain base class this raises AttributeError, which the
        `except RecursionError` guard cannot catch."""
        g16 = Gaussian16Input(filename=gaussian_qmmm_inputfile_2layer)
        with pytest.raises(AttributeError, match="partition"):
            g16.oniom_charge
        with pytest.raises(AttributeError, match="partition"):
            g16.oniom_multiplicity

    def test_has_frozen_coordinates_and_indices(
        self, gaussian_frozen_opt_inputfile
    ):
        g16 = Gaussian16Input(filename=gaussian_frozen_opt_inputfile)
        assert g16.has_frozen_coordinates
        assert g16.frozen_coordinate_indices == list(range(1, 11))
        assert g16.free_coordinate_indices == [11, 12, 13, 14]

    def test_no_frozen_coordinates_returns_none(
        self, gaussian_opt_genecp_inputfile
    ):
        g16 = Gaussian16Input(filename=gaussian_opt_genecp_inputfile)
        assert not g16.has_frozen_coordinates
        assert g16.frozen_coordinate_indices is None
        assert g16.free_coordinate_indices is None

    def test_gen_genecp_group_and_light_heavy_elements_direct(
        self, gaussian_opt_genecp_inputfile
    ):
        g16 = Gaussian16Input(filename=gaussian_opt_genecp_inputfile)
        assert g16.gen_genecp_group is not None
        assert g16.light_elements == ["H", "C", "O"]
        assert g16.light_elements_basis == "def2svp"
        assert g16.heavy_elements == ["Pd"]
        assert g16.heavy_elements_basis == "def2-tzvppd"

    def test_genecp_derived_properties_none_without_gen_basis(
        self, gaussian_modred_inputfile
    ):
        g16 = Gaussian16Input(filename=gaussian_modred_inputfile)
        assert g16.genecp_section is None
        assert g16.light_elements is None
        assert g16.light_elements_basis is None
        assert g16.heavy_elements is None
        assert g16.heavy_elements_basis is None
        assert g16.custom_solvent is None
        assert g16.custom_solvent_group is None

    def test_gen_genecp_group_modred_and_solvent_present(
        self, modred_genecp_custom_solvent_inputfile
    ):
        g16 = Gaussian16Input(filename=modred_genecp_custom_solvent_inputfile)
        assert "modred" in g16.route_string
        assert "solvent=generic" in g16.route_string
        assert g16.gen_genecp_group == g16.content_groups[4:-1]
        assert g16.custom_solvent_group == g16.content_groups[-1]

    def test_gen_genecp_group_modred_only(self, modred_gen_inputfile):
        g16 = Gaussian16Input(filename=modred_gen_inputfile)
        assert "modred" in g16.route_string
        assert "solvent=generic" not in g16.route_string
        assert g16.gen_genecp_group == g16.content_groups[4:]

    def test_gen_genecp_group_solvent_only(self, tmp_path):
        """Neither of the existing genecp fixtures combines a custom
        solvent without modred, so a minimal synthetic input covers
        this branch of _get_gen_genecp_group."""
        path = tmp_path / "solvent_only_genecp.com"
        path.write_text(
            "\n".join(
                [
                    "%chk=t.chk",
                    "%mem=4GB",
                    "# opt mn15/genecp scrf=(smd,solvent=generic,read)",
                    "",
                    "title",
                    "",
                    "0 1",
                    "C 0.0 0.0 0.0",
                    "H 0.0 0.0 1.0",
                    "",
                    "C 0",
                    "def2svp",
                    "****",
                    "",
                    "stoichiometry=CH4",
                    "solventname=customsolvent",
                    "eps=10.0",
                    "",
                ]
            )
        )
        g16 = Gaussian16Input(filename=str(path))
        assert "modred" not in g16.route_string
        assert "solvent=generic" in g16.route_string
        assert g16.gen_genecp_group == g16.content_groups[3:-1]
        assert g16.custom_solvent is not None
        assert "solventname=customsolvent" in g16.custom_solvent

    def test_constrained_atoms_getter(self, gaussian_frozen_opt_inputfile):
        g16 = Gaussian16Input(filename=gaussian_frozen_opt_inputfile)
        assert g16.constrained_atoms == g16.coordinate_block.constrained_atoms

    def test_constrained_atoms_setter_recurses_infinitely(
        self, gaussian_opt_genecp_inputfile
    ):
        """Documents BUGS_FOUND.md #43: the setter reassigns
        self.constrained_atoms, recursing into itself forever instead
        of storing the value anywhere."""
        g16 = Gaussian16Input(filename=gaussian_opt_genecp_inputfile)
        with pytest.raises(RecursionError):
            g16.constrained_atoms = [1, 2]

    def test_qmmm_2layer_oniom_charge_and_real_charge(
        self, gaussian_qmmm_inputfile_2layer
    ):
        g16 = Gaussian16QMMMInput(filename=gaussian_qmmm_inputfile_2layer)
        assert g16.partition == {
            "high level atoms": ["2-5"],
            "low level atoms": ["6-9"],
        }
        assert g16.oniom_charge == {
            "charge_total": "0",
            "int_charge": "0",
        }
        assert g16.oniom_multiplicity == {"real_multiplicity": "1"}
        assert g16.real_charge == 0
        assert g16.int_charge == 0
        assert g16.real_multiplicity == 1

    def test_qmmm_3layer_partition_includes_medium_atoms(
        self, gaussian_qmmm_inputfile_3layer
    ):
        g16 = Gaussian16QMMMInput(filename=gaussian_qmmm_inputfile_3layer)
        assert "medium level atoms" in g16.partition

    def test_gen_genecp_group_none_for_semiempirical_without_basis(
        self, tmp_path
    ):
        path = tmp_path / "semiempirical_no_basis.com"
        path.write_text(
            "\n".join(
                [
                    "%chk=t.chk",
                    "%mem=4GB",
                    "# opt freq pm6",
                    "",
                    "title",
                    "",
                    "0 1",
                    "C 0.0 0.0 0.0",
                    "H 0.0 0.0 1.0",
                    "",
                ]
            )
        )
        g16 = Gaussian16Input(filename=str(path))
        assert g16.basis is None
        assert g16.gen_genecp_group is None

    def test_qmmm_2layer_model_charge_crashes_with_keyerror(
        self, gaussian_qmmm_inputfile_2layer
    ):
        """Documents BUGS_FOUND.md #45: for a 2-layer ONIOM system, the
        QMMM subclass's own _get_oniom_charge_and_multiplicity keeps
        charge_total/real_multiplicity/int_charge (not model_charge),
        so model_charge/model_multiplicity crash with KeyError."""
        g16 = Gaussian16QMMMInput(filename=gaussian_qmmm_inputfile_2layer)
        with pytest.raises(KeyError, match="model_charge"):
            g16.model_charge
        with pytest.raises(KeyError, match="model_multiplicity"):
            g16.model_multiplicity


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
        assert len(g16_output.contributions) == 50
        assert g16_output.transitions[0] == [
            "104A -> 108A",
            "105A -> 107A",
            "106A -> 107A",
            "106A -> 108A",
            "105B -> 106B",
            "106A <- 107A",
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
        assert g16_output.contributions[0] == [
            2.4,
            1.5,
            87.5,
            1.1,
            6.8,
            1.5,
        ]
        assert g16_output.contributions[-1] == [
            3.0,
            2.2,
            1.9,
            9.7,
            62.5,
            3.2,
        ]

        assert (
            g16_output.total_core_hours
            == g16_output.total_service_unit
            == 361.7
        )
        assert g16_output.total_elapsed_walltime == 6.4
        mol = g16_output.molecule
        assert not mol.has_vibrations

    def test_contribution_percentage_spin_scaling(self):
        output = type("Output", (), {})()
        output.contribution_coefficients = [[0.5, -0.3]]

        output.spin = "restricted"
        assert Gaussian16Output.contributions.func(output) == [[50.0, 18.0]]

        output.spin = "unrestricted"
        assert Gaussian16Output.contributions.func(output) == [[25.0, 9.0]]

        output.spin = None
        with pytest.raises(ValueError, match="Unknown spin type"):
            Gaussian16Output.contributions.func(output)

    def test_singlet_opt_output(self, gaussian_singlet_opt_outfile):
        assert os.path.exists(gaussian_singlet_opt_outfile)
        g16_output = Gaussian16Output(filename=gaussian_singlet_opt_outfile)
        assert g16_output.version == "G16RevB.01"
        assert g16_output.file_date == "2024-06-20 18:09:26"
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
        assert g16_output.temperature_in_K == 298.15
        assert g16_output.pressure_in_atm == 1.0
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
        assert g16_output.zero_point_energy == 0.284336
        assert np.isclose(
            g16_output.thermal_vibration_correction,
            190.931 * kcal_per_mol_to_hartree - 0.284336,
            atol=1e-6,
        )
        assert np.isclose(
            g16_output.thermal_rotation_correction,
            0.889 * kcal_per_mol_to_hartree,
            atol=1e-6,
        )
        assert np.isclose(
            g16_output.thermal_translation_correction,
            0.889 * kcal_per_mol_to_hartree,
            atol=1e-6,
        )
        assert g16_output.thermal_energy_correction == 0.307101
        assert g16_output.thermal_enthalpy_correction == 0.308045
        assert g16_output.thermal_gibbs_free_energy_correction == 0.225790
        assert g16_output.internal_energy == -1863.733079
        assert g16_output.enthalpy == -1863.732135
        assert g16_output.gibbs_free_energy == -1863.814390
        assert np.isclose(
            g16_output.electronic_entropy,
            0.000 * 1e-3 * kcal_per_mol_to_hartree,
            atol=1e-3,
        )
        assert np.isclose(
            g16_output.vibrational_entropy,
            90.556 * 1e-3 * kcal_per_mol_to_hartree,
            atol=1e-3,
        )
        assert np.isclose(
            g16_output.rotational_entropy,
            37.462 * 1e-3 * kcal_per_mol_to_hartree,
            atol=1e-3,
        )
        assert np.isclose(
            g16_output.translational_entropy,
            45.103 * 1e-3 * kcal_per_mol_to_hartree,
            atol=1e-3,
        )
        assert g16_output.has_dipole_moment
        assert np.allclose(
            g16_output.all_dipole_moments[-1],
            np.array([4.7915, -0.1097, 0.4554]),
            rtol=1e-4,
        )
        assert np.isclose(
            g16_output.all_dipole_moment_magnitudes[-1], 4.8143, rtol=1e-4
        )
        assert g16_output.all_point_groups[-1] == "C1"
        assert np.allclose(
            g16_output.all_rotational_constants()[-1],
            np.array([0.16245 * 1e9, 0.07382 * 1e9, 0.05332 * 1e9]),
            rtol=1e-4,
        )

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
        assert g16_output.has_dipole_moment
        assert g16_output.num_dipole_moments == 3
        assert np.allclose(
            g16_output.all_dipole_moments[-1],
            np.array([-1.6500, -5.4954, -2.3627]),
            rtol=1e-4,
        )
        assert np.isclose(
            g16_output.all_dipole_moment_magnitudes[-1], 6.2052, rtol=1e-4
        )
        assert g16_output.all_point_groups[-1] == "C1"
        assert np.allclose(
            g16_output.all_rotational_constants()[-1],
            np.array([0.06229 * 1e9, 0.05513 * 1e9, 0.04690 * 1e9]),
            rtol=1e-4,
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
        assert g16_link_opt.spin == "unrestricted"
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
        assert g16_link_ts.spin == "unrestricted"
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
        assert g16_link_modred.spin == "unrestricted"
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
        optimized_flags = [
            mol.is_optimized_structure for mol in g16_genecp.all_structures
        ]
        assert optimized_flags == [False] * 9 + [True, True]
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

    def test_read_full_gen_outputfile(self, gaussian_full_gen_outfile):
        assert os.path.exists(gaussian_full_gen_outfile)
        g16 = Gaussian16Output(filename=gaussian_full_gen_outfile)
        assert g16.normal_termination
        assert g16.gen_genecp == "gen"
        # Light elements use named basis
        assert g16.light_elements == ["H", "C"]
        assert g16.light_elements_basis == "6-31g(d)"
        # Heavy elements has explicit orbital basis
        assert g16.heavy_elements == ["Cl", "Br"]
        heavy_basis = g16.heavy_elements_basis

        cl_shells = heavy_basis["Cl"]
        assert [shell["shell"] for shell in cl_shells] == [
            "S",
            "S",
            "S",
            "S",
            "S",
            "S",
            "P",
            "P",
            "P",
            "P",
            "P",
        ]
        cl_first_shell = cl_shells[0]
        assert cl_first_shell["shell"] == "S"
        assert len(cl_first_shell["primitives"]) == 6
        cl_first_exp, cl_first_coeff = cl_first_shell["primitives"][0]
        assert np.isclose(cl_first_exp, 1.0581900000e05)
        assert np.isclose(cl_first_coeff, 7.3800000000e-04)
        cl_last_shell = cl_shells[-1]
        assert cl_last_shell["shell"] == "P"
        assert len(cl_last_shell["primitives"]) == 1
        cl_last_exp, cl_last_coeff = cl_last_shell["primitives"][0]
        assert np.isclose(cl_last_exp, 1.0943700000e-01)
        assert np.isclose(cl_last_coeff, 1.0000000000e00)

        br_shells = heavy_basis["Br"]
        assert [shell["shell"] for shell in br_shells] == [
            "S",
            "S",
            "S",
            "S",
            "S",
            "S",
            "S",
            "S",
            "P",
            "P",
            "P",
            "P",
            "P",
            "P",
            "P",
            "D",
            "D",
        ]
        br_first_shell = br_shells[0]
        assert br_first_shell["shell"] == "S"
        assert len(br_first_shell["primitives"]) == 6
        br_first_exp, br_first_coeff = br_first_shell["primitives"][0]
        assert np.isclose(br_first_exp, 4.3970000000e05)
        assert np.isclose(br_first_coeff, 8.1300000000e-04)
        br_last_shell = br_shells[-1]
        assert br_last_shell["shell"] == "D"
        assert len(br_last_shell["primitives"]) == 1
        br_last_exp, br_last_coeff = br_last_shell["primitives"][0]
        assert np.isclose(br_last_exp, 1.5350000000e00)
        assert np.isclose(br_last_coeff, 1.0000000000e00)

    def test_read_full_genecp_outputfile(self, gaussian_full_genecp_outfile):
        assert os.path.exists(gaussian_full_genecp_outfile)
        g16 = Gaussian16Output(filename=gaussian_full_genecp_outfile)
        assert g16.normal_termination
        assert g16.gen_genecp == "genecp"
        # Light element (Cl) uses named basis
        assert g16.light_elements == ["Cl"]
        assert g16.light_elements_basis == "def2svp"
        # Heavy element (Ag) has explicit orbital basis
        assert g16.heavy_elements == ["Ag"]

        ag_shells = g16.heavy_elements_basis["Ag"]
        assert [s["shell"] for s in ag_shells] == [
            "S",
            "S",
            "S",
            "S",
            "S",
            "S",
            "P",
            "P",
            "P",
            "P",
            "D",
            "D",
            "D",
            "F",
        ]
        ag_first_shell = ag_shells[0]
        assert ag_first_shell["shell"] == "S"
        assert len(ag_first_shell["primitives"]) == 2
        ag_first_exp, ag_first_coef = ag_first_shell["primitives"][0]
        assert np.isclose(ag_first_exp, 1.9000000000e01)
        assert np.isclose(ag_first_coef, -1.6600104141e-01)
        ag_last_shell = ag_shells[-1]
        assert ag_last_shell["shell"] == "F"
        assert len(ag_last_shell["primitives"]) == 1
        ag_last_exp, ag_last_coeff = ag_last_shell["primitives"][0]
        assert np.isclose(ag_last_exp, 1.3971100000e00)
        assert np.isclose(ag_last_coeff, 1.0000000000e00)

        ecp = g16.heavy_elements_ecp
        ag_ecp = ecp["Ag"]
        assert ag_ecp["n_valence_electrons"] == 19
        channel_names = [ch["name"] for ch in ag_ecp["channels"]]
        assert channel_names == ["F and up", "S - F", "P - F", "D - F"]
        ag_first_channel = ag_ecp["channels"][0]
        assert ag_first_channel["name"] == "F and up"
        assert len(ag_first_channel["terms"]) == 2
        r_pow, ag_first_exp, ag_first_coef, spin_orbit_coef = ag_first_channel[
            "terms"
        ][0]
        assert r_pow == 2
        assert np.isclose(ag_first_exp, 14.22)
        assert np.isclose(ag_first_coef, -33.68992012)
        assert np.isclose(spin_orbit_coef, 0.0)
        ag_last_channel = ag_ecp["channels"][-1]
        assert ag_last_channel["name"] == "D - F"
        assert len(ag_last_channel["terms"]) == 4
        r_pow, ag_last_exp, ag_last_coef, spin_orbit_coef = ag_last_channel[
            "terms"
        ][0]
        assert r_pow == 2
        assert np.isclose(ag_last_exp, 10.21)
        assert np.isclose(ag_last_coef, 73.71926087)
        assert np.isclose(spin_orbit_coef, 0.0)

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

        # since use_frozen is False, this is
        # not included in the output structure
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

    def test_normal_termination_with_trailing_blank_lines(
        self, gaussian_ts_genecp_outfile, tmp_path
    ):
        with open(gaussian_ts_genecp_outfile, "r") as f:
            contents = f.read()

        output_with_trailing_blanks = tmp_path / "pd_genecp_ts_trailing.log"
        with open(output_with_trailing_blanks, "w") as f:
            f.write(contents + "\n\n")

        g16_output = Gaussian16Output(
            filename=str(output_with_trailing_blanks)
        )
        assert g16_output.normal_termination

    def test_oldform_redundant_coordinates_atomic_numbers(self, tmp_path):
        outputfile = tmp_path / "old_form_numeric_coords.log"
        outputfile.write_text(
            "\n".join(
                [
                    " ----------------------------------------------------------------------",
                    " # opt b3lyp/gen",
                    " ----------------------------------------------------------------------",
                    ' Structure from the checkpoint file:  "Pd_insertion_ts_r.chk"',
                    " Charge =  0 Multiplicity = 1",
                    " Redundant internal coordinates found in file.  (old form).",
                    " 46.0,0,0.000000,0.000000,0.000000",
                    " H,0,0.000000,0.000000,1.000000",
                    " Recover connectivity data from disk.",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16_output = Gaussian16Output(filename=str(outputfile))
        assert g16_output.symbols == ["Pd", "H"]
        assert g16_output.all_structures == []
        assert g16_output.last_structure.chemical_symbols == ["Pd", "H"]
        assert g16_output.molecule.chemical_symbols == ["Pd", "H"]

    def test_pd_insertion_ts_r_logfile(
        self, gaussian_pd_insertion_ts_r_outfile
    ):
        g16_output = Gaussian16Output(
            filename=gaussian_pd_insertion_ts_r_outfile
        )
        assert g16_output.normal_termination
        assert g16_output.charge == 0
        assert g16_output.multiplicity == 1
        assert len(g16_output.symbols) == g16_output.molecule.num_atoms
        assert "Pd" in g16_output.symbols
        assert "Pd" in g16_output.molecule.chemical_symbols

    def test_energy_extraction_from_gaussian_output_file(
        self, gaussian_quintet_opt_outfile
    ):
        g16_out = Gaussian16Output(filename=gaussian_quintet_opt_outfile)
        h_from_file = g16_out.enthalpy
        assert np.isclose(h_from_file, -7521.416016, rtol=1e-4)
        g_from_file = g16_out.gibbs_free_energy
        assert np.isclose(g_from_file, -7521.548653, rtol=1e-4)

    def test_custom_solvent_smd_generic(self, gaussian_smd_generic_outfile):
        g16_generic = Gaussian16Output(filename=gaussian_smd_generic_outfile)
        assert g16_generic.normal_termination
        custom_solvent = g16_generic.custom_solvent
        assert (
            custom_solvent["SolventName"]
            == "1,1,1,3,3,3-HEXAFLUOROPROPAN-2-OL"
        )
        assert custom_solvent["Eps"] == 16.7
        assert custom_solvent["EpsInf"] == 1.625625
        assert custom_solvent["HbondAcidity"] == 0.77
        assert custom_solvent["HbondBasicity"] == 0.10
        assert custom_solvent["SurfaceTensionAtInterface"] == 23.23
        assert custom_solvent["CarbonAromaticity"] == 0.0
        assert custom_solvent["ElectronegativeHalogenicity"] == 0.60


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

    def test_molecules(self, gaussian_benzene_opt_outfile):
        mol = Molecule.from_filepath(
            gaussian_benzene_opt_outfile
        )  # benzene molecule
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


class TestGaussian16pKaOutput:
    """Tests for Gaussian16pKaOutput class for pKa thermochemistry calculations.

    Reference values are from .dat files in tests/data/GaussianTests/outputs/
    Generated at T=373.15K, c=1.0 mol/L, csg=100 cm^-1, ch=100 cm^-1

    5PQ_Me_ts1_no_pd_opt.dat values (HA - protonated acid):
        E = -345.741944 hartree
        ZPE = 0.133804 hartree
        H = -345.595097 hartree
        qh-H = -345.596472 hartree
        T.S = 0.053929 hartree
        T.qh-S = 0.052278 hartree
        G(T) = -345.649026 hartree
        qh-G(T) = -345.648751 hartree

    5PQ_Me_ts1_b_no_pd_opt.dat values (A- - conjugate base):
        E = -344.915399 hartree
        ZPE = 0.115987 hartree
        H = -344.786580 hartree
        qh-H = -344.787739 hartree
        T.S = 0.052926 hartree
        T.qh-S = 0.051766 hartree
        G(T) = -344.839506 hartree
        qh-G(T) = -344.839505 hartree
    """

    # Reference values from 5PQ_Me_ts1_no_pd_opt.dat at T=373.15K
    HA_E = -345.741944
    HA_ZPE = 0.133804
    HA_H = -345.595097
    HA_QH_H = -345.596472
    HA_TS = 0.053929
    HA_QH_TS = 0.052278
    HA_G = -345.649026
    HA_QH_G = -345.648751

    # Reference values from 5PQ_Me_ts1_b_no_pd_opt.dat at T=373.15K
    A_E = -344.915399
    A_ZPE = 0.115987
    A_H = -344.786580
    A_QH_H = -344.787739
    A_TS = 0.052926
    A_QH_TS = 0.051766
    A_G = -344.839506
    A_QH_G = -344.839505

    # Reference pKa for collidine (2,4,6-trimethylpyridine)
    PKA_COLLIDINE_REFERENCE = 6.75

    def test_init_with_default_settings(
        self, gaussian_pKa_HA_optimization_outputfile
    ):
        """Test initialization with default thermochemistry settings."""
        output = Gaussian16pKaOutput(
            filename=gaussian_pKa_HA_optimization_outputfile
        )
        assert output.filename == gaussian_pKa_HA_optimization_outputfile
        assert output.temperature == 298.15
        assert output.concentration == 1.0
        assert output.pressure == 1.0
        assert output.cutoff_entropy_grimme == 100.0
        assert output.cutoff_enthalpy == 100.0
        assert output.energy_units == "hartree"

    def test_init_with_custom_settings(
        self, gaussian_pKa_HA_optimization_outputfile
    ):
        """Test initialization with custom thermochemistry settings."""
        output = Gaussian16pKaOutput(
            filename=gaussian_pKa_HA_optimization_outputfile,
            temperature=373.15,
            concentration=1.0,
            pressure=1.0,
            cutoff_entropy_grimme=100.0,
            cutoff_enthalpy=100.0,
            energy_units="hartree",
        )
        assert output.temperature == 373.15
        assert output.concentration == 1.0
        assert output.cutoff_entropy_grimme == 100.0
        assert output.cutoff_enthalpy == 100.0
        assert output.energy_units == "hartree"

    def test_electronic_energy_ha(
        self, gaussian_pKa_HA_optimization_outputfile
    ):
        """Test electronic energy for HA matches reference value."""
        output = Gaussian16pKaOutput(
            filename=gaussian_pKa_HA_optimization_outputfile,
            temperature=373.15,
            concentration=1.0,
            cutoff_entropy_grimme=100.0,
            cutoff_enthalpy=100.0,
            energy_units="hartree",
        )
        E = output.electronic_energy_in_units
        assert np.isclose(E, self.HA_E, rtol=1e-6)

    def test_zero_point_energy_ha(
        self, gaussian_pKa_HA_optimization_outputfile
    ):
        """Test ZPE for HA matches reference value."""
        output = Gaussian16pKaOutput(
            filename=gaussian_pKa_HA_optimization_outputfile,
            temperature=373.15,
            concentration=1.0,
            cutoff_entropy_grimme=100.0,
            cutoff_enthalpy=100.0,
            energy_units="hartree",
        )
        zpe = output.zero_point_energy_in_units
        assert np.isclose(zpe, self.HA_ZPE, rtol=1e-4)

    def test_enthalpy_ha(self, gaussian_pKa_HA_optimization_outputfile):
        """Test enthalpy for HA matches reference value."""
        output = Gaussian16pKaOutput(
            filename=gaussian_pKa_HA_optimization_outputfile,
            temperature=373.15,
            concentration=1.0,
            cutoff_entropy_grimme=100.0,
            cutoff_enthalpy=100.0,
            energy_units="hartree",
        )
        H = output.enthalpy_in_units
        assert np.isclose(H, self.HA_H, rtol=1e-5)

    def test_qh_enthalpy_ha(self, gaussian_pKa_HA_optimization_outputfile):
        """Test qh-H for HA matches reference value."""
        output = Gaussian16pKaOutput(
            filename=gaussian_pKa_HA_optimization_outputfile,
            temperature=373.15,
            concentration=1.0,
            cutoff_entropy_grimme=100.0,
            cutoff_enthalpy=100.0,
            energy_units="hartree",
        )
        qh_H = output.qh_enthalpy_in_units
        assert np.isclose(qh_H, self.HA_QH_H, rtol=1e-5)

    def test_gibbs_free_energy_ha(
        self, gaussian_pKa_HA_optimization_outputfile
    ):
        """Test G(T) for HA matches reference value."""
        output = Gaussian16pKaOutput(
            filename=gaussian_pKa_HA_optimization_outputfile,
            temperature=373.15,
            concentration=1.0,
            cutoff_entropy_grimme=100.0,
            cutoff_enthalpy=100.0,
            energy_units="hartree",
        )
        G = output.gibbs_free_energy_in_units
        assert np.isclose(G, self.HA_G, rtol=1e-5)

    def test_qh_gibbs_free_energy_ha(
        self, gaussian_pKa_HA_optimization_outputfile
    ):
        """Test qh-G(T) for HA matches reference value."""
        output = Gaussian16pKaOutput(
            filename=gaussian_pKa_HA_optimization_outputfile,
            temperature=373.15,
            concentration=1.0,
            cutoff_entropy_grimme=100.0,
            cutoff_enthalpy=100.0,
            energy_units="hartree",
        )
        qh_G = output.qh_gibbs_free_energy
        assert np.isclose(qh_G, self.HA_QH_G, rtol=1e-5)

    def test_electronic_energy_a(self, gaussian_pKa_A_optimization_outputfile):
        """Test electronic energy for A- matches reference value."""
        output = Gaussian16pKaOutput(
            filename=gaussian_pKa_A_optimization_outputfile,
            temperature=373.15,
            concentration=1.0,
            cutoff_entropy_grimme=100.0,
            cutoff_enthalpy=100.0,
            energy_units="hartree",
        )
        E = output.electronic_energy_in_units
        assert np.isclose(E, self.A_E, rtol=1e-6)

    def test_qh_gibbs_free_energy_a(
        self, gaussian_pKa_A_optimization_outputfile
    ):
        """Test qh-G(T) for A- matches reference value."""
        output = Gaussian16pKaOutput(
            filename=gaussian_pKa_A_optimization_outputfile,
            temperature=373.15,
            concentration=1.0,
            cutoff_entropy_grimme=100.0,
            cutoff_enthalpy=100.0,
            energy_units="hartree",
        )
        qh_G = output.qh_gibbs_free_energy
        assert np.isclose(qh_G, self.A_QH_G, rtol=1e-5)

    def test_compute_thermochemistry_ha(
        self, gaussian_pKa_HA_optimization_outputfile
    ):
        """Test compute_thermochemistry returns all values for HA."""
        output = Gaussian16pKaOutput(
            filename=gaussian_pKa_HA_optimization_outputfile,
            temperature=373.15,
            concentration=1.0,
            cutoff_entropy_grimme=100.0,
            cutoff_enthalpy=100.0,
            energy_units="hartree",
        )
        result = output.compute_thermochemistry()

        # Check structure name
        assert result["structure"] == "5PQ_Me_ts1_no_pd_opt"

        # Check all values match reference
        assert np.isclose(result["electronic_energy"], self.HA_E, rtol=1e-6)
        assert np.isclose(result["zero_point_energy"], self.HA_ZPE, rtol=1e-4)
        assert np.isclose(result["enthalpy"], self.HA_H, rtol=1e-5)
        assert np.isclose(result["qh_enthalpy"], self.HA_QH_H, rtol=1e-5)
        assert np.isclose(result["gibbs_free_energy"], self.HA_G, rtol=1e-5)
        assert np.isclose(
            result["qh_gibbs_free_energy"], self.HA_QH_G, rtol=1e-5
        )

    def test_compute_pka_thermochemistry_ha_and_a(
        self,
        gaussian_pKa_HA_optimization_outputfile,
        gaussian_pKa_A_optimization_outputfile,
    ):
        """Test compute_pka_thermochemistry with exact reference values."""
        results = Gaussian16pKaOutput.compute_pka_thermochemistry(
            ha_file=gaussian_pKa_HA_optimization_outputfile,
            a_file=gaussian_pKa_A_optimization_outputfile,
            temperature=373.15,
            concentration=1.0,
            cutoff_entropy_grimme=100.0,
            cutoff_enthalpy=100.0,
            energy_units="hartree",
        )

        # Check settings
        assert results["settings"]["temperature"] == 373.15
        assert results["settings"]["concentration"] == 1.0
        assert results["settings"]["cutoff_entropy_grimme"] == 100.0
        assert results["settings"]["cutoff_enthalpy"] == 100.0
        assert results["settings"]["energy_units"] == "hartree"

        # Check HA values
        assert results["HA"]["name"] == "HA"
        assert np.isclose(results["HA"]["E"], self.HA_E, rtol=1e-6)
        assert np.isclose(results["HA"]["qh_G"], self.HA_QH_G, rtol=1e-5)
        assert np.isclose(results["HA"]["ZPE"], self.HA_ZPE, rtol=1e-4)
        assert np.isclose(results["HA"]["H"], self.HA_H, rtol=1e-5)
        assert np.isclose(results["HA"]["qh_H"], self.HA_QH_H, rtol=1e-5)
        assert np.isclose(results["HA"]["G"], self.HA_G, rtol=1e-5)

        # Check A- values
        assert results["A"]["name"] == "A-"
        assert np.isclose(results["A"]["E"], self.A_E, rtol=1e-6)
        assert np.isclose(results["A"]["qh_G"], self.A_QH_G, rtol=1e-5)
        assert np.isclose(results["A"]["ZPE"], self.A_ZPE, rtol=1e-4)
        assert np.isclose(results["A"]["H"], self.A_H, rtol=1e-5)
        assert np.isclose(results["A"]["qh_H"], self.A_QH_H, rtol=1e-5)
        assert np.isclose(results["A"]["G"], self.A_G, rtol=1e-5)

    def test_deprotonation_energy_difference(
        self,
        gaussian_pKa_HA_optimization_outputfile,
        gaussian_pKa_A_optimization_outputfile,
    ):
        """Test that deprotonation energy difference is calculated correctly.

        ΔE = E(A-) - E(HA) should be positive (deprotonation is endothermic)
        Δqh-G = qh-G(A-) - qh-G(HA) should also be positive
        """
        results = Gaussian16pKaOutput.compute_pka_thermochemistry(
            ha_file=gaussian_pKa_HA_optimization_outputfile,
            a_file=gaussian_pKa_A_optimization_outputfile,
            temperature=373.15,
            concentration=1.0,
            cutoff_entropy_grimme=100.0,
            cutoff_enthalpy=100.0,
        )

        delta_E = results["A"]["E"] - results["HA"]["E"]
        delta_qh_G = results["A"]["qh_G"] - results["HA"]["qh_G"]

        # Expected values from .dat files
        expected_delta_E = self.A_E - self.HA_E
        expected_delta_qh_G = self.A_QH_G - self.HA_QH_G

        assert np.isclose(delta_E, expected_delta_E, rtol=1e-6)
        assert np.isclose(delta_qh_G, expected_delta_qh_G, rtol=1e-5)

        # Deprotonation should be endothermic (ΔE > 0)
        assert delta_E > 0
        assert delta_qh_G > 0

    def test_print_pka_thermochemistry_summary(
        self,
        gaussian_pKa_HA_optimization_outputfile,
        gaussian_pKa_A_optimization_outputfile,
        capsys,
    ):
        """Test thermochemistry output for pKa calculation shows correct values.

        This test verifies that individual thermochemistry values can be
        extracted from the output objects for HA and A- species.
        """
        # Get thermochemistry for HA and A-
        results = Gaussian16pKaOutput.compute_pka_thermochemistry(
            ha_file=gaussian_pKa_HA_optimization_outputfile,
            a_file=gaussian_pKa_A_optimization_outputfile,
            temperature=373.15,
            concentration=1.0,
            cutoff_entropy_grimme=100.0,
            cutoff_enthalpy=100.0,
        )

        # Verify HA values match reference
        assert np.isclose(results["HA"]["E"], self.HA_E, rtol=1e-6)
        assert np.isclose(results["HA"]["qh_G"], self.HA_QH_G, rtol=1e-5)

        # Verify A- values match reference
        assert np.isclose(results["A"]["E"], self.A_E, rtol=1e-6)
        assert np.isclose(results["A"]["qh_G"], self.A_QH_G, rtol=1e-5)

    def test_thermochemistry_property_caching(
        self, gaussian_pKa_HA_optimization_outputfile
    ):
        """Test that thermochemistry object is cached."""
        output = Gaussian16pKaOutput(
            filename=gaussian_pKa_HA_optimization_outputfile,
            temperature=373.15,
        )

        # Access thermochemistry twice
        thermo1 = output.thermochemistry
        thermo2 = output.thermochemistry

        # Should be the same object (cached)
        assert thermo1 is thermo2

    def test_energy_units_conversion_kcal_mol(
        self, gaussian_pKa_HA_optimization_outputfile
    ):
        """Test energy conversion to kcal/mol."""
        output_hartree = Gaussian16pKaOutput(
            filename=gaussian_pKa_HA_optimization_outputfile,
            temperature=373.15,
            energy_units="hartree",
        )
        output_kcal = Gaussian16pKaOutput(
            filename=gaussian_pKa_HA_optimization_outputfile,
            temperature=373.15,
            energy_units="kcal/mol",
        )

        E_hartree = output_hartree.electronic_energy_in_units
        E_kcal = output_kcal.electronic_energy_in_units

        # 1 hartree ≈ 627.5094740631 kcal/mol
        assert np.isclose(E_kcal / E_hartree, 627.5094740631, rtol=0.001)

    def test_energy_units_conversion_kj_mol(
        self, gaussian_pKa_HA_optimization_outputfile
    ):
        """Test energy conversion to kJ/mol."""
        output_hartree = Gaussian16pKaOutput(
            filename=gaussian_pKa_HA_optimization_outputfile,
            temperature=373.15,
            energy_units="hartree",
        )
        output_kj = Gaussian16pKaOutput(
            filename=gaussian_pKa_HA_optimization_outputfile,
            temperature=373.15,
            energy_units="kJ/mol",
        )

        E_hartree = output_hartree.electronic_energy_in_units
        E_kj = output_kj.electronic_energy_in_units

        # 1 hartree ≈ 2625.5002 kJ/mol
        assert np.isclose(E_kj / E_hartree, 2625.5002, rtol=0.001)

    # =========================================================================
    # pKa Calculation Tests - Dual-level Proton Exchange Scheme
    # =========================================================================

    def test_compute_pka(
        self,
        gaussian_pKa_HA_optimization_outputfile,
        gaussian_pKa_A_optimization_outputfile,
        gaussian_pKa_HB_optimization_outputfile,
        gaussian_pKa_B_optimization_outputfile,
        gaussian_pKa_HA_single_point_outputfile,
        gaussian_pKa_A_single_point_outputfile,
        gaussian_pKa_HB_single_point_outputfile,
        gaussian_pKa_B_single_point_outputfile,
    ):
        """Test pKa calculation using Dual-level Proton Exchange scheme.

        Uses 5PQ_Me_ts1 as target acid (HA/A-) and collidine as reference (HB/B-).
        Reference pKa of collidine = 6.75

        The dual-level approach uses:
        1. Gas-phase frequency calculations for thermal corrections (G_corr)
        2. Solvent single-point calculations for E_solv
        3. G_soln = E_solv + G_corr for solution free energy (in Hartree/au)
        4. Proton exchange scheme: HA + B⁻ → A⁻ + HB
        """
        result = Gaussian16pKaOutput.compute_pka(
            ha_gas_file=gaussian_pKa_HA_optimization_outputfile,
            a_gas_file=gaussian_pKa_A_optimization_outputfile,
            href_gas_file=gaussian_pKa_HB_optimization_outputfile,
            ref_gas_file=gaussian_pKa_B_optimization_outputfile,
            ha_solv_file=gaussian_pKa_HA_single_point_outputfile,
            a_solv_file=gaussian_pKa_A_single_point_outputfile,
            href_solv_file=gaussian_pKa_HB_single_point_outputfile,
            ref_solv_file=gaussian_pKa_B_single_point_outputfile,
            pka_reference=self.PKA_COLLIDINE_REFERENCE,
            temperature=373.15,
        )

        # Check that result contains expected keys
        assert "pKa" in result
        assert "pKa_reference" in result
        assert "delta_G_soln_kcal_mol" in result
        assert "delta_G_soln_au" in result
        assert "temperature" in result

        # Check solution free energies are present (in Hartree/au)
        assert "G_soln_HA_au" in result
        assert "G_soln_A_au" in result
        assert "G_soln_HRef_au" in result
        assert "G_soln_Ref_au" in result

        # Check solvent SP energies are present (in Hartree/au)
        assert "E_solv_HA_au" in result
        assert "E_solv_A_au" in result
        assert "E_solv_HRef_au" in result
        assert "E_solv_Ref_au" in result

        # Check thermal corrections are present (in Hartree/au)
        assert "G_corr_HA_au" in result
        assert "G_corr_A_au" in result
        assert "G_corr_HRef_au" in result
        assert "G_corr_Ref_au" in result

        # Check gas-phase electronic energies are present (in Hartree/au)
        assert "E_gas_HA_au" in result
        assert "E_gas_A_au" in result
        assert "E_gas_HRef_au" in result
        assert "E_gas_Ref_au" in result

        # Verify reference pKa is stored correctly
        assert result["pKa_reference"] == self.PKA_COLLIDINE_REFERENCE

        # Verify temperature is stored correctly
        assert result["temperature"] == 373.15

        assert np.isclose(result["pKa"], 52.7025859, rtol=1e-6)

    def test_compute_pka_direct_scheme(
        self,
        gaussian_pKa_HA_optimization_outputfile,
        gaussian_pKa_A_optimization_outputfile,
        gaussian_pKa_HA_single_point_outputfile,
        gaussian_pKa_A_single_point_outputfile,
    ):
        """Test direct dissociation via unified compute_pka(scheme='direct')."""
        delta_G_proton = -265.9
        temperature = 373.15
        result = Gaussian16pKaOutput.compute_pka(
            ha_gas_file=gaussian_pKa_HA_optimization_outputfile,
            a_gas_file=gaussian_pKa_A_optimization_outputfile,
            ha_solv_file=gaussian_pKa_HA_single_point_outputfile,
            a_solv_file=gaussian_pKa_A_single_point_outputfile,
            scheme="direct",
            delta_G_proton=delta_G_proton,
            temperature=temperature,
        )

        HARTREE_TO_KCAL = 627.5094740631
        G_soln_HA_kcal = result["G_soln_HA_au"] * HARTREE_TO_KCAL
        G_soln_A_kcal = result["G_soln_A_au"] * HARTREE_TO_KCAL
        expected_delta_G_diss = G_soln_A_kcal + delta_G_proton - G_soln_HA_kcal

        assert result["scheme"] == "direct"
        assert np.isclose(
            result["delta_G_diss_kcal_mol"], expected_delta_G_diss, rtol=1e-6
        )

    def test_compute_pka_energy_values(
        self,
        gaussian_pKa_HA_optimization_outputfile,
        gaussian_pKa_A_optimization_outputfile,
        gaussian_pKa_HB_optimization_outputfile,
        gaussian_pKa_B_optimization_outputfile,
        gaussian_pKa_HA_single_point_outputfile,
        gaussian_pKa_A_single_point_outputfile,
        gaussian_pKa_HB_single_point_outputfile,
        gaussian_pKa_B_single_point_outputfile,
    ):
        """Test that dual-level calculation uses correct energy values.

        All energies are in Hartree (au) except ΔG_soln which is also
        provided in kcal/mol for the pKa formula.

        Verifies:
        - E_solv values from solvent SP files (Hartree)
        - G_corr = qh-G(T) - E_gas from gas-phase files (Hartree)
        - G_soln = E_solv + G_corr (Hartree)
        - ΔG_soln in both au and kcal/mol
        """
        result = Gaussian16pKaOutput.compute_pka(
            ha_gas_file=gaussian_pKa_HA_optimization_outputfile,
            a_gas_file=gaussian_pKa_A_optimization_outputfile,
            href_gas_file=gaussian_pKa_HB_optimization_outputfile,
            ref_gas_file=gaussian_pKa_B_optimization_outputfile,
            ha_solv_file=gaussian_pKa_HA_single_point_outputfile,
            a_solv_file=gaussian_pKa_A_single_point_outputfile,
            href_solv_file=gaussian_pKa_HB_single_point_outputfile,
            ref_solv_file=gaussian_pKa_B_single_point_outputfile,
            pka_reference=self.PKA_COLLIDINE_REFERENCE,
            temperature=373.15,
        )

        HARTREE_TO_KCAL = 627.5094740631

        # Verify G_soln = E_solv + G_corr for each species (all in Hartree/au)
        for species in ["HA", "A", "HRef", "Ref"]:
            E_solv_au = result[f"E_solv_{species}_au"]
            G_corr_au = result[f"G_corr_{species}_au"]
            G_soln_au = result[f"G_soln_{species}_au"]
            expected_G_soln_au = E_solv_au + G_corr_au
            assert np.isclose(
                G_soln_au, expected_G_soln_au, rtol=1e-10
            ), f"G_soln_{species}_au mismatch: {G_soln_au} vs {expected_G_soln_au}"

        # Verify ΔG_soln in Hartree (au)
        # ΔG_soln = [G(A⁻)_soln + G(HRef)_soln] - [G(HA)_soln + G(Ref⁻)_soln]
        expected_delta_G_au = (
            result["G_soln_A_au"] + result["G_soln_HRef_au"]
        ) - (result["G_soln_HA_au"] + result["G_soln_Ref_au"])
        assert np.isclose(
            result["delta_G_soln_au"], expected_delta_G_au, rtol=1e-10
        )

        # Verify ΔG_soln conversion to kcal/mol
        expected_delta_G_kcal = expected_delta_G_au * HARTREE_TO_KCAL
        assert np.isclose(
            result["delta_G_soln_kcal_mol"], expected_delta_G_kcal, rtol=1e-6
        )

    def test_print_pka_summary(
        self,
        gaussian_pKa_HA_optimization_outputfile,
        gaussian_pKa_A_optimization_outputfile,
        gaussian_pKa_HB_optimization_outputfile,
        gaussian_pKa_B_optimization_outputfile,
        gaussian_pKa_HA_single_point_outputfile,
        gaussian_pKa_A_single_point_outputfile,
        gaussian_pKa_HB_single_point_outputfile,
        gaussian_pKa_B_single_point_outputfile,
        capsys,
    ):
        """Test that print_pka_summary outputs correct format.

        All energies should be displayed in Hartree (au) except ΔG_soln
        which is shown in both au and kcal/mol.
        """
        Gaussian16pKaOutput.print_pka_summary(
            ha_gas_file=gaussian_pKa_HA_optimization_outputfile,
            a_gas_file=gaussian_pKa_A_optimization_outputfile,
            href_gas_file=gaussian_pKa_HB_optimization_outputfile,
            ref_gas_file=gaussian_pKa_B_optimization_outputfile,
            ha_solv_file=gaussian_pKa_HA_single_point_outputfile,
            a_solv_file=gaussian_pKa_A_single_point_outputfile,
            href_solv_file=gaussian_pKa_HB_single_point_outputfile,
            ref_solv_file=gaussian_pKa_B_single_point_outputfile,
            pka_reference=self.PKA_COLLIDINE_REFERENCE,
            temperature=373.15,
        )

        captured = capsys.readouterr()
        output = captured.out

        # Check header
        assert "Dual-level Proton Exchange Scheme" in output
        assert "HA + Ref⁻ → A⁻ + HRef" in output

        # Check method description is present
        assert "G_corr = qh-G(T) - E_gas" in output
        assert "G_soln = E_solv + G_corr" in output

        # Check sections are present with correct units (au)
        assert "Gas-Phase Electronic Energies (E_gas, au)" in output
        assert "Thermal Corrections" in output
        assert "Solvent Single-Point Energies (E_solv, au)" in output
        assert (
            "Solution Free Energies (G_soln = E_solv + G_corr, au)" in output
        )

        # Check ΔG_soln is shown in both units
        assert "ΔG_soln" in output
        assert "kcal/mol" in output

        # Check computed pKa is displayed
        assert "Computed pKa(HA)" in output


class TestGaussian16OutputAdditionalCoverage:
    """Additional direct-property coverage for Gaussian16Output/subclasses,
    targeting edge cases and rarely-exercised branches not reached by the
    higher-level fixture-driven tests above: empty/blank files, missing
    thermochemistry sections (SP-only jobs), link-job structure assembly,
    ONIOM helpers, WBI helpers, PBC helpers, and pKa error paths. Real
    fixtures are used wherever a natural one exists; small synthetic
    outputs are used only for edge cases no real fixture covers."""

    def test_normal_termination_empty_file(self, tmp_path):
        outputfile = tmp_path / "empty.log"
        outputfile.write_text("")
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.contents == []
        assert g16.normal_termination is False

    def test_normal_termination_blank_lines_only(self, tmp_path):
        outputfile = tmp_path / "blank_only.log"
        outputfile.write_text("\n\n\n")
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.contents == ["", "", ""]
        assert g16.normal_termination is False

    def test_gen_genecp_none_for_semiempirical(
        self, gaussian_semiempirical_pm6_output_file
    ):
        g16 = Gaussian16Output(filename=gaussian_semiempirical_pm6_output_file)
        assert g16.basis is None
        assert g16.gen_genecp is None

    def test_genecp_info_gen_route_but_no_basis_block_found(self, tmp_path):
        """gen_genecp is not None (route says genecp) but the output never
        actually prints a 'General basis read from cards:' block, so the
        parsing loop in _genecp_info runs to completion without matching."""
        outputfile = tmp_path / "no_basis_block.log"
        outputfile.write_text(
            "\n".join(
                [
                    " ----------------------------------------------------------------------",
                    " # opt mn15/genecp",
                    " ----------------------------------------------------------------------",
                    " Symbolic Z-matrix:",
                    " Charge =  0 Multiplicity = 1",
                    " C                     0.0000    0.0000    0.0000",
                    " H                     0.0000    0.0000    1.0000",
                    "",
                    " NAtoms=      2 NQM=        2 NQMF=       0",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.gen_genecp == "genecp"
        assert g16.heavy_elements is None
        assert g16.heavy_elements_basis is None
        assert g16.light_elements is None
        assert g16.light_elements_basis is None
        assert g16.heavy_elements_ecp is None

    def test_genecp_info_empty_symbols_raises_and_is_caught(self, tmp_path):
        """gen_genecp is not None but the Symbolic Z-matrix coordinate
        block is empty, so self.symbols raises ValueError (no symbols
        found), which _genecp_info's try/except catches, bailing out
        early with the all-None defaults."""
        outputfile = tmp_path / "no_symbols.log"
        outputfile.write_text(
            "\n".join(
                [
                    " ----------------------------------------------------------------------",
                    " # opt mn15/genecp",
                    " ----------------------------------------------------------------------",
                    " Symbolic Z-matrix:",
                    " Charge =  0 Multiplicity = 1",
                    "",
                    " General basis read from cards:  (5D, 7F)",
                    " Centers:       1",
                    " def2svp",
                    " ****",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        with pytest.raises(ValueError, match="No symbols found"):
            g16.symbols
        assert g16.heavy_elements is None
        assert g16.light_elements is None

    def test_genecp_info_multi_center_and_blank_line_handling(self, tmp_path):
        """Exercises: blank line preceding a 'Centers:' continuation line,
        multiple center numbers on a single light-element block (loop runs
        more than once), a second light-element block reusing the already
        -set light_elements_basis, and a heavy-element block with more
        than one center number."""
        outputfile = tmp_path / "genecp_multi_center.log"
        outputfile.write_text(
            "\n".join(
                [
                    " ----------------------------------------------------------------------",
                    " # opt mn15/genecp",
                    " ----------------------------------------------------------------------",
                    " Symbolic Z-matrix:",
                    " Charge =  0 Multiplicity = 1",
                    " C                     0.0000    0.0000    0.0000",
                    " H                     0.0000    0.0000    1.0000",
                    " O                     0.0000    0.0000    2.0000",
                    " N                     0.0000    0.0000    3.0000",
                    " Cl                    0.0000    0.0000    4.0000",
                    " Br                    0.0000    0.0000    5.0000",
                    "",
                    " NAtoms=      6 NQM=        6 NQMF=       0",
                    " General basis read from cards:  (5D, 7F)",
                    " Centers:       1      2",
                    " def2svp",
                    " ****",
                    "",
                    " Centers:       3      4",
                    " def2svp",
                    " ****",
                    " Centers:       5      6",
                    " S   1 1.00",
                    "     Exponent=  1.0000000000D+01 Coefficients=  1.0000000000D+00",
                    " ****",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert set(g16.light_elements) == {"N", "O", "C", "H"}
        assert g16.light_elements_basis == "def2svp"
        assert set(g16.heavy_elements) == {"Br", "Cl"}
        assert "Cl" in g16.heavy_elements_basis
        assert "Br" in g16.heavy_elements_basis

    def test_genecp_info_centers_line_with_nondigit_token(self, tmp_path):
        """A 'Centers:' line with a non-numeric token is tolerated: the
        token is silently skipped (ValueError caught) while the digit
        tokens are still parsed."""
        outputfile = tmp_path / "genecp_bad_center_token.log"
        outputfile.write_text(
            "\n".join(
                [
                    " ----------------------------------------------------------------------",
                    " # opt mn15/genecp",
                    " ----------------------------------------------------------------------",
                    " Symbolic Z-matrix:",
                    " Charge =  0 Multiplicity = 1",
                    " C                     0.0000    0.0000    0.0000",
                    "",
                    " NAtoms=      1 NQM=        1 NQMF=       0",
                    " General basis read from cards:  (5D, 7F)",
                    " Centers:     1 x",
                    " def2svp",
                    " ****",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.light_elements == ["C"]
        assert g16.light_elements_basis == "def2svp"

    def test_genecp_info_centers_line_at_end_of_file(self, tmp_path):
        """A 'Centers:' line with nothing following it before EOF exercises
        the 'if j >= len(self.contents): break' guard."""
        outputfile = tmp_path / "genecp_centers_eof.log"
        outputfile.write_text(
            "\n".join(
                [
                    " ----------------------------------------------------------------------",
                    " # opt mn15/genecp",
                    " ----------------------------------------------------------------------",
                    " Symbolic Z-matrix:",
                    " Charge =  0 Multiplicity = 1",
                    " C                     0.0000    0.0000    0.0000",
                    "",
                    " NAtoms=      1 NQM=        1 NQMF=       0",
                    " General basis read from cards:  (5D, 7F)",
                    " Centers:     1",
                ]
            )
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.light_elements is None
        assert g16.heavy_elements is None

    def test_custom_solvent_none_when_no_marker(
        self, gaussian_singlet_opt_outfile
    ):
        g16 = Gaussian16Output(filename=gaussian_singlet_opt_outfile)
        assert g16.custom_solvent is None

    def test_num_steps_none_without_scan(self, gaussian_link_sp_outputfile):
        g16 = Gaussian16Output(filename=gaussian_link_sp_outputfile)
        assert g16.num_steps is None

    def test_thermochemistry_none_fields_for_sp_job(
        self, gaussian_link_sp_outputfile
    ):
        """An SP-only link output has no frequency/thermochemistry
        section at all, so all of these thermal/entropy correction
        properties fall through their loops to the implicit/explicit
        None return."""
        g16 = Gaussian16Output(filename=gaussian_link_sp_outputfile)
        assert g16.zero_point_energy is None
        assert g16.thermal_vibration_correction is None
        assert g16.thermal_rotation_correction is None
        assert g16.thermal_translation_correction is None
        assert g16.thermal_energy_correction is None
        assert g16.thermal_enthalpy_correction is None
        assert g16.thermal_gibbs_free_energy_correction is None
        assert g16.internal_energy is None
        assert g16.enthalpy is None
        assert g16.electronic_entropy_no_temperature_in_SI is None
        assert g16.electronic_entropy is None
        assert g16.vibrational_entropy_no_temperature_in_SI is None
        assert g16.vibrational_entropy is None
        assert g16.rotational_entropy_no_temperature_in_SI is None
        assert g16.rotational_entropy is None
        assert g16.translational_entropy_no_temperature_in_SI is None
        assert g16.translational_entropy is None
        assert g16.entropy_in_J_per_mol_per_K is None
        assert g16.entropy is None
        assert g16.entropy_times_temperature is None
        assert g16.gibbs_free_energy is None
        assert g16.convergence_criterion_not_met is False
        assert g16.has_forces is False
        assert g16.forces is None
        assert g16.temperature_in_K is None
        assert g16.rotational_symmetry_number is None
        assert g16.service_units_by_jobs == g16.cpu_runtime_by_jobs_core_hours

    def test_spin_none_when_no_scf_done_line(self, tmp_path):
        outputfile = tmp_path / "no_scf_done.log"
        outputfile.write_text(
            "\n".join(
                [
                    " ----------------------------------------------------------------------",
                    " # opt mn15/def2svp",
                    " ----------------------------------------------------------------------",
                    " Symbolic Z-matrix:",
                    " Charge =  0 Multiplicity = 1",
                    " C                     0.0000    0.0000    0.0000",
                    "",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.spin is None

    def test_parse_explicit_basis_block_empty_returns_no_shells(self):
        assert Gaussian16Output._parse_explicit_basis_block([]) == []

    def test_parse_explicit_basis_block_primitive_before_any_shell_header(
        self,
    ):
        """A primitive line appearing before any shell header is ignored
        (current_shell is still None), so no shell dict is produced."""
        block_lines = [
            "    Exponent=  1.0000000000D+01 Coefficients=  1.0000000000D+00",
        ]
        assert Gaussian16Output._parse_explicit_basis_block(block_lines) == []

    def test_parse_pseudopotential_section_only_reachable_via_genecp_info(
        self, gaussian_full_genecp_outfile
    ):
        """_parse_pseudopotential_section is only ever invoked from
        _genecp_info, after self.symbols has already succeeded and been
        cached, so calling it directly still exercises the normal
        (non-exception) path."""
        g16 = Gaussian16Output(filename=gaussian_full_genecp_outfile)
        result = g16._parse_pseudopotential_section()
        assert "Ag" in result

    def test_spin_none_for_method_without_r_or_u_prefix(self, tmp_path):
        """Some composite/theory labels (e.g. printed for CBS-type or
        other composite methods) do not begin with 'R' or 'U', so spin
        falls through to the else branch and returns None."""
        outputfile = tmp_path / "no_ru_spin.log"
        outputfile.write_text(
            "\n".join(
                [
                    " ----------------------------------------------------------------------",
                    " # cbs-qb3",
                    " ----------------------------------------------------------------------",
                    " SCF Done:  E(CBS-QB3) =  -1.234567890     A.U. after   10 cycles",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.spin is None

    def test_custom_solvent_marker_present_but_no_solvent_line(self, tmp_path):
        """The non-standard PCM marker is present, but no 'Solvent...:'
        line ever follows, so the parsing loop runs to completion without
        ever setting `inside = True`, and params stays empty -> None."""
        outputfile = tmp_path / "custom_solvent_no_name.log"
        outputfile.write_text(
            "\n".join(
                [
                    " Using the following non-standard input for PCM:",
                    " Some other unrelated line.",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.custom_solvent is None

    def test_input_coordinates_block_no_markers_present(self, tmp_path):
        outputfile = tmp_path / "no_coord_markers.log"
        outputfile.write_text(
            "\n".join(
                [
                    " Just some text.",
                    " Nothing relevant here.",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.input_coordinates_block.coordinate_block == []

    def test_input_coordinates_block_symbolic_zmatrix_runs_to_eof(
        self, tmp_path
    ):
        """No trailing blank line after the coordinates: the inner loop
        exhausts self.contents[i+2:] normally instead of breaking on a
        blank line."""
        outputfile = tmp_path / "zmatrix_eof.log"
        outputfile.write_text(
            "\n".join(
                [
                    " Symbolic Z-matrix:",
                    " Charge =  0 Multiplicity = 1",
                    " C                     0.0000    0.0000    0.0000",
                    " H                     0.0000    0.0000    1.0000",
                ]
            )
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.symbols == ["C", "H"]

    def test_input_coordinates_block_symbolic_zmatrix_skips_extra_charge_line(
        self, tmp_path
    ):
        """A second 'Charge =' line inside the coordinate block (as seen
        in QM/MM output) is skipped rather than treated as an atom."""
        outputfile = tmp_path / "zmatrix_extra_charge.log"
        outputfile.write_text(
            "\n".join(
                [
                    " Symbolic Z-matrix:",
                    " Charge =  0 Multiplicity = 1",
                    " Charge =  0 Multiplicity = 1 for low level calculation on real system.",
                    " C                     0.0000    0.0000    0.0000",
                    " H                     0.0000    0.0000    1.0000",
                    "",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.symbols == ["C", "H"]

    def test_input_coordinates_block_first_zmatrix_block_empty_second_valid(
        self, tmp_path
    ):
        """The first 'Symbolic Z-matrix:' occurrence has nothing after it
        (blank line right away), so the outer loop must continue past it
        and pick up the coordinates from the second occurrence."""
        outputfile = tmp_path / "zmatrix_two_blocks.log"
        outputfile.write_text(
            "\n".join(
                [
                    " Symbolic Z-matrix:",
                    " Charge =  0 Multiplicity = 1",
                    "",
                    " Symbolic Z-matrix:",
                    " Charge =  0 Multiplicity = 1",
                    " C                     0.0000    0.0000    0.0000",
                    " H                     0.0000    0.0000    1.0000",
                    "",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.symbols == ["C", "H"]

    def test_input_coordinates_block_redundant_form_with_junk_and_charge(
        self, tmp_path
    ):
        """Covers several branches of the 'Redundant internal coordinates'
        old-form parsing path in a single file: a leading blank line
        (continue, not break, since nothing collected yet), a malformed
        non-numeric atom token, a non-integer atomic number token, a
        non-numeric coordinate token, a junk line, and a 'Charge =' line
        -- all skipped -- followed by valid old-form coordinate lines that
        run straight to EOF (no trailing blank line)."""
        outputfile = tmp_path / "redundant_old_form_junk.log"
        outputfile.write_text(
            "\n".join(
                [
                    " Redundant internal coordinates found in file.  (old form).",
                    "",
                    " not,a,valid,coordinate,line,at,all",
                    " abc,0,0.000000,0.000000,0.000000",
                    " 45.5,0,0.000000,0.000000,0.000000",
                    " 46.0,0,abc,0.000000,0.000000",
                    " Charge =  0 Multiplicity = 1",
                    " 46.0,0,0.000000,0.000000,0.000000",
                    " 1.0,0,0.000000,0.000000,1.000000",
                ]
            )
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.symbols == ["Pd", "H"]

    def test_input_coordinates_block_redundant_form_blank_after_data(
        self, tmp_path
    ):
        """A blank line appearing AFTER coordinates have already been
        collected terminates the inner loop via the break at line 566,
        distinct from the 'blank line before any data' continue case."""
        outputfile = tmp_path / "redundant_blank_after_data.log"
        outputfile.write_text(
            "\n".join(
                [
                    " Redundant internal coordinates found in file.  (old form).",
                    " 46.0,0,0.000000,0.000000,0.000000",
                    " 1.0,0,0.000000,0.000000,1.000000",
                    "",
                    " Recover connectivity data from disk.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.symbols == ["Pd", "H"]

    def test_thermal_corrections_none_when_component_keyword_missing(
        self, tmp_path
    ):
        """The 'E (Thermal) ... CV ...' header line is found (and
        zero_point_energy is available, satisfying
        thermal_vibration_correction's extra guard), but none of the
        Electronic/Vibrational/Rotational/Translational/Total component
        lines that should follow it are present, so each of these
        properties' inner search loop exhausts without a match and the
        outer loop simply keeps scanning (eventually returning None)."""
        outputfile = tmp_path / "thermal_missing_components.log"
        outputfile.write_text(
            "\n".join(
                [
                    " Zero-point correction=                          0.284336",
                    " E (Thermal)             CV                       S",
                    "                          KCal/Mol        Cal/Mol-Kelvin",
                    " Nothing relevant follows here at all.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.zero_point_energy == 0.284336
        assert g16.thermal_vibration_correction is None
        assert g16.thermal_rotation_correction is None
        assert g16.thermal_translation_correction is None
        assert g16.electronic_entropy_no_temperature_in_SI is None
        assert g16.vibrational_entropy_no_temperature_in_SI is None
        assert g16.rotational_entropy_no_temperature_in_SI is None
        assert g16.translational_entropy_no_temperature_in_SI is None
        assert g16.entropy_in_J_per_mol_per_K is None

    def test_hirshfeld_heavy_atoms_three_token_line_and_eof_no_terminator(
        self, tmp_path
    ):
        """Covers a heavy-atom Hirshfeld data line with neither 4 nor 5
        tokens (charge only, no CM5/spin -- the elif's False branch),
        and the block running to EOF without a 'Tot'/blank terminator."""
        outputfile = tmp_path / "hirshfeld_heavy_edge.log"
        outputfile.write_text(
            "\n".join(
                [
                    " Hirshfeld charges with hydrogens summed into heavy atoms:",
                    "       Q-H",
                    "  1 C    0.100000",
                ]
            )
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.hirshfeld_charges_heavy_atoms == {"C1": 0.1}

    def test_input_coordinates_block_redundant_form_never_succeeds(
        self, tmp_path
    ):
        """No valid old-form coordinate line ever appears, so the inner
        loop exhausts with an empty list, the outer break is skipped, and
        the outer loop continues (and ultimately exhausts too, since no
        further marker exists)."""
        outputfile = tmp_path / "redundant_old_form_all_junk.log"
        outputfile.write_text(
            "\n".join(
                [
                    " Redundant internal coordinates found in file.  (old form).",
                    " junk line one, not coordinates",
                    " junk line two, still not coordinates",
                ]
            )
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.input_coordinates_block.coordinate_block == []

    def test_num_atoms_charge_multiplicity_none_when_absent(self, tmp_path):
        outputfile = tmp_path / "no_natoms_charge.log"
        outputfile.write_text(
            "\n".join(
                [
                    " ----------------------------------------------------------------------",
                    " # opt mn15/def2svp",
                    " ----------------------------------------------------------------------",
                    " Just some unrelated text.",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.num_atoms is None
        assert g16.charge is None
        assert g16.multiplicity is None
        assert g16.num_basis_functions is None
        assert g16.num_primitive_gaussians is None
        assert g16.num_cartesian_basis_functions is None
        assert g16.all_dipole_moments == []
        assert g16.all_dipole_moment_magnitudes == []
        assert g16.has_dipole_moment is False
        assert g16.route_string == "# opt mn15/def2svp"
        assert g16.pressure_in_atm is None
        assert g16.mass is None
        assert g16._get_moments_of_inertia_and_principal_axes() is None

    def test_mulliken_and_hirshfeld_loop_exhaustion_and_none_returns(
        self, tmp_path
    ):
        """Covers the 'section header found but no terminating marker
        line before EOF' branches for both the plain and heavy-atom
        Mulliken parsers, plus the 'section never found at all' None
        -returning branch for the heavy-atom Mulliken parser."""
        outputfile = tmp_path / "mulliken_no_terminator.log"
        outputfile.write_text(
            "\n".join(
                [
                    " Mulliken charges:",
                    "               1",
                    "     1  C    0.100000",
                ]
            )
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.mulliken_atomic_charges == {"C1": 0.1}
        assert g16.mulliken_atomic_charges_heavy_atoms is None
        assert g16.mulliken_spin_densities_heavy_atoms is None

        outputfile2 = tmp_path / "mulliken_heavy_no_terminator.log"
        outputfile2.write_text(
            "\n".join(
                [
                    " Mulliken charges with hydrogens summed into heavy atoms:",
                    "               1",
                    "     1  C    0.100000",
                ]
            )
        )
        g16_heavy = Gaussian16Output(filename=str(outputfile2))
        assert g16_heavy.mulliken_atomic_charges_heavy_atoms == {"C1": 0.1}

    def test_hirshfeld_charges_raises_indexerror_when_section_absent(
        self, gaussian_singlet_opt_outfile
    ):
        """BUG: unlike hirshfeld_charges_heavy_atoms (which gracefully
        returns None when the Hirshfeld section is absent),
        _get_hirshfeld_charges_spins_dipoles_cm5 unconditionally indexes
        all_hirshfeld_charges[-1] etc. without checking for emptiness, so
        hirshfeld_charges/hirshfeld_spin_densities/hirshfeld_dipoles/
        hirshfeld_cm5_charges crash with IndexError instead of returning
        None for a file with no Hirshfeld analysis section at all."""
        g16 = Gaussian16Output(filename=gaussian_singlet_opt_outfile)
        with pytest.raises(IndexError):
            g16.hirshfeld_charges
        # the heavy-atom counterpart handles the same "absent" case
        # gracefully by returning None instead of crashing
        assert g16.hirshfeld_charges_heavy_atoms is None
        assert g16.hirshfeld_spin_densities_heavy_atoms is None
        assert g16.hirshfeld_cm5_charges_heavy_atoms is None

    def test_hirshfeld_cm5_charges_heavy_atoms_wrong_type_when_spin_present(
        self, gaussian_rc_hirshfeld_outfile
    ):
        """BUG: when both Hirshfeld charges and spin densities are
        present with hydrogens summed into heavy atoms (open-shell
        case), _get_hirshfeld_charges_spin_densities_cm5_charges_heavy_atoms
        returns the raw list `all_cm5_charges_heavy_atoms` for the CM5
        component instead of `all_cm5_charges_heavy_atoms[-1]` like the
        other two return values and like the closed-shell branch below
        it. hirshfeld_cm5_charges_heavy_atoms therefore returns a
        one-element list-of-dicts instead of a dict for any open-shell
        Hirshfeld calculation, unlike its closed-shell counterpart."""
        g16 = Gaussian16Output(filename=gaussian_rc_hirshfeld_outfile)
        # spin densities ARE present for this fixture (open-shell)
        assert g16.hirshfeld_spin_densities_heavy_atoms is not None
        result = g16.hirshfeld_cm5_charges_heavy_atoms
        assert isinstance(result, list)  # should be a dict, like the
        # closed-shell branch (see test_read_hirshfeld_charges_outputfile)
        assert isinstance(result[0], dict)

    def test_hirshfeld_no_terminator_before_eof(self, tmp_path):
        """The non-heavy-atom Hirshfeld block's inner loop runs to EOF
        without ever hitting the 'Tot' or blank-line terminator."""
        outputfile = tmp_path / "hirshfeld_no_terminator.log"
        outputfile.write_text(
            "\n".join(
                [
                    " Hirshfeld charges, spin densities, dipoles, and CM5 charges",
                    "       Q-H        Spin       Dipole X   Dipole Y   Dipole Z    Q-CM5",
                    "  1 C    0.100000   0.000000   0.010000   0.020000   0.030000   0.150000",
                ]
            )
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.hirshfeld_charges == {"C1": 0.1}

    def test_oniom_partition_alt_format_with_medium_layer(self, tmp_path):
        """CH3COOH_qmmm.log uses the 'alternative' ONIOM coordinate-block
        format (coordinates start 7 lines after 'Symbolic Z-matrix:'
        because multiple 'Charge =' lines are echoed for the different
        ONIOM sub-systems), and includes atoms in all three H/M/L
        layers, plus 'med'/'low' on 'mid' and 'med' on 'model' charge
        /multiplicity lines not exercised by the 2-layer ONIOM fixture
        used elsewhere."""
        g16 = Gaussian16Output(
            filename=os.path.join(
                "tests", "data", "GaussianTests", "outputs", "CH3COOH_qmmm.log"
            )
        )
        partition = g16.oniom_partition
        assert "high level atoms" in partition
        assert "medium level atoms" in partition
        assert "low level atoms" in partition

        charge_mult = g16.oniom_get_charge_and_multiplicity
        assert charge_mult["medium-level, mid system"] == (0, 1)
        assert charge_mult["low-level, mid system"] == (0, 1)
        assert charge_mult["medium-level, model system"] == (0, 1)

    def test_to_dataset_is_a_noop(self, gaussian_singlet_opt_outfile):
        g16 = Gaussian16Output(filename=gaussian_singlet_opt_outfile)
        assert g16.to_dataset() is None

    def test_moments_of_inertia_full_parse(
        self, gaussian_pKa_HA_optimization_outputfile
    ):
        """This real fixture happens to demonstrate a known Gaussian
        formatting quirk: when three eigenvalues run together with no
        separating whitespace ('229.315721660.916151828.89264'), the
        combined token fails float() parsing, so the code's own
        exception handler substitutes one 3-element inf array in place
        of the (unparseable) 3 separate eigenvalues -- a graceful,
        already-handled degradation, not a bug."""
        g16 = Gaussian16Output(
            filename=gaussian_pKa_HA_optimization_outputfile
        )
        moments, axes = g16._get_moments_of_inertia_and_principal_axes()
        assert len(moments) == 1
        assert np.all(np.isinf(moments[0]))
        assert axes.shape[0] == 3

    def test_moments_of_inertia_crashes_when_section_absent(
        self, gaussian_link_sp_outputfile
    ):
        """BUG: when 'Principal axes and moments of inertia' is not found
        at all (e.g. an SP-only job, which never prints that banner),
        _get_moments_of_inertia_and_principal_axes falls off the end of
        the function and implicitly returns a single None (not a
        (None, None) tuple). moments_of_inertia and
        moments_of_inertia_principal_axes both unconditionally unpack
        this return value as a 2-tuple, so both crash with TypeError
        instead of gracefully returning None like almost every other
        'section not found' property in this class does."""
        g16 = Gaussian16Output(filename=gaussian_link_sp_outputfile)
        assert "Principal axes and moments of inertia" not in "\n".join(
            g16.contents
        )
        with pytest.raises(TypeError):
            g16.moments_of_inertia
        with pytest.raises(TypeError):
            g16.moments_of_inertia_principal_axes

    def test_wbi_sections_run_to_eof_without_terminator(self, tmp_path):
        """Each of natural_atomic_orbitals, natural_population_analysis,
        and electronic_configuration searches for its own header line and
        then scans forward for a terminator line ('WARNING'/'Summary of
        Natural Population Analysis', a '===' divider, or 'Wiberg bond
        index matrix' respectively). When the relevant section is the
        last thing in the file, that inner loop exhausts self.contents
        without ever finding the terminator -- tested here one section
        per minimal file so the sections don't bleed into each other."""
        nao_file = tmp_path / "nao_no_terminator.log"
        nao_file.write_text(
            "\n".join(
                [
                    " NAO  Atom  No  lang   Type(AO)    Occupancy      Energy",
                    " ---------------------------------------------------",
                    "    1    Ni    1  S      Cor( 1S)     1.99858       -2.68937",
                ]
            )
        )
        g16_nao = Gaussian16WBIOutput(filename=str(nao_file))
        assert (
            g16_nao.natural_atomic_orbitals["Ni1"]["NAO_Ni1"]["occupancy"]
            == 1.99858
        )

        npa_file = tmp_path / "npa_no_terminator.log"
        npa_file.write_text(
            "\n".join(
                [
                    " Atom  No    Charge         Core      Valence    Rydberg      Total",
                    " ---------------------------------------------------",
                    " Ni     1     0.52827        10.0      15.0        1.0         27.47173",
                    " Ni     1     0.52827        10.0      15.0        1.0         27.47173",
                ]
            )
        )
        g16_npa = Gaussian16WBIOutput(filename=str(npa_file))
        assert g16_npa.natural_charges["Ni1"] == 0.52827

        econf_file = tmp_path / "econf_no_terminator.log"
        econf_file.write_text(
            "\n".join(
                [
                    " Natural Electron Configuration",
                    " ---------------------------------------------------",
                    " Ni    1     [core]4S(0.27)3d(8.70)4p(0.51)",
                ]
            )
        )
        g16_econf = Gaussian16WBIOutput(filename=str(econf_file))
        assert (
            g16_econf.electronic_configuration["Ni1"]
            == "[core]4S(0.27)3d(8.70)4p(0.51)"
        )

    def test_wbi_properties_none_or_empty_for_non_wbi_file(
        self, gaussian_singlet_opt_outfile
    ):
        g16 = Gaussian16WBIOutput(filename=gaussian_singlet_opt_outfile)
        assert g16.nbo_version is None
        assert g16.natural_atomic_orbitals == {}
        assert g16.natural_population_analysis == {}
        assert g16.natural_charges == {}
        assert g16.total_electrons == {}
        assert g16.electronic_configuration == {}

    def test_pbc_properties_none_for_non_pbc_file(
        self, gaussian_singlet_opt_outfile
    ):
        g16 = Gaussian16OutputWithPBC(filename=gaussian_singlet_opt_outfile)
        assert g16._parse("anything") is None
        assert g16.pbc is None
        assert g16.input_translation_vectors is None
        assert g16.final_translation_vector is None

    def test_pka_output_raises_valueerror_without_frequency_data(
        self, gaussian_pKa_HA_single_point_outputfile
    ):
        """The SP-only file has no frequency section, so
        Thermochemistry's derived quantities are all None, and each of
        these *_in_units properties (except electronic_energy_in_units,
        which doesn't need frequency data) raises a descriptive
        ValueError instead of silently returning None."""
        output = Gaussian16pKaOutput(
            filename=gaussian_pKa_HA_single_point_outputfile
        )
        # electronic energy doesn't require frequency data
        assert output.electronic_energy_in_units is not None
        with pytest.raises(ValueError, match="zero-point energy"):
            output.zero_point_energy_in_units
        with pytest.raises(ValueError, match="enthalpy"):
            output.enthalpy_in_units
        with pytest.raises(ValueError, match="qh-enthalpy"):
            output.qh_enthalpy_in_units
        with pytest.raises(ValueError, match="Gibbs free energy"):
            output.gibbs_free_energy_in_units
        with pytest.raises(ValueError, match="qh-Gibbs free energy"):
            output.qh_gibbs_free_energy

    def test_pka_output_thermochemical_properties_full(
        self, gaussian_pKa_HA_optimization_outputfile
    ):
        output = Gaussian16pKaOutput(
            filename=gaussian_pKa_HA_optimization_outputfile,
            temperature=373.15,
        )
        props = output.thermochemical_properties
        assert set(props.keys()) == {
            "electronic_energy",
            "zero_point_energy",
            "enthalpy",
            "qh_enthalpy",
            "gibbs_free_energy",
            "qh_gibbs_free_energy",
        }
        assert props["electronic_energy"] == pytest.approx(
            -345.741944, abs=1e-4
        )

    def test_route_string_none_when_no_hash_line(self, tmp_path):
        outputfile = tmp_path / "no_route.log"
        outputfile.write_text(
            "\n".join(
                [
                    " Just some unrelated text with no route line at all.",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.route_string is None

    def test_route_string_spanning_two_lines(self, tmp_path):
        outputfile = tmp_path / "route_two_lines.log"
        outputfile.write_text(
            "\n".join(
                [
                    " ----------------------------------------------------------------------",
                    " # opt freq mn15 def2svp scrf=(smd,solvent=generic,read)",
                    "  additional continued keyword",
                    " ----------------------------------------------------------------------",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert "additional continued keyword" in g16.route_string

    def test_route_string_spanning_three_lines(self, tmp_path):
        outputfile = tmp_path / "route_three_lines.log"
        outputfile.write_text(
            "\n".join(
                [
                    " ----------------------------------------------------------------------",
                    " # opt freq mn15 def2svp",
                    "  scrf=(smd,solvent=generic,read)",
                    "  additional continued keyword",
                    " ----------------------------------------------------------------------",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert "additional continued keyword" in g16.route_string
        assert "generic" in g16.route_string

    def test_modredundant_group_on_output_for_failed_modred_and_scan(
        self, gaussian_failed_modred_outfile, gaussian_failed_scan_outfile
    ):
        g16_modred = Gaussian16Output(filename=gaussian_failed_modred_outfile)
        assert g16_modred.modredundant_group is not None
        assert len(g16_modred.modredundant_group) > 0

        g16_scan = Gaussian16Output(filename=gaussian_failed_scan_outfile)
        assert g16_scan.modredundant_group is not None

    def test_frozen_and_free_coordinate_indices_none_without_frozen(
        self, gaussian_singlet_opt_outfile
    ):
        g16 = Gaussian16Output(filename=gaussian_singlet_opt_outfile)
        assert g16.has_frozen_coordinates is False
        assert g16.frozen_coordinate_indices is None
        assert g16.free_coordinate_indices is None
        assert g16.frozen_elements == []
        assert g16.free_elements == []

    def test_num_forces_and_optimized_structure_none_for_abnormal_termination(
        self, gaussian_ts_genecp_outfile, gaussian_failed_modred_outfile
    ):
        g16 = Gaussian16Output(filename=gaussian_ts_genecp_outfile)
        assert g16.num_forces == len(g16.forces)

        g16_failed = Gaussian16Output(filename=gaussian_failed_modred_outfile)
        assert not g16_failed.normal_termination
        assert g16_failed.optimized_structure is None

    def test_link_job_structure_assembly_sp_and_ts(
        self, gaussian_link_sp_outputfile, gaussian_link_ts_outputfile
    ):
        """Exercises the is_link branches of _get_all_molecular_structures:
        normal-termination SP link job (keep_last_only after drop_first),
        and abnormal-termination link job with multiple carried-over
        frames (drop_first then safe_min_lengths truncation)."""
        g16_sp = Gaussian16Output(filename=gaussian_link_sp_outputfile)
        assert g16_sp.is_link
        assert g16_sp.normal_termination
        assert g16_sp.jobtype == "sp"
        structures_sp = g16_sp.all_structures
        assert len(structures_sp) == 1

        g16_ts = Gaussian16Output(filename=gaussian_link_ts_outputfile)
        assert g16_ts.is_link
        assert not g16_ts.normal_termination
        structures_ts = g16_ts.all_structures
        assert len(structures_ts) >= 1

    def test_all_structures_no_mulliken_charges_attached(self, tmp_path):
        """When no Mulliken section is printed at all, the final structure
        does not get mulliken_atomic_charges/mulliken_spin_densities
        attached (the 'is not None' guards take their False branch)."""
        outputfile = tmp_path / "no_mulliken.log"
        outputfile.write_text(
            "\n".join(
                [
                    " ----------------------------------------------------------------------",
                    " # opt hf/sto-3g",
                    " ----------------------------------------------------------------------",
                    " Symbolic Z-matrix:",
                    " Charge =  0 Multiplicity = 1",
                    " C                     0.0000    0.0000    0.0000",
                    " H                     0.0000    0.0000    1.0000",
                    "",
                    " NAtoms=      2 NQM=        2 NQMF=       0",
                    "                         Standard orientation:",
                    " ---------------------------------------------------------------------",
                    " Center     Atomic      Atomic             Coordinates (Angstroms)",
                    " Number     Number       Type             X           Y           Z",
                    " ---------------------------------------------------------------------",
                    "      1          6           0        0.000000    0.000000    0.000000",
                    "      2          1           0        0.000000    0.000000    1.000000",
                    " ---------------------------------------------------------------------",
                    " SCF Done:  E(RHF) =  -38.0000000     A.U. after   10 cycles",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.mulliken_atomic_charges is None
        assert g16.mulliken_spin_densities is None
        structures = g16.all_structures
        assert len(structures) == 1
        assert not hasattr(structures[-1], "mulliken_atomic_charges") or (
            structures[-1].mulliken_atomic_charges is None
        )

    def test_absorptions_in_nm_and_oscillatory_strengths(self, td_outputfile):
        g16 = Gaussian16Output(filename=td_outputfile)
        assert len(g16.absorptions_in_nm) == 50
        assert g16.absorptions_in_nm[0] == 1601.13
        assert len(g16.oscillatory_strengths) == 50
        assert g16.oscillatory_strengths[0] == 0.0084

    def test_alpha_eigenvalues_none_when_absent(self, tmp_path):
        outputfile = tmp_path / "no_eigenvalues.log"
        outputfile.write_text(
            "\n".join(
                [
                    " ----------------------------------------------------------------------",
                    " # opt hf/sto-3g",
                    " ----------------------------------------------------------------------",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.alpha_occ_eigenvalues == []
        assert g16.alpha_virtual_eigenvalues is None

    def test_read_transitions_edge_cases(self, tmp_path):
        """Covers two edge branches of
        _read_transitions_and_contribution_coefficients: a non-blank,
        non-transition-matching line appearing before any transition line
        has been found for a state (falls through via plain increment),
        and an 'Excited State' header that is the very last line in the
        file (the inner while loop runs zero iterations)."""
        outputfile = tmp_path / "td_edge_cases.log"
        outputfile.write_text(
            "\n".join(
                [
                    "Excited State   1:  Singlet-A  1.0 eV  100.0 nm  f=0.1",
                    " This state for optimization and/or second-order correction.",
                    "   104A -> 108A        0.15573",
                    "",
                    "Excited State   2:  Singlet-A  2.0 eV  200.0 nm  f=0.2",
                ]
            )
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        transitions = g16.transitions
        assert transitions[0] == ["104A -> 108A"]
        assert transitions[1] == []

    def test_genecp_info_heavy_elements_none_when_gen_genecp_none(
        self, gaussian_semiempirical_pm6_output_file
    ):
        """When gen_genecp is None (e.g. semiempirical calc), _genecp_info
        returns its all-None defaults immediately without ever touching
        self.symbols."""
        g16 = Gaussian16Output(filename=gaussian_semiempirical_pm6_output_file)
        assert g16.gen_genecp is None
        assert g16.heavy_elements is None
        assert g16.heavy_elements_basis is None
        assert g16.heavy_elements_ecp is None
        assert g16.light_elements is None
        assert g16.light_elements_basis is None

    def test_genecp_info_not_atom_symbols_branch_is_dead_code(
        self, monkeypatch, tmp_path
    ):
        """BUG (#105): _genecp_info's `if not atom_symbols: return result`
        can never fire through any real call path, since self.symbols
        either raises (caught above) or returns a non-empty list. Forcing
        it via monkeypatch is the only way to reach it directly."""
        outputfile = tmp_path / "genecp_dead_branch.log"
        outputfile.write_text(
            "\n".join(
                [
                    " ----------------------------------------------------------------------",
                    " # opt mn15/genecp",
                    " ----------------------------------------------------------------------",
                    " Symbolic Z-matrix:",
                    " Charge =  0 Multiplicity = 1",
                    " C                     0.0000    0.0000    0.0000",
                    "",
                    " NAtoms=      1 NQM=        1 NQMF=       0",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        monkeypatch.setattr(type(g16), "symbols", property(lambda self: []))
        assert g16.heavy_elements is None
        assert g16.light_elements is None

    def test_genecp_info_general_basis_marker_at_eof(self, tmp_path):
        """'General basis read from cards:' is the very last line in the
        file, so the while loop's condition is False on its very first
        check (j >= len(self.contents) immediately, before the loop body
        ever runs)."""
        outputfile = tmp_path / "genecp_marker_eof.log"
        outputfile.write_text(
            "\n".join(
                [
                    " ----------------------------------------------------------------------",
                    " # opt mn15/genecp",
                    " ----------------------------------------------------------------------",
                    " Symbolic Z-matrix:",
                    " Charge =  0 Multiplicity = 1",
                    " C                     0.0000    0.0000    0.0000",
                    "",
                    " NAtoms=      1 NQM=        1 NQMF=       0",
                    " General basis read from cards:  (5D, 7F)",
                ]
            )
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.light_elements is None
        assert g16.heavy_elements is None

    def test_genecp_info_blank_line_between_centers_and_basis_name(
        self, tmp_path
    ):
        """A blank line between 'Centers:' and its basis-name content
        line exercises the 'skip blank lines' inner loop."""
        outputfile = tmp_path / "genecp_centers_blank.log"
        outputfile.write_text(
            "\n".join(
                [
                    " ----------------------------------------------------------------------",
                    " # opt mn15/genecp",
                    " ----------------------------------------------------------------------",
                    " Symbolic Z-matrix:",
                    " Charge =  0 Multiplicity = 1",
                    " C                     0.0000    0.0000    0.0000",
                    "",
                    " NAtoms=      1 NQM=        1 NQMF=       0",
                    " General basis read from cards:  (5D, 7F)",
                    " Centers:       1",
                    "",
                    " def2svp",
                    " ****",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.light_elements == ["C"]
        assert g16.light_elements_basis == "def2svp"

    def test_genecp_info_out_of_range_center_numbers(self, tmp_path):
        """A 'Centers:' line listing a center number outside
        1..len(atom_symbols) is silently skipped for both the heavy
        (explicit-orbital) and light (named-basis) branches, without
        crashing or being recorded."""
        outputfile = tmp_path / "genecp_out_of_range_centers.log"
        outputfile.write_text(
            "\n".join(
                [
                    " ----------------------------------------------------------------------",
                    " # opt mn15/genecp",
                    " ----------------------------------------------------------------------",
                    " Symbolic Z-matrix:",
                    " Charge =  0 Multiplicity = 1",
                    " C                     0.0000    0.0000    0.0000",
                    "",
                    " NAtoms=      1 NQM=        1 NQMF=       0",
                    " General basis read from cards:  (5D, 7F)",
                    " Centers:      99",
                    " S   1 1.00",
                    "     Exponent=  1.0000000000D+01 Coefficients=  1.0000000000D+00",
                    " ****",
                    " Centers:      99",
                    " def2svp",
                    " ****",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.heavy_elements is None
        assert g16.light_elements is None

    def test_parse_pseudopotential_section_returns_empty_when_symbols_fail(
        self, tmp_path
    ):
        """Calling _parse_pseudopotential_section directly on a file whose
        self.symbols raises exercises its own try/except (independent of
        _genecp_info's identical guard, which is never reached since this
        method is called directly here)."""
        outputfile = tmp_path / "no_symbols_for_ecp.log"
        outputfile.write_text(
            "\n".join(
                [
                    " Just some text with no coordinate block at all.",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        with pytest.raises(ValueError):
            g16.symbols
        assert g16._parse_pseudopotential_section() == {}

    def test_parse_pseudopotential_section_line_matching_term_shape_but_not_is_term(
        self, tmp_path
    ):
        """A line with exactly 4 tokens whose first token is a digit but
        whose second token has no '.' fails the stricter is_term token
        check yet still matches the looser ecp_term_pattern regex, so
        the final 'channel name' elif's `not term_re.match(line)` is
        False and the line falls through untouched back to the next
        loop iteration instead of being treated as a channel name."""
        outputfile = tmp_path / "ecp_channel_fallthrough.log"
        outputfile.write_text(
            "\n".join(
                [
                    " Symbolic Z-matrix:",
                    " Charge =  0 Multiplicity = 1",
                    " Ag                    0.0000    0.0000    0.0000",
                    "",
                    " NAtoms=      1 NQM=        1 NQMF=       0",
                    " Pseudopotential Parameters",
                    " ======================================================================",
                    " ======================================================================",
                    " Center     Number     Number of atoms",
                    " ----------------------------------------------------------------------",
                    "   1     19",
                    " 1 abc def ghi",
                    " ======================================================================",
                ]
            )
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16._parse_pseudopotential_section() == {}

    def test_link_job_drop_first_and_keep_last_only_falsy_branches(
        self, tmp_path
    ):
        """A minimal link-sp job (normal termination, jobtype 'sp') with
        two Standard orientation frames and a Forces block, but no SCF
        energies, no rotational constants, and no point group data,
        exercises the falsy (data-absent) branches of both drop_first()
        and keep_last_only() inside _get_all_molecular_structures --
        except for their `if orientations:`/`if orientations_pbc:`
        guards, which are dead code (see BUGS_FOUND.md #107): orientations
        is always non-empty when drop_first/keep_last_only run (guarded
        by their call sites), and orientations_pbc always mirrors
        orientations' length 1:1. Also covers keep_last_only's forces
        truthy branch, since Forces data (but not energies/rot_consts/
        point_groups) is present here."""
        outputfile = tmp_path / "link_sp_falsy_branches.log"
        outputfile.write_text(
            "\n".join(
                [
                    " ----------------------------------------------------------------------",
                    " # um062x def2tzvp stable=opt guess=mix",
                    " ----------------------------------------------------------------------",
                    " # um062x def2tzvp",
                    " ----------------------------------------------------------------------",
                    " Symbolic Z-matrix:",
                    " Charge =  0 Multiplicity = 1",
                    " C                     0.0000    0.0000    0.0000",
                    "",
                    " NAtoms=      1 NQM=        1 NQMF=       0",
                    "                         Standard orientation:",
                    " ---------------------------------------------------------------------",
                    " Center     Atomic      Atomic             Coordinates (Angstroms)",
                    " Number     Number       Type             X           Y           Z",
                    " ---------------------------------------------------------------------",
                    "      1          6           0        0.000000    0.000000    0.000000",
                    " ---------------------------------------------------------------------",
                    "                         Standard orientation:",
                    " ---------------------------------------------------------------------",
                    " Center     Atomic      Atomic             Coordinates (Angstroms)",
                    " Number     Number       Type             X           Y           Z",
                    " ---------------------------------------------------------------------",
                    "      1          6           0        0.000000    0.000000    1.000000",
                    " ---------------------------------------------------------------------",
                    " Center     Atomic                   Forces (Hartrees/Bohr)",
                    " Number     Number              X              Y              Z",
                    " -------------------------------------------------------------------",
                    "      1          6           0.000046905   -0.000110437   -0.000107477",
                    " -------------------------------------------------------------------",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.is_link
        assert g16.jobtype == "sp"
        assert g16.normal_termination
        assert g16.energies == []
        assert g16.all_rotational_constants(mode="physical") == []
        assert g16.all_point_groups == []
        structures = g16.all_structures
        assert len(structures) == 1
        assert structures[-1].positions.tolist() == [[0.0, 0.0, 1.0]]

    def test_forces_no_terminator_at_true_eof(self, tmp_path):
        """The Forces block's inner loop exhausts self.contents cleanly
        (no closing divider) when the block is the literal last content
        in the file, with nothing after it to trip up the parser."""
        outputfile = tmp_path / "forces_eof.log"
        outputfile.write_text(
            "\n".join(
                [
                    " Center     Atomic                   Forces (Hartrees/Bohr)",
                    " Number     Number              X              Y              Z",
                    " -------------------------------------------------------------------",
                    "      1          6           0.000046905   -0.000110437   -0.000107477",
                ]
            )
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        forces = g16.forces
        assert len(forces) == 1
        assert forces[0].shape == (1, 3)

    def test_forces_no_terminator_followed_by_trailing_content_crashes(
        self, tmp_path
    ):
        """BUG (#106): when a Forces block has no closing divider AND is
        followed by further non-blank content later in the file (e.g.
        the standard termination line), the parser keeps scanning and
        tries to parse that trailing content as force data, crashing
        with ValueError instead of stopping at the table's natural
        end."""
        outputfile = tmp_path / "forces_no_terminator_trailing.log"
        outputfile.write_text(
            "\n".join(
                [
                    " Center     Atomic                   Forces (Hartrees/Bohr)",
                    " Number     Number              X              Y              Z",
                    " -------------------------------------------------------------------",
                    "      1          6           0.000046905   -0.000110437   -0.000107477",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        with pytest.raises(ValueError):
            g16.forces

    def test_align_lengths_to_orientations_trims_longer_energies(
        self, tmp_path
    ):
        """Two identical Standard orientation frames get deduplicated
        down to one by clean_duplicate_structure, but the two SCF Done
        energies printed alongside them are untouched by dedup, so
        align_lengths_to_orientations must right-trim energies (now
        longer than the deduplicated orientations list) back down to
        match."""
        std_block = [
            "                         Standard orientation:",
            " ---------------------------------------------------------------------",
            " Center     Atomic      Atomic             Coordinates (Angstroms)",
            " Number     Number       Type             X           Y           Z",
            " ---------------------------------------------------------------------",
            "      1          6           0        0.000000    0.000000    0.000000",
            " ---------------------------------------------------------------------",
        ]
        lines = [
            " ----------------------------------------------------------------------",
            " # opt mn15/def2svp",
            " ----------------------------------------------------------------------",
            " Symbolic Z-matrix:",
            " Charge =  0 Multiplicity = 1",
            " C                     0.0000    0.0000    0.0000",
            "",
            " NAtoms=      1 NQM=        1 NQMF=       0",
        ]
        lines += std_block
        lines.append(
            " SCF Done:  E(RHF) =  -38.0000000     A.U. after   10 cycles"
        )
        lines += std_block
        lines.append(
            " SCF Done:  E(RHF) =  -38.0000001     A.U. after   10 cycles"
        )
        lines.append(
            " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023."
        )
        outputfile = tmp_path / "align_lengths_trim.log"
        outputfile.write_text("\n".join(lines) + "\n")
        g16 = Gaussian16Output(filename=str(outputfile))
        assert len(g16.standard_orientations) == 2
        assert g16.energies == [-38.0, -38.0000001]
        structures = g16.all_structures
        assert len(structures) == 1
        assert structures[0].energy == -38.0

    def test_include_intermediate_scan_multiple_optimized_indices(
        self, tmp_path
    ):
        """A synthetic multi-scan-point job with two fully-completed scan
        points (each with its own final optimized step) gives
        optimized_steps_indices more than one entry, so the
        is_optimized-tagging loop in _get_all_molecular_structures
        actually iterates more than once. A third scan point is recorded
        (via a 'Step number ... on scan point 3 out of 3' line) without a
        corresponding fourth orientation frame, so its mapped index (3)
        falls outside the valid range for the 3-frame is_optimized list,
        exercising the loop's `0 <= idx < len(is_optimized)` False branch
        (which loops back without setting anything) in addition to the
        True branch."""
        std_block_lines = [
            "                         Standard orientation:",
            " ---------------------------------------------------------------------",
            " Center     Atomic      Atomic             Coordinates (Angstroms)",
            " Number     Number       Type             X           Y           Z",
            " ---------------------------------------------------------------------",
            "      1          6           0        0.000000    0.000000    0.000000",
            " ---------------------------------------------------------------------",
        ]
        lines = [
            " ----------------------------------------------------------------------",
            " # opt modredundant mn15/def2svp",
            " ----------------------------------------------------------------------",
            " Symbolic Z-matrix:",
            " Charge =  0 Multiplicity = 1",
            " C                     0.0000    0.0000    0.0000",
            "",
            " NAtoms=      1 NQM=        1 NQMF=       0",
        ]
        lines += std_block_lines
        lines.append(
            " Step number   1 out of a maximum of  100 on scan point"
            "     1 out of     2"
        )
        lines += std_block_lines
        lines.append(
            " Step number   2 out of a maximum of  100 on scan point"
            "     1 out of     2"
        )
        lines += std_block_lines
        lines.append(
            " Step number   1 out of a maximum of  100 on scan point"
            "     2 out of     2"
        )
        lines.append(
            " Step number   1 out of a maximum of  100 on scan point"
            "     3 out of     3"
        )
        outputfile = tmp_path / "multi_scan_point.log"
        outputfile.write_text("\n".join(lines) + "\n")
        g16 = Gaussian16Output(
            filename=str(outputfile), include_intermediate=True
        )
        assert g16.optimized_steps_indices == [1, 2, 3]
        structures = g16.all_structures
        assert len(structures) == 3
        assert [s.is_optimized_structure for s in structures] == [
            False,
            True,
            True,
        ]

    def test_energies_uses_mp2_energies_when_present(self, tmp_path):
        """When EUMP2 lines are present, `energies` returns mp2_energies
        rather than falling back to scf_energies."""
        outputfile = tmp_path / "mp2_energies.log"
        outputfile.write_text(
            "\n".join(
                [
                    " SCF Done:  E(RHF) =  -76.0000000     A.U. after   10 cycles",
                    " EUMP2 =    -0.7635026712D+02",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.scf_energies == [-76.0]
        assert g16.mp2_energies == [-76.35026712]
        assert g16.energies == [-76.35026712]

    def test_modredundant_group_none_when_route_has_no_modred(self, tmp_path):
        outputfile = tmp_path / "no_modred_route.log"
        outputfile.write_text(
            "\n".join(
                [
                    " ----------------------------------------------------------------------",
                    " # opt mn15/def2svp",
                    " ----------------------------------------------------------------------",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.modredundant_group is None

    def test_modredundant_group_no_terminator_runs_to_eof(self, tmp_path):
        """The ModRedundant section is the last content in the file, so
        the inner blank-line-terminator search loop exhausts
        self.contents instead of hitting `break`."""
        outputfile = tmp_path / "modred_no_terminator.log"
        outputfile.write_text(
            "\n".join(
                [
                    " ----------------------------------------------------------------------",
                    " # opt modredundant mn15/def2svp",
                    " ----------------------------------------------------------------------",
                    " The following ModRedundant input section has been read:",
                    " B 1 2 F",
                ]
            )
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.modredundant_group == ["B 1 2 F"]

    def test_gen_genecp_none_for_non_gen_basis(self, tmp_path):
        """A basis set string that does not contain 'gen' (e.g. a
        standard Pople/Karlsruhe basis) makes _get_gen_genecp fall
        through to its final `return None`."""
        outputfile = tmp_path / "non_gen_basis.log"
        outputfile.write_text(
            "\n".join(
                [
                    " ----------------------------------------------------------------------",
                    " # opt mn15/def2svp",
                    " ----------------------------------------------------------------------",
                    " Normal termination of Gaussian 16 at Wed Nov  8 08:36:34 2023.",
                ]
            )
            + "\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.basis is not None
        assert "gen" not in g16.basis
        assert g16.gen_genecp is None

    def test_basis_function_counts_parsed_from_header_line(self, tmp_path):
        """A single real-format 'N basis functions, M primitive
        gaussians, K cartesian basis functions' line satisfies all
        three properties' search conditions at once."""
        outputfile = tmp_path / "basis_counts.log"
        outputfile.write_text(
            "    43 basis functions,    85 primitive gaussians,"
            "    46 cartesian basis functions\n"
        )
        g16 = Gaussian16Output(filename=str(outputfile))
        assert g16.num_basis_functions == 43
        assert g16.num_primitive_gaussians == 85
        assert g16.num_cartesian_basis_functions == 46

    def test_vibrational_modes_empty_when_frequencies_line_near_eof(
        self, tmp_path
    ):
        """A 'Frequencies --' line within the last 4 lines of the file
        makes the mode-row inner loop's iterable (self.contents[i+5:])
        empty, so it runs zero iterations instead of collecting rows or
        breaking on a mismatch."""
        outputfile = tmp_path / "vib_modes_near_eof.log"
        outputfile.write_text(" Frequencies --   123.4")
        g16 = Gaussian16Output(filename=str(outputfile))
        modes = g16.vibrational_modes
        assert len(modes) == 3
        assert all(m.size == 0 for m in modes)

    def test_entropy_in_j_per_mol_per_k_with_real_thermochemistry_data(
        self, gaussian_koh_linear_opt_outfile
    ):
        """A real fixture with a full 'E (Thermal) ... CV ... S' table
        (including its 'Total' row) exercises the successful-parse
        return path of entropy_in_J_per_mol_per_K, and by extension
        entropy and entropy_times_temperature."""
        g16 = Gaussian16Output(filename=gaussian_koh_linear_opt_outfile)
        assert g16.entropy_in_J_per_mol_per_K is not None
        assert g16.entropy is not None
        assert g16.entropy > 0
        if g16.temperature_in_K:
            assert g16.entropy_times_temperature is not None
