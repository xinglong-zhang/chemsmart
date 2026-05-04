import pytest

from chemsmart.io.gaussian.route import GaussianRoute
from chemsmart.io.molecules.structure import Molecule, QMMMMolecule
from chemsmart.jobs.gaussian.settings import (
    GaussianJobSettings,
    GaussianLinkJobSettings,
    GaussianQMMMJobSettings,
)
from chemsmart.jobs.settings import read_molecular_job_yaml


class TestGaussianJobSettings:
    def test_merge_dict(self):
        settings = GaussianJobSettings.default()
        settings.charge = None
        settings.multiplicity = None

        merged_settings = settings.merge({"charge": -1, "multiplicity": 3})
        assert merged_settings.charge == -1
        assert merged_settings.multiplicity == 3

    def test_merge_other_settings(self):
        settings1 = GaussianJobSettings.default()
        settings1.charge = None
        settings1.multiplicity = None

        settings2 = GaussianJobSettings.default()
        settings2.charge = 0
        settings2.multiplicity = 1

        merged_settings = settings1.merge(
            settings2, keywords=("charge", "multiplicity")
        )
        assert merged_settings.charge == 0
        assert merged_settings.multiplicity == 1

    def test_get_settings_from_yaml_gas_solv(
        self, gaussian_yaml_settings_gas_solv
    ):
        all_project_settings = read_molecular_job_yaml(
            gaussian_yaml_settings_gas_solv, program="gaussian"
        )
        opt_settings_dict = all_project_settings["opt"]
        opt_settings = GaussianJobSettings.from_dict(opt_settings_dict)
        assert isinstance(opt_settings, GaussianJobSettings)
        assert opt_settings.chk is True  # from defaults file
        assert opt_settings.functional == "b3lyp empiricaldispersion=gd3bj"
        assert opt_settings.basis == "def2svp"
        assert opt_settings.solvent_model is None
        assert opt_settings.solvent_id is None

        modred_settings_dict = all_project_settings["modred"]
        modred_settings = GaussianJobSettings.from_dict(modred_settings_dict)
        assert isinstance(modred_settings, GaussianJobSettings)
        assert modred_settings.chk is True  # from defaults file
        assert modred_settings.functional == "b3lyp empiricaldispersion=gd3bj"
        assert modred_settings.basis == "def2svp"
        assert modred_settings.solvent_model is None
        assert modred_settings.solvent_id is None

        sp_settings_dict = all_project_settings["sp"]
        sp_settings = GaussianJobSettings.from_dict(sp_settings_dict)
        assert isinstance(sp_settings, GaussianJobSettings)
        assert sp_settings.chk is True  # from defaults file
        assert sp_settings.functional == "b3lyp empiricaldispersion=gd3bj"
        assert sp_settings.basis == "def2tzvp"
        assert sp_settings.solvent_model == "smd"
        assert sp_settings.solvent_id == "toluene"

    def test_get_settings_from_yaml_solv(self, gaussian_yaml_settings_solv):
        all_project_settings = read_molecular_job_yaml(
            gaussian_yaml_settings_solv, program="gaussian"
        )
        opt_settings_dict = all_project_settings["opt"]
        opt_settings = GaussianJobSettings.from_dict(opt_settings_dict)
        assert isinstance(opt_settings, GaussianJobSettings)
        assert opt_settings.chk is True  # from defaults file
        assert opt_settings.functional == "m062x"
        assert opt_settings.basis == "def2tzvp"
        assert opt_settings.solvent_model == "smd"
        assert opt_settings.solvent_id == "toluene"

        modred_settings_dict = all_project_settings["modred"]
        modred_settings = GaussianJobSettings.from_dict(modred_settings_dict)
        assert isinstance(modred_settings, GaussianJobSettings)
        assert modred_settings.chk is True  # from defaults file
        assert modred_settings.functional == "m062x"
        assert modred_settings.basis == "def2tzvp"
        assert modred_settings.solvent_model == "smd"
        assert modred_settings.solvent_id == "toluene"

        sp_settings_dict = all_project_settings["sp"]
        sp_settings = GaussianJobSettings.from_dict(sp_settings_dict)
        assert isinstance(sp_settings, GaussianJobSettings)
        assert sp_settings.chk is True  # from defaults file
        assert sp_settings.functional == "m062x"
        assert sp_settings.basis == "def2tzvp"
        assert sp_settings.solvent_model == "smd"
        assert sp_settings.solvent_id == "toluene"

    def test_read_gaussian_hf_comfile(self, hf_com_filepath):
        settings = GaussianJobSettings.from_comfile(hf_com_filepath)
        assert isinstance(settings, GaussianJobSettings)
        assert settings.ab_initio == "hf"
        assert settings.functional is None
        assert settings.basis == "6-31g"
        assert settings.solvent_model is None
        assert settings.solvent_id is None

    def test_read_gaussian_settings_from_orca_inp(self, water_sp_input_path):
        settings = GaussianJobSettings.from_inpfile(water_sp_input_path)
        assert isinstance(settings, GaussianJobSettings)
        assert settings.ab_initio == "hf"
        assert settings.functional is None
        assert settings.basis == "def2-svp"
        assert settings.solvent_model is None
        assert settings.solvent_id is None


class TestGaussianQMMMJobSettings:
    def test_qmmm_settings(self):
        settings1 = GaussianQMMMJobSettings(
            high_level_functional="b3lyp",
            high_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_charge=0,
            real_multiplicity=1,
            high_level_atoms=[1, 2, 3],
            parent_jobtype="opt",
            freq=True,
        )
        assert settings1.route_string == "# opt freq oniom(b3lyp/6-31g(d):uff)"

        settings2 = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_charge=0,
            real_multiplicity=1,
            high_level_atoms=[1, 2, 3],
            parent_jobtype="sp",
        )
        assert (
            settings2.route_string
            == "# oniom(mn15/def2svp:b3lyp/6-31g(d):uff)"
        )

        settings3 = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            high_level_force_field="uff",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_charge=0,
            real_multiplicity=1,
            high_level_atoms=[1, 2, 3],
            parent_jobtype="sp",
        )
        # assert settings3.route_string == "#
        # oniom(mn15/def2svp:uff:b3lyp/6-31g(d):uff)"
        # ValueError: For high level of theory, one should
        # specify only functional/basis or force field!
        with pytest.raises(ValueError):
            settings3.route_string

        # settings with solvent specification
        settings4 = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_charge=0,
            real_multiplicity=1,
            high_level_atoms=[1, 2, 3],
            solvent_model="smd",
            solvent_id="toluene",
            parent_jobtype="sp",
        )
        assert (
            settings4.route_string
            == "# oniom(mn15/def2svp:b3lyp/6-31g(d):uff) scrf=(smd,solvent=toluene)"
        )

        # settings with solvent specification for opt job
        settings5 = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_charge=0,
            real_multiplicity=1,
            high_level_atoms=[1, 2, 3],
            parent_jobtype="opt",
            freq=True,
            solvent_model="smd",
            solvent_id="toluene",
        )
        assert (
            settings5.route_string
            == "# opt freq oniom(mn15/def2svp:b3lyp/6-31g(d):uff) scrf=(smd,solvent=toluene)"
        )

        # settings with solvent specification for ts job
        settings5 = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_charge=0,
            real_multiplicity=1,
            high_level_atoms=[1, 2, 3],
            parent_jobtype="ts",
            freq=True,
            solvent_model="smd",
            solvent_id="toluene",
        )
        assert (
            settings5.route_string
            == "# opt=(ts,calcfc,noeigentest) freq oniom(mn15/def2svp:b3lyp/6-31g(d):uff) scrf=(smd,solvent=toluene)"
        )

        # settings with solvent specification for ts job
        settings6 = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_charge=0,
            real_multiplicity=1,
            high_level_atoms=[1, 2, 3],
            parent_jobtype="ts",
            freq=False,
            numfreq=True,
            solvent_model="smd",
            solvent_id="toluene",
        )
        assert (
            settings6.route_string
            == "# opt=(ts,calcfc,noeigentest) freq oniom(mn15/def2svp:b3lyp/6-31g(d):uff) scrf=(smd,solvent=toluene)"
        )

    def test_qmmm_additional_route_parameters(self):
        """Regression test: -r/additional_route_parameters must appear in the
        QMMM route string.  Previously _get_route_string_from_jobtype() in
        GaussianQMMMJobSettings never appended this field, so keywords like
        'scf=xqc' were silently dropped."""

        # Plain opt+freq QMMM job with extra route keyword
        settings = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="genecp",
            low_level_force_field="PM6",
            charge_total=0,
            mult_total=2,
            charge_high=2,
            mult_high=2,
            high_level_atoms=list(range(1, 66)),
            parent_jobtype="opt",
            freq=True,
            additional_route_parameters="scf=xqc",
        )
        assert settings.route_string == (
            "# opt freq oniom(mn15/genecp:PM6) scf=xqc"
        ), (
            "additional_route_parameters were not appended to the QMMM "
            "route string"
        )

        # Verify it also works without solvation on a plain sp job
        settings_sp = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            low_level_force_field="PM6",
            charge_total=0,
            mult_total=1,
            high_level_atoms=[1, 2, 3],
            parent_jobtype="sp",
            additional_route_parameters="scf=xqc opt",
        )
        assert settings_sp.route_string == (
            "# oniom(mn15/def2svp:PM6) scf=xqc opt"
        )

        # Verify it works alongside solvation
        settings_solv = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
            high_level_atoms=[1, 2, 3],
            parent_jobtype="opt",
            freq=True,
            solvent_model="smd",
            solvent_id="water",
            additional_route_parameters="scf=xqc",
        )
        assert settings_solv.route_string == (
            "# opt freq oniom(mn15/def2svp:b3lyp/6-31g(d):uff) "
            "scrf=(smd,solvent=water) scf=xqc"
        )

        # Duplicate guard: if the parameter is already in the string
        # it should not appear twice
        settings_dup = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
            high_level_atoms=[1, 2, 3],
            parent_jobtype="opt",
            additional_route_parameters="opt",  # 'opt' is already in string
        )
        route = settings_dup.route_string
        assert (
            route.count("opt") == 1
        ), "Duplicate opt keyword should not be appended when already present"

    def test_qmmm_additional_opt_options_in_route(self):
        """Regression test: -o/additional_opt_options_in_route must be merged
        into the opt/ts/modred keyword in the QMMM route string.
        Previously this field was silently ignored."""

        # opt parent with extra opt option
        s_opt = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="genecp",
            low_level_force_field="PM6",
            charge_total=0,
            mult_total=1,
            high_level_atoms=[1, 2, 3],
            parent_jobtype="opt",
            additional_opt_options_in_route="maxstep=8",
        )
        assert s_opt.route_string == "# opt=(maxstep=8) oniom(mn15/genecp:PM6)"

        # ts parent with extra opt option (no calcall)
        s_ts = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            low_level_force_field="UFF",
            charge_total=0,
            mult_total=1,
            high_level_atoms=[1, 2, 3],
            parent_jobtype="ts",
            additional_opt_options_in_route="maxstep=5",
        )
        assert s_ts.route_string == (
            "# opt=(ts,calcfc,noeigentest,maxstep=5) oniom(mn15/def2svp:UFF)"
        )

        # ts parent with calcall replaces calcfc
        s_ts_calcall = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            low_level_force_field="UFF",
            charge_total=0,
            mult_total=1,
            high_level_atoms=[1, 2, 3],
            parent_jobtype="ts",
            additional_opt_options_in_route="calcall",
        )
        assert s_ts_calcall.route_string == (
            "# opt=(ts,noeigentest,calcall) oniom(mn15/def2svp:UFF)"
        ), "calcall should replace calcfc in ts QMMM route"

        # modred parent with extra opt option
        s_modred = GaussianQMMMJobSettings(
            high_level_functional="b3lyp",
            high_level_basis="6-31g(d)",
            low_level_force_field="UFF",
            charge_total=0,
            mult_total=1,
            high_level_atoms=[1, 2, 3],
            parent_jobtype="modred",
            additional_opt_options_in_route="maxstep=10",
        )
        assert s_modred.route_string == (
            "# opt=(modredundant,maxstep=10) oniom(b3lyp/6-31g(d):UFF)"
        )

        # without additional_opt_options_in_route, opt keyword is plain
        s_plain = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            low_level_force_field="PM6",
            charge_total=0,
            mult_total=1,
            high_level_atoms=[1, 2, 3],
            parent_jobtype="opt",
        )
        assert s_plain.route_string == "# opt oniom(mn15/def2svp:PM6)"

        # empty string and whitespace-only must not produce opt=() or opt=(  )
        for blank in ("", "   ", "\t"):
            s_blank = GaussianQMMMJobSettings(
                high_level_functional="mn15",
                high_level_basis="def2svp",
                low_level_force_field="PM6",
                charge_total=0,
                mult_total=1,
                high_level_atoms=[1, 2, 3],
                parent_jobtype="opt",
                additional_opt_options_in_route=blank,
            )
            assert s_blank.route_string == "# opt oniom(mn15/def2svp:PM6)", (
                f"blank opt option {blank!r} should produce plain 'opt', "
                f"got: {s_blank.route_string}"
            )

        # same guard for ts parent
        s_ts_blank = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            low_level_force_field="UFF",
            charge_total=0,
            mult_total=1,
            high_level_atoms=[1, 2, 3],
            parent_jobtype="ts",
            additional_opt_options_in_route="  ",
        )
        assert s_ts_blank.route_string == (
            "# opt=(ts,calcfc,noeigentest) oniom(mn15/def2svp:UFF)"
        ), "whitespace-only opt option for ts should fall back to plain ts keyword"

        # same guard for modred parent
        s_modred_blank = GaussianQMMMJobSettings(
            high_level_functional="b3lyp",
            high_level_basis="6-31g(d)",
            low_level_force_field="UFF",
            charge_total=0,
            mult_total=1,
            high_level_atoms=[1, 2, 3],
            parent_jobtype="modred",
            additional_opt_options_in_route="",
        )
        assert s_modred_blank.route_string == (
            "# opt=modredundant oniom(b3lyp/6-31g(d):UFF)"
        ), "empty opt option for modred should fall back to plain opt=modredundant"

    def test_qmmm_settings_for_atoms(
        self,
        gaussian_inputs_test_directory,
        gaussian_semiempirical_pm6_output_file,
    ):
        mol1 = QMMMMolecule(
            molecule=Molecule.from_filepath(
                gaussian_semiempirical_pm6_output_file
            )
        )

        settings1 = QMMMMolecule(
            symbols=mol1.symbols,
            positions=mol1.positions,
            high_level_atoms=[1, 2, 3, 8, 9, 10],
            bonded_atoms=[[1, 2], [2, 3], [8, 9], [9, 10]],
        )
        assert settings1.high_level_atoms == [1, 2, 3, 8, 9, 10]
        assert settings1.partition_level_strings == [
            "H",
            "H",
            "H",
            "L",
            "L",
            "L",
            "L",
            "H",
            "H",
            "H",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
        ]

        # Test for 3-layer ONIOM calculation with example
        # input dppeFeCl2_phenyldioxazolone_qmmm.com
        # mol2 = QMMM(molecule=Molecule._read_gaussian_inputfile(
        #     os.path.join(
        #         gaussian_inputs_test_directory,
        #         "qmmm/dppeFeCl2_phenyldioxazolone_qmmm.com",
        #     ),
        # ))
        # settings2 = QMMM(
        #     symbols=mol2.symbols,
        #     positions=mol2.positions,
        #     high_level_atoms="[18-28,29-39,40-50,51-61,62-72]",
        #     medium_level_atoms=[1, 2, 3, 16],
        #     bonded_atoms=[[2, 18], [2, 29], [1, 40], [1, 51], [16, 62]],
        # )
        # assert settings2.high_level_atoms == list(range(18, 29)) + list(
        #     range(29, 40)
        # ) + list(range(40, 51)) + list(range(51, 62)) + list(range(62, 73))
        # assert settings2.medium_level_atoms == [1, 2, 3, 16]
        # assert settings2.low_level_atoms == list(range(4, 16)) + [17]

    def test_qmmm_settings_for_charge_and_multiplicity(self):
        # test cases for 3-layer ONIOM model
        settings1 = GaussianQMMMJobSettings(
            jobtype="opt",
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_charge=0,
            real_multiplicity=1,
        )
        assert (
            settings1.charge_and_multiplicity_string
            == "0 1 0 1 0 1 0 1 0 1 0 1"
        )

        settings2 = GaussianQMMMJobSettings(
            jobtype="sp",
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_charge=0,
            real_multiplicity=1,
            int_charge=1,
            int_multiplicity=3,
        )
        assert (
            settings2.charge_and_multiplicity_string
            == "0 1 1 3 1 3 1 3 1 3 1 3"
        )

        settings3 = GaussianQMMMJobSettings(
            jobtype="ts",
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_charge=0,
            real_multiplicity=1,
            int_charge=-1,
            int_multiplicity=2,
            model_charge=0,
            model_multiplicity=1,
        )
        assert (
            settings3.charge_and_multiplicity_string
            == "0 1 -1 2 -1 2 0 1 0 1 0 1"
        )

        # test cases for 2-layer ONIOM model

        settings4 = GaussianQMMMJobSettings(
            jobtype="opt",
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            real_charge=0,
            real_multiplicity=1,
        )
        assert settings4.charge_and_multiplicity_string == "0 1 0 1 0 1"

        settings5 = GaussianQMMMJobSettings(
            jobtype="sp",
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            real_charge=0,
            real_multiplicity=1,
            model_charge=-1,
            model_multiplicity=2,
        )
        assert settings5.charge_and_multiplicity_string == "0 1 -1 2 -1 2"


class TestGaussianRoute:
    gas_route = "opt=(ts,calcfc,noeigentest) freq m062x def2svp"
    gas_route2 = "opt b3lyp 6-31G(d) empiricaldispersion=gd3bj"
    solv_route = "mn15 def2svp scrf=(smd,solvent=DiChloroMethane)"
    solv_route2 = "# opt=(calcfc,ts,noeigentest,maxstep=5) freq\nscrf=(smd,solvent=chlorobenzene) nosymm scf=qc def2svp m06"

    def test_settings_from_route(self):
        route_object1 = GaussianRoute(route_string=self.gas_route)
        assert isinstance(route_object1, object)
        assert route_object1.jobtype == "ts"
        assert route_object1.freq is True
        assert route_object1.functional == "m062x"
        assert route_object1.basis == "def2svp"
        assert route_object1.solvent_model is None
        assert route_object1.solvent_id is None

    def test_settings_from_route2(self):
        route_object2 = GaussianRoute(route_string=self.gas_route2)
        assert isinstance(route_object2, object)
        assert route_object2.jobtype == "opt"
        assert route_object2.freq is False
        assert route_object2.functional == "b3lyp empiricaldispersion=gd3bj"
        assert route_object2.basis == "6-31G(d)".lower()
        assert route_object2.solvent_model is None
        assert route_object2.solvent_id is None

    def test_settings_from_route3(self):
        route_object3 = GaussianRoute(route_string=self.solv_route)
        assert isinstance(route_object3, object)
        assert route_object3.jobtype == "sp"
        assert route_object3.freq is False
        assert route_object3.functional == "mn15"
        assert route_object3.basis == "def2svp"
        assert route_object3.solvent_model == "smd"
        assert route_object3.solvent_id == "DiChloroMethane".lower()

    def test_settings_from_route4(self):
        route_object4 = GaussianRoute(route_string=self.solv_route2)
        assert isinstance(route_object4, object)
        assert route_object4.jobtype == "ts"
        assert route_object4.freq is True
        assert route_object4.functional == "m06"
        assert route_object4.basis == "def2svp"
        assert route_object4.solvent_model == "smd"
        assert route_object4.solvent_id == "chlorobenzene"


class TestGaussianJobFromComFile:
    def test_reads_com_file(self, gaussian_opt_inputfile):
        com_settings = GaussianJobSettings.from_comfile(gaussian_opt_inputfile)
        assert com_settings.chk is True
        assert com_settings.jobtype == "opt"
        assert com_settings.freq is True
        assert com_settings.functional == "m062x"
        assert com_settings.basis == "def2svp"
        assert com_settings.charge == 0
        assert com_settings.multiplicity == 1
        assert com_settings.solvent_model is None
        assert com_settings.solvent_id is None

    def test_update_solvent(self, gaussian_opt_inputfile, tmpdir):
        com_settings = GaussianJobSettings.from_comfile(gaussian_opt_inputfile)
        assert com_settings.solvent_model is None
        assert com_settings.solvent_id is None
        com_settings.update_solvent(solvent_model="smd", solvent_id="toluene")
        assert com_settings.solvent_model == "smd"
        assert com_settings.solvent_id == "toluene"

    def test_include_solvent(self, gaussian_opt_inputfile, tmpdir):
        com_settings = GaussianJobSettings.from_filepath(
            gaussian_opt_inputfile
        )
        assert com_settings.solvent_model is None
        assert com_settings.solvent_id is None

        com_settings.modify_solvent(
            remove_solvent=False, solvent_model="cpcm", solvent_id="water"
        )
        assert com_settings.solvent_model == "cpcm"
        assert com_settings.solvent_id == "water"

        com_settings.modify_solvent(remove_solvent=True)
        assert com_settings.solvent_model is None
        assert com_settings.solvent_id is None

    def test_cli_group_solvent_options_propagate_to_opt(
        self, gaussian_yaml_settings_gas_solv_project_name
    ):
        """Solvent options given at the gaussian group level propagate to opt settings.

        This simulates: ``gaussian -sm smd -si water -so iterative opt``
        where the solvent options live on the *group* command and are merged
        into the subcommand settings via :meth:`GaussianJobSettings.merge`.
        """
        from chemsmart.settings.gaussian import GaussianProjectSettings

        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        opt_settings = project_settings.opt_settings()
        # Project has no solvent for the gas/opt path
        assert opt_settings.solvent_model is None
        assert opt_settings.solvent_id is None

        # Simulate what the gaussian group CLI callback does when
        # -sm smd -si water -so iterative are supplied
        job_settings = GaussianJobSettings.default()
        job_settings.solvent_model = "smd"
        job_settings.solvent_id = "water"
        job_settings.additional_solvent_options = "iterative"
        keywords = (
            "charge",
            "multiplicity",
            "solvent_model",
            "solvent_id",
            "additional_solvent_options",
        )

        # Simulate what the opt subcommand does
        opt_settings = opt_settings.merge(job_settings, keywords=keywords)

        assert opt_settings.solvent_model == "smd"
        assert opt_settings.solvent_id == "water"
        assert opt_settings.additional_solvent_options == "iterative"
        assert (
            "scrf=(smd,solvent=water,iterative)" in opt_settings.route_string
        )

    def test_cli_group_solvent_options_propagate_to_td(self):
        """Solvent options given at the gaussian group level propagate to td settings.

        This simulates: ``gaussian -sm smd -si water -so iterative td``
        """
        from chemsmart.jobs.gaussian.settings import GaussianTDDFTJobSettings

        base_settings = GaussianJobSettings.default()
        base_settings.functional = "cam-b3lyp"
        base_settings.basis = "def2svp"

        # Simulate what the gaussian group CLI callback does
        job_settings = GaussianJobSettings.default()
        job_settings.solvent_model = "smd"
        job_settings.solvent_id = "water"
        job_settings.additional_solvent_options = "iterative"
        keywords = (
            "charge",
            "multiplicity",
            "solvent_model",
            "solvent_id",
            "additional_solvent_options",
        )

        # Simulate what the td subcommand does
        td_settings = base_settings.merge(job_settings, keywords=keywords)
        td_settings = GaussianTDDFTJobSettings(**td_settings.__dict__)
        td_settings.states = "singlets"
        td_settings.root = 1
        td_settings.nstates = 3

        assert td_settings.solvent_model == "smd"
        assert td_settings.solvent_id == "water"
        assert td_settings.additional_solvent_options == "iterative"
        assert "scrf=(smd,solvent=water,iterative)" in td_settings.route_string
        assert "TD(" in td_settings.route_string

    def test_cli_group_remove_solvent_overrides_project_settings(
        self, gaussian_yaml_settings_gas_solv_project_name
    ):
        """``--remove-solvent`` at the gaussian group level clears project solvent.

        This simulates: ``gaussian --remove-solvent sp`` when the project's
        solvent-phase sp settings carry a solvent model.
        """
        from chemsmart.settings.gaussian import GaussianProjectSettings

        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        sp_settings = project_settings.sp_settings()
        # Project has solvent for the sp (solv) path
        assert sp_settings.solvent_model == "smd"
        assert sp_settings.solvent_id == "toluene"

        # Simulate what the gaussian group CLI callback does for --remove-solvent
        job_settings = GaussianJobSettings.default()
        job_settings.solvent_model = None
        job_settings.solvent_id = None
        job_settings.custom_solvent = None
        keywords = (
            "charge",
            "multiplicity",
            "solvent_model",
            "solvent_id",
            "custom_solvent",
        )

        # Simulate what the sp subcommand does (merge picks up the None values)
        sp_settings = sp_settings.merge(job_settings, keywords=keywords)

        assert sp_settings.solvent_model is None
        assert sp_settings.solvent_id is None
        assert "scrf" not in sp_settings.route_string


class TestGaussianJobFromLogFile:
    def test_accumulates_settings(self, tmpdir, gaussian_ts_genecp_outfile):
        settings = GaussianJobSettings.from_logfile(gaussian_ts_genecp_outfile)
        assert settings.functional == "mn15"
        assert settings.basis == "genecp"
        assert settings.solvent_model is None
        assert settings.solvent_id is None

        settings.functional = "b3lyp"
        assert settings.solvent_model is None
        assert settings.functional == "b3lyp"

        settings.basis = "def2tzvp"
        assert settings.functional == "b3lyp"
        assert settings.basis == "def2tzvp"

        settings.update_solvent(solvent_model="smd", solvent_id="toluene")
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "toluene"

    def test_reads_gaussian_ts_genecp_outfile(
        self, tmpdir, gaussian_ts_genecp_outfile
    ):
        settings = GaussianJobSettings.from_logfile(gaussian_ts_genecp_outfile)
        assert settings.jobtype == "ts"
        assert settings.functional == "mn15"
        assert settings.basis == "genecp"
        assert settings.solvent_model is None
        assert settings.solvent_id is None

    def test_reads_gaussian_pm6_outfile(
        self, tmpdir, gaussian_semiempirical_pm6_output_file
    ):
        settings = GaussianJobSettings.from_logfile(
            gaussian_semiempirical_pm6_output_file
        )
        assert settings.jobtype == "opt"
        assert settings.ab_initio is None
        assert settings.functional is None
        assert settings.basis is None
        assert settings.semiempirical == "PM6"
        assert settings.solvent_model is None
        assert settings.solvent_id is None


class TestGaussianPBCJob:
    def test_writes_gaussian_input_from_pbc_comfile(
        self, tmpdir, gaussian_pbc_3d_outputfile
    ):
        settings = GaussianJobSettings.from_filepath(
            filepath=gaussian_pbc_3d_outputfile
        )
        assert settings.jobtype == "sp"
        assert settings.functional.lower() == "pbepbe"
        assert settings.basis.lower() == "6-31g(d,p)/auto"
        assert settings.additional_route_parameters.lower() == "scf=tight"


class TestGaussianLinkJobSettingsGuess:
    """Tests for guess= formatting in GaussianLinkJobSettings route strings."""

    _COMMON = dict(functional="um062x", basis="def2svp", charge=0, multiplicity=1)

    def _route(self, guess):
        s = GaussianLinkJobSettings(guess=guess, **self._COMMON)
        return s._get_route_string_from_jobtype()

    def test_single_guess_option_no_parentheses(self):
        """Single option must appear without parentheses: guess=mix"""
        assert "guess=mix" in self._route("mix")
        assert "guess=(mix)" not in self._route("mix")

    def test_multiple_guess_options_with_parentheses(self):
        """Multiple comma-separated options must be wrapped: guess=(mix,always)"""
        assert "guess=(mix,always)" in self._route("mix,always")

    def test_pre_parenthesized_input_no_double_wrapping(self):
        """Already-parenthesized input must not produce double parentheses."""
        route = self._route("(mix,always)")
        assert "guess=(mix,always)" in route
        assert "guess=((mix,always))" not in route

    def test_pre_parenthesized_input_with_whitespace(self):
        """Whitespace around parenthesized input must be handled."""
        route = self._route(" (mix,always) ")
        assert "guess=(mix,always)" in route
        assert "guess=((mix,always))" not in route
