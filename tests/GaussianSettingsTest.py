from chemsmart.jobs.settings import read_molecular_job_yaml
from chemsmart.io.gaussian.route import GaussianRoute
from chemsmart.jobs.gaussian.settings import GaussianJobSettings


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
            gaussian_yaml_settings_gas_solv
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
            gaussian_yaml_settings_solv
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


class TestGaussianRoute:
    gas_route = "opt=(ts,calcfc,noeigentest) freq m062x def2svp"
    gas_route2 = "opt b3lyp 6-31G(d) empiricaldispersion=gd3bj"
    solv_route = "mn15 def2svp scrf=(smd,solvent=DiChloroMethane)"
    solv_route2 = "# opt=(calcfc,ts,noeigentest,maxstep=5) freq\nscrf=(smd,solvent=chlorobenzene) nosymm scf=qc def2svp m06"

    def test_settings_from_route(self):
        route_object1 = GaussianRoute(route_string=self.gas_route)
        assert isinstance(route_object1, object)
        assert route_object1.job_type == "ts"
        assert route_object1.freq is True
        assert route_object1.functional == "m062x"
        assert route_object1.basis == "def2svp"
        assert route_object1.solvent_model is None
        assert route_object1.solvent_id is None

    def test_settings_from_route2(self):
        route_object2 = GaussianRoute(route_string=self.gas_route2)
        assert isinstance(route_object2, object)
        assert route_object2.job_type == "opt"
        assert route_object2.freq is False
        assert route_object2.functional == "b3lyp empiricaldispersion=gd3bj"
        assert route_object2.basis == "6-31G(d)".lower()
        assert route_object2.solvent_model is None
        assert route_object2.solvent_id is None

    def test_settings_from_route3(self):
        route_object3 = GaussianRoute(route_string=self.solv_route)
        assert isinstance(route_object3, object)
        assert route_object3.job_type == "sp"
        assert route_object3.freq is False
        assert route_object3.functional == "mn15"
        assert route_object3.basis == "def2svp"
        assert route_object3.solvent_model == "smd"
        assert route_object3.solvent_id == "DiChloroMethane".lower()

    def test_settings_from_route4(self):
        route_object4 = GaussianRoute(route_string=self.solv_route2)
        assert isinstance(route_object4, object)
        assert route_object4.job_type == "ts"
        assert route_object4.freq is True
        assert route_object4.functional == "m06"
        assert route_object4.basis == "def2svp"
        assert route_object4.solvent_model == "smd"
        assert route_object4.solvent_id == "chlorobenzene"


class TestGaussianJobFromComFile:
    def test_reads_com_file(self, gaussian_opt_inputfile):
        com_settings = GaussianJobSettings.from_comfile(gaussian_opt_inputfile)
        assert com_settings.chk is True
        assert com_settings.job_type == "opt"
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


class TestGaussianPBCJob:
    def test_writes_gaussian_input_from_pbc_comfile(
        self, tmpdir, gaussian_pbc_3d_outputfile
    ):
        settings = GaussianJobSettings.from_filepath(
            filepath=gaussian_pbc_3d_outputfile
        )
        assert settings.job_type == "sp"
        assert settings.functional.lower() == "pbepbe"
        assert settings.basis.lower() == "6-31g(d,p)/auto"
        assert settings.additional_route_parameters.lower() == "scf=tight"
