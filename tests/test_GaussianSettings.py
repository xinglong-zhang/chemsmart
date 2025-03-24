import os

import pytest

from chemsmart.io.gaussian.route import GaussianRoute
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.settings import (
    GaussianJobSettings,
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
            functional_high="b3lyp",
            basis_high="6-31g(d)",
            force_field_low="uff",
            high_level_charge=0,
            high_level_multiplicity=1,
            high_level_atoms=[1, 2, 3],
            job_type="opt",
            freq=True,
        )
        assert settings1.route_string == "# opt freq oniom(b3lyp/6-31g(d):uff)"

        settings2 = GaussianQMMMJobSettings(
            functional_high="mn15",
            basis_high="def2svp",
            functional_medium="b3lyp",
            basis_medium="6-31g(d)",
            force_field_low="uff",
            high_level_charge=0,
            high_level_multiplicity=1,
            high_level_atoms=[1, 2, 3],
        )
        assert (
            settings2.route_string
            == "# oniom(mn15/def2svp:b3lyp/6-31g(d):uff)"
        )

        settings3 = GaussianQMMMJobSettings(
            functional_high="mn15",
            basis_high="def2svp",
            force_field_high="uff",
            functional_medium="b3lyp",
            basis_medium="6-31g(d)",
            force_field_low="uff",
            high_level_charge=0,
            high_level_multiplicity=1,
            high_level_atoms=[1, 2, 3],
        )
        # assert settings3.route_string == "# oniom(mn15/def2svp:uff:b3lyp/6-31g(d):uff)"
        # ValueError: For high level of theory, one should specify only functional/basis or force field!
        with pytest.raises(ValueError):
            settings3.route_string

        # settings with solvent specification
        settings4 = GaussianQMMMJobSettings(
            functional_high="mn15",
            basis_high="def2svp",
            functional_medium="b3lyp",
            basis_medium="6-31g(d)",
            force_field_low="uff",
            high_level_charge=0,
            high_level_multiplicity=1,
            high_level_atoms=[1, 2, 3],
            solvent_model="smd",
            solvent_id="toluene",
        )
        assert (
            settings4.route_string
            == "# oniom(mn15/def2svp:b3lyp/6-31g(d):uff) scrf=(smd,solvent=toluene)"
        )

        # settings with solvent specification for opt job
        settings5 = GaussianQMMMJobSettings(
            functional_high="mn15",
            basis_high="def2svp",
            functional_medium="b3lyp",
            basis_medium="6-31g(d)",
            force_field_low="uff",
            high_level_charge=0,
            high_level_multiplicity=1,
            high_level_atoms=[1, 2, 3],
            job_type="opt",
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
            functional_high="mn15",
            basis_high="def2svp",
            functional_medium="b3lyp",
            basis_medium="6-31g(d)",
            force_field_low="uff",
            high_level_charge=0,
            high_level_multiplicity=1,
            high_level_atoms=[1, 2, 3],
            job_type="ts",
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
            functional_high="mn15",
            basis_high="def2svp",
            functional_medium="b3lyp",
            basis_medium="6-31g(d)",
            force_field_low="uff",
            high_level_charge=0,
            high_level_multiplicity=1,
            high_level_atoms=[1, 2, 3],
            job_type="ts",
            freq=False,
            numfreq=True,
            solvent_model="smd",
            solvent_id="toluene",
        )
        assert (
            settings6.route_string
            == "# opt=(ts,calcfc,noeigentest) freq oniom(mn15/def2svp:b3lyp/6-31g(d):uff) scrf=(smd,solvent=toluene)"
        )

    def test_qmmm_settings_for_atoms(self, gaussian_inputs_test_directory):
        mol1 = Molecule.from_pubchem("81184")

        settings1 = Molecule(
            symbols=mol1.symbols,
            positions=mol1.positions,
            high_level_atoms=[1, 2, 3, 8, 9, 10],
            bonded_atoms=[[1, 2], [2, 3], [8, 9], [9, 10]],
        )
        assert settings1.high_level_atoms == [1, 2, 3, 8, 9, 10]
        assert settings1.low_level_atoms == [
            4,
            5,
            6,
            7,
            11,
            12,
            13,
            14,
            15,
            16,
            17,
            18,
            19,
            20,
            21,
            22,
            23,
            24,
            25,
            26,
            27,
        ]
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

        # Test for 3-layer ONIOM calculation with example input dppeFeCl2_phenyldioxazolone_qmmm.com
        mol2 = Molecule._read_gaussian_inputfile(
            os.path.join(
                gaussian_inputs_test_directory,
                "qmmm/dppeFeCl2_phenyldioxazolone_qmmm.com",
            )
        )
        settings2 = Molecule(
            symbols=mol2.symbols,
            positions=mol2.positions,
            high_level_atoms="[18-28,29-39,40-50,51-61,62-72]",
            medium_level_atoms=[1, 2, 3, 16],
            bonded_atoms=[[2, 18], [2, 29], [1, 40], [1, 51], [16, 62]],
        )
        assert settings2.high_level_atoms == list(range(18, 29)) + list(
            range(29, 40)
        ) + list(range(40, 51)) + list(range(51, 62)) + list(range(62, 73))
        assert settings2.medium_level_atoms == [1, 2, 3, 16]
        assert settings2.low_level_atoms == list(range(4, 16)) + [17]

    def test_qmmm_settings_for_charge_and_multiplicity(self):
        # test cases for 3-layer ONIOM model
        settings1 = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_low_charge=0,
            real_low_multiplicity=1,
        )
        assert (
            settings1.charge_and_multiplicity_string
            == "0 1 0 1 0 1 0 1 0 1 0 1"
        )

        settings2 = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_low_charge=0,
            real_low_multiplicity=1,
            int_med_charge=1,
            int_med_multiplicity=3,
        )
        assert (
            settings2.charge_and_multiplicity_string
            == "0 1 1 3 1 3 1 3 1 3 1 3"
        )

        settings3 = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_low_charge=0,
            real_low_multiplicity=1,
            int_med_charge=-1,
            int_med_multiplicity=2,
            int_low_charge=0,
            int_low_multiplicity=1,
        )
        assert (
            settings3.charge_and_multiplicity_string
            == "0 1 -1 2 0 1 -1 2 -1 2 -1 2"
        )

        settings4 = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_low_charge=0,
            real_low_multiplicity=1,
            int_med_charge=0,
            int_med_multiplicity=3,
            int_low_charge=0,
            int_low_multiplicity=1,
            model_high_charge=1,
            model_high_multiplicity=2,
        )
        assert (
            settings4.charge_and_multiplicity_string
            == "0 1 0 3 0 1 1 2 1 2 1 2"
        )

        settings5 = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_low_charge=0,
            real_low_multiplicity=1,
            int_med_charge=1,
            int_med_multiplicity=3,
            int_low_charge=0,
            int_low_multiplicity=1,
            model_high_charge=1,
            model_high_multiplicity=2,
            model_med_charge=1,
            model_med_multiplicity=2,
        )
        assert (
            settings5.charge_and_multiplicity_string
            == "0 1 1 3 0 1 1 2 1 2 1 2"
        )

        settings6 = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_low_charge=0,
            real_low_multiplicity=1,
            int_med_charge=-1,
            int_med_multiplicity=2,
            int_low_charge=0,
            int_low_multiplicity=1,
            model_high_charge=1,
            model_high_multiplicity=2,
            model_med_charge=1,
            model_med_multiplicity=2,
            model_low_charge=1,
            model_low_multiplicity=2,
        )
        assert (
            settings6.charge_and_multiplicity_string
            == "0 1 -1 2 0 1 1 2 1 2 1 2"
        )
        # test cases for 2-layer ONIOM model
        settings7 = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            low_level_force_field="uff",
            real_low_charge=0,
            real_low_multiplicity=1,
        )
        assert settings7.charge_and_multiplicity_string == "0 1 0 1 0 1"

        settings8 = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            real_low_charge=0,
            real_low_multiplicity=1,
            model_high_charge=-1,
            model_high_multiplicity=2,
        )
        assert settings8.charge_and_multiplicity_string == "0 1 -1 2 -1 2"

        settings9 = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            real_low_charge=0,
            real_low_multiplicity=1,
            model_high_charge=-1,
            model_high_multiplicity=2,
            model_low_charge=0,
            model_low_multiplicity=1,
        )
        assert settings9.charge_and_multiplicity_string == "0 1 -1 2 0 1"


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
