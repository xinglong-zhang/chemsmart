import os
from unittest.mock import MagicMock, patch

import pytest
import yaml

from chemsmart.io.gaussian.route import GaussianRoute
from chemsmart.io.molecules.structure import Molecule, QMMMMolecule
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.pka import GaussianpKaJob
from chemsmart.jobs.gaussian.settings import (
    GaussianIRCJobSettings,
    GaussianJobSettings,
    GaussianLinkJobSettings,
    GaussianpKaJobSettings,
    GaussianQMMMJobSettings,
)
from chemsmart.jobs.settings import (
    MolecularJobSettings,
    read_molecular_job_yaml,
)


def _mol_with_h(charge=None, multiplicity=None, frozen_atoms=None):
    mol = Molecule(
        symbols=["C", "H", "H", "H", "H"],
        positions=[
            [0.0, 0.0, 0.0],
            [0.63, 0.63, 0.63],
            [0.63, -0.63, -0.63],
            [-0.63, 0.63, -0.63],
            [-0.63, -0.63, 0.63],
        ],
        frozen_atoms=frozen_atoms,
    )
    mol.charge = charge
    mol.multiplicity = multiplicity
    return mol


def _mol(symbols=("H", "H"), positions=None):
    if positions is None:
        positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]
    return Molecule(symbols=list(symbols), positions=positions)


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

    def test_get_settings_from_yaml_no_defaults_file_gaussian(self, tmp_path):
        project_yaml = tmp_path / "project.yaml"
        project_yaml.write_text(
            yaml.safe_dump({"gas": {"functional": "m062x"}, "solv": {}})
        )
        all_project_settings = read_molecular_job_yaml(
            str(project_yaml), program="gaussian"
        )
        assert all_project_settings["opt"]["functional"] == "m062x"

    def test_get_settings_from_yaml_no_defaults_file_orca(self, tmp_path):
        project_yaml = tmp_path / "project.yaml"
        project_yaml.write_text(
            yaml.safe_dump({"gas": {"functional": "m062x"}, "solv": {}})
        )
        all_project_settings = read_molecular_job_yaml(
            str(project_yaml), program="orca"
        )
        assert all_project_settings["opt"]["functional"] == "m062x"

    def test_get_settings_from_yaml_no_defaults_file_xtb(self, tmp_path):
        project_yaml = tmp_path / "project.yaml"
        project_yaml.write_text(yaml.safe_dump({"gas": {}, "solv": {}}))
        all_project_settings = read_molecular_job_yaml(
            str(project_yaml), program="xtb"
        )
        # Falls back to XTBJobSettings.default() when no defaults.yaml.
        assert all_project_settings["opt"]["gfn_version"] == "gfn2"
        assert all_project_settings["opt"]["jobtype"] == "opt"

    def test_get_settings_from_yaml_no_defaults_file_other_program(
        self, tmp_path
    ):
        project_yaml = tmp_path / "project.yaml"
        project_yaml.write_text(yaml.safe_dump({"solv": {}}))
        all_project_settings = read_molecular_job_yaml(
            str(project_yaml), program="other"
        )
        # default_config stays {} for an unrecognized program; jobtype is
        # still populated per job.
        assert all_project_settings["opt"]["jobtype"] == "opt"

    def test_get_settings_from_yaml_gas_config_bad_key_raises(self, tmp_path):
        # See BUGS_FOUND.md: read_molecular_job_yaml has an unguarded
        # duplicate call to update_dict_with_existing_keys(..., gas_config)
        # immediately before the try/except that's meant to catch exactly
        # this failure and fall back to qmmm_config. The unguarded call
        # always raises first, so the except branch is unreachable dead
        # code -- a bad key in gas_config always propagates as a ValueError
        # instead of triggering the documented qmmm_config fallback.
        defaults_yaml = tmp_path / "defaults.yaml"
        defaults_yaml.write_text(
            yaml.safe_dump({"functional": None, "basis": None})
        )
        project_yaml = tmp_path / "project.yaml"
        project_yaml.write_text(
            yaml.safe_dump(
                {
                    "gas": {"high_level_functional": "MN15"},
                    "solv": {},
                    "qmmm": {"functional": "pbe"},
                }
            )
        )
        with pytest.raises(ValueError, match="high_level_functional"):
            read_molecular_job_yaml(str(project_yaml), program="gaussian")

    def test_get_settings_from_yaml_td_config(self, tmp_path):
        defaults_yaml = tmp_path / "defaults.yaml"
        defaults_yaml.write_text(
            yaml.safe_dump({"functional": None, "basis": None})
        )
        project_yaml = tmp_path / "project.yaml"
        project_yaml.write_text(
            yaml.safe_dump({"solv": {}, "td": {"functional": "cam-b3lyp"}})
        )
        all_project_settings = read_molecular_job_yaml(
            str(project_yaml), program="gaussian"
        )
        assert all_project_settings["td"]["functional"] == "cam-b3lyp"
        assert all_project_settings["td"]["jobtype"] == "td"

    def test_get_settings_from_yaml_qmmm_config(
        self, gaussian_yaml_settings_qmmm
    ):
        all_project_settings = read_molecular_job_yaml(
            gaussian_yaml_settings_qmmm, program="gaussian"
        )
        qmmm_settings = all_project_settings["qmmm"]
        assert qmmm_settings["jobtype"] == "qmmm"
        assert qmmm_settings["high_level_functional"] == "MN15"
        assert qmmm_settings["low_level_force_field"] == "UFF"

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
    def test_builtin_mm_helpers(self, tmpdir):
        amber = GaussianQMMMJobSettings(
            low_level_force_field="AMBER=HardFirst"
        )
        uff = GaussianQMMMJobSettings(low_level_force_field="UFF")
        dreiding = GaussianQMMMJobSettings(low_level_force_field="dreiding")
        pm6 = GaussianQMMMJobSettings(low_level_force_field="PM6")
        empty = GaussianQMMMJobSettings()

        assert amber.uses_builtin_mm()
        assert uff.uses_builtin_mm()
        assert dreiding.uses_builtin_mm()
        assert not pm6.uses_builtin_mm()
        assert not empty.uses_builtin_mm()

        assert amber.requires_mm_atom_info()
        assert not uff.requires_mm_atom_info()
        assert not pm6.requires_mm_atom_info()

        assert GaussianQMMMJobSettings.format_mm_atom_label("C", None) == "C"
        assert (
            GaussianQMMMJobSettings.format_mm_atom_label(
                "C", ("CT", 0.03, None, None)
            )
            == "C-CT-0.03"
        )
        assert (
            GaussianQMMMJobSettings.format_mm_atom_label(
                "O", ("OH", -0.65, None, None)
            )
            == "O-OH--0.65"
        )
        assert GaussianQMMMJobSettings.format_mm_link_atom(3, None) == "H 3"
        assert (
            GaussianQMMMJobSettings.format_mm_link_atom(
                3, ("CT", 0.03, "HC", 0.09)
            )
            == "H-HC-0.09 3"
        )

        indexed = os.path.join(tmpdir, "indexed.dat")
        with open(indexed, "w") as handle:
            handle.write("# comment\n")
            handle.write("2 CT 0.03 HC 0.09\n")
            handle.write("1 OH -0.65\n")
        records = GaussianQMMMJobSettings.load_mm_atom_info(
            indexed, num_atoms=2
        )
        assert records[0][0] == "OH"
        assert records[0][1] == -0.65
        assert records[1][2] == "HC"
        assert records[1][3] == 0.09

        ordered = os.path.join(tmpdir, "ordered.dat")
        with open(ordered, "w") as handle:
            handle.write("OH -0.65\n")
            handle.write("CT 0.03 HC 0.09\n")
        records = GaussianQMMMJobSettings.load_mm_atom_info(
            ordered, num_atoms=2
        )
        assert [r[0] for r in records] == ["OH", "CT"]

        bad = os.path.join(tmpdir, "bad.dat")
        with open(bad, "w") as handle:
            handle.write("OH -0.65\n")
        with pytest.raises(ValueError, match="expected 2"):
            GaussianQMMMJobSettings.load_mm_atom_info(bad, num_atoms=2)

        mixed = os.path.join(tmpdir, "mixed.dat")
        with open(mixed, "w") as handle:
            handle.write("1 CT 0.03\n")
            handle.write("OH -0.65\n")
        with pytest.raises(ValueError, match="mixes indexed and unindexed"):
            GaussianQMMMJobSettings.load_mm_atom_info(mixed, num_atoms=2)

        with pytest.raises(FileNotFoundError, match="not found"):
            GaussianQMMMJobSettings.load_mm_atom_info(
                os.path.join(tmpdir, "missing.dat"), num_atoms=1
            )

        short = os.path.join(tmpdir, "short.dat")
        with open(short, "w") as handle:
            handle.write("OH\n")
        with pytest.raises(ValueError, match="Invalid MM atom info line"):
            GaussianQMMMJobSettings.load_mm_atom_info(short, num_atoms=1)

        bad_charge = os.path.join(tmpdir, "bad_charge.dat")
        with open(bad_charge, "w") as handle:
            handle.write("CT not-a-float\n")
        with pytest.raises(ValueError, match="Invalid MM atom info line"):
            GaussianQMMMJobSettings.load_mm_atom_info(bad_charge, num_atoms=1)

        out_of_range = os.path.join(tmpdir, "oor.dat")
        with open(out_of_range, "w") as handle:
            handle.write("1 CT 0.03\n")
            handle.write("3 OH -0.65\n")
        with pytest.raises(ValueError, match="out-of-range"):
            GaussianQMMMJobSettings.load_mm_atom_info(
                out_of_range, num_atoms=2
            )

        missing_index = os.path.join(tmpdir, "missing_index.dat")
        with open(missing_index, "w") as handle:
            handle.write("1 CT 0.03\n")
            handle.write("3 OH -0.65\n")
        with pytest.raises(ValueError, match="missing atom indices"):
            GaussianQMMMJobSettings.load_mm_atom_info(
                missing_index, num_atoms=3
            )

        medium_amber = GaussianQMMMJobSettings(
            medium_level_force_field="amber"
        )
        assert medium_amber.uses_builtin_mm()
        assert medium_amber.requires_mm_atom_info()

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
        assert (
            settings1.route_string
            == "# opt freq oniom(b3lyp/6-31g(d):uff) geom=connectivity"
        )

        settings_medium_theory_no_atoms = GaussianQMMMJobSettings(
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
            settings_medium_theory_no_atoms.route_string
            == "# oniom(mn15/def2svp:uff) geom=connectivity"
        )
        assert (
            settings_medium_theory_no_atoms.charge_and_multiplicity_string
            == "0 1 0 1 0 1"
        )

        settings2 = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_charge=0,
            real_multiplicity=1,
            high_level_atoms=[1, 2, 3],
            medium_level_atoms=[4],
            parent_jobtype="sp",
        )
        assert (
            settings2.route_string
            == "# oniom(mn15/def2svp:b3lyp/6-31g(d):uff) geom=connectivity"
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
            medium_level_atoms=[4],
            parent_jobtype="sp",
        )
        # assert settings3.route_string == "#
        # oniom(mn15/def2svp:uff:b3lyp/6-31g(d):uff) geom=connectivity"
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
            medium_level_atoms=[4],
            solvent_model="smd",
            solvent_id="toluene",
            parent_jobtype="sp",
        )
        assert (
            settings4.route_string
            == "# oniom(mn15/def2svp:b3lyp/6-31g(d):uff) geom=connectivity scrf=(smd,solvent=toluene)"
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
            medium_level_atoms=[4],
            parent_jobtype="opt",
            freq=True,
            solvent_model="smd",
            solvent_id="toluene",
        )
        assert (
            settings5.route_string
            == "# opt freq oniom(mn15/def2svp:b3lyp/6-31g(d):uff) geom=connectivity scrf=(smd,solvent=toluene)"
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
            medium_level_atoms=[4],
            parent_jobtype="ts",
            freq=True,
            solvent_model="smd",
            solvent_id="toluene",
        )
        assert (
            settings5.route_string
            == "# opt=(ts,calcfc,noeigentest) freq oniom(mn15/def2svp:b3lyp/6-31g(d):uff) geom=connectivity scrf=(smd,solvent=toluene)"
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
            medium_level_atoms=[4],
            parent_jobtype="ts",
            freq=False,
            numfreq=True,
            solvent_model="smd",
            solvent_id="toluene",
        )
        assert (
            settings6.route_string
            == "# opt=(ts,calcfc,noeigentest) freq oniom(mn15/def2svp:b3lyp/6-31g(d):uff) geom=connectivity scrf=(smd,solvent=toluene)"
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
            medium_level_atoms=[4],
            parent_jobtype="opt",
            freq=True,
            solvent_model="smd",
            solvent_id="water",
            additional_route_parameters="scf=xqc",
        )
        assert settings_solv.route_string == (
            "# opt freq oniom(mn15/def2svp:b3lyp/6-31g(d):uff) geom=connectivity "
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
            "# opt=(ts,calcfc,noeigentest,maxstep=5) oniom(mn15/def2svp:UFF) geom=connectivity"
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
            "# opt=(ts,noeigentest,calcall) oniom(mn15/def2svp:UFF) geom=connectivity"
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
            "# opt=(modredundant,maxstep=10) oniom(b3lyp/6-31g(d):UFF) geom=connectivity"
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
            "# opt=(ts,calcfc,noeigentest) oniom(mn15/def2svp:UFF) geom=connectivity"
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
            "# opt=modredundant oniom(b3lyp/6-31g(d):UFF) geom=connectivity"
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
            medium_level_atoms=[4],
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
            medium_level_atoms=[4],
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
            medium_level_atoms=[4],
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
        assert route_object2.functional == "b3lyp-d3bj"
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
    def test_reads_custom_solvent_from_logfile_as_settings_text(
        self, gaussian_smd_generic_outfile
    ):
        """Parsed output parameters can be reused or safely overridden."""
        settings = GaussianJobSettings.from_logfile(
            gaussian_smd_generic_outfile
        )

        assert isinstance(settings.custom_solvent, str)
        assert (
            "SolventName=1,1,1,3,3,3-HEXAFLUOROPROPAN-2-OL"
            in settings.custom_solvent
        )
        assert "Eps=16.7" in settings.custom_solvent
        assert "ElectronegativeHalogenicity=0.6" in settings.custom_solvent

        replacement = "SolventName=replacement\nEps=2.5\n"
        project_settings = GaussianJobSettings.default()
        project_settings.custom_solvent = replacement
        merged = project_settings.merge(settings)

        assert merged.custom_solvent == replacement

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

    def test_reads_oldform_redundant_coordinates_with_atomic_numbers(
        self, tmp_path
    ):
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
        settings = GaussianJobSettings.from_logfile(str(outputfile))
        assert settings.jobtype == "opt"
        assert settings.functional == "b3lyp"
        assert settings.basis == "gen"
        assert settings.charge == 0
        assert settings.multiplicity == 1

    def test_reads_pd_insertion_ts_r_logfile(
        self, gaussian_pd_insertion_ts_r_outfile
    ):
        settings = GaussianJobSettings.from_logfile(
            gaussian_pd_insertion_ts_r_outfile
        )
        assert settings.charge == 0
        assert settings.multiplicity == 1
        assert settings.functional == "b3lyp-d3"


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


class TestGaussianpKaJobSettings:
    """Tests for GaussianpKaJobSettings and GaussianpKaJob classes."""

    def test_init_custom_values(self):
        """Test initialization with custom values."""
        settings = GaussianpKaJobSettings(
            proton_index=10,
            scheme="proton exchange",
            solvent_model="PCM",
            solvent_id="water",
            charge=0,  # Protonated form charge (inherited from parent)
            multiplicity=1,  # Protonated form multiplicity (inherited from parent)
            conjugate_base_charge=-1,
            conjugate_base_multiplicity=1,
            functional="B3LYP",
            basis="6-311+G(d,p)",
        )
        assert settings.proton_index == 10
        assert settings.scheme == "proton exchange"
        assert settings.solvent_model == "PCM"
        assert settings.solvent_id == "water"
        assert settings.charge == 0
        assert settings.multiplicity == 1
        assert settings.protonated_charge == 0
        assert settings.protonated_multiplicity == 1
        assert settings.conjugate_base_charge == -1
        assert settings.conjugate_base_multiplicity == 1
        assert settings.functional == "B3LYP"
        assert settings.basis == "6-311+G(d,p)"

    def test_direct_cycle_no_reference(self):
        """Test that direct cycle does not require reference acid."""
        settings = GaussianpKaJobSettings(
            proton_index=10,
            scheme="direct",
            charge=0,
            multiplicity=1,
        )
        assert settings.scheme == "direct"
        assert settings.reference_file is None  # Not needed for direct cycle
        assert settings.delta_G_proton == -265.9  # Default value

    def test_direct_cycle_custom_delta_g(self):
        """Test direct cycle with custom delta_G_proton."""
        settings = GaussianpKaJobSettings(
            proton_index=10,
            scheme="direct",
            delta_G_proton=-270.0,
            charge=0,
            multiplicity=1,
        )
        assert settings.delta_G_proton == -270.0

    def test_proton_exchange_with_reference_file(
        self, single_molecule_xyz_file
    ):
        """Test proton exchange cycle with reference acid file."""
        settings = GaussianpKaJobSettings(
            proton_index=10,
            scheme="proton exchange",
            reference_file=single_molecule_xyz_file,
            reference_proton_index=1,
            reference_charge=0,
            reference_multiplicity=1,
            charge=0,
            multiplicity=1,
            functional="B3LYP",
            basis="6-31G*",
        )
        assert settings.scheme == "proton exchange"
        assert settings.reference_file == single_molecule_xyz_file
        assert settings.reference_proton_index == 1
        assert settings.reference_charge == 0
        assert settings.reference_multiplicity == 1
        assert settings.has_reference_file is True

    def test_proton_exchange_without_reference_file(self):
        """Test proton exchange cycle without reference file."""
        settings = GaussianpKaJobSettings(
            proton_index=10,
            scheme="proton exchange",
            charge=0,
            multiplicity=1,
        )
        assert settings.scheme == "proton exchange"
        assert settings.reference_file is None
        assert settings.has_reference_file is False

    def test_default_scheme(self):
        """Test that default thermodynamic cycle is proton exchange."""
        settings = GaussianpKaJobSettings(
            proton_index=10,
            charge=0,
            multiplicity=1,
        )
        assert settings.scheme == "proton exchange"
        assert settings.has_reference_file is False

    def test_reference_validation(self, single_molecule_xyz_file):
        """Test that reference settings validation works."""
        settings = GaussianpKaJobSettings(
            proton_index=10,
            scheme="proton exchange",
            reference_file=single_molecule_xyz_file,
            reference_proton_index=1,
            reference_charge=0,
            reference_multiplicity=1,
            charge=0,
            multiplicity=1,
        )
        # Should not raise
        settings.validate_reference_settings()

    def test_reference_validation_missing_proton_index(
        self, single_molecule_xyz_file
    ):
        """Test validation fails when reference_proton_index is missing."""
        settings = GaussianpKaJobSettings(
            proton_index=10,
            scheme="proton exchange",
            reference_file=single_molecule_xyz_file,
            reference_charge=0,
            reference_multiplicity=1,
            charge=0,
            multiplicity=1,
        )
        import pytest

        with pytest.raises(ValueError, match="reference_proton_index"):
            settings.validate_reference_settings()

    def test_gas_phase_optimization_settings(self, single_molecule_xyz_file):
        """Test that gas phase optimization has no solvent."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(
            proton_index=proton_index,
            functional="B3LYP",
            basis="6-31G*",
            solvent_model="SMD",
            solvent_id="water",
        )

        prot_settings, conj_base_settings = (
            settings._create_gas_phase_job_settings(mol)
        )

        # Gas phase should have no solvent
        assert prot_settings.solvent_model is None
        assert prot_settings.solvent_id is None
        assert conj_base_settings.solvent_model is None
        assert conj_base_settings.solvent_id is None
        # Should use same functional/basis
        assert prot_settings.functional == "B3LYP"
        assert prot_settings.basis == "6-31G*"

    def test_solution_phase_sp_settings(self, single_molecule_xyz_file):
        """Test that solution phase SP uses same level of theory with solvent."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(
            proton_index=proton_index,
            functional="B3LYP",
            basis="6-31G*",
            solvent_model="SMD",
            solvent_id="water",
        )

        prot_sp_settings, conj_base_sp_settings = (
            settings._create_solution_phase_sp_settings(mol)
        )

        # Solution phase should have solvent
        assert prot_sp_settings.solvent_model == "SMD"
        assert prot_sp_settings.solvent_id == "water"
        assert conj_base_sp_settings.solvent_model == "SMD"
        assert conj_base_sp_settings.solvent_id == "water"
        # Should use SAME functional/basis as gas phase for error cancellation
        assert prot_sp_settings.functional == "B3LYP"
        assert prot_sp_settings.basis == "6-31G*"
        assert conj_base_sp_settings.functional == "B3LYP"
        assert conj_base_sp_settings.basis == "6-31G*"

    def test_protonated_charge_multiplicity_properties(self):
        """Test that protonated_charge/multiplicity are aliases for charge/multiplicity."""
        settings = GaussianpKaJobSettings(
            proton_index=10,
            charge=2,
            multiplicity=3,
        )
        # Properties should return the same values
        assert settings.protonated_charge == settings.charge
        assert settings.protonated_multiplicity == settings.multiplicity

        # Setting via property should update the underlying attribute
        settings.protonated_charge = 5
        assert settings.charge == 5
        settings.protonated_multiplicity = 4
        assert settings.multiplicity == 4

    def test_reference_gas_phase_job_settings(self, single_molecule_xyz_file):
        """Test reference acid gas phase optimization settings."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        ref_proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(
            proton_index=1,
            scheme="proton exchange",
            reference_file=single_molecule_xyz_file,
            reference_proton_index=ref_proton_index,
            reference_charge=0,
            reference_multiplicity=1,
            charge=0,
            multiplicity=1,
            functional="B3LYP",
            basis="6-31G*",
            solvent_model="SMD",
            solvent_id="water",
        )

        ref_acid_settings, ref_cb_settings = (
            settings._create_reference_gas_phase_job_settings()
        )

        # Gas phase should have no solvent
        assert ref_acid_settings.solvent_model is None
        assert ref_acid_settings.solvent_id is None
        assert ref_cb_settings.solvent_model is None
        assert ref_cb_settings.solvent_id is None
        # Should use same functional/basis
        assert ref_acid_settings.functional == "B3LYP"
        assert ref_acid_settings.basis == "6-31G*"
        # Check charge/multiplicity
        assert ref_acid_settings.charge == 0
        assert ref_acid_settings.multiplicity == 1
        assert ref_cb_settings.charge == -1  # Default: reference_charge - 1
        assert ref_cb_settings.multiplicity == 1

    def test_reference_solution_phase_sp_settings(
        self, single_molecule_xyz_file
    ):
        """Test reference acid solution phase SP settings."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        ref_proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(
            proton_index=1,
            scheme="proton exchange",
            reference_file=single_molecule_xyz_file,
            reference_proton_index=ref_proton_index,
            reference_charge=0,
            reference_multiplicity=1,
            charge=0,
            multiplicity=1,
            functional="B3LYP",
            basis="6-31G*",
            solvent_model="SMD",
            solvent_id="water",
        )

        ref_acid_sp_settings, ref_cb_sp_settings = (
            settings._create_reference_solution_phase_sp_settings()
        )

        # Solution phase should have solvent
        assert ref_acid_sp_settings.solvent_model == "SMD"
        assert ref_acid_sp_settings.solvent_id == "water"
        assert ref_cb_sp_settings.solvent_model == "SMD"
        assert ref_cb_sp_settings.solvent_id == "water"
        # Same functional/basis for error cancellation
        assert ref_acid_sp_settings.functional == "B3LYP"
        assert ref_acid_sp_settings.basis == "6-31G*"

    def test_get_reference_molecule(self, single_molecule_xyz_file):
        """Test loading reference molecule from file."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        ref_proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(
            proton_index=1,
            scheme="proton exchange",
            reference_file=single_molecule_xyz_file,
            reference_proton_index=ref_proton_index,
            reference_charge=0,
            reference_multiplicity=1,
            charge=0,
            multiplicity=1,
        )

        ref_mol = settings.get_reference_molecule()
        assert ref_mol is not None
        assert ref_mol.charge == 0
        assert ref_mol.multiplicity == 1
        assert len(ref_mol) == len(mol)

    def test_get_reference_conjugate_base_molecule(
        self, single_molecule_xyz_file
    ):
        """Test creating reference conjugate base by removing proton."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        ref_proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(
            proton_index=1,
            scheme="proton exchange",
            reference_file=single_molecule_xyz_file,
            reference_proton_index=ref_proton_index,
            reference_charge=0,
            reference_multiplicity=1,
            charge=0,
            multiplicity=1,
        )

        ref_cb_mol = settings.get_reference_conjugate_base_molecule()
        assert ref_cb_mol is not None
        assert ref_cb_mol.charge == -1  # Default: reference_charge - 1
        assert ref_cb_mol.multiplicity == 1
        assert len(ref_cb_mol) == len(mol) - 1  # One H removed

    def test_create_conjugate_base_molecule(self, single_molecule_xyz_file):
        """Test creating conjugate base molecule by removing a proton."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1
        original_num_atoms = len(mol)

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        assert (
            len(h_indices) > 0
        ), "Test molecule must have at least one hydrogen"

        proton_index = h_indices[0]
        settings = GaussianpKaJobSettings(proton_index=proton_index)

        conjugate_base = settings._create_conjugate_base_molecule(mol)

        # Check that one atom was removed
        assert len(conjugate_base) == original_num_atoms - 1
        # Check that charge decreased by 1
        assert conjugate_base.charge == -1
        # Check that multiplicity is preserved
        assert conjugate_base.multiplicity == 1

    def test_create_conjugate_base_molecule_custom_charge(
        self, single_molecule_xyz_file
    ):
        """Test creating conjugate base with custom charge/multiplicity."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 1
        mol.multiplicity = 2

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(
            proton_index=proton_index,
            conjugate_base_charge=0,
            conjugate_base_multiplicity=1,
        )

        conjugate_base = settings._create_conjugate_base_molecule(mol)

        # Custom values should override defaults
        assert conjugate_base.charge == 0
        assert conjugate_base.multiplicity == 1

    def test_create_conjugate_base_molecule_no_proton_index(
        self, single_molecule_xyz_file
    ):
        """Test that error is raised when proton_index is not specified."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        settings = GaussianpKaJobSettings()

        with pytest.raises(ValueError, match="proton_index must be specified"):
            settings._create_conjugate_base_molecule(mol)

    def test_create_conjugate_base_molecule_invalid_index(
        self, single_molecule_xyz_file
    ):
        """Test that error is raised for out-of-range proton index."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        settings = GaussianpKaJobSettings(proton_index=999)

        with pytest.raises(ValueError, match="out of range"):
            settings._create_conjugate_base_molecule(mol)

    def test_create_conjugate_base_molecule_not_hydrogen(
        self, single_molecule_xyz_file
    ):
        """Test that error is raised when index is not a hydrogen."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)

        # Find a non-hydrogen atom index
        non_h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s != "H"]
        assert len(non_h_indices) > 0

        settings = GaussianpKaJobSettings(proton_index=non_h_indices[0])

        with pytest.raises(ValueError, match="not hydrogen"):
            settings._create_conjugate_base_molecule(mol)

    def test_create_job_settings(self, single_molecule_xyz_file):
        """Test creating gas phase job settings for both forms."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(
            proton_index=proton_index,
            functional="B3LYP",
            basis="6-31G*",
            solvent_model="SMD",
            solvent_id="water",
        )

        prot_settings, conj_base_settings = (
            settings._create_gas_phase_job_settings(mol)
        )

        # Check protonated settings - GAS PHASE (no solvent)
        assert isinstance(prot_settings, GaussianJobSettings)
        assert prot_settings.charge == 0
        assert prot_settings.multiplicity == 1
        assert prot_settings.functional == "B3LYP"
        assert prot_settings.basis == "6-31G*"
        assert prot_settings.jobtype == "opt"
        assert prot_settings.freq is True
        assert prot_settings.solvent_model is None  # Gas phase
        assert prot_settings.solvent_id is None

        # Check conjugate base settings - GAS PHASE (no solvent)
        assert isinstance(conj_base_settings, GaussianJobSettings)
        assert conj_base_settings.charge == -1
        assert conj_base_settings.multiplicity == 1
        assert conj_base_settings.functional == "B3LYP"
        assert conj_base_settings.basis == "6-31G*"
        assert conj_base_settings.jobtype == "opt"
        assert conj_base_settings.freq is True
        assert conj_base_settings.solvent_model is None  # Gas phase
        assert conj_base_settings.solvent_id is None

    def test_create_molecules(self, single_molecule_xyz_file):
        """Test creating both protonated and conjugate base molecules."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1
        original_num_atoms = len(mol)

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(proton_index=proton_index)

        prot_mol, conj_base_mol = settings._create_molecules(mol)

        # Check protonated molecule
        assert len(prot_mol) == original_num_atoms
        assert prot_mol.charge == 0
        assert prot_mol.multiplicity == 1

        # Check conjugate base molecule
        assert len(conj_base_mol) == original_num_atoms - 1
        assert conj_base_mol.charge == -1
        assert conj_base_mol.multiplicity == 1

    def test_conjugate_base_molecule_method(self, single_molecule_xyz_file):
        """Test the public conjugate_base_molecule method."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(proton_index=proton_index)

        conjugate_base = settings.conjugate_base_molecule(mol)

        assert len(conjugate_base) == len(mol) - 1
        assert conjugate_base.charge == -1

    def test_conjugate_pair_molecules_method(self, single_molecule_xyz_file):
        """Test the public conjugate_pair_molecules method."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(proton_index=proton_index)

        prot_mol, conj_base_mol = settings.conjugate_pair_molecules(mol)

        assert len(prot_mol) == len(mol)
        assert len(conj_base_mol) == len(mol) - 1
        assert prot_mol.charge == 0
        assert conj_base_mol.charge == -1

    def test_conjugate_pair_job_settings_method(
        self, single_molecule_xyz_file
    ):
        """Test the public conjugate_pair_job_settings method returns gas phase settings."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(
            proton_index=proton_index,
            functional="B3LYP",
            basis="6-31G*",
            solvent_model="SMD",
            solvent_id="water",
        )

        prot_settings, conj_base_settings = (
            settings.conjugate_pair_job_settings(mol)
        )

        assert isinstance(prot_settings, GaussianJobSettings)
        assert isinstance(conj_base_settings, GaussianJobSettings)
        assert prot_settings.charge == 0
        assert conj_base_settings.charge == -1
        # Should be gas phase (no solvent for optimization)
        assert prot_settings.solvent_model is None
        assert conj_base_settings.solvent_model is None

    # ------------------------------------------------------------------
    # GaussianpKaJob integration tests
    # ------------------------------------------------------------------
    def test_job_init_valid_settings(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        """Test initialization with valid pKa settings."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(
            proton_index=proton_index,
            functional="B3LYP",
            basis="6-31G*",
        )

        job = GaussianpKaJob(
            molecule=mol,
            settings=settings,
            label="test_pka",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        assert isinstance(job, GaussianpKaJob)
        assert job.TYPE == "g16pka"
        assert job.label == "test_pka"

    def test_init_invalid_settings_type(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        """Test that error is raised for non-pKa settings."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)

        settings = GaussianJobSettings(functional="B3LYP", basis="6-31G*")

        with pytest.raises(
            ValueError, match="must be instance of GaussianpKaJobSettings"
        ):
            GaussianpKaJob(
                molecule=mol,
                settings=settings,
                label="test_pka",
                jobrunner=gaussian_jobrunner_no_scratch,
            )

    def test_init_no_proton_index(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        """Test that error is raised when proton_index is not specified."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)

        settings = GaussianpKaJobSettings(functional="B3LYP", basis="6-31G*")

        with pytest.raises(ValueError, match="proton_index must be specified"):
            GaussianpKaJob(
                molecule=mol,
                settings=settings,
                label="test_pka",
                jobrunner=gaussian_jobrunner_no_scratch,
            )

    def test_pka_jobs_property(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        """Test that pka_jobs returns both jobs."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(
            proton_index=proton_index,
            functional="B3LYP",
            basis="6-31G*",
        )

        job = GaussianpKaJob(
            molecule=mol,
            settings=settings,
            label="test_pka",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        pka_jobs = job.opt_jobs
        assert len(pka_jobs) == 2
        assert isinstance(pka_jobs[0], GaussianOptJob)
        assert isinstance(pka_jobs[1], GaussianOptJob)

    def test_protonated_job_property(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        """Test protonated_job property."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(
            proton_index=proton_index,
            functional="B3LYP",
            basis="6-31G*",
        )

        job = GaussianpKaJob(
            molecule=mol,
            settings=settings,
            label="test_pka",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        protonated_job = job.protonated_job
        assert isinstance(protonated_job, GaussianOptJob)
        assert protonated_job.label == "test_pka_HA_opt"
        assert protonated_job.settings.charge == 0

    def test_conjugate_base_job_property(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        """Test conjugate_base_job property."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(
            proton_index=proton_index,
            functional="B3LYP",
            basis="6-31G*",
        )

        job = GaussianpKaJob(
            molecule=mol,
            settings=settings,
            label="test_pka",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        conjugate_base_job = job.conjugate_base_job
        assert isinstance(conjugate_base_job, GaussianOptJob)
        assert conjugate_base_job.label == "test_pka_A_opt"
        assert conjugate_base_job.settings.charge == -1

    def test_protonated_molecule_property(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        """Test protonated_molecule method on settings."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1
        original_num_atoms = len(mol)

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(proton_index=proton_index)

        protonated_mol = settings.protonated_molecule(mol)
        assert len(protonated_mol) == original_num_atoms
        assert protonated_mol.charge == 0

    def test_conjugate_base_molecule_property(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        """Test conjugate_base_molecule method on settings."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1
        original_num_atoms = len(mol)

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(proton_index=proton_index)

        conjugate_base_mol = settings.conjugate_base_molecule(mol)
        # print(len(conjugate_base_mol))
        # print(conjugate_base_mol)
        # conjugate_base_mol = job.conjugate_base_molecule
        assert len(conjugate_base_mol) == original_num_atoms - 1
        assert conjugate_base_mol.charge == -1

    def test_job_labels(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        """Test that job labels are correctly generated."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(proton_index=proton_index)

        job = GaussianpKaJob(
            molecule=mol,
            settings=settings,
            label="acetic_acid_pka",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        protonated_job, conjugate_base_job = job.opt_jobs

        assert protonated_job.label == "acetic_acid_pka_HA_opt"
        assert conjugate_base_job.label == "acetic_acid_pka_A_opt"


class TestGaussianLinkJobSettingsGuess:
    """Tests for guess= formatting in GaussianLinkJobSettings route strings."""

    _COMMON = dict(
        functional="um062x", basis="def2svp", charge=0, multiplicity=1
    )

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


class TestInitValidation:
    def test_expands_user_home_in_gen_genecp_file(self):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            gen_genecp_file="~/genecp.txt",
        )
        assert "~" not in settings.gen_genecp_file

    def test_forces_and_freq_both_true_raises(self):
        with pytest.raises(ValueError, match="Frequency and Force"):
            GaussianJobSettings(
                functional="b3lyp", basis="sto-3g", forces=True, freq=True
            )


class TestMerge:
    def test_merge_with_none_keywords_merges_everything(self):
        base = GaussianJobSettings(functional="b3lyp", basis="sto-3g")
        other = GaussianJobSettings(functional="wb97xd", basis="def2-tzvp")
        merged = base.merge(other, keywords=None)
        assert merged.functional == "wb97xd"
        assert merged.basis == "def2-tzvp"


class TestDunderMethods:
    def test_getitem_returns_attribute(self):
        settings = GaussianJobSettings(functional="b3lyp", basis="sto-3g")
        assert settings["functional"] == "b3lyp"

    def test_eq_different_types_returns_not_implemented(self):
        settings = GaussianJobSettings(functional="b3lyp", basis="sto-3g")
        assert settings.__eq__(object()) is NotImplemented

    def test_eq_equal_settings(self):
        s1 = GaussianJobSettings(functional="b3lyp", basis="sto-3g")
        s2 = GaussianJobSettings(functional="b3lyp", basis="sto-3g")
        assert s1 == s2

    def test_eq_unequal_settings_logs_diff(self, caplog):
        s1 = GaussianJobSettings(functional="b3lyp", basis="sto-3g")
        s2 = GaussianJobSettings(functional="wb97xd", basis="def2-tzvp")
        with caplog.at_level("INFO"):
            assert not (s1 == s2)
        assert "not equal" in caplog.text


class TestAlternateConstructors:
    def test_from_outfile_merges_orca_settings(self):
        mock_orca_settings = GaussianJobSettings(
            functional="wb97xd", basis="def2-tzvp"
        )
        with patch(
            "chemsmart.io.orca.output.ORCAOutput.read_settings",
            return_value=mock_orca_settings,
        ):
            settings = GaussianJobSettings.from_outfile("fake.out")
        assert settings.functional == "wb97xd"
        assert settings.basis == "def2-tzvp"

    def test_from_logfile_falls_back_to_pbc_on_value_error(self):
        mock_pbc_settings = GaussianJobSettings(
            functional="b3lyp", basis="sto-3g"
        )
        with (
            patch(
                "chemsmart.io.gaussian.output.Gaussian16Output.read_settings",
                side_effect=ValueError("not a valid non-PBC log"),
            ),
            patch(
                "chemsmart.io.gaussian.output.Gaussian16OutputWithPBC.read_settings",
                return_value=mock_pbc_settings,
            ),
        ):
            settings = GaussianJobSettings.from_logfile("fake.log")
        assert settings.functional == "b3lyp"

    def test_from_filepath_inp_extension(self):
        from unittest.mock import MagicMock

        mock_orca_settings = GaussianJobSettings(
            functional="wb97xd", basis="def2-tzvp"
        )
        mock_orca_input = MagicMock()
        mock_orca_input.read_settings.return_value = mock_orca_settings
        with patch(
            "chemsmart.io.orca.input.ORCAInput",
            return_value=mock_orca_input,
        ):
            settings = GaussianJobSettings.from_filepath("fake.inp")
        assert settings.functional == "wb97xd"

    def test_from_filepath_unsupported_extension_raises(self):
        with pytest.raises(ValueError, match="Could not create"):
            GaussianJobSettings.from_filepath("fake.xyz")


class TestRouteStringSetterAndUserInput:
    def test_route_string_setter(self):
        settings = GaussianJobSettings(functional="b3lyp", basis="sto-3g")
        settings.route_string = "# custom route"
        assert settings._route_string == "# custom route"

    def test_user_input_route_without_dieze_tag(self):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            route_to_be_written="opt freq",
        )
        assert settings.route_string == "# opt freq"

    def test_user_input_route_with_dieze_tag(self):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            dieze_tag="p",
            route_to_be_written="opt freq",
        )
        assert settings.route_string == "#p opt freq"

    def test_user_input_route_already_prefixed(self):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            route_to_be_written="#p opt freq",
        )
        assert settings.route_string == "#p opt freq"


class TestGetLightElements:
    def test_returns_none_without_heavy_elements(self):
        settings = GaussianJobSettings(functional="b3lyp", basis="sto-3g")
        assert settings.get_light_elements(_mol()) is None


class TestDiezeTagAndOptFreq:
    def test_custom_dieze_tag(self):
        settings = GaussianJobSettings(
            functional="b3lyp", basis="sto-3g", jobtype="sp", dieze_tag="p"
        )
        assert settings._get_dieze_tag().startswith("#p")

    @pytest.mark.parametrize(
        "jobtype,expected",
        [
            ("opt", "opt=(maxstep=5)"),
            ("modred", "opt=(modredundant,maxstep=5)"),
            ("scan", "opt=(modredundant,maxstep=5)"),
            ("sp", ""),
        ],
    )
    def test_additional_opt_options_by_jobtype(self, jobtype, expected):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype=jobtype,
            additional_opt_options_in_route="maxstep=5",
        )
        route = settings._get_dieze_tag()
        if expected:
            assert expected in route
        if jobtype in ("modred",):
            assert settings.freq is True
        if jobtype in ("scan", "sp"):
            assert settings.freq is False

    def test_ts_additional_opt_options_without_calcall(self):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ts",
            additional_opt_options_in_route="maxstep=5",
        )
        route = settings._get_dieze_tag()
        assert "opt=(ts,calcfc,noeigentest,maxstep=5)" in route

    def test_ts_additional_opt_options_with_calcall(self):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ts",
            additional_opt_options_in_route="calcall",
        )
        route = settings._get_dieze_tag()
        assert "opt=(ts,noeigentest,calcall)" in route

    def test_jobtype_not_matching_any_branch_falls_through(self):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="nci",
            additional_opt_options_in_route=None,
        )
        route = settings._get_dieze_tag()
        assert "opt" not in route

    def test_jobtype_not_matching_any_branch_with_additional_opt_options(
        self,
    ):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="nci",
            additional_opt_options_in_route="maxstep=5",
        )
        route = settings._get_dieze_tag()
        assert "opt" not in route

    def test_numfreq_without_freq(self):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="sp",
            freq=False,
            numfreq=True,
        )
        route = settings._get_dieze_tag()
        assert "freq=numer" in route

    def test_freq_and_numfreq_both_true_raises(self):
        # jobtype="sp" would itself force freq=False as a side effect of
        # _get_dieze_tag(), so use a jobtype that doesn't touch freq/numfreq
        # to let the manually-set combination reach the final check.
        settings = GaussianJobSettings(
            functional="b3lyp", basis="sto-3g", jobtype="nci"
        )
        settings.freq = True
        settings.numfreq = True
        with pytest.raises(ValueError, match="cannot be True at the same"):
            settings._get_dieze_tag()


class TestLevelOfTheoryString:
    def test_semiempirical_with_basis_warns(self, caplog):
        settings = GaussianJobSettings(
            semiempirical="PM6", basis="sto-3g", jobtype="sp"
        )
        with caplog.at_level("WARNING"):
            route = settings._get_level_of_theory_string()
        assert "PM6" in route
        assert "not required" in caplog.text

    def test_semiempirical_without_basis_no_warning(self, caplog):
        settings = GaussianJobSettings(semiempirical="PM6", jobtype="sp")
        with caplog.at_level("WARNING"):
            route = settings._get_level_of_theory_string()
        assert "PM6" in route
        assert "not required" not in caplog.text

    def test_ab_initio_with_basis_succeeds(self):
        settings = GaussianJobSettings(
            ab_initio="MP2", basis="sto-3g", jobtype="sp"
        )
        route = settings._get_level_of_theory_string()
        assert "MP2 sto-3g" in route

    def test_ab_initio_without_basis_raises(self):
        settings = GaussianJobSettings(ab_initio="MP2", jobtype="sp")
        with pytest.raises(ValueError, match="ab initio methods"):
            settings._get_level_of_theory_string()

    def test_functional_without_basis_raises(self):
        settings = GaussianJobSettings(functional="b3lyp", jobtype="sp")
        with pytest.raises(ValueError, match="DFT methods"):
            settings._get_level_of_theory_string()

    def test_both_ab_initio_and_functional_raises(self):
        settings = GaussianJobSettings(
            ab_initio="MP2",
            functional="b3lyp",
            basis="sto-3g",
            jobtype="sp",
        )
        with pytest.raises(ValueError, match="Both ab initio and DFT"):
            settings._get_level_of_theory_string()

    def test_no_method_provided_raises(self):
        settings = GaussianJobSettings(jobtype="sp")
        with pytest.raises(ValueError, match="No computational method"):
            settings._get_level_of_theory_string()

    def test_forces_keyword_added(self):
        settings = GaussianJobSettings(
            functional="b3lyp", basis="sto-3g", jobtype="sp", forces=True
        )
        route = settings._get_level_of_theory_string()
        assert " force" in route

    def test_custom_solvent_without_model_or_id_uses_pcm_default(self):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="sp",
            custom_solvent="SolventName=custom",
        )
        route = settings._get_level_of_theory_string()
        assert "scrf=(pcm,read" in route

    def test_incomplete_solvation_spec_raises(self):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="sp",
            solvent_model="SMD",
        )
        with pytest.raises(ValueError, match="Both solvent model"):
            settings._get_level_of_theory_string()

    def test_nci_jobtype_adds_output_wfn(self):
        settings = GaussianJobSettings(
            functional="b3lyp", basis="sto-3g", jobtype="nci"
        )
        route = settings._get_route_string_from_jobtype()
        assert "output=wfn" in route

    def test_wbi_jobtype_adds_pop_nboread(self):
        settings = GaussianJobSettings(
            functional="b3lyp", basis="sto-3g", jobtype="wbi"
        )
        route = settings._get_route_string_from_jobtype()
        assert "pop=nboread" in route


class TestGenecpSection:
    def test_raises_when_neither_elements_nor_file_specified(self):
        settings = GaussianJobSettings(functional="b3lyp", basis="genecp")
        with pytest.raises(ValueError, match="Could not get GenECPSection"):
            settings.get_genecp_section(_mol())


class TestDetermineBasisKeyword:
    def test_no_heavy_elements_in_structure_without_light_basis_returns_basis(
        self,
    ):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="genecp",
            heavy_elements=["Fe"],
            light_elements_basis=None,
        )
        # Molecule contains no Fe, so no heavy elements are present.
        result = settings.determine_basis_keyword(_mol(symbols=["H", "H"]))
        assert result == "genecp"


class TestCheckSolvent:
    def test_unsupported_solvent_model_raises(self):
        settings = GaussianJobSettings(functional="b3lyp", basis="sto-3g")
        with pytest.raises(ValueError, match="not in"):
            settings._check_solvent("totally_bogus_solvent_model")


class TestBackwardsCompatThermodynamicCycle:
    def test_thermodynamic_cycle_kwarg_maps_to_scheme(self, caplog):
        with caplog.at_level("WARNING"):
            settings = GaussianpKaJobSettings(
                proton_index=2, thermodynamic_cycle="direct"
            )
        assert settings.scheme == "direct"
        assert "deprecated" in caplog.text


class TestTitleAndDelegatingMethods:
    def test_explicit_title_is_not_overridden(self):
        settings = GaussianpKaJobSettings(proton_index=2, title="My Title")
        assert settings.title == "My Title"

    def test_default_title_applied_when_falsy(self):
        settings = GaussianpKaJobSettings(proton_index=2, title=None)
        assert settings.title == "Gaussian pKa calculation job"

    def test_reference_pair_molecules_delegates(self):
        settings = GaussianpKaJobSettings(
            proton_index=2,
            reference_file="ref.xyz",
            reference_proton_index=2,
            reference_charge=0,
            reference_multiplicity=1,
        )
        from unittest.mock import patch

        with patch.object(
            settings, "get_reference_molecule", return_value=_mol_with_h()
        ):
            acid, base = settings.reference_pair_molecules()
        assert acid is not None
        assert base.charge == -1

    def test_reference_pair_job_settings_delegates(self):
        settings = GaussianpKaJobSettings(
            proton_index=2,
            reference_file="ref.xyz",
            reference_proton_index=2,
            reference_charge=0,
            reference_multiplicity=1,
            functional="b3lyp",
            basis="sto-3g",
        )
        acid, base = settings.reference_pair_job_settings()
        assert acid.jobtype == "opt"
        assert base.jobtype == "opt"

    def test_reference_pair_sp_job_settings_delegates(self):
        settings = GaussianpKaJobSettings(
            proton_index=2,
            reference_file="ref.xyz",
            reference_proton_index=2,
            reference_charge=0,
            reference_multiplicity=1,
            functional="b3lyp",
            basis="sto-3g",
        )
        acid, base = settings.reference_pair_sp_job_settings()
        assert acid.jobtype == "sp"
        assert base.jobtype == "sp"

    def test_conjugate_pair_sp_job_settings_delegates(self):
        settings = GaussianpKaJobSettings(
            proton_index=2, functional="b3lyp", basis="sto-3g"
        )
        mol = _mol_with_h(charge=0, multiplicity=1)
        prot, base = settings.conjugate_pair_sp_job_settings(mol)
        assert prot.jobtype == "sp"
        assert base.jobtype == "sp"


class TestValidateReferenceSettings:
    def test_no_reference_file_is_a_no_op(self):
        settings = GaussianpKaJobSettings(proton_index=2)
        settings.validate_reference_settings()  # should not raise

    def test_missing_reference_charge_raises(self):
        settings = GaussianpKaJobSettings(
            proton_index=2,
            reference_file="ref.xyz",
            reference_proton_index=1,
            reference_multiplicity=1,
        )
        with pytest.raises(ValueError, match="reference_charge"):
            settings.validate_reference_settings()

    def test_missing_reference_multiplicity_raises(self):
        settings = GaussianpKaJobSettings(
            proton_index=2,
            reference_file="ref.xyz",
            reference_proton_index=1,
            reference_charge=0,
        )
        with pytest.raises(ValueError, match="reference_multiplicity"):
            settings.validate_reference_settings()


class TestGetReferenceMolecule:
    def test_raises_without_reference_file(self):
        settings = GaussianpKaJobSettings(proton_index=2)
        with pytest.raises(ValueError, match="Reference file not provided"):
            settings.get_reference_molecule()


class TestCreateReferenceConjugateBaseMolecule:
    def _settings(self, **overrides):
        kwargs = dict(
            proton_index=2,
            reference_file="ref.xyz",
            reference_proton_index=2,
            reference_charge=0,
            reference_multiplicity=1,
        )
        kwargs.update(overrides)
        return GaussianpKaJobSettings(**kwargs)

    def test_missing_reference_proton_index_raises(self):
        settings = self._settings(reference_proton_index=None)
        with pytest.raises(
            ValueError, match="reference_proton_index must be specified"
        ):
            settings._create_reference_conjugate_base_molecule(_mol_with_h())

    def test_out_of_range_reference_proton_index_raises(self):
        settings = self._settings(reference_proton_index=999)
        with pytest.raises(ValueError, match="out of range"):
            settings._create_reference_conjugate_base_molecule(_mol_with_h())

    def test_non_hydrogen_reference_proton_index_raises(self):
        settings = self._settings(reference_proton_index=1)  # index 1 = C
        with pytest.raises(ValueError, match="not hydrogen"):
            settings._create_reference_conjugate_base_molecule(_mol_with_h())

    def test_frozen_atoms_are_filtered(self):
        settings = self._settings(reference_proton_index=2)
        mol = _mol_with_h(frozen_atoms=[0, 0, 0, 0, 0])
        result = settings._create_reference_conjugate_base_molecule(mol)
        assert len(result.frozen_atoms) == 4

    def test_explicit_reference_conjugate_base_charge_and_multiplicity(self):
        settings = self._settings(
            reference_proton_index=2,
            reference_conjugate_base_charge=-5,
            reference_conjugate_base_multiplicity=3,
        )
        result = settings._create_reference_conjugate_base_molecule(
            _mol_with_h()
        )
        assert result.charge == -5
        assert result.multiplicity == 3

    def test_default_reference_conjugate_base_charge_and_multiplicity(self):
        settings = self._settings(reference_proton_index=2)
        result = settings._create_reference_conjugate_base_molecule(
            _mol_with_h()
        )
        assert result.charge == -1  # reference_charge - 1
        assert result.multiplicity == 1  # reference_multiplicity


class TestReferenceGasPhaseAndSolutionPhaseSettings:
    def _settings(self, **overrides):
        kwargs = dict(
            proton_index=2,
            reference_file="ref.xyz",
            reference_proton_index=2,
            reference_charge=0,
            reference_multiplicity=1,
            functional="b3lyp",
            basis="sto-3g",
        )
        kwargs.update(overrides)
        return GaussianpKaJobSettings(**kwargs)

    def test_gas_phase_default_conjugate_base_charge_and_mult(self):
        settings = self._settings()
        acid, base = settings._create_reference_gas_phase_job_settings()
        assert acid.charge == 0
        assert base.charge == -1
        assert base.multiplicity == 1

    def test_gas_phase_explicit_conjugate_base_charge_and_mult(self):
        settings = self._settings(
            reference_conjugate_base_charge=-7,
            reference_conjugate_base_multiplicity=2,
        )
        _, base = settings._create_reference_gas_phase_job_settings()
        assert base.charge == -7
        assert base.multiplicity == 2

    def test_solution_phase_default_conjugate_base_charge_and_mult(self):
        settings = self._settings()
        acid, base = settings._create_reference_solution_phase_sp_settings()
        assert acid.charge == 0
        assert base.charge == -1
        assert base.multiplicity == 1

    def test_solution_phase_explicit_conjugate_base_charge_and_mult(self):
        settings = self._settings(
            reference_conjugate_base_charge=-7,
            reference_conjugate_base_multiplicity=2,
        )
        _, base = settings._create_reference_solution_phase_sp_settings()
        assert base.charge == -7
        assert base.multiplicity == 2


class TestProtonatedMolecule:
    def test_explicit_settings_charge_and_multiplicity_used(self):
        settings = GaussianpKaJobSettings(
            proton_index=2, charge=2, multiplicity=3
        )
        mol = _mol_with_h(charge=None, multiplicity=None)
        result = settings.protonated_molecule(mol)
        assert result.charge == 2
        assert result.multiplicity == 3

    def test_falls_back_to_molecule_charge_and_multiplicity(self):
        settings = GaussianpKaJobSettings(proton_index=2)
        settings.charge = None
        settings.multiplicity = None
        mol = _mol_with_h(charge=None, multiplicity=None)
        result = settings.protonated_molecule(mol)
        assert result.charge == 0
        assert result.multiplicity == 1


class TestCreateGasPhaseJobSettings:
    def test_explicit_settings_charge_and_multiplicity_used(self):
        settings = GaussianpKaJobSettings(
            proton_index=2,
            functional="b3lyp",
            basis="sto-3g",
            charge=1,
            multiplicity=2,
        )
        mol = _mol_with_h(charge=None, multiplicity=None)
        prot, _ = settings._create_gas_phase_job_settings(mol)
        assert prot.charge == 1
        assert prot.multiplicity == 2

    def test_explicit_conjugate_base_charge_and_multiplicity(self):
        settings = GaussianpKaJobSettings(
            proton_index=2,
            functional="b3lyp",
            basis="sto-3g",
            charge=0,
            multiplicity=1,
            conjugate_base_charge=-9,
            conjugate_base_multiplicity=4,
        )
        mol = _mol_with_h(charge=None, multiplicity=None)
        _, base = settings._create_gas_phase_job_settings(mol)
        assert base.charge == -9
        assert base.multiplicity == 4


class TestCreateMolecules:
    def test_explicit_settings_charge_and_multiplicity_used(self):
        settings = GaussianpKaJobSettings(
            proton_index=2, charge=3, multiplicity=2
        )
        mol = _mol_with_h(charge=None, multiplicity=None)
        prot, _ = settings._create_molecules(mol)
        assert prot.charge == 3
        assert prot.multiplicity == 2

    def test_falls_back_to_defaults_when_molecule_charge_none(self):
        settings = GaussianpKaJobSettings(proton_index=2)
        settings.charge = None
        settings.multiplicity = None
        mol = _mol_with_h(charge=None, multiplicity=None)
        prot, _ = settings._create_molecules(mol)
        assert prot.charge == 0
        assert prot.multiplicity == 1


class TestCreateSolutionPhaseSpSettings:
    def test_explicit_settings_charge_and_multiplicity_used(self):
        settings = GaussianpKaJobSettings(
            proton_index=2,
            functional="b3lyp",
            basis="sto-3g",
            charge=1,
            multiplicity=2,
        )
        mol = _mol_with_h(charge=None, multiplicity=None)
        prot, _ = settings._create_solution_phase_sp_settings(mol)
        assert prot.charge == 1
        assert prot.multiplicity == 2

    def test_explicit_conjugate_base_charge_and_multiplicity(self):
        settings = GaussianpKaJobSettings(
            proton_index=2,
            functional="b3lyp",
            basis="sto-3g",
            charge=0,
            multiplicity=1,
            conjugate_base_charge=-9,
            conjugate_base_multiplicity=4,
        )
        mol = _mol_with_h(charge=None, multiplicity=None)
        _, base = settings._create_solution_phase_sp_settings(mol)
        assert base.charge == -9
        assert base.multiplicity == 4

    def test_falls_back_to_defaults_when_molecule_charge_none(self):
        settings = GaussianpKaJobSettings(
            proton_index=2, functional="b3lyp", basis="sto-3g"
        )
        settings.charge = None
        settings.multiplicity = None
        mol = _mol_with_h(charge=None, multiplicity=None)
        prot, _ = settings._create_solution_phase_sp_settings(mol)
        assert prot.charge == 0
        assert prot.multiplicity == 1


class TestLegacyChargeMultiplicityFallback:
    def test_canonical_charge_intermediate_used_directly(self):
        settings = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
            charge_intermediate=5,
            mult_intermediate=6,
        )
        assert settings.charge_intermediate == 5
        assert settings.mult_intermediate == 6
        assert settings.int_charge == 5
        assert settings.int_multiplicity == 6


class TestRouteStringDiezeTagAndIrc:
    def test_dieze_tag_prepended(self):
        settings = GaussianQMMMJobSettings(
            high_level_functional="b3lyp",
            high_level_basis="6-31g(d)",
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
            parent_jobtype="sp",
            dieze_tag="p",
        )
        assert settings.route_string.startswith("#p")

    def test_irc_parent_jobtype_adds_irc_keyword(self):
        settings = GaussianQMMMJobSettings(
            high_level_functional="b3lyp",
            high_level_basis="6-31g(d)",
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
            parent_jobtype="irc",
        )
        assert " irc" in settings.route_string


class TestOniomStringBranches:
    def test_medium_only_level_specified(self):
        # medium_level_atoms must be set for the medium layer to count --
        # a medium level of theory with no atoms assigned to it is treated
        # as absent (see _get_oniom_string's `and self.medium_level_atoms`
        # guard).
        settings = GaussianQMMMJobSettings(
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            medium_level_atoms=[2],
            charge_total=0,
            mult_total=1,
        )
        oniom_string = settings._get_oniom_string()
        assert oniom_string == " oniom(b3lyp/6-31g(d))"

    def test_low_only_level_specified(self):
        settings = GaussianQMMMJobSettings(
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
        )
        oniom_string = settings._get_oniom_string()
        assert oniom_string == " oniom(uff)"

    def test_no_levels_specified_bare_oniom(self):
        settings = GaussianQMMMJobSettings(charge_total=0, mult_total=1)
        assert settings._get_oniom_string() == " oniom"


class TestChargeAndMultiplicityErrors:
    def test_two_layer_inconsistent_specification_raises(self):
        # Only model_charge given (not model_multiplicity) for a 2-layer
        # job (no medium level) -- doesn't match any of the "fill from
        # a single reference level" patterns, so it must raise.
        settings = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
            model_charge=-1,
        )
        with pytest.raises(ValueError, match="cannot override"):
            settings.charge_and_multiplicity_string

    def test_three_layer_fills_from_intermediate_when_only_real_and_int_set(
        self,
    ):
        # This is the only reachable "partial fill" branch other than
        # "only real specified": see BUGS_FOUND.md for why the
        # "int_med only" and "model_high only" fill branches are dead
        # code (int_med/int_low and model_high/model_med/model_low are
        # always equal, since they all derive from the same single
        # int_charge/int_multiplicity and model_charge/model_multiplicity
        # constructor parameters respectively).
        # medium_level_atoms must be set, otherwise this is treated as a
        # 2-layer job regardless of int_charge/int_multiplicity (see the
        # `not self.medium_level_atoms` guard in _get_charge_and_multiplicity).
        settings = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            medium_level_atoms=[2],
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
            int_charge=-1,
            int_multiplicity=2,
        )
        assert (
            settings.charge_and_multiplicity_string
            == "0 1 -1 2 -1 2 -1 2 -1 2 -1 2"
        )

    def test_three_layer_inconsistent_specification_raises(self):
        # int_charge set without int_multiplicity, and model_charge set
        # without model_multiplicity: doesn't match any single-level
        # fill pattern, so it must raise.
        settings = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
            int_charge=1,
            model_charge=0,
        )
        with pytest.raises(ValueError, match="cannot override"):
            settings.charge_and_multiplicity_string


class TestEquality:
    def test_different_types_returns_not_implemented(self):
        settings = GaussianQMMMJobSettings(
            high_level_functional="b3lyp",
            high_level_basis="6-31g(d)",
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
        )
        assert settings.__eq__(object()) is NotImplemented

    def test_equal_settings(self):
        kwargs = dict(
            high_level_functional="b3lyp",
            high_level_basis="6-31g(d)",
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
        )
        s1 = GaussianQMMMJobSettings(**kwargs)
        s2 = GaussianQMMMJobSettings(**kwargs)
        assert s1 == s2

    def test_unequal_settings_logs_diff(self, caplog):
        s1 = GaussianQMMMJobSettings(
            high_level_functional="b3lyp",
            high_level_basis="6-31g(d)",
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
        )
        s2 = GaussianQMMMJobSettings(
            high_level_functional="wb97xd",
            high_level_basis="def2-tzvp",
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
        )
        with caplog.at_level("INFO"):
            assert not (s1 == s2)
        assert "not equal" in caplog.text


class TestIRCRouteString:
    def test_jobtype_neither_ircf_nor_ircr_leaves_direction_unset(self):
        settings = GaussianIRCJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="irc",
            predictor="LQA",
            recorrect="never",
            direction="forward",
        )
        settings._get_route_string_from_jobtype()
        # direction is untouched by the ircf/ircr-specific branches
        assert settings.direction == "forward"

    def test_predictor_and_recorrect_both_specified(self):
        settings = GaussianIRCJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ircf",
            predictor="LQA",
            recorrect="never",
        )
        route = settings._get_route_string_from_jobtype()
        assert "irc(LQA,calcfc,recorrect=never" in route
        assert "forward" in route
        assert settings.direction == "forward"

    def test_ircr_sets_reverse_direction(self):
        settings = GaussianIRCJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ircr",
            predictor="LQA",
            recorrect="never",
        )
        route = settings._get_route_string_from_jobtype()
        assert "reverse" in route
        assert settings.direction == "reverse"

    def test_neither_predictor_nor_recorrect_uses_basic_route(self):
        settings = GaussianIRCJobSettings(
            functional="b3lyp", basis="sto-3g", jobtype="ircf"
        )
        route = settings._get_route_string_from_jobtype()
        assert "irc(calcfc,recalc=" in route
        assert "LQA" not in route

    def test_only_predictor_specified_raises(self):
        settings = GaussianIRCJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ircf",
            predictor="LQA",
        )
        with pytest.raises(ValueError, match="Only one of predictor"):
            settings._get_route_string_from_jobtype()

    def test_only_recorrect_specified_raises(self):
        settings = GaussianIRCJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ircf",
            recorrect="never",
        )
        with pytest.raises(ValueError, match="Only one of predictor"):
            settings._get_route_string_from_jobtype()

    def test_flat_irc_applies_defaults(self):
        settings = GaussianIRCJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ircf",
            flat_irc=True,
        )
        route = settings._get_route_string_from_jobtype()
        assert settings.predictor == "LQA"
        assert settings.recorrect == "never"
        assert settings.recalc_step == -5
        assert "irc(LQA,calcfc,recorrect=never" in route

    def test_flat_irc_does_not_override_explicit_recalc_step(self):
        settings = GaussianIRCJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ircf",
            flat_irc=True,
            recalc_step=3,
        )
        settings._get_route_string_from_jobtype()
        assert settings.recalc_step == 3

    def test_additional_route_parameters_appended_once(self):
        settings = GaussianIRCJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ircf",
            additional_route_parameters="scf=qc",
        )
        route = settings._get_route_string_from_jobtype()
        assert route.count("scf=qc") == 1

    def test_additional_route_parameters_not_duplicated_if_present(self):
        # The base class (via super()._get_route_string_from_jobtype())
        # unconditionally appends additional_route_parameters once
        # already; IRC's own duplicate-check at the end of the method
        # then correctly recognizes it's already present and skips
        # re-appending it, so it still only shows up once overall.
        settings = GaussianIRCJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ircf",
            additional_route_parameters="iop(1/8=3)",
        )
        route = settings._get_route_string_from_jobtype()
        assert route.count("iop(1/8=3)") == 1


class TestLinkRouteString:
    def test_link_route_string_uses_custom_link_route(self):
        settings = GaussianLinkJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="opt",
            link_route="#p opt",
        )
        result = settings.link_route_string
        assert "b3lyp" in result
        assert "sto-3g" in result
        assert "geom=check" in result
        assert "guess=read" in result

    def test_link_route_string_skips_already_present_tokens(self):
        settings = GaussianLinkJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="opt",
            link_route="#p opt b3lyp sto-3g geom=check guess=read",
        )
        result = settings.link_route_string
        assert result.count("b3lyp") == 1
        assert result.count("sto-3g") == 1
        assert result.count("geom=check") == 1
        assert result.count("guess=read") == 1

    def test_link_route_string_falls_back_to_jobtype_logic(self):
        settings = GaussianLinkJobSettings(
            functional="b3lyp", basis="sto-3g", jobtype="opt"
        )
        result = settings.link_route_string
        assert "geom=check" in result
        assert "guess=read" in result

    def test_stability_route_removes_opt_and_freq_keywords(self):
        settings = GaussianLinkJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="opt",
            freq=True,
            stable="opt",
            guess="mix",
        )
        route = settings._get_route_string_from_jobtype()
        assert "opt" not in route.replace("stable=opt", "")
        assert "freq" not in route
        assert "stable=opt" in route
        assert "guess=mix" in route

    def test_guess_with_comma_gets_parenthesized(self):
        settings = GaussianLinkJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="opt",
            guess="mix,always",
        )
        route = settings._get_route_string_from_jobtype()
        assert "guess=(mix,always)" in route

    def test_guess_with_surrounding_parens_normalized(self):
        settings = GaussianLinkJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="opt",
            guess="(mix,always)",
        )
        route = settings._get_route_string_from_jobtype()
        assert "guess=(mix,always)" in route
        assert "guess=((mix,always))" not in route

    def test_no_stable_or_guess_omits_both(self):
        settings = GaussianLinkJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="opt",
            stable=None,
            guess=None,
        )
        route = settings._get_route_string_from_jobtype()
        assert "stable=" not in route
        assert "guess=" not in route

    def test_link_route_for_non_irc_jobtype(self):
        settings = GaussianLinkJobSettings(
            functional="b3lyp", basis="sto-3g", jobtype="opt"
        )
        route = settings._get_link_route_string_from_jobtype()
        assert "geom=check" in route
        assert "guess=read" in route

    def test_link_route_for_irc_jobtype_delegates_to_irc_settings(self):
        settings = GaussianLinkJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ircf",
            predictor="LQA",
            recorrect="never",
        )
        route = settings._get_link_route_string_from_jobtype()
        assert "irc(LQA" in route
        assert "geom=check" in route
        assert "guess=read" in route

    def test_link_route_for_irc_does_not_duplicate_existing_tokens(self):
        # If the IRC-generated route already ends up containing
        # "geom=check"/"guess=read" (unlikely in practice, but the guard
        # exists), the method must not duplicate them. We simulate this
        # by checking the guard directly via the additional_route_parameters
        # channel, which is appended by the underlying IRC settings.
        settings = GaussianLinkJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ircf",
            predictor="LQA",
            recorrect="never",
            additional_route_parameters="geom=check guess=read",
        )
        route = settings._get_link_route_string_from_jobtype()
        assert route.count("geom=check") == 1
        assert route.count("guess=read") == 1


class TestBuildGaussianPkaSettings:
    def test_builds_settings_with_defaults(self):
        opt_settings = GaussianJobSettings(
            functional="b3lyp", basis="sto-3g", charge=0, multiplicity=1
        )
        shared = {
            "reference": "ref.log",
            "delta_g_proton": -270.3,
            "reference_color_code": "red",
            "skip_completed": True,
        }
        settings = GaussianpKaJobSettings.build_gaussian_pka_settings(
            proton_index=1, shared=shared, opt_settings=opt_settings
        )
        assert settings.proton_index == 1
        assert settings.reference_file == "ref.log"
        assert settings.delta_G_proton == -270.3
        assert settings.functional == "b3lyp"
        assert settings.basis == "sto-3g"
        # cli_only keys must not leak into the constructed settings kwargs
        assert not hasattr(settings, "reference_color_code")
        # Defaults applied since neither pka_kwargs nor opt/sp settings
        # specify solvent info.
        assert settings.solvent_model == "SMD"
        assert settings.solvent_id == "water"

    def test_solvent_settings_prefer_explicit_shared_values(self):
        opt_settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            solvent_model="PCM",
            solvent_id="dmso",
        )
        shared = {"solvent_model": "SMD", "solvent_id": "acetonitrile"}
        settings = GaussianpKaJobSettings.build_gaussian_pka_settings(
            proton_index=2, shared=shared, opt_settings=opt_settings
        )
        assert settings.solvent_model == "SMD"
        assert settings.solvent_id == "acetonitrile"

    def test_solvent_settings_from_opt_settings_crashes(self):
        # Documents a real bug: when opt_settings carries a non-None
        # solvent_model/solvent_id and `shared` doesn't override it, the
        # resolved value ends up in both `pka_kwargs` (via the
        # solvent_model/solvent_id fallback resolution) and `opt_kwargs`
        # (via the vars(opt_settings) passthrough, since at the time
        # opt_kwargs is built pka_kwargs doesn't contain solvent_model/
        # solvent_id yet). Passing both dicts as **kwargs to the same
        # call then raises "got multiple values for keyword argument".
        # See BUGS_FOUND.md.
        opt_settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            solvent_model="PCM",
            solvent_id="dmso",
        )
        with pytest.raises(TypeError, match="multiple values"):
            GaussianpKaJobSettings.build_gaussian_pka_settings(
                proton_index=2, shared={}, opt_settings=opt_settings
            )

    def test_solvent_settings_fall_back_to_sp_settings(self):
        opt_settings = GaussianJobSettings(functional="b3lyp", basis="sto-3g")
        sp_settings = GaussianJobSettings(
            functional="wb97xd",
            basis="def2-tzvp",
            solvent_model="SMD",
            solvent_id="methanol",
        )
        settings = GaussianpKaJobSettings.build_gaussian_pka_settings(
            proton_index=3,
            shared={},
            opt_settings=opt_settings,
            sp_settings=sp_settings,
        )
        assert settings.solvent_model == "SMD"
        assert settings.solvent_id == "methanol"

    def test_opt_kwargs_do_not_override_explicit_pka_kwargs(self):
        opt_settings = GaussianJobSettings(
            functional="b3lyp", basis="sto-3g", charge=5
        )
        shared = {"charge": 0}
        settings = GaussianpKaJobSettings.build_gaussian_pka_settings(
            proton_index=1, shared=shared, opt_settings=opt_settings
        )
        assert settings.charge == 0


class TestTDDFTRouteString:
    def test_eqsolv_none_omits_eqsolv_option(self):
        from chemsmart.jobs.gaussian.settings import GaussianTDDFTJobSettings

        settings = GaussianTDDFTJobSettings(
            functional="cam-b3lyp",
            basis="def2svp",
            jobtype="sp",
            eqsolv=None,
        )
        route = settings._get_route_string_from_jobtype()
        assert "TD(singlets,nstates=3,root=1)" in route

    def test_valid_eqsolv_option_included(self):
        from chemsmart.jobs.gaussian.settings import GaussianTDDFTJobSettings

        settings = GaussianTDDFTJobSettings(
            functional="cam-b3lyp",
            basis="def2svp",
            jobtype="sp",
            eqsolv="noneqsolv",
        )
        route = settings._get_route_string_from_jobtype()
        assert "TD(singlets,nstates=3,root=1,noneqsolv)" in route

    def test_invalid_eqsolv_option_raises_assertion(self):
        from chemsmart.jobs.gaussian.settings import GaussianTDDFTJobSettings

        settings = GaussianTDDFTJobSettings(
            functional="cam-b3lyp",
            basis="def2svp",
            jobtype="sp",
            eqsolv="bogus",
        )
        with pytest.raises(AssertionError, match="equilibrium solvation"):
            settings._get_route_string_from_jobtype()


class TestBaseCheckSolvent:
    def test_check_solvent_is_a_no_op_on_the_base_class(self):
        # GaussianJobSettings/ORCAJobSettings override _check_solvent with
        # real validation; the base MolecularJobSettings implementation is
        # a deliberate no-op hook for subclasses that don't need it.
        settings = MolecularJobSettings()
        assert settings._check_solvent("anything") is None
        settings.update_solvent(solvent_model="pcm")
        assert settings.solvent_model == "pcm"


class TestCustomSolventValidation:
    def test_rejects_non_string_custom_solvent(self):
        with pytest.raises(ValueError, match="must be a string"):
            GaussianJobSettings(custom_solvent=123)

    def test_custom_solvent_free_text_stored_with_trailing_newline(self):
        settings = GaussianJobSettings(custom_solvent="eps=78.4")
        assert settings.custom_solvent == "eps=78.4\n"

    def test_custom_solvent_free_text_preserves_existing_newline(self):
        settings = GaussianJobSettings(custom_solvent="eps=78.4\n")
        assert settings.custom_solvent == "eps=78.4\n"

    def test_custom_solvent_from_file(self, tmp_path):
        solvent_file = tmp_path / "solvent.txt"
        solvent_file.write_text("eps=78.4\nrsolv=1.0\n")
        settings = GaussianJobSettings(custom_solvent=str(solvent_file))
        assert settings.custom_solvent == "eps=78.4\nrsolv=1.0\n"

    def test_set_custom_solvent_via_file_missing_raises(self):
        settings = GaussianJobSettings()
        with pytest.raises(ValueError, match="does not exist"):
            settings.set_custom_solvent_via_file("/no/such/file.txt")

    def test_custom_solvent_defaults_to_none(self):
        settings = GaussianJobSettings()
        assert settings.custom_solvent is None


class TestSolventHelpers:
    def test_remove_solvent_clears_all_solvent_fields(self):
        settings = GaussianJobSettings(solvent_model="smd", solvent_id="water")
        settings.custom_solvent = "eps=78.4\n"
        settings.remove_solvent()
        assert settings.solvent_model is None
        assert settings.solvent_id is None
        assert settings.custom_solvent is None

    def test_update_solvent_only_updates_given_fields(self):
        settings = GaussianJobSettings(solvent_model="smd", solvent_id="water")
        settings.update_solvent(solvent_id="toluene")
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "toluene"

    def test_update_solvent_no_op_when_nothing_given(self):
        settings = GaussianJobSettings(solvent_model="smd", solvent_id="water")
        settings.update_solvent()
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "water"

    def test_modify_solvent_remove_flag_removes(self):
        settings = GaussianJobSettings(solvent_model="smd", solvent_id="water")
        settings.modify_solvent(remove_solvent=True)
        assert settings.solvent_model is None

    def test_modify_solvent_without_remove_updates(self):
        settings = GaussianJobSettings()
        settings.modify_solvent(
            remove_solvent=False, solvent_model="cpcm", solvent_id="dmso"
        )
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id == "dmso"


class TestFromDict:
    def test_from_dict_constructs_settings(self):
        settings = GaussianJobSettings.from_dict(
            {"functional": "b3lyp", "basis": "6-31g(d)"}
        )
        assert settings.functional == "b3lyp"
        assert settings.basis == "6-31g(d)"


class TestFromDatabase:
    def test_raises_when_file_missing(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="Database file not"):
            GaussianJobSettings.from_database(
                filepath=str(tmp_path / "nope.db")
            )

    def test_structure_id_and_record_selector_mutually_exclusive(
        self, tmp_path
    ):
        db_file = tmp_path / "test.db"
        db_file.write_text("")
        with patch("chemsmart.database.database.Database"):
            with pytest.raises(ValueError, match="not both"):
                GaussianJobSettings.from_database(
                    filepath=str(db_file),
                    structure_id="abc123",
                    record_index=1,
                )

    def test_structure_id_not_found_raises(self, tmp_path):
        db_file = tmp_path / "test.db"
        db_file.write_text("")
        mock_db = MagicMock()
        mock_db.get_structure_by_partial_id.return_value = "full-sid"
        mock_db.get_structure.return_value = None
        with patch(
            "chemsmart.database.database.Database", return_value=mock_db
        ):
            with pytest.raises(ValueError, match="No structure found"):
                GaussianJobSettings.from_database(
                    filepath=str(db_file), structure_id="abc"
                )

    def test_structure_id_populates_charge_and_multiplicity(self, tmp_path):
        db_file = tmp_path / "test.db"
        db_file.write_text("")
        mock_db = MagicMock()
        mock_db.get_structure_by_partial_id.return_value = "full-sid"
        mock_db.get_structure.return_value = {
            "charge": 1,
            "multiplicity": 2,
        }
        with patch(
            "chemsmart.database.database.Database", return_value=mock_db
        ):
            settings = GaussianJobSettings.from_database(
                filepath=str(db_file), structure_id="abc"
            )
        assert settings.charge == 1
        assert settings.multiplicity == 2
        assert "chemsmart database" in settings.title

    def test_record_not_found_returns_none(self, tmp_path):
        db_file = tmp_path / "test.db"
        db_file.write_text("")
        with (
            patch("chemsmart.database.database.Database"),
            patch(
                "chemsmart.database.utils.resolve_record", return_value=None
            ),
        ):
            result = GaussianJobSettings.from_database(
                filepath=str(db_file), record_index=1
            )
        assert result is None

    def test_record_selector_populates_settings_from_meta_and_structure(
        self, tmp_path
    ):
        db_file = tmp_path / "test.db"
        db_file.write_text("")
        fake_record = {
            "meta": {
                "method": "b3lyp",
                "basis": "6-31g(d)",
                "jobtype": "opt",
                "solvent_model": "smd",
                "solvent_id": "water",
                "custom_solvent": None,
            },
            "molecules": [{"charge": 0, "multiplicity": 1}],
        }
        with (
            patch("chemsmart.database.database.Database"),
            patch(
                "chemsmart.database.utils.resolve_record",
                return_value=fake_record,
            ),
        ):
            settings = GaussianJobSettings.from_database(
                filepath=str(db_file), record_index=1, structure_index="-1"
            )

        assert settings.charge == 0
        assert settings.multiplicity == 1
        assert settings.functional == "b3lyp"
        assert settings.basis == "6-31g(d)"
        assert settings.jobtype == "opt"
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "water"

    def test_record_with_slice_structure_index_raises(self, tmp_path):
        db_file = tmp_path / "test.db"
        db_file.write_text("")
        fake_record = {"meta": {}, "molecules": [{}]}
        with (
            patch("chemsmart.database.database.Database"),
            patch(
                "chemsmart.database.utils.resolve_record",
                return_value=fake_record,
            ),
        ):
            with pytest.raises(ValueError, match="one structure at a time"):
                GaussianJobSettings.from_database(
                    filepath=str(db_file),
                    record_index=1,
                    structure_index=":",
                )

    def test_record_structure_index_out_of_range_raises(self, tmp_path):
        db_file = tmp_path / "test.db"
        db_file.write_text("")
        fake_record = {"meta": {}, "molecules": [{"charge": 0}]}
        with (
            patch("chemsmart.database.database.Database"),
            patch(
                "chemsmart.database.utils.resolve_record",
                return_value=fake_record,
            ),
        ):
            with pytest.raises(ValueError, match="out of range"):
                GaussianJobSettings.from_database(
                    filepath=str(db_file),
                    record_index=1,
                    structure_index="5",
                )
