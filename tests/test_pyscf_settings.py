"""Tests for the PySCF project and job-settings contracts."""

import hashlib
import os
import shutil
import sys
from pathlib import Path
from unittest.mock import PropertyMock, patch

import pytest
import yaml

from chemsmart.jobs.pyscf.settings import PySCFJobSettings
from chemsmart.settings.executable import PySCFExecutable
from chemsmart.settings.pyscf import (
    PySCFProjectSettings,
    PySCFProjectSettingsManager,
    YamlPySCFProjectSettings,
)
from chemsmart.settings.user import CHEMSMARTUserSettings


class TestPySCFJobSettings:
    @pytest.mark.parametrize(
        ("multiplicity", "spin", "restricted"),
        [(1, 0, True), (2, 1, False), (3, 2, False)],
    )
    def test_spin_is_two_s_not_multiplicity(
        self, multiplicity, spin, restricted
    ):
        settings = PySCFJobSettings(
            functional="b3lyp",
            basis="def2-svp",
            multiplicity=multiplicity,
        )

        assert settings.spin == spin
        assert settings.is_restricted is restricted

    def test_merge_whitelist_preserves_project_receipt(self):
        project = PySCFJobSettings(
            functional="b3lyp",
            basis="def2-svp",
            charge=0,
            multiplicity=1,
            _project_yaml_digest="a" * 64,
        )
        overrides = PySCFJobSettings(
            functional="pbe0",
            basis="def2-tzvp",
            charge=1,
            multiplicity=2,
        )

        merged = project.merge(overrides, keywords=("charge", "multiplicity"))

        assert merged.functional == "b3lyp"
        assert merged.basis == "def2-svp"
        assert merged.charge == 1
        assert merged.multiplicity == 2
        assert merged.project_yaml_digest == "a" * 64

    def test_unsupported_genecp_request_is_rejected(self):
        settings = PySCFJobSettings(
            functional="b3lyp",
            basis="def2-svp",
            gen_genecp_file="heavy.ecp",
        )

        with pytest.raises(ValueError, match="gen_genecp_file"):
            settings.validate()

    def test_only_hf_is_supported_as_an_ab_initio_method(self):
        settings = PySCFJobSettings(ab_initio="mp2", basis="def2-svp")

        with pytest.raises(ValueError, match="supports only 'hf'"):
            settings.validate()

    def test_ab_initio_and_dft_functional_are_mutually_exclusive(self):
        settings = PySCFJobSettings(
            ab_initio="hf", functional="pbe0", basis="def2-svp"
        )

        with pytest.raises(ValueError, match="either 'ab_initio: hf'"):
            settings.validate()

    @pytest.mark.parametrize(
        ("kwargs", "message"),
        [
            *[
                (
                    {"functional": functional, "basis": "def2-svp"},
                    "perturbative correlation term",
                )
                for functional in ("b2plyp", "pbe0-2", "wb97x-2")
            ],
            (
                {
                    "ab_initio": "hf",
                    "basis": "def2-svp",
                    "defgrid": "defgrid2",
                },
                "DFT-only",
            ),
            (
                {
                    "functional": "pbe0",
                    "basis": "def2-svp",
                    "density_fit": False,
                    "aux_basis": "def2-universal-jkfit",
                },
                "cannot be applied",
            ),
            (
                {
                    "functional": "pbe0",
                    "basis": "def2-svp",
                    "solvent_model": "cpcm",
                },
                "solvent_id is required",
            ),
        ],
    )
    def test_settings_that_would_overclaim_applied_science_are_rejected(
        self, kwargs, message
    ):
        with pytest.raises(ValueError, match=message):
            PySCFJobSettings(**kwargs).validate()

    @pytest.mark.parametrize(
        ("kwargs", "message"),
        [
            ({"basis": "def2-svp"}, "No method specified"),
            ({"functional": "b3lyp"}, "No basis set specified"),
        ],
    )
    def test_method_and_basis_are_required(self, kwargs, message):
        with pytest.raises(ValueError, match=message):
            PySCFJobSettings(**kwargs).validate()


class TestPySCFProjectSettings:
    def test_yaml_maps_gas_and_solv_to_three_supported_jobs(
        self, pyscf_test_project_yaml
    ):
        project = PySCFProjectSettings.from_project("test")

        sp = project.sp_settings()
        opt = project.opt_settings()
        hess = project.hess_settings()

        assert project.PROJECT_NAME == "test"
        assert (sp.jobtype, opt.jobtype, hess.jobtype) == (
            "sp",
            "opt",
            "hess",
        )
        assert (opt.solvent_model, opt.solvent_id) == (None, None)
        assert (hess.solvent_model, hess.solvent_id) == (None, None)
        assert (sp.solvent_model, sp.solvent_id) == ("cpcm", "water")
        assert opt.freq is True
        assert hess.freq is True
        assert sp.freq is False

        expected_digest = hashlib.sha256(
            pyscf_test_project_yaml.read_bytes()
        ).hexdigest()
        assert {
            sp.project_yaml_digest,
            opt.project_yaml_digest,
            hess.project_yaml_digest,
        } == {expected_digest}

    def test_settings_copies_do_not_mutate_project(self):
        project = PySCFProjectSettings.from_project("test")
        first = project.opt_settings()
        first.functional = "pbe0"

        assert project.opt_settings().functional == "b3lyp"


def test_all_packaged_server_templates_declare_pyscf_environment():
    server_dir = (
        Path(__file__).parent.parent
        / "chemsmart"
        / "settings"
        / "templates"
        / ".chemsmart"
        / "server"
    )

    files = sorted(server_dir.glob("*.yaml"))
    assert files
    for filename in files:
        config = yaml.safe_load(filename.read_text())
        assert config["PYSCF"]["LOCAL_RUN"] is True
        assert config["PYSCF"]["SCRATCH"] is False
        assert "conda activate" in config["PYSCF"]["CONDA_ENV"]


@pytest.mark.parametrize(
    "request_field",
    [{"semiempirical": "pm6"}, {"modred": [[1, 2]]}],
)
def test_settings_the_driver_cannot_apply_are_refused(request_field):
    settings = PySCFJobSettings(
        functional="b3lyp", basis="def2-svp", **request_field
    )
    with pytest.raises(ValueError, match="not supported"):
        settings.validate()


def test_project_yaml_digest_survives_merge():
    settings = PySCFJobSettings(functional="b3lyp", basis="def2-svp")
    assert settings.project_yaml_digest is None
    settings.project_yaml_digest = "abc"
    assert settings.merge({}).project_yaml_digest == "abc"


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"ab_initio": "mp2"}, "Unknown ab_initio"),
        ({"ab_initio": "hf"}, "either 'ab_initio: hf' or 'functional'"),
        ({"functional": "b2plyp"}, "Double-hybrid"),
        ({"solvent_id": "water"}, "without a solvent_model"),
        ({"solvent_model": "cosmo-rs", "solvent_id": "w"}, "not available"),
        ({"solvent_model": "cpcm"}, "solvent_id is required"),
        ({"point_group_tol": 0.0}, "must be positive"),
        ({"opt_solver": "bfgs"}, "Unknown opt_solver"),
        ({"engine": "tpu"}, "Unknown engine"),
        ({"defgrid": "defgrid9"}, "Unknown defgrid"),
        ({"density_fit": False, "aux_basis": "x"}, "aux_basis cannot"),
        ({"multiplicity": 0}, "Multiplicity must be"),
    ],
)
def test_impossible_requests_fail_validation(kwargs, message):
    settings = {"functional": "b3lyp", "basis": "def2-svp", **kwargs}
    with pytest.raises(ValueError, match=message):
        PySCFJobSettings(**settings).validate()


def test_hf_cannot_use_a_dft_grid():
    settings = PySCFJobSettings(
        ab_initio="hf", basis="def2-svp", defgrid="defgrid2"
    )
    with pytest.raises(ValueError, match="DFT-only"):
        settings.validate()


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"basis": "def2-svp"}, "No method specified"),
        ({"functional": "b3lyp"}, "No basis set"),
    ],
)
def test_method_and_basis_are_required(kwargs, message):
    with pytest.raises(ValueError, match=message):
        PySCFJobSettings(**kwargs).validate()


def test_method_name_and_spin_of_an_unset_method():
    settings = PySCFJobSettings()
    assert settings.method_name is None
    assert settings.spin is None
    assert settings.xc is None


def test_hf_is_named_and_has_no_functional():
    settings = PySCFJobSettings(ab_initio="HF", basis="def2-svp")
    assert settings.method_name == "hf"
    assert settings.xc is None


def test_merge_all_takes_every_field_of_another_settings_object():
    base = PySCFJobSettings(functional="b3lyp", basis="def2-svp")
    other = PySCFJobSettings(functional="pbe", basis="sto-3g")
    merged = base.merge(other, merge_all=True)
    assert (merged.functional, merged.basis) == ("pbe", "sto-3g")


def test_merge_keywords_limit_what_crosses_over():
    base = PySCFJobSettings(functional="b3lyp", basis="def2-svp")
    other = PySCFJobSettings(functional="pbe", basis="sto-3g")
    merged = base.merge(other, keywords=("basis",))
    assert (merged.functional, merged.basis) == ("b3lyp", "sto-3g")


class TestPySCFProjectResolution:
    def test_unknown_project_lists_where_to_put_one(self):
        with pytest.raises(FileNotFoundError, match="No PySCF project"):
            PySCFProjectSettings.from_project("no_such_project")

    def test_user_project_takes_precedence(
        self, pyscf_yaml_settings_test, tmp_path
    ):
        shutil.copy(pyscf_yaml_settings_test, tmp_path / "mine.yaml")
        with patch.object(
            CHEMSMARTUserSettings,
            "user_pyscf_settings_dir",
            new_callable=PropertyMock,
            return_value=str(tmp_path),
        ):
            project = PySCFProjectSettings.from_project("mine")
        assert project.PROJECT_NAME == "mine"
        assert project.hess_settings().jobtype == "hess"

    def test_yaml_records_its_digest(self, pyscf_yaml_settings_test):
        project = YamlPySCFProjectSettings.from_yaml(pyscf_yaml_settings_test)
        assert len(project.sp_settings().project_yaml_digest) == 64

    def test_missing_jobtype_is_reported(self, pyscf_yaml_settings_test):
        with patch(
            "chemsmart.jobs.settings.read_molecular_job_yaml", return_value={}
        ):
            with pytest.raises(RuntimeError, match="cannot be found"):
                YamlPySCFProjectSettings.from_yaml(pyscf_yaml_settings_test)

    def test_manager_requires_a_filename(self):
        with pytest.raises(ValueError, match="filename is not specified"):
            PySCFProjectSettingsManager(filename=None)

    @pytest.mark.parametrize(
        ("method", "jobtype"),
        [
            ("sp_settings", "sp"),
            ("opt_settings", "opt"),
            ("hess_settings", "hess"),
        ],
    )
    def test_base_project_defines_each_jobtype(self, method, jobtype):
        settings = getattr(PySCFProjectSettings(), method)()
        assert settings.jobtype == jobtype
        assert settings.functional is None


class TestPySCFUserSettings:
    def test_project_files_are_discovered(
        self, pyscf_yaml_settings_test, tmp_path
    ):
        shutil.copy(pyscf_yaml_settings_test, tmp_path / "mine.yaml")
        user = CHEMSMARTUserSettings()
        with patch.object(
            CHEMSMARTUserSettings,
            "user_pyscf_settings_dir",
            new_callable=PropertyMock,
            return_value=str(tmp_path),
        ):
            assert user.all_available_pyscf_projects == ["mine"]


class TestPySCFExecutable:
    def test_configured_folder_names_the_interpreter(self, tmp_path):
        executable = PySCFExecutable(executable_folder=str(tmp_path))
        path = executable.get_executable()
        assert os.path.dirname(path) == str(tmp_path)
        assert os.path.basename(path) in ("python", "python.exe")

    def test_no_folder_uses_the_running_interpreter(self):
        executable = PySCFExecutable(executable_folder=None)
        assert executable.get_executable() == sys.executable
