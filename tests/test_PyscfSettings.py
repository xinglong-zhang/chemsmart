"""Tests for the PySCF project and job-settings contracts."""

import hashlib
from pathlib import Path

import pytest
import yaml

from chemsmart.jobs.pyscf.settings import PySCFJobSettings
from chemsmart.settings.pyscf import PySCFProjectSettings


@pytest.fixture()
def pyscf_test_project_yaml(pyscf_yaml_settings_test):
    return Path(pyscf_yaml_settings_test)


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
        assert (sp.solvent_model, sp.solvent_id) == ("smd", "water")
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
