"""Tests for the PySCF project and job-settings contracts."""

import hashlib
from pathlib import Path

import pytest
import yaml

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.pyscf.singlepoint import PySCFSinglePointJob
from chemsmart.jobs.pyscf.settings import PySCFJobSettings
from chemsmart.settings.pyscf import PySCFProjectSettings


@pytest.fixture()
def pyscf_test_project_yaml():
    return (
        Path(__file__).parent
        / "data"
        / "PySCFTests"
        / "project_yaml"
        / "test.yaml"
    )


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
            PySCFJobSettings(jobtype="sp", **kwargs).validate()

    @pytest.mark.parametrize(
        ("kwargs", "message"),
        [
            ({"basis": "def2-svp"}, "No method specified"),
            ({"functional": "b3lyp"}, "No basis set specified"),
        ],
    )
    def test_method_and_basis_are_required(self, kwargs, message):
        with pytest.raises(ValueError, match=message):
            PySCFJobSettings(jobtype="sp", **kwargs).validate()

    @pytest.mark.parametrize(
        ("kwargs", "message"),
        [
            ({"jobtype": "ts"}, "jobtype"),
            ({"density_fit": 1}, "strict boolean"),
            ({"freq": 1}, "strict boolean"),
            ({"scf_tol": 0.0}, "finite and > 0"),
            ({"scf_tol": float("nan")}, "finite and > 0"),
            ({"scf_maxiter": 0}, "positive integer"),
            ({"scf_maxiter": 1.5}, "positive integer"),
            ({"opt_maxsteps": True}, "positive integer"),
        ],
    )
    def test_scalar_boolean_and_jobtype_contracts(self, kwargs, message):
        values = {
            "functional": "b3lyp",
            "basis": "def2-svp",
            "jobtype": "sp",
        }
        values.update(kwargs)

        with pytest.raises(ValueError, match=message):
            PySCFJobSettings(**values).validate()

    def test_sp_rejects_frequency_flag(self):
        with pytest.raises(ValueError, match="sp rejects freq=True"):
            PySCFJobSettings(
                functional="b3lyp",
                basis="def2-svp",
                jobtype="sp",
                freq=True,
            ).validate()

    @pytest.mark.parametrize("jobtype", ["sp", "opt", "hess"])
    def test_gpu_settings_cover_sp_opt_and_hess(self, jobtype):
        settings = PySCFJobSettings(
            functional="b3lyp",
            basis="def2-svp",
            jobtype=jobtype,
            freq=jobtype == "hess",
            engine="gpu",
        )

        assert settings.validate() is settings

    def test_job_preserves_source_state_unless_settings_are_explicit(self):
        source = Molecule(
            symbols=["O", "H", "H"],
            positions=[
                [0.0, 0.0, 0.0],
                [0.8, 0.0, 0.5],
                [-0.8, 0.0, 0.5],
            ],
            charge=-1,
            multiplicity=2,
        )
        inherited = PySCFSinglePointJob(
            molecule=source,
            settings=PySCFJobSettings(
                functional="b3lyp",
                basis="def2-svp",
                jobtype="sp",
            ),
            label="water_anion",
        )
        explicit = PySCFSinglePointJob(
            molecule=source,
            settings=PySCFJobSettings(
                functional="b3lyp",
                basis="def2-svp",
                jobtype="sp",
                charge=0,
                multiplicity=1,
            ),
            label="water_neutral",
        )

        assert (inherited.settings.charge, inherited.settings.multiplicity) == (
            -1,
            2,
        )
        assert (explicit.settings.charge, explicit.settings.multiplicity) == (
            0,
            1,
        )
        assert (source.charge, source.multiplicity) == (-1, 2)


class TestPySCFProjectSettings:
    def test_yaml_loads_three_explicit_stage_sections(
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

    def test_gas_solv_migration_input_resolves_to_stage_settings(self, tmp_path):
        path = tmp_path / "migration.yaml"
        path.write_text(
            yaml.safe_dump(
                {
                    "gas": {
                        "functional": "b3lyp",
                        "basis": "def2-svp",
                        "freq": False,
                    },
                    "solv": {
                        "functional": "b3lyp",
                        "basis": "def2-svp",
                        "solvent_model": "smd",
                        "solvent_id": "water",
                    },
                }
            ),
            encoding="utf-8",
        )

        project = PySCFProjectSettings.from_yaml(path)

        assert project.sp_settings().solvent_id == "water"
        assert project.opt_settings().freq is False
        assert project.hess_settings().freq is True
        canonical = yaml.safe_load(project.render_canonical_yaml())
        assert tuple(canonical) == ("hess", "opt", "sp")
        assert "gas" not in canonical
        assert "solv" not in canonical
        assert canonical["hess"]["freq"] is True

    def test_only_invoked_stage_is_materialized(self, tmp_path):
        path = tmp_path / "lazy.yaml"
        path.write_text(
            yaml.safe_dump(
                {
                    "sp": {"functional": "b3lyp", "basis": "def2-svp"},
                    "opt": {
                        "functional": "b3lyp",
                        "basis": "def2-svp",
                        "invented_option": "future-stage-error",
                    },
                }
            ),
            encoding="utf-8",
        )

        project = PySCFProjectSettings.from_yaml(path)

        assert project.sp_settings().jobtype == "sp"
        with pytest.raises(ValueError, match="invented_option"):
            project.opt_settings()

    def test_inconsistent_duplicate_stage_and_migration_definition_fails(
        self, tmp_path
    ):
        path = tmp_path / "conflict.yaml"
        path.write_text(
            yaml.safe_dump(
                {
                    "sp": {"functional": "pbe0", "basis": "def2-svp"},
                    "solv": {"functional": "b3lyp", "basis": "def2-svp"},
                }
            ),
            encoding="utf-8",
        )

        project = PySCFProjectSettings.from_yaml(path)

        with pytest.raises(ValueError, match="Conflicting PySCF definitions"):
            project.sp_settings()

    def test_consistent_duplicate_definition_is_accepted(self, tmp_path):
        settings = {"functional": "b3lyp", "basis": "def2-svp"}
        path = tmp_path / "consistent.yaml"
        path.write_text(
            yaml.safe_dump({"sp": settings, "solv": settings}),
            encoding="utf-8",
        )

        project = PySCFProjectSettings.from_yaml(path)

        assert project.sp_settings().functional == "b3lyp"

    def test_unknown_top_level_section_is_rejected(self, tmp_path):
        path = tmp_path / "unknown-section.yaml"
        path.write_text(
            yaml.safe_dump(
                {
                    "sp": {"functional": "b3lyp", "basis": "def2-svp"},
                    "native_input": {"route": "forbidden"},
                }
            ),
            encoding="utf-8",
        )

        with pytest.raises(ValueError, match="native_input"):
            PySCFProjectSettings.from_yaml(path)

    def test_explicit_yaml_filepath_is_resolved_as_exact_project(
        self, tmp_path, pyscf_test_project_yaml
    ):
        path = tmp_path / "exact-project.yml"
        path.write_bytes(pyscf_test_project_yaml.read_bytes())

        project = PySCFProjectSettings.from_project(path)

        assert project.PROJECT_NAME == "exact-project"
        assert project.sp_settings().project_yaml_digest == hashlib.sha256(
            path.read_bytes()
        ).hexdigest()

    def test_missing_named_project_fails_closed(self, tmp_path, monkeypatch):
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(tmp_path))

        with pytest.raises(FileNotFoundError, match="not-a-project"):
            PySCFProjectSettings.from_project("not-a-project")


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
