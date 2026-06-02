
import importlib

import pytest
from click.testing import CliRunner

from chemsmart.cli.gromacs.gromacs import gromacs


em_module = importlib.import_module("chemsmart.cli.gromacs.em")


class DummyGromacsEMJob:
    """
    Dummy job class used to test CLI argument flow without creating a real job
    or running real GROMACS.
    """

    captured_from_settings = {}
    captured_direct_init = {}

    @classmethod
    def from_project_settings(
        cls,
        settings,
        molecule=None,
        jobrunner=None,
        skip_completed=False,
        **kwargs,
    ):
        cls.captured_from_settings = {
            "settings": settings,
            "molecule": molecule,
            "jobrunner": jobrunner,
            "skip_completed": skip_completed,
            "kwargs": kwargs,
        }
        return cls()

    def __init__(
        self,
        molecule=None,
        label=None,
        jobrunner=None,
        mdp_file=None,
        ions_mdp_file=None,
        structure_file=None,
        input_pdb=None,
        top_file=None,
        itp_files=None,
        index_file=None,
        workflow=None,
        skip_completed=False,
        **kwargs,
    ):
        self.__class__.captured_direct_init = {
            "molecule": molecule,
            "label": label,
            "jobrunner": jobrunner,
            "mdp_file": mdp_file,
            "ions_mdp_file": ions_mdp_file,
            "structure_file": structure_file,
            "input_pdb": input_pdb,
            "top_file": top_file,
            "itp_files": itp_files,
            "index_file": index_file,
            "workflow": workflow,
            "skip_completed": skip_completed,
            "kwargs": kwargs,
        }


@pytest.fixture(autouse=True)
def patch_gromacs_em_job(monkeypatch):
    """
    Replace the real GromacsEMJob used by the CLI with a dummy class.

    This keeps the test focused on CLI argument parsing and job creation flow.
    """
    DummyGromacsEMJob.captured_from_settings = {}
    DummyGromacsEMJob.captured_direct_init = {}

    monkeypatch.setattr(
        em_module,
        "GromacsEMJob",
        DummyGromacsEMJob,
    )


@pytest.fixture
def demo_project_yaml(tmp_path):
    """
    Create a minimal prepared GROMACS project YAML and required input files.
    """
    yaml_file = tmp_path / "project.yaml"

    yaml_file.write_text(
        """
project:
  name: prepared_em
  type: gromacs
  job_type: em
  mode: prepared

inputs:
  mdp_file: em.mdp
  structure_file: input.gro
  topology_file: topol.top
""",
        encoding="utf-8",
    )

    for filename in ["em.mdp", "input.gro", "topol.top"]:
        (tmp_path / filename).write_text(
            "dummy content",
            encoding="utf-8",
        )

    return yaml_file


def test_cli_project_yaml_creates_job(demo_project_yaml):
    """
    Test YAML-driven CLI mode.

    Expected flow:
        gromacs -p project.yaml em
        -> GromacsProjectSettings
        -> GromacsEMJob.from_project_settings()
    """
    runner = CliRunner()

    result = runner.invoke(
        gromacs,
        [
            "-p",
            str(demo_project_yaml),
            "em",
        ],
        obj={"molecule": None},
    )

    assert result.exit_code == 0, result.output
    assert result.exception is None

    captured = DummyGromacsEMJob.captured_from_settings
    assert captured

    settings = captured["settings"]

    assert settings.project_name == "prepared_em"
    assert settings.workflow == "prepared"
    assert settings.job_type == "em"

    assert settings.mdp_file == demo_project_yaml.parent / "em.mdp"
    assert settings.structure_file == demo_project_yaml.parent / "input.gro"
    assert settings.top_file == demo_project_yaml.parent / "topol.top"

    assert captured["molecule"] is None
    assert captured["jobrunner"] is None


def test_cli_project_yaml_option_creates_job(demo_project_yaml):
    """
    Test explicit --project-yaml mode.

    This should behave the same as -p project.yaml.
    """
    runner = CliRunner()

    result = runner.invoke(
        gromacs,
        [
            "--project-yaml",
            str(demo_project_yaml),
            "em",
        ],
        obj={"molecule": None},
    )

    assert result.exit_code == 0, result.output
    assert result.exception is None

    captured = DummyGromacsEMJob.captured_from_settings
    assert captured

    settings = captured["settings"]

    assert settings.project_name == "prepared_em"
    assert settings.workflow == "prepared"
    assert settings.job_type == "em"

    assert settings.mdp_file == demo_project_yaml.parent / "em.mdp"
    assert settings.structure_file == demo_project_yaml.parent / "input.gro"
    assert settings.top_file == demo_project_yaml.parent / "topol.top"


def test_cli_direct_options_creates_job(tmp_path):
    """
    Test direct CLI option mode.

    Expected flow:
        gromacs em --mdp em.mdp --structure input.gro --top topol.top
        -> GromacsEMJob(...)
    """
    runner = CliRunner()

    mdp_file = tmp_path / "em.mdp"
    structure_file = tmp_path / "input.gro"
    top_file = tmp_path / "topol.top"

    for file_path in [mdp_file, structure_file, top_file]:
        file_path.write_text(
            "dummy content",
            encoding="utf-8",
        )

    result = runner.invoke(
        gromacs,
        [
            "em",
            "--mdp",
            str(mdp_file),
            "--structure",
            str(structure_file),
            "--top",
            str(top_file),
        ],
        obj={"molecule": None},
    )

    assert result.exit_code == 0, result.output
    assert result.exception is None

    captured = DummyGromacsEMJob.captured_direct_init
    assert captured

    assert captured["label"] == "input_em"
    assert captured["workflow"] == "prepared"

    assert captured["mdp_file"] == mdp_file
    assert captured["structure_file"] == structure_file
    assert captured["top_file"] == top_file

    assert captured["itp_files"] == []
    assert captured["index_file"] is None
    assert captured["molecule"] is None
    assert captured["jobrunner"] is None


def test_cli_direct_options_accepts_itp_files(tmp_path):
    """
    Test that --itp can be provided multiple times.
    """
    runner = CliRunner()

    mdp_file = tmp_path / "em.mdp"
    structure_file = tmp_path / "input.gro"
    top_file = tmp_path / "topol.top"
    ligand_itp = tmp_path / "ligand.itp"
    forcefield_itp = tmp_path / "forcefield.itp"

    for file_path in [
        mdp_file,
        structure_file,
        top_file,
        ligand_itp,
        forcefield_itp,
    ]:
        file_path.write_text(
            "dummy content",
            encoding="utf-8",
        )

    result = runner.invoke(
        gromacs,
        [
            "em",
            "--mdp",
            str(mdp_file),
            "--structure",
            str(structure_file),
            "--top",
            str(top_file),
            "--itp",
            str(ligand_itp),
            "--itp",
            str(forcefield_itp),
        ],
        obj={"molecule": None},
    )

    assert result.exit_code == 0, result.output
    assert result.exception is None

    captured = DummyGromacsEMJob.captured_direct_init
    assert captured

    assert captured["itp_files"] == [
        ligand_itp,
        forcefield_itp,
    ]


def test_cli_rejects_missing_project_yaml(tmp_path):
    """
    Test that an invalid -p path fails early.
    """
    runner = CliRunner()

    missing_yaml = tmp_path / "missing_project.yaml"

    result = runner.invoke(
        gromacs,
        [
            "-p",
            str(missing_yaml),
            "em",
        ],
        obj={"molecule": None},
    )

    assert result.exit_code != 0
    assert result.exception is not None