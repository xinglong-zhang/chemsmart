from pathlib import Path

import pytest

from chemsmart.jobs.gromacs.job import GromacsEMJob, GromacsJob
from chemsmart.settings.gromacs import GromacsProjectSettings


def test_gromacs_project_settings_from_dict_converts_paths():
    settings = GromacsProjectSettings.from_dict(
        {
            "project_name": "prepared_em",
            "workflow": "prepared",
            "mdp_file": "em.mdp",
            "structure_file": "input.gro",
            "top_file": "topol.top",
            "tpr_file": "em.tpr",
            "index_file": "index.ndx",
            "itp_files": ["forcefield.itp", "ligand.itp"],
            "temperature": 300.0,
            "timestep": 0.002,
            "thermostat": "V-rescale",
            "constraint_algorithm": "LINCS",
        }
    )

    assert settings.project_name == "prepared_em"
    assert settings.workflow == "prepared"

    assert settings.mdp_file == Path("em.mdp")
    assert settings.structure_file == Path("input.gro")
    assert settings.top_file == Path("topol.top")
    assert settings.tpr_file == Path("em.tpr")
    assert settings.index_file == Path("index.ndx")
    assert settings.itp_files == [
        Path("forcefield.itp"),
        Path("ligand.itp"),
    ]

    assert settings.temperature == 300.0
    assert settings.timestep == 0.002
    assert settings.thermostat == "V-rescale"
    assert settings.constraint_algorithm == "LINCS"


def test_gromacs_project_settings_validate_prepared_inputs_passes():
    settings = GromacsProjectSettings.from_dict(
        {
            "mdp_file": "em.mdp",
            "structure_file": "input.gro",
            "top_file": "topol.top",
        }
    )

    settings.validate_prepared_inputs()


def test_gromacs_project_settings_validate_prepared_inputs_raises():
    settings = GromacsProjectSettings.from_dict(
        {
            "mdp_file": "em.mdp",
            "structure_file": "input.gro",
        }
    )

    with pytest.raises(
        ValueError,
        match="Missing required GROMACS prepared input settings",
    ):
        settings.validate_prepared_inputs()


def test_gromacs_project_settings_to_job_kwargs():
    settings = GromacsProjectSettings.from_dict(
        {
            "workflow": "prepared",
            "mdp_file": "em.mdp",
            "structure_file": "input.gro",
            "top_file": "topol.top",
            "tpr_file": "em.tpr",
            "index_file": "index.ndx",
            "itp_files": ["forcefield.itp"],
        }
    )

    job_kwargs = settings.to_job_kwargs()

    assert job_kwargs["mdp_file"] == Path("em.mdp")
    assert job_kwargs["structure_file"] == Path("input.gro")
    assert job_kwargs["top_file"] == Path("topol.top")
    assert job_kwargs["tpr_file"] == Path("em.tpr")
    assert job_kwargs["index_file"] == Path("index.ndx")
    assert job_kwargs["itp_files"] == [Path("forcefield.itp")]
    assert job_kwargs["workflow"] == "prepared"


def test_gromacs_project_settings_defaults_to_prepared_workflow():
    settings = GromacsProjectSettings.from_dict(
        {
            "mdp_file": "em.mdp",
            "structure_file": "input.gro",
            "top_file": "topol.top",
        }
    )

    assert settings.workflow == "prepared"


def test_gromacs_project_settings_accepts_project_level_parameters():
    settings = GromacsProjectSettings.from_dict(
        {
            "mdp_file": "em.mdp",
            "structure_file": "input.gro",
            "top_file": "topol.top",
            "force_field": "amber99sb-ildn",
            "water_model": "tip3p",
            "temperature": 300.0,
            "pressure": 1.0,
            "timestep": 0.002,
            "thermostat": "V-rescale",
            "barostat": "Parrinello-Rahman",
            "constraints": "h-bonds",
            "constraint_algorithm": "LINCS",
        }
    )

    assert settings.force_field == "amber99sb-ildn"
    assert settings.water_model == "tip3p"
    assert settings.temperature == 300.0
    assert settings.pressure == 1.0
    assert settings.timestep == 0.002
    assert settings.thermostat == "V-rescale"
    assert settings.barostat == "Parrinello-Rahman"
    assert settings.constraints == "h-bonds"
    assert settings.constraint_algorithm == "LINCS"


def test_gromacs_project_settings_validate_full_setup_inputs_passes():
    settings = GromacsProjectSettings.from_dict(
        {
            "workflow": "full_setup",
            "input_pdb": "input.pdb",
            "mdp_file": "em.mdp",
            "force_field": "amber99sb-ildn",
            "water_model": "tip3p",
        }
    )

    settings.validate_full_setup_inputs()


def test_gromacs_project_settings_validate_full_setup_inputs_raises():
    settings = GromacsProjectSettings.from_dict(
        {
            "workflow": "full_setup",
            "input_pdb": "input.pdb",
            "mdp_file": "em.mdp",
        }
    )

    with pytest.raises(
        ValueError,
        match="Missing required GROMACS full_setup settings",
    ):
        settings.validate_full_setup_inputs()


def test_gromacs_job_can_be_created_from_project_settings():
    settings = GromacsProjectSettings.from_dict(
        {
            "project_name": "prepared_em",
            "workflow": "prepared",
            "mdp_file": "em.mdp",
            "structure_file": "input.gro",
            "top_file": "topol.top",
            "tpr_file": "em.tpr",
            "index_file": "index.ndx",
            "itp_files": ["forcefield.itp"],
        }
    )

    job = GromacsJob.from_project_settings(
        settings=settings,
        molecule=None,
        jobrunner=None,
    )

    assert job.label == "prepared_em"
    assert job.workflow == "prepared"
    assert job.mdp_file == Path("em.mdp")
    assert job.structure_file == Path("input.gro")
    assert job.top_file == Path("topol.top")
    assert job.tpr_file == Path("em.tpr")
    assert job.index_file == Path("index.ndx")
    assert job.itp_files == [Path("forcefield.itp")]


def test_gromacs_em_job_can_be_created_from_project_settings():
    settings = GromacsProjectSettings.from_dict(
        {
            "project_name": "prepared_em",
            "workflow": "prepared",
            "mdp_file": "em.mdp",
            "structure_file": "input.gro",
            "top_file": "topol.top",
            "tpr_file": "em.tpr",
            "index_file": "index.ndx",
            "itp_files": ["forcefield.itp"],
        }
    )

    job = GromacsEMJob.from_project_settings(
        settings=settings,
        molecule=None,
        jobrunner=None,
    )

    assert isinstance(job, GromacsEMJob)
    assert job.label == "prepared_em"
    assert job.workflow == "prepared"
    assert job.mdp_file == Path("em.mdp")
    assert job.structure_file == Path("input.gro")
    assert job.top_file == Path("topol.top")
    assert job.tpr_file == Path("em.tpr")
    assert job.index_file == Path("index.ndx")
    assert job.itp_files == [Path("forcefield.itp")]


def test_gromacs_project_settings_from_yaml(tmp_path):
    yaml_file = tmp_path / "project.yaml"
    yaml_file.write_text(
        """
project:
  name: prepared_em
  type: gromacs
  mode: prepared

gromacs_settings:
  force_field: user_provided
  water_model: user_provided
  timestep: 0.002
  temperature: 300.0
  pressure: 1.0
  thermostat: V-rescale
  barostat: Parrinello-Rahman
  constraints: h-bonds
  constraint_algorithm: LINCS

inputs:
  mdp_file: em.mdp
  structure_file: input.gro
  topology_file: topol.top
  tpr_file: em.tpr
  index_file: index.ndx
  itp_files:
    - forcefield.itp
    - ligand.itp
""",
        encoding="utf-8",
    )

    settings = GromacsProjectSettings.from_yaml(yaml_file)

    assert settings.project_name == "prepared_em"
    assert settings.workflow == "prepared"

    assert settings.force_field == "user_provided"
    assert settings.water_model == "user_provided"
    assert settings.timestep == 0.002
    assert settings.temperature == 300.0
    assert settings.pressure == 1.0
    assert settings.thermostat == "V-rescale"
    assert settings.barostat == "Parrinello-Rahman"
    assert settings.constraints == "h-bonds"
    assert settings.constraint_algorithm == "LINCS"

    assert settings.project_dir == tmp_path
    assert settings.mdp_file == tmp_path / "em.mdp"
    assert settings.structure_file == tmp_path / "input.gro"
    assert settings.top_file == tmp_path / "topol.top"
    assert settings.tpr_file == tmp_path / "em.tpr"
    assert settings.index_file == tmp_path / "index.ndx"
    assert settings.itp_files == [
        tmp_path / "forcefield.itp",
        tmp_path / "ligand.itp",
    ]


def test_gromacs_em_job_can_be_created_from_yaml_settings(tmp_path):
    yaml_file = tmp_path / "project.yaml"
    yaml_file.write_text(
        """
project:
  name: prepared_em
  type: gromacs
  mode: prepared

inputs:
  mdp_file: em.mdp
  structure_file: input.gro
  topology_file: topol.top
  tpr_file: em.tpr
  index_file: index.ndx
  itp_files:
    - forcefield.itp
""",
        encoding="utf-8",
    )

    settings = GromacsProjectSettings.from_yaml(yaml_file)

    job = GromacsEMJob.from_project_settings(
        settings=settings,
        molecule=None,
        jobrunner=None,
    )

    assert job.label == "prepared_em"
    assert job.workflow == "prepared"
    assert job.mdp_file == tmp_path / "em.mdp"
    assert job.structure_file == tmp_path / "input.gro"
    assert job.top_file == tmp_path / "topol.top"
    assert job.tpr_file == tmp_path / "em.tpr"
    assert job.index_file == tmp_path / "index.ndx"
    assert job.itp_files == [tmp_path / "forcefield.itp"]


def test_gromacs_project_settings_keeps_absolute_paths(tmp_path):
    mdp_file = tmp_path / "em.mdp"

    settings = GromacsProjectSettings.from_dict(
        {
            "project_dir": tmp_path / "project",
            "mdp_file": mdp_file,
            "structure_file": "input.gro",
            "top_file": "topol.top",
        }
    )

    assert settings.mdp_file == mdp_file
    assert settings.structure_file == tmp_path / "project" / "input.gro"
    assert settings.top_file == tmp_path / "project" / "topol.top"


def test_gromacs_project_settings_with_overrides_ignores_none(tmp_path):
    yaml_file = tmp_path / "project.yaml"
    yaml_file.write_text(
        """
project:
  name: prepared_em
  mode: prepared

inputs:
  mdp_file: em.mdp
  structure_file: input.gro
  topology_file: topol.top
""",
        encoding="utf-8",
    )

    settings = GromacsProjectSettings.from_yaml(yaml_file)
    updated = settings.with_overrides(
        mdp_file=None,
        workflow="full_setup",
        input_pdb="input.pdb",
        force_field="amber99sb-ildn",
        water_model="tip3p",
    )

    assert updated.mdp_file == tmp_path / "em.mdp"
    assert updated.workflow == "full_setup"
    assert updated.input_pdb == tmp_path / "input.pdb"
    assert updated.force_field == "amber99sb-ildn"
    assert updated.water_model == "tip3p"
