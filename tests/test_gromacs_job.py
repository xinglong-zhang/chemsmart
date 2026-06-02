# tests/test_gromacs_job.py

from pathlib import Path
import pytest

from chemsmart.jobs.gromacs.job import GromacsEMJob, GromacsJob
from chemsmart.settings.gromacs import GromacsProjectSettings

def test_gromacs_job_can_be_created_from_project_settings(tmp_path):
    # Prepare virtual project settings
    settings = GromacsProjectSettings.from_dict({
        "project_name": "prepared_em",
        "workflow": "prepared",
        "mdp_file": tmp_path / "em.mdp",
        "structure_file": tmp_path / "input.gro",
        "top_file": tmp_path / "topol.top",
        "tpr_file": tmp_path / "em.tpr",
        "index_file": tmp_path / "index.ndx",
        "itp_files": [tmp_path / "forcefield.itp"],
    })

    # Create a Job through settings
    job = GromacsEMJob.from_project_settings(
        settings=settings,
        molecule=None,
        jobrunner=None
    )

    assert isinstance(job, GromacsEMJob)
    assert job.label == "prepared_em"
    assert job.workflow == "prepared"
    assert job.mdp_file == tmp_path / "em.mdp"
    assert job.structure_file == tmp_path / "input.gro"
    assert job.top_file == tmp_path / "topol.top"
    assert job.tpr_file == tmp_path / "em.tpr"
    assert job.index_file == tmp_path / "index.ndx"
    assert job.itp_files == [tmp_path / "forcefield.itp"]

def test_gromacs_job_default_tpr(tmp_path):
    settings = GromacsProjectSettings.from_dict({
        "project_name": "prepared_em",
        "workflow": "prepared",
        "mdp_file": tmp_path / "em.mdp",
        "structure_file": tmp_path / "input.gro",
        "top_file": tmp_path / "topol.top",
    })

    job = GromacsEMJob.from_project_settings(
        settings=settings,
        molecule=None,
        jobrunner=None
    )

    # Set folder
    job.set_folder(tmp_path)

    # Automatic tpr file generation
    expected_tpr = tmp_path / f"{job.label}.tpr"
    assert job.tpr_file == expected_tpr

def test_gromacs_job_has_required_prepared_inputs(tmp_path):
    settings = GromacsProjectSettings.from_dict({
        "mdp_file": tmp_path / "em.mdp",
        "structure_file": tmp_path / "input.gro",
        "top_file": tmp_path / "topol.top",
    })

    job = GromacsEMJob.from_project_settings(
        settings=settings,
        molecule=None,
        jobrunner=None
    )

    job.set_folder(tmp_path)

    # Return ture when all necessary files exist
    assert job.has_required_prepared_inputs() is False
    # Because the file does not exist on the disk

    # Create a file simulation existence
    for file_attr in ["mdp_file", "structure_file", "top_file"]:
        f = getattr(job, file_attr)
        f.write_text("dummy content", encoding="utf-8")

    assert job.has_required_prepared_inputs() is True