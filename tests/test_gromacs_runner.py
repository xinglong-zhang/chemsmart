import pytest

from chemsmart.jobs.gromacs.job import GromacsEMJob
from chemsmart.jobs.gromacs.runner import GromacsJobRunner


def test_gromacs_em_job_type_matches_runner_jobtypes():
    assert GromacsEMJob.TYPE in GromacsJobRunner.JOBTYPES


def test_gromacs_em_job_stores_file_attributes(tmp_path):
    mdp_file = tmp_path / "em.mdp"
    structure_file = tmp_path / "input.gro"
    top_file = tmp_path / "topol.top"
    tpr_file = tmp_path / "em.tpr"

    mdp_file.write_text("integrator = steep\n", encoding="utf-8")
    structure_file.write_text("dummy structure\n", encoding="utf-8")
    top_file.write_text("dummy topology\n", encoding="utf-8")

    job = GromacsEMJob(
        molecule=None,
        label="em",
        jobrunner=None,
        mdp_file=mdp_file,
        structure_file=structure_file,
        top_file=top_file,
    )

    job.set_folder(str(tmp_path))

    assert job.mdp_file == mdp_file
    assert job.structure_file == structure_file
    assert job.top_file == top_file
    assert job.tpr_file == tpr_file
    assert job.has_required_prepared_inputs()


def test_gromacs_runner_get_command_uses_tpr_file(tmp_path):
    mdp_file = tmp_path / "em.mdp"
    structure_file = tmp_path / "input.gro"
    top_file = tmp_path / "topol.top"
    tpr_file = tmp_path / "em.tpr"

    mdp_file.write_text("integrator = steep\n", encoding="utf-8")
    structure_file.write_text("dummy structure\n", encoding="utf-8")
    top_file.write_text("dummy topology\n", encoding="utf-8")
    tpr_file.write_text("dummy tpr\n", encoding="utf-8")

    job = GromacsEMJob(
        molecule=None,
        label="em",
        jobrunner=None,
        mdp_file=mdp_file,
        structure_file=structure_file,
        top_file=top_file,
        tpr_file=tmp_path / "em.tpr",
    )

    job.set_folder(str(tmp_path))

    runner = GromacsJobRunner.__new__(GromacsJobRunner)
    command = runner._get_command(job)

    assert command == f"gmx mdrun -deffnm {tmp_path / 'em'}"


def test_gromacs_em_job_keeps_user_provided_tpr_file(tmp_path):
    mdp_file = tmp_path / "em.mdp"
    structure_file = tmp_path / "input.gro"
    top_file = tmp_path / "topol.top"
    custom_tpr_file = tmp_path / "custom_em.tpr"

    mdp_file.write_text("integrator = steep\n", encoding="utf-8")
    structure_file.write_text("dummy structure\n", encoding="utf-8")
    top_file.write_text("dummy topology\n", encoding="utf-8")
    custom_tpr_file.write_text("dummy tpr\n", encoding="utf-8")

    job = GromacsEMJob(
        molecule=None,
        label="em",
        jobrunner=None,
        mdp_file=mdp_file,
        structure_file=structure_file,
        top_file=top_file,
        tpr_file=custom_tpr_file,
    )

    new_folder = tmp_path / "new_folder"
    new_folder.mkdir()
    job.set_folder(str(new_folder))

    assert job.tpr_file == custom_tpr_file


def test_gromacs_runner_validate_inputs_passes_with_required_files(tmp_path):
    mdp_file = tmp_path / "em.mdp"
    structure_file = tmp_path / "input.gro"
    top_file = tmp_path / "topol.top"

    mdp_file.write_text("integrator = steep\n", encoding="utf-8")
    structure_file.write_text("dummy structure\n", encoding="utf-8")
    top_file.write_text("dummy topology\n", encoding="utf-8")

    job = GromacsEMJob(
        molecule=None,
        label="em",
        jobrunner=None,
        mdp_file=mdp_file,
        structure_file=structure_file,
        top_file=top_file,
    )

    runner = GromacsJobRunner.__new__(GromacsJobRunner)

    runner._validate_gromacs_inputs(job)


def test_gromacs_runner_validate_inputs_raises_for_missing_files(tmp_path):
    mdp_file = tmp_path / "em.mdp"
    structure_file = tmp_path / "input.gro"

    mdp_file.write_text("integrator = steep\n", encoding="utf-8")
    structure_file.write_text("dummy structure\n", encoding="utf-8")

    job = GromacsEMJob(
        molecule=None,
        label="em",
        jobrunner=None,
        mdp_file=mdp_file,
        structure_file=structure_file,
        top_file=None,
    )

    runner = GromacsJobRunner.__new__(GromacsJobRunner)

    with pytest.raises(FileNotFoundError, match="top_file"):
        runner._validate_gromacs_inputs(job)
