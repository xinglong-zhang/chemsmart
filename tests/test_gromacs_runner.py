import os
from pathlib import Path

import pytest

from chemsmart.jobs.gromacs.job import GromacsEMJob
from chemsmart.jobs.gromacs.runner import GromacsJobRunner
from chemsmart.settings.executable import GromacsExecutable


def _make_runner(gmx_executable="gmx"):
    """
    Create a minimal GromacsJobRunner instance for unit tests.

    __new__ is used to avoid requiring a real server object.
    """
    runner = GromacsJobRunner.__new__(GromacsJobRunner)
    runner.gmx_executable = gmx_executable
    return runner


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
    assert job.workflow == "prepared"
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
        tpr_file=tpr_file,
    )

    job.set_folder(str(tmp_path))

    runner = _make_runner()
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

    runner = _make_runner()

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

    runner = _make_runner()

    with pytest.raises(FileNotFoundError, match="top_file"):
        runner._validate_gromacs_inputs(job)


def test_gromacs_runner_get_grompp_command_without_index(tmp_path):
    mdp_file = tmp_path / "em.mdp"
    structure_file = tmp_path / "input.gro"
    top_file = tmp_path / "topol.top"
    tpr_file = tmp_path / "em.tpr"

    job = GromacsEMJob(
        molecule=None,
        label="em",
        jobrunner=None,
        mdp_file=mdp_file,
        structure_file=structure_file,
        top_file=top_file,
        tpr_file=tpr_file,
    )

    runner = _make_runner()

    command = runner._get_grompp_command(job)

    assert command == [
        "gmx",
        "grompp",
        "-f",
        str(mdp_file),
        "-c",
        str(structure_file),
        "-p",
        str(top_file),
        "-o",
        str(tpr_file),
    ]


def test_gromacs_runner_get_grompp_command_with_index(tmp_path):
    mdp_file = tmp_path / "em.mdp"
    structure_file = tmp_path / "input.gro"
    top_file = tmp_path / "topol.top"
    tpr_file = tmp_path / "em.tpr"
    index_file = tmp_path / "index.ndx"

    job = GromacsEMJob(
        molecule=None,
        label="em",
        jobrunner=None,
        mdp_file=mdp_file,
        structure_file=structure_file,
        top_file=top_file,
        tpr_file=tpr_file,
        index_file=index_file,
    )

    runner = _make_runner()

    command = runner._get_grompp_command(job)

    assert command == [
        "gmx",
        "grompp",
        "-f",
        str(mdp_file),
        "-c",
        str(structure_file),
        "-p",
        str(top_file),
        "-o",
        str(tpr_file),
        "-n",
        str(index_file),
    ]


def test_gromacs_runner_get_mdrun_command(tmp_path):
    tpr_file = tmp_path / "em.tpr"

    job = GromacsEMJob(
        molecule=None,
        label="em",
        jobrunner=None,
        mdp_file=tmp_path / "em.mdp",
        structure_file=tmp_path / "input.gro",
        top_file=tmp_path / "topol.top",
        tpr_file=tpr_file,
    )

    runner = _make_runner()

    command = runner._get_mdrun_command(job)

    assert command == [
        "gmx",
        "mdrun",
        "-deffnm",
        str(tmp_path / "em"),
    ]


def test_gromacs_runner_get_commands_for_prepared_workflow(tmp_path):
    mdp_file = tmp_path / "em.mdp"
    structure_file = tmp_path / "input.gro"
    top_file = tmp_path / "topol.top"
    tpr_file = tmp_path / "em.tpr"

    job = GromacsEMJob(
        molecule=None,
        label="em",
        jobrunner=None,
        mdp_file=mdp_file,
        structure_file=structure_file,
        top_file=top_file,
        tpr_file=tpr_file,
        workflow="prepared",
    )

    runner = _make_runner()

    commands = runner._get_commands(job)

    assert commands == [
        [
            "gmx",
            "grompp",
            "-f",
            str(mdp_file),
            "-c",
            str(structure_file),
            "-p",
            str(top_file),
            "-o",
            str(tpr_file),
        ],
        [
            "gmx",
            "mdrun",
            "-deffnm",
            str(tmp_path / "em"),
        ],
    ]


def test_gromacs_runner_accepts_custom_gmx_executable(tmp_path):
    tpr_file = tmp_path / "em.tpr"

    job = GromacsEMJob(
        molecule=None,
        label="em",
        jobrunner=None,
        mdp_file=tmp_path / "em.mdp",
        structure_file=tmp_path / "input.gro",
        top_file=tmp_path / "topol.top",
        tpr_file=tpr_file,
    )

    runner = _make_runner(gmx_executable="/opt/gromacs/bin/gmx")

    command = runner._get_mdrun_command(job)

    assert command == [
        "/opt/gromacs/bin/gmx",
        "mdrun",
        "-deffnm",
        str(tmp_path / "em"),
    ]


def test_gromacs_runner_get_pdb2gmx_command(tmp_path):
    input_pdb = tmp_path / "input.pdb"
    output_gro = tmp_path / "processed.gro"
    top_file = tmp_path / "topol.top"

    job = GromacsEMJob(
        molecule=None,
        label="em",
        jobrunner=None,
        mdp_file=tmp_path / "em.mdp",
        structure_file=output_gro,
        top_file=top_file,
    )
    job.input_pdb = input_pdb
    job.output_gro = output_gro
    job.force_field = "amber99sb-ildn"
    job.water_model = "tip3p"

    runner = _make_runner()

    command = runner._get_pdb2gmx_command(job)

    assert command == [
        "gmx",
        "pdb2gmx",
        "-f",
        str(input_pdb),
        "-o",
        str(output_gro),
        "-p",
        str(top_file),
        "-ff",
        "amber99sb-ildn",
        "-water",
        "tip3p",
    ]


def test_gromacs_runner_get_editconf_command(tmp_path):
    input_structure = tmp_path / "processed.gro"
    output_structure = tmp_path / "boxed.gro"

    job = GromacsEMJob(
        molecule=None,
        label="em",
        jobrunner=None,
        mdp_file=tmp_path / "em.mdp",
        structure_file=input_structure,
        top_file=tmp_path / "topol.top",
    )
    job.editconf_input_file = input_structure
    job.editconf_output_file = output_structure
    job.box_type = "cubic"
    job.box_distance = 1.0

    runner = _make_runner()

    command = runner._get_editconf_command(job)

    assert command == [
        "gmx",
        "editconf",
        "-f",
        str(input_structure),
        "-o",
        str(output_structure),
        "-bt",
        "cubic",
        "-d",
        "1.0",
    ]


def test_gromacs_runner_get_solvate_command(tmp_path):
    input_structure = tmp_path / "boxed.gro"
    output_structure = tmp_path / "solvated.gro"
    top_file = tmp_path / "topol.top"

    job = GromacsEMJob(
        molecule=None,
        label="em",
        jobrunner=None,
        mdp_file=tmp_path / "em.mdp",
        structure_file=tmp_path / "input.gro",
        top_file=top_file,
    )
    job.boxed_structure_file = input_structure
    job.solvated_structure_file = output_structure

    runner = _make_runner()

    command = runner._get_solvate_command(job)

    assert command == [
        "gmx",
        "solvate",
        "-cp",
        str(input_structure),
        "-o",
        str(output_structure),
        "-p",
        str(top_file),
    ]


def test_gromacs_runner_get_genion_command(tmp_path):
    input_tpr = tmp_path / "ions.tpr"
    output_structure = tmp_path / "ionized.gro"
    top_file = tmp_path / "topol.top"

    job = GromacsEMJob(
        molecule=None,
        label="em",
        jobrunner=None,
        mdp_file=tmp_path / "em.mdp",
        structure_file=tmp_path / "input.gro",
        top_file=top_file,
    )
    job.ions_tpr_file = input_tpr
    job.ions_output_file = output_structure
    job.positive_ion = "NA"
    job.negative_ion = "CL"
    job.neutral = True

    runner = _make_runner()

    command = runner._get_genion_command(job)

    assert command == [
        "gmx",
        "genion",
        "-s",
        str(input_tpr),
        "-o",
        str(output_structure),
        "-p",
        str(top_file),
        "-pname",
        "NA",
        "-nname",
        "CL",
        "-neutral",
    ]


def test_gromacs_runner_full_setup_is_not_implemented():
    job = GromacsEMJob(
        molecule=None,
        label="em",
        jobrunner=None,
        workflow="full_setup",
    )

    runner = _make_runner()

    with pytest.raises(NotImplementedError):
        runner._get_commands(job)


def test_gromacs_runner_raises_for_unsupported_workflow():
    job = GromacsEMJob(
        molecule=None,
        label="em",
        jobrunner=None,
        workflow="unsupported",
    )

    runner = _make_runner()

    with pytest.raises(ValueError, match="Unsupported GROMACS workflow"):
        runner._get_commands(job)


def test_gromacs_executable_defaults_to_gmx():
    executable = GromacsExecutable()

    assert executable.get_executable() == "gmx"


def test_gromacs_executable_uses_executable_folder(tmp_path):
    executable = GromacsExecutable(executable_folder=tmp_path)

    assert executable.get_executable() == os.path.join(str(tmp_path), "gmx")