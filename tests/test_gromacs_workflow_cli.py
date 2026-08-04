"""
Tests for the staged GROMACS workflow CLI command.
"""

import importlib

from click.testing import CliRunner

from chemsmart.cli.gromacs.gromacs import gromacs

workflow_module = importlib.import_module(
    "chemsmart.cli.gromacs.workflow"
)


class DummyGromacsWorkflow:
    """
    Capture workflow constructor arguments without running GROMACS.
    """

    captured = {}

    def __init__(self, **kwargs):
        self.__class__.captured = dict(kwargs)


def test_gromacs_workflow_command_is_registered():
    runner = CliRunner()

    result = runner.invoke(
        gromacs,
        ["--help"],
        obj={"molecule": None},
    )

    assert result.exit_code == 0, result.output
    assert "workflow" in result.output


def test_gromacs_workflow_cli_creates_workflow(
    monkeypatch,
    tmp_path,
):
    DummyGromacsWorkflow.captured = {}

    monkeypatch.setattr(
        workflow_module,
        "GromacsWorkflow",
        DummyGromacsWorkflow,
    )

    structure_file = tmp_path / "input.gro"
    top_file = tmp_path / "topol.top"
    index_file = tmp_path / "index.ndx"
    ligand_itp = tmp_path / "ligand.itp"
    npt_mdp = tmp_path / "npt.mdp"
    output_folder = tmp_path / "workflow"

    for path in (
        structure_file,
        top_file,
        index_file,
        ligand_itp,
        npt_mdp,
    ):
        path.write_text(
            "dummy content\n",
            encoding="utf-8",
        )

    runner = CliRunner()

    result = runner.invoke(
        gromacs,
        [
            "workflow",
            "--structure",
            str(structure_file),
            "--top",
            str(top_file),
            "--index",
            str(index_file),
            "--itp",
            str(ligand_itp),
            "--npt-mdp",
            str(npt_mdp),
            "--output-folder",
            str(output_folder),
            "--skip-completed",
        ],
        obj={
            "molecule": None,
            "jobrunner": None,
        },
    )

    assert result.exit_code == 0, result.output
    assert result.exception is None

    captured = DummyGromacsWorkflow.captured

    assert captured
    assert captured["structure_file"] == structure_file
    assert captured["top_file"] == top_file
    assert captured["index_file"] == index_file
    assert captured["itp_files"] == [ligand_itp]
    assert captured["folder"] == output_folder
    assert captured["skip_completed"] is True

    assert captured["mdp_files"] == {
        "em": None,
        "nvt": None,
        "npt": npt_mdp,
        "md": None,
    }


def test_gromacs_workflow_cli_requires_structure_and_top():
    runner = CliRunner()

    result = runner.invoke(
        gromacs,
        ["workflow"],
        obj={"molecule": None},
    )

    assert result.exit_code != 0
    assert "Missing required workflow options" in result.output