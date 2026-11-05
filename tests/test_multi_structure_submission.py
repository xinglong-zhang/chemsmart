from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest
from click.testing import CliRunner

from chemsmart.cli.run import run
from chemsmart.cli.sub import sub
from chemsmart.io.molecules.structure import Molecule
from chemsmart.settings.server import Server


@pytest.mark.parametrize("kind", ["npz", "xyz"])
@pytest.mark.parametrize("program", ["xtb", "gaussian", "orca"])
def test_multi_structure_submission_replays_one_structure(
    tmp_path, monkeypatch, kind, program
):
    source = tmp_path / ("molecules." + kind)
    if kind == "npz":
        np.savez(
            source,
            coords=np.array([[[0.0, 0.0, 0.0]], [[1.0, 0.0, 0.0]]]),
            charges=np.array([[2], [2]]),
        )
    else:
        source.write_text("1\nfirst\nHe 0 0 0\n1\nsecond\nHe 1 0 0\n")
    server = Server(
        "test", SCHEDULER="SLURM", NUM_CORES=8, MEM_GB=24, NUM_GPUS=0
    )
    monkeypatch.chdir(tmp_path)
    with (
        patch("chemsmart.cli.sub.Server.current", return_value=server),
        patch.object(server, "submit") as submit,
    ):
        result = CliRunner().invoke(
            sub,
            [
                "--test",
                program,
                "-p",
                "test",
                "-f",
                str(source),
                "-i",
                "1,2",
                "-c",
                "0",
                "-m",
                "1",
                "sp",
            ],
        )
    assert result.exit_code == 0, (result.output, result.exception)
    assert submit.call_count == 2
    for index, call in enumerate(submit.call_args_list):
        args = call.kwargs["cli_args"]
        prepared = Path(args[args.index("--filename") + 1])
        assert prepared != source and prepared.is_file()
        assert "--index" not in args
        molecule = Molecule.from_filepath(str(prepared))
        np.testing.assert_allclose(molecule.positions[0], [float(index), 0, 0])
        captured = []

        def capture(job, _):
            captured.append(job)
            job.is_complete = lambda: True

        with (
            patch("chemsmart.cli.run.JobRunner") as runner,
            patch(
                "chemsmart.cli.run._run_single_job",
                side_effect=capture,
            ),
        ):
            runner.return_value = (
                submit.call_args_list[index].kwargs["job"].jobrunner
            )
            replay = CliRunner().invoke(run, args)
        assert replay.exit_code == 0, (replay.output, replay.exception)
        assert len(captured) == 1
        assert captured[0].label == call.kwargs["job"].label


def test_unsupported_multistructure_input_is_not_submitted(
    tmp_path, monkeypatch
):
    source = tmp_path / "structures.trj"
    source.write_text("1\nfirst\nHe 0 0 0\n1\nsecond\nHe 1 0 0\n")
    server = Server(
        "test", SCHEDULER="SLURM", NUM_CORES=8, MEM_GB=24, NUM_GPUS=0
    )
    monkeypatch.chdir(tmp_path)
    with (
        patch("chemsmart.cli.sub.Server.current", return_value=server),
        patch.object(server, "submit") as submit,
    ):
        result = CliRunner().invoke(
            sub,
            [
                "--test",
                "xtb",
                "-p",
                "test",
                "-f",
                str(source),
                "-i",
                ":",
                "-c",
                "0",
                "-m",
                "1",
                "sp",
            ],
        )
    assert result.exit_code != 0
    assert "convert the input to XYZ" in result.output
    submit.assert_not_called()


@pytest.mark.parametrize("count", [1, 2])
def test_batch_label_requires_one_replayed_job(count):
    import click

    from chemsmart.jobs.job import Job
    from chemsmart.utils.cli import MyCommand

    jobs = [
        Job(molecule=None, label=f"job_{i}", jobrunner=None)
        for i in range(count)
    ]
    for job in jobs:
        job.is_complete = lambda: True

    @click.command("replay-test", cls=MyCommand)
    def replay_test():
        return jobs

    run.add_command(replay_test)
    try:
        with (
            patch("chemsmart.cli.run.JobRunner"),
            patch("chemsmart.cli.run._run_single_job") as execute,
        ):
            result = CliRunner().invoke(
                run, ["--batch-label", "selected", "replay-test"]
            )
        if count == 1:
            assert result.exit_code == 0, (result.output, result.exception)
            assert execute.call_args.args[0].label == "selected"
        else:
            assert isinstance(result.exception, ValueError)
            assert "exactly one job" in str(result.exception)
            execute.assert_not_called()
    finally:
        run.commands.pop("replay-test", None)


def test_batch_label_rejects_multiple_molecules():
    from types import SimpleNamespace

    from chemsmart.cli.utils import build_jobs

    context = SimpleNamespace(
        obj={
            "jobrunner": None,
            "molecules": [object(), object()],
            "molecule_indices": [1, 2],
            "label": "original",
            "batch_label": "selected",
        }
    )
    with pytest.raises(ValueError, match="exactly one molecule"):
        build_jobs(context, object, None, True, {})


def test_prepared_structure_arguments_require_input_filename():
    from chemsmart.cli.utils import prepared_structure_arguments

    with pytest.raises(ValueError, match="input filename"):
        prepared_structure_arguments(
            ["xtb", "sp"], "/tmp/structure.xyz", "molecule"
        )


@pytest.mark.parametrize("program", ["xtb", "gaussian", "orca"])
def test_node_batch_prepares_multiframe_xyz(tmp_path, monkeypatch, program):
    source = tmp_path / "structures.xyz"
    source.write_text("1\nfirst\nHe 0 0 0\n1\nsecond\nHe 1 0 0\n")
    server = Server(
        "test", SCHEDULER="SLURM", NUM_CORES=8, MEM_GB=24, NUM_GPUS=0
    )
    monkeypatch.chdir(tmp_path)
    with (
        patch("chemsmart.cli.sub.Server.current", return_value=server),
        patch.object(server, "submit_batch_job") as submit,
    ):
        result = CliRunner().invoke(
            sub,
            [
                "-N",
                "2",
                "--test",
                program,
                "-p",
                "test",
                "-f",
                str(source),
                "-i",
                ":",
                "-c",
                "0",
                "-m",
                "1",
                "sp",
            ],
        )
    assert result.exit_code == 0, (result.output, result.exception)
    args_list = submit.call_args.kwargs["cli_args"]
    assert len(args_list) == 2
    for index, args in enumerate(args_list):
        path = args[args.index("--filename") + 1]
        molecule = Molecule.from_filepath(path)
        np.testing.assert_allclose(molecule.positions[0], [index, 0, 0])
        assert "--index" not in args
