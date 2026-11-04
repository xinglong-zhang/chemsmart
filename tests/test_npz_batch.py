from pathlib import Path
from unittest.mock import patch

import numpy as np
from click.testing import CliRunner

from chemsmart.cli.sub import sub
from chemsmart.io.molecules.structure import Molecule
from chemsmart.settings.server import Server


def test_npz_batch_prepares_selected_structures(tmp_path, monkeypatch):
    source = Path(__file__).parent / "data/IOTests/revQM9_subset.npz"
    original = Molecule.from_filepath(
        str(source), index="1,2", return_list=True
    )
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
                "-n",
                "1",
                "--test",
                "xtb",
                "-p",
                "test",
                "-f",
                str(source.resolve()),
                "-i",
                "1,2",
                "-c",
                "0",
                "-m",
                "1",
                "sp",
            ],
        )
    assert result.exit_code == 0, result.output
    calls = submit.call_args.kwargs
    assert calls["num_nodes"] == 2
    assert len(calls["jobs"]) == 2
    assert all(not Path(job.folder).exists() for job in calls["jobs"])
    for args, molecule in zip(calls["cli_args"], original):
        assert args[args.index("--num-cores") + 1] == "1"
        assert "--index" not in args
        prepared = Molecule.from_filepath(args[args.index("--filename") + 1])
        np.testing.assert_array_equal(prepared.symbols, molecule.symbols)
        np.testing.assert_allclose(
            prepared.positions, molecule.positions, atol=1e-5
        )


def test_batch_replay_uses_original_folder_without_temporary_directories(
    tmp_path, monkeypatch
):
    from chemsmart.cli.run import run

    source = Path(__file__).parent / "data/IOTests/revQM9_subset.npz"
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
                "-n",
                "1",
                "--test",
                "xtb",
                "-p",
                "test",
                "-f",
                str(source.resolve()),
                "-i",
                "1,2",
                "-c",
                "0",
                "-m",
                "1",
                "sp",
            ],
        )
    assert result.exit_code == 0, result.output
    jobs = submit.call_args.kwargs["jobs"]
    arguments = submit.call_args.kwargs["cli_args"]
    assert not any(Path(job.folder).exists() for job in jobs)

    executed = []

    def execute(job, runner):
        executed.append((job.label, Path(job.folder)))
        Path(job.folder).mkdir(parents=True, exist_ok=True)
        Path(job.outputfile).write_text("test execution output")

    with (
        patch.object(Server, "current", return_value=server),
        patch("chemsmart.cli.run._run_single_job", side_effect=execute),
        patch("chemsmart.jobs.xtb.job.XTBJob.is_complete", return_value=True),
    ):
        for args in arguments:
            result = CliRunner().invoke(run, args)
            assert result.exit_code == 0, result.output
    for original, (label, folder) in zip(jobs, executed):
        assert label == original.label
        assert folder == tmp_path / label
        assert (folder / f"{label}.out").is_file()
    assert not list(tmp_path.glob("structure_*_sp"))
