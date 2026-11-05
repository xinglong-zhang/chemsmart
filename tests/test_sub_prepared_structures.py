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
        with (
            patch("chemsmart.cli.run.JobRunner") as runner,
            patch(
                "chemsmart.cli.run._run_single_job",
                side_effect=lambda job, _: captured.append(job),
            ),
        ):
            runner.return_value = (
                submit.call_args_list[index].kwargs["job"].jobrunner
            )
            replay = CliRunner().invoke(run, args)
        assert replay.exit_code == 0, (replay.output, replay.exception)
        assert len(captured) == 1
        assert captured[0].label == call.kwargs["job"].label
