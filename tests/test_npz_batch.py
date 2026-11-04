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
    for args, molecule in zip(calls["cli_args"], original):
        assert args[args.index("--num-cores") + 1] == "1"
        assert "--index" not in args
        prepared = Molecule.from_filepath(args[args.index("--filename") + 1])
        np.testing.assert_array_equal(prepared.symbols, molecule.symbols)
        np.testing.assert_allclose(
            prepared.positions, molecule.positions, atol=1e-5
        )
