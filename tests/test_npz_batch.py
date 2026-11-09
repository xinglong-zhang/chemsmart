import json
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest
from click.testing import CliRunner

from chemsmart.cli.main import maybe_submit_deferred_npz
from chemsmart.cli.sub import sub
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.xyz.xyzfile import XYZFile
from chemsmart.jobs.batch import (
    _replace_task_structure,
    run_streaming_npz_batch,
)
from chemsmart.settings.server import Server
from tests.helpers import fake_executable_configuration


def test_multiframe_xyz_streams_selected_frames(tmp_path):
    source = tmp_path / "frames.xyz"
    source.write_text(
        "1\nfirst\nH 0 0 0\n1\nsecond\nH 0 0 1\n1\nthird\nH 0 0 2\n"
    )
    results = list(XYZFile(str(source)).iter_molecule_results("2:"))
    assert [index for index, _, error in results if error is None] == [2, 3]
    assert [molecule.positions[0, 2] for _, molecule, _ in results] == [1, 2]
    assert [
        index
        for index, _, _ in XYZFile(str(source)).iter_molecule_results("-1")
    ] == [3]


def test_multiframe_xyz_uses_deferred_submission(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    server = Server(
        "test", SCHEDULER="SLURM", NUM_CORES=8, MEM_GB=24, NUM_GPUS=0
    )
    with (
        patch("chemsmart.settings.server.Server.current", return_value=server),
        patch(
            "chemsmart.settings.executable.XTBExecutable.from_servername",
            return_value=fake_executable_configuration(),
        ),
    ):
        assert maybe_submit_deferred_npz(
            ["sub", "-N", "1", "--test", "xtb", "-f", "frames.xyz", "sp"]
        )
    manifest = json.loads(Path("frames_sp_deferred_npz.json").read_text())
    assert manifest["format"] == ".xyz"


def test_npz_node_batch_is_submitted_before_archive_is_read(
    tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    server = Server(
        "test", SCHEDULER="SLURM", NUM_CORES=8, MEM_GB=24, NUM_GPUS=0
    )
    arguments = [
        "sub",
        "-N",
        "2",
        "-n",
        "1",
        "xtb",
        "-f",
        "/data/large.npz",
        "-i",
        ":",
        "sp",
    ]

    with (
        patch("chemsmart.settings.server.Server.current", return_value=server),
        patch(
            "chemsmart.settings.executable.XTBExecutable.from_servername",
            return_value=fake_executable_configuration(),
        ),
        patch(
            "chemsmart.utils.cluster.ClusterHelper.get_running_jobs",
            return_value=([], []),
        ),
        patch("chemsmart.settings.server.subprocess.run") as submit,
    ):
        submit.return_value.returncode = 0
        assert maybe_submit_deferred_npz(arguments) is True

    submit.assert_called_once_with(
        ["sbatch", "large_sp_deferred_npz.sh"], check=False
    )
    script = Path("large_sp_deferred_npz.sh").read_text()
    assert "#SBATCH --nodes=2" in script
    assert "run_streaming_npz_batch" in script
    assert "python -c" in script
    assert sys.executable not in script
    manifest = Path("large_sp_deferred_npz.json").read_text()
    assert '"queue_capacity": 100' in manifest
    assert '"filename": "/data/large.npz"' in manifest
    assert '"base_label": "large_sp"' in manifest


def test_deferred_npz_test_mode_writes_without_submitting(
    tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    server = Server(
        "test", SCHEDULER="SLURM", NUM_CORES=8, MEM_GB=24, NUM_GPUS=0
    )
    with (
        patch("chemsmart.settings.server.Server.current", return_value=server),
        patch(
            "chemsmart.settings.executable.XTBExecutable.from_servername",
            return_value=fake_executable_configuration(),
        ),
        patch("chemsmart.settings.server.subprocess.run") as submit,
    ):
        assert maybe_submit_deferred_npz(
            ["sub", "-N", "1", "--test", "xtb", "-f", "large.npz", "sp"]
        )

    submit.assert_not_called()
    assert Path("large_sp_deferred_npz.sh").is_file()


def test_equals_options_use_deferred_submission(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    server = Server(
        "test", SCHEDULER="SLURM", NUM_CORES=8, MEM_GB=24, NUM_GPUS=0
    )
    with (
        patch("chemsmart.settings.server.Server.current", return_value=server),
        patch(
            "chemsmart.settings.executable.XTBExecutable.from_servername",
            return_value=fake_executable_configuration(),
        ),
    ):
        assert maybe_submit_deferred_npz(
            [
                "sub",
                "--num-nodes=2",
                "--test",
                "xtb",
                "--filename=large.npz",
                "--index=1,2",
                "sp",
            ]
        )

    manifest = json.loads(Path("large_sp_deferred_npz.json").read_text())
    assert manifest["num_nodes"] == 2
    assert manifest["index"] == "1,2"
    assert "--num-nodes=2" not in manifest["task_arguments"]


def test_equals_filename_is_replaced_for_worker():
    arguments = _replace_task_structure(
        ["xtb", "--filename=data.npz", "--index=4", "sp"],
        "prepared.xyz",
        "data_idx4",
    )
    assert "--filename=prepared.xyz" in arguments
    assert not any(argument.startswith("--index") for argument in arguments)


@pytest.mark.parametrize(
    ("label_arguments", "expected_label"),
    [
        (["--append-label", "trial"], "large_trial"),
        (["--label", "custom"], "custom"),
    ],
)
def test_streaming_batch_preserves_program_label_rules(
    tmp_path, monkeypatch, label_arguments, expected_label
):
    monkeypatch.chdir(tmp_path)
    server = Server(
        "test", SCHEDULER="SLURM", NUM_CORES=8, MEM_GB=24, NUM_GPUS=0
    )
    with (
        patch("chemsmart.settings.server.Server.current", return_value=server),
        patch(
            "chemsmart.settings.executable.XTBExecutable.from_servername",
            return_value=fake_executable_configuration(),
        ),
    ):
        assert maybe_submit_deferred_npz(
            [
                "sub",
                "-N",
                "1",
                "--test",
                "xtb",
                "-f",
                "large.npz",
                *label_arguments,
                "sp",
            ]
        )

    manifest = json.loads(
        Path(f"{expected_label}_deferred_npz.json").read_text()
    )
    assert manifest["base_label"] == expected_label


@pytest.mark.parametrize(
    ("program", "executable_class"),
    [
        ("gaussian", "GaussianExecutable"),
        ("orca", "ORCAExecutable"),
    ],
)
def test_program_title_is_not_consumed_as_slurm_time(
    tmp_path, monkeypatch, program, executable_class
):
    monkeypatch.chdir(tmp_path)
    server = Server(
        "test",
        SCHEDULER="SLURM",
        NUM_CORES=8,
        MEM_GB=24,
        NUM_GPUS=0,
        NUM_HOURS=12,
    )
    with (
        patch("chemsmart.settings.server.Server.current", return_value=server),
        patch(
            f"chemsmart.settings.executable.{executable_class}.from_servername",
            return_value=fake_executable_configuration(),
        ),
    ):
        assert maybe_submit_deferred_npz(
            [
                "sub",
                "-N",
                "1",
                "--test",
                program,
                "-f",
                "large.npz",
                "-t",
                "molecular title",
                "sp",
            ]
        )

    manifest = json.loads(Path("large_sp_deferred_npz.json").read_text())
    assert ["-t", "molecular title"] == manifest["task_arguments"][
        manifest["task_arguments"]
        .index("-t") : manifest["task_arguments"]
        .index("-t")
        + 2
    ]
    script = Path("large_sp_deferred_npz.sh").read_text()
    assert "#SBATCH --time=12:00:00" in script
    assert "#SBATCH --time=molecular title" not in script


def test_streaming_script_preserves_gpu_account_and_email(
    tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    server = Server(
        "test", SCHEDULER="SLURM", NUM_CORES=8, MEM_GB=24, NUM_GPUS=2
    )
    with (
        patch("chemsmart.settings.server.Server.current", return_value=server),
        patch(
            "chemsmart.settings.executable.XTBExecutable.from_servername",
            return_value=fake_executable_configuration(),
        ),
        patch.dict(
            "chemsmart.settings.submitters.user_settings.data",
            {"PROJECT": "chem", "EMAIL": "chemist@example.com"},
            clear=True,
        ),
    ):
        assert maybe_submit_deferred_npz(
            ["sub", "-N", "1", "--test", "xtb", "-f", "large.npz", "sp"]
        )

    script = Path("large_sp_deferred_npz.sh").read_text()
    assert "#SBATCH --gres=gpu:2" in script
    assert "#SBATCH --account=chem" in script
    assert "#SBATCH --mail-user=chemist@example.com" in script


def test_duplicate_batch_does_not_overwrite_manifest(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    server = Server(
        "test", SCHEDULER="SLURM", NUM_CORES=8, MEM_GB=24, NUM_GPUS=0
    )
    with (
        patch("chemsmart.settings.server.Server.current", return_value=server),
        patch(
            "chemsmart.utils.cluster.ClusterHelper.get_running_jobs",
            return_value=([], ["large_sp"]),
        ),
        pytest.raises(
            SystemExit, match="Duplicate job NOT submitted: large_sp"
        ),
    ):
        maybe_submit_deferred_npz(
            ["sub", "-N", "1", "xtb", "-f", "large.npz", "sp"]
        )

    assert not Path("large_sp_deferred_npz.json").exists()
    assert not Path("large_sp_deferred_npz.sh").exists()


def test_streaming_batch_removes_each_prepared_structure(
    tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    source = Path(__file__).parent / "data/IOTests/revQM9_subset.npz"
    manifest = {
        "filename": str(source),
        "index": "1,2",
        "base_label": "sample",
        "queue_capacity": 1,
        "num_nodes": 1,
        "cores_per_node": 2,
        "memory_per_node_gb": 4,
        "gpus_per_node": 0,
        "task_cores": 1,
        "task_memory_gb": 1,
        "task_gpus": 0,
        "task_arguments": ["-n", "1", "xtb", "-f", str(source), "sp"],
    }
    manifest_path = tmp_path / "sample_deferred_npz.json"
    manifest_path.write_text(json.dumps(manifest))

    with patch("chemsmart.jobs.batch.subprocess.run") as run:
        run.return_value.returncode = 0
        assert run_streaming_npz_batch(manifest_path) == 0

    assert run.call_count == 2
    assert not (tmp_path / "sample_deferred_npz_inputs").exists()
    report = json.loads(
        (tmp_path / "sample_deferred_npz.report.json").read_text()
    )
    assert report == {
        "successful": 2,
        "failed": 0,
        "unfinished": 0,
        "total": 2,
        "status_file": str(tmp_path / "sample_deferred_npz.status.jsonl"),
    }
    commands = [call.args[0] for call in run.call_args_list]
    assert all(command[0] == "srun" for command in commands)
    assert all("--batch-label" in command for command in commands)
    assert all(
        call.kwargs["env"]["OMP_NUM_THREADS"] == "1"
        for call in run.call_args_list
    )


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
