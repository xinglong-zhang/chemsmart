import json
from types import SimpleNamespace

from chemsmart.jobs.batch import (
    _batch_capacity,
    run_batch_chunk,
    run_slurm_batch,
)


def test_batch_capacity_respects_core_and_memory_limits():
    manifest = {
        "num_nodes": 5,
        "cores_per_node": 64,
        "memory_per_node_gb": 375,
        "gpus_per_node": 0,
        "task_cores": 1,
        "task_memory_gb": 10,
        "task_gpus": 0,
    }

    assert _batch_capacity(manifest) == 185


def test_batch_runner_continues_after_task_failure(tmp_path, mocker, capsys):
    manifest_path = tmp_path / "batch.json"
    manifest_path.write_text(
        json.dumps(
            {
                "num_nodes": 2,
                "cores_per_node": 4,
                "memory_per_node_gb": 8,
                "gpus_per_node": 0,
                "task_cores": 1,
                "task_memory_gb": 2,
                "task_gpus": 0,
                "tasks": [
                    {"label": f"structure_{i}", "runscript": f"run_{i}.py"}
                    for i in range(3)
                ],
            }
        )
    )
    commands = []

    def run_task(command, check, stderr, stdout=None):
        commands.append(command)
        if command[0] == "srun":
            run_batch_chunk(*command[-3:])
            return SimpleNamespace(returncode=0)
        return SimpleNamespace(
            returncode=1 if command[-1] == "run_1.py" else 0
        )

    mocker.patch("chemsmart.jobs.batch.subprocess.run", side_effect=run_task)

    assert run_slurm_batch(manifest_path) == 1
    assert len(commands) == 4
    assert sum(command[0] == "srun" for command in commands) == 1
    assert "--cpus-per-task=4" in commands[0]

    report = json.loads((tmp_path / "batch.report.json").read_text())
    assert report["successful"] == 2
    assert report["failed"] == 1
    assert report["unfinished"] == 0
    assert [task["status"] for task in report["tasks"]] == [
        "successful",
        "failed",
        "successful",
    ]
    assert "successful=2, failed=1, unfinished=0" in capsys.readouterr().out


def test_batch_resume_skips_successful_tasks(tmp_path, mocker):
    manifest_path = tmp_path / "batch.json"
    manifest_path.write_text(
        json.dumps(
            {
                "num_nodes": 1,
                "cores_per_node": 2,
                "memory_per_node_gb": 4,
                "gpus_per_node": 0,
                "task_cores": 1,
                "task_memory_gb": 1,
                "task_gpus": 0,
                "tasks": [
                    {"label": f"structure_{i}", "runscript": f"run_{i}.py"}
                    for i in range(3)
                ],
            }
        )
    )
    (tmp_path / "batch.status.jsonl").write_text(
        json.dumps(
            {"index": 1, "label": "structure_0", "status": "successful"}
        )
        + "\n"
    )
    commands = []

    def run_task(command, check, stderr, stdout=None):
        commands.append(command)
        if command[0] == "srun":
            run_batch_chunk(*command[-3:])
        return SimpleNamespace(returncode=0)

    mocker.patch("chemsmart.jobs.batch.subprocess.run", side_effect=run_task)
    assert run_slurm_batch(manifest_path) == 0
    assert [command[-1] for command in commands if command[0] != "srun"] == [
        "run_1.py",
        "run_2.py",
    ]


def test_step_launch_failure_leaves_tasks_unfinished(tmp_path, mocker):
    manifest_path = tmp_path / "batch.json"
    manifest_path.write_text(
        json.dumps(
            {
                "num_nodes": 1,
                "cores_per_node": 2,
                "memory_per_node_gb": 4,
                "gpus_per_node": 0,
                "task_cores": 1,
                "task_memory_gb": 1,
                "task_gpus": 0,
                "tasks": [
                    {"label": "structure_1", "runscript": "run_1.py"},
                    {"label": "structure_2", "runscript": "run_2.py"},
                ],
            }
        )
    )
    mocker.patch(
        "chemsmart.jobs.batch.subprocess.run",
        return_value=SimpleNamespace(returncode=1),
    )
    assert run_slurm_batch(manifest_path) == 1
    report = json.loads((tmp_path / "batch.report.json").read_text())
    assert report["successful"] == 0
    assert report["failed"] == 0
    assert report["unfinished"] == 2
