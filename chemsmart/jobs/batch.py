import json
import logging
import os
import re
import subprocess
import sys
import threading
from concurrent.futures import (
    FIRST_COMPLETED,
    ThreadPoolExecutor,
    as_completed,
    wait,
)
from pathlib import Path

logger = logging.getLogger(__name__)


def _replace_task_structure(arguments, filename, label):
    """Point one replayed CLI invocation at a prepared XYZ structure."""
    arguments = list(arguments)
    for option in ("-f", "--filename"):
        if option in arguments:
            arguments[arguments.index(option) + 1] = str(filename)
            break
        equals_prefix = f"{option}="
        for position, argument in enumerate(arguments):
            if option.startswith("--") and argument.startswith(equals_prefix):
                arguments[position] = f"{option}={filename}"
                break
    for option in ("-i", "--index", "--si", "--structure-index"):
        if option in arguments:
            position = arguments.index(option)
            del arguments[position : position + 2]
            break
        equals_prefix = f"{option}="
        for position, argument in enumerate(arguments):
            if option.startswith("--") and argument.startswith(equals_prefix):
                del arguments[position]
                break
    return ["--batch-label", label, *arguments]


def run_streaming_npz_batch(manifest_path):
    """Prepare NPZ or multi-frame XYZ lazily in one SLURM allocation."""
    from chemsmart.io.npz import NPZFile
    from chemsmart.io.xyz.xyzfile import XYZFile

    manifest_path = Path(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    workers = _batch_capacity(manifest)
    queue_capacity = int(manifest.get("queue_capacity", 100))
    if queue_capacity <= 0:
        raise ValueError("queue capacity must be a positive integer.")
    max_in_flight = workers + queue_capacity
    prefix = manifest_path.with_suffix("")
    input_dir = prefix.parent / f"{prefix.name}_inputs"
    input_dir.mkdir(exist_ok=True)
    status_file = Path(f"{prefix}.status.jsonl")
    report_file = Path(f"{prefix}.report.json")
    status_file.write_text("")
    status_lock = threading.Lock()

    def record(index, label, status, **details):
        _append_status(
            status_file,
            {"index": index, "label": label, "status": status, **details},
            status_lock,
        )

    def run_task(index, molecule, xyz_path):
        label = f"{manifest['base_label']}_idx{index}"
        stdout_path = Path(f"{prefix}_{index:06d}.out")
        stderr_path = Path(f"{prefix}_{index:06d}.err")
        try:
            molecule.write_xyz(str(xyz_path), mode="w")
            record(index, label, "running")
            args = _replace_task_structure(
                manifest["task_arguments"], xyz_path, label
            )
            python = (
                "import sys; from chemsmart.cli.run import run; "
                "run(sys.argv[1:])"
            )
            command = [
                "srun",
                "--exclusive",
                "--nodes=1",
                "--ntasks=1",
                f"--cpus-per-task={manifest['task_cores']}",
                f"--mem={manifest['task_memory_gb']}G",
            ]
            if manifest.get("task_gpus", 0):
                command.append(f"--gres=gpu:{manifest['task_gpus']}")
            command.extend([sys.executable, "-c", python, *args])
            environment = os.environ.copy()
            environment["OMP_NUM_THREADS"] = str(manifest["task_cores"])
            with (
                stdout_path.open("w") as stdout,
                stderr_path.open("w") as stderr,
            ):
                result = subprocess.run(
                    command,
                    stdout=stdout,
                    stderr=stderr,
                    env=environment,
                    check=False,
                )
            status = "successful" if result.returncode == 0 else "failed"
            record(
                index,
                label,
                status,
                returncode=result.returncode,
                stdout=str(stdout_path),
                stderr=str(stderr_path),
            )
            return result.returncode
        except Exception as error:
            record(index, label, "failed", error=str(error))
            return 1
        finally:
            try:
                xyz_path.unlink()
            except FileNotFoundError:
                pass

    failed = 0
    total = 0
    pending = set()
    with ThreadPoolExecutor(max_workers=workers) as executor:
        source = (
            XYZFile(manifest["filename"])
            if manifest.get("format") == ".xyz"
            else NPZFile(manifest["filename"])
        )
        for index, molecule, error in source.iter_molecule_results(
            manifest["index"]
        ):
            total += 1
            label = f"{manifest['base_label']}_idx{index}"
            if error is not None:
                record(index, label, "failed", error=str(error))
                failed += 1
                continue
            while len(pending) >= max_in_flight:
                done, pending = wait(pending, return_when=FIRST_COMPLETED)
                failed += sum(future.result() != 0 for future in done)
            xyz_path = input_dir / f"structure_{index:06d}.xyz"
            record(index, label, "pending")
            pending.add(executor.submit(run_task, index, molecule, xyz_path))
        for future in as_completed(pending):
            failed += future.result() != 0

    try:
        input_dir.rmdir()
    except OSError:
        pass
    successful = total - failed
    report = {
        "successful": successful,
        "failed": failed,
        "unfinished": 0,
        "total": total,
        "status_file": str(status_file),
    }
    temporary_report = Path(f"{report_file}.tmp")
    temporary_report.write_text(json.dumps(report, indent=2) + "\n")
    os.replace(temporary_report, report_file)
    print(
        "Streaming NPZ batch complete: "
        f"successful={successful}, failed={failed}, unfinished=0"
    )
    return int(failed > 0)


def _append_status(status_file, record, lock):
    with lock, open(status_file, "a") as stream:
        stream.write(json.dumps(record) + "\n")
        stream.flush()


def _safe_label(label):
    label = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(label or "task"))
    return label[:100] or "task"


def _batch_capacity(manifest):
    nodes = manifest["num_nodes"]
    cores_per_node = manifest["cores_per_node"]
    memory_per_node = manifest["memory_per_node_gb"]
    gpus_per_node = manifest.get("gpus_per_node", 0)
    task_cores = manifest["task_cores"]
    task_memory = manifest["task_memory_gb"]
    task_gpus = manifest.get("task_gpus", 0)

    values = {
        "node count": nodes,
        "cores per node": cores_per_node,
        "memory per node": memory_per_node,
        "task core count": task_cores,
        "task memory": task_memory,
    }
    if task_gpus:
        values["task GPU count"] = task_gpus
    for name, value in values.items():
        if not isinstance(value, int) or value <= 0:
            raise ValueError(f"{name} must be a positive integer.")
    if task_cores > cores_per_node:
        raise ValueError("A subtask cannot request more cores than one node.")
    if task_memory > memory_per_node:
        raise ValueError("A subtask cannot request more memory than one node.")
    if task_gpus > gpus_per_node:
        raise ValueError("A subtask cannot request more GPUs than one node.")

    per_node = min(
        cores_per_node // task_cores,
        memory_per_node // task_memory,
    )
    if task_gpus:
        per_node = min(per_node, gpus_per_node // task_gpus)
    return nodes * per_node


def run_slurm_batch(manifest_path):
    """Run manifest tasks concurrently as exclusive Slurm steps."""
    manifest_path = Path(manifest_path)
    with manifest_path.open() as stream:
        manifest = json.load(stream)

    tasks = manifest["tasks"]
    if not tasks:
        print("Batch complete: successful=0, failed=0, unfinished=0")
        return 0

    workers = min(len(tasks), _batch_capacity(manifest))
    prefix = manifest_path.with_suffix("")
    status_file = Path(f"{prefix}.status.jsonl")
    report_file = Path(f"{prefix}.report.json")
    status_file.write_text("")
    status_lock = threading.Lock()

    def update_status(index, task, status, **details):
        record = {
            "index": index,
            "label": task["label"],
            "status": status,
            **details,
        }
        _append_status(status_file, record, status_lock)

    for index, task in enumerate(tasks, start=1):
        label = _safe_label(task["label"])
        task["stdout"] = str(
            prefix.parent / f"{prefix.name}_{index:05d}_{label}.out"
        )
        task["stderr"] = str(
            prefix.parent / f"{prefix.name}_{index:05d}_{label}.err"
        )
        update_status(index, task, "pending")

    def run_task(index, task):
        update_status(index, task, "running")
        command = [
            "srun",
            "--exclusive",
            "--nodes=1",
            "--ntasks=1",
            f"--cpus-per-task={manifest['task_cores']}",
            f"--mem={manifest['task_memory_gb']}G",
        ]
        if manifest.get("task_gpus", 0):
            command.append(f"--gres=gpu:{manifest['task_gpus']}")
        command.extend([sys.executable, task["runscript"]])

        try:
            with (
                open(task["stdout"], "w") as stdout,
                open(task["stderr"], "w") as stderr,
            ):
                result = subprocess.run(
                    command,
                    check=False,
                    stdout=stdout,
                    stderr=stderr,
                )
            status = "successful" if result.returncode == 0 else "failed"
            details = {}
            if result.returncode != 0:
                details["error"] = (
                    f"srun exited with status {result.returncode}"
                )
        except Exception as error:
            status = "failed"
            details = {"error": str(error)}

        update_status(
            index,
            task,
            status,
            stdout=task["stdout"],
            stderr=task["stderr"],
            **details,
        )

    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(run_task, index, task)
            for index, task in enumerate(tasks, start=1)
        ]
        for future in as_completed(futures):
            future.result()

    latest = {}
    with status_file.open() as stream:
        for line in stream:
            record = json.loads(line)
            latest[record["index"]] = record

    records = [latest[index] for index in sorted(latest)]
    successful = sum(r["status"] == "successful" for r in records)
    failed = sum(r["status"] == "failed" for r in records)
    unfinished = len(tasks) - successful - failed
    report = {
        "successful": successful,
        "failed": failed,
        "unfinished": unfinished,
        "tasks": records,
    }
    temporary_report = Path(f"{report_file}.tmp")
    temporary_report.write_text(json.dumps(report, indent=2) + "\n")
    os.replace(temporary_report, report_file)
    print(
        "Batch complete: "
        f"successful={successful}, failed={failed}, unfinished={unfinished}"
    )
    logger.info(
        "Batch complete: successful=%s, failed=%s, unfinished=%s",
        successful,
        failed,
        unfinished,
    )
    return int(failed > 0 or unfinished > 0)
