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
    _batch_capacity(manifest)
    queue_capacity = int(manifest.get("queue_capacity", 100))
    if queue_capacity <= 0:
        raise ValueError("queue capacity must be a positive integer.")
    active_nodes = min(manifest["num_nodes"], queue_capacity)
    chunk_size = max(1, queue_capacity // active_nodes)
    prefix = manifest_path.with_suffix("")
    input_dir = prefix.parent / f"{prefix.name}_inputs"
    input_dir.mkdir(exist_ok=True)
    status_file = Path(f"{prefix}.status.jsonl")
    report_file = Path(f"{prefix}.report.json")
    previous_successes = set()
    if status_file.exists():
        latest = {}
        for line in status_file.read_text().splitlines():
            try:
                item = json.loads(line)
            except json.JSONDecodeError:
                continue
            latest[item["index"]] = item["status"]
        previous_successes = {
            index for index, status in latest.items() if status == "successful"
        }
    status_lock = threading.Lock()

    def record(index, label, status, **details):
        _append_status(
            status_file,
            {"index": index, "label": label, "status": status, **details},
            status_lock,
        )

    def run_chunk(number, entries):
        result_path = Path(f"{prefix}.chunk_{number:05d}.json")
        result_path.unlink(missing_ok=True)
        step_error = Path(f"{prefix}.chunk_{number:05d}.err")
        command = [
            "srun",
            "--exclusive",
            "--nodes=1",
            "--ntasks=1",
            f"--cpus-per-task={manifest['cores_per_node']}",
            f"--mem={manifest['memory_per_node_gb']}G",
        ]
        if manifest.get("gpus_per_node", 0):
            command.append(f"--gres=gpu:{manifest['gpus_per_node']}")
        command.extend(
            [
                sys.executable,
                "-c",
                "import sys; from chemsmart.jobs.batch import run_streaming_chunk; "
                "raise SystemExit(run_streaming_chunk(*sys.argv[1:]))",
                str(manifest_path),
                json.dumps(entries),
                str(result_path),
            ]
        )
        try:
            with step_error.open("w") as stderr:
                result = subprocess.run(command, stderr=stderr, check=False)
            if result.returncode != 0 or not result_path.exists():
                raise RuntimeError(f"Slurm step failed; see {step_error}")
            return json.loads(result_path.read_text())
        except Exception as error:
            logger.error("Streaming batch stopped: %s", error)
            return None

    failed = 0
    successful = 0
    total = 0
    pending = {}
    stopped = False
    chunk = []
    chunk_number = 0
    with ThreadPoolExecutor(max_workers=active_nodes) as executor:
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
            if index in previous_successes:
                successful += 1
                continue
            if error is not None:
                record(index, label, "failed", error=str(error))
                failed += 1
                continue
            while len(pending) >= active_nodes:
                done, _ = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    pending.pop(future)
                    results = future.result()
                    if results is None:
                        stopped = True
                        break
                    for item in results:
                        record(
                            item.pop("index"),
                            item.pop("label"),
                            item.pop("status"),
                            **item,
                        )
                        failed += item.get("returncode", 0) != 0
                        successful += item.get("returncode", 0) == 0
                if stopped:
                    break
            if stopped:
                break
            xyz_path = input_dir / f"structure_{index:06d}.xyz"
            molecule.write_xyz(str(xyz_path), mode="w")
            record(index, label, "pending")
            chunk.append([index, str(xyz_path)])
            if len(chunk) < chunk_size:
                continue
            chunk_number += 1
            pending[executor.submit(run_chunk, chunk_number, chunk)] = chunk
            chunk = []
        if chunk and not stopped:
            chunk_number += 1
            pending[executor.submit(run_chunk, chunk_number, chunk)] = chunk
        for future in as_completed(pending):
            results = future.result()
            if results is None:
                stopped = True
                continue
            for item in results:
                record(
                    item.pop("index"),
                    item.pop("label"),
                    item.pop("status"),
                    **item,
                )
                failed += item.get("returncode", 0) != 0
                successful += item.get("returncode", 0) == 0

    try:
        input_dir.rmdir()
    except OSError:
        pass
    unfinished = total - successful - failed
    report = {
        "successful": successful,
        "failed": failed,
        "unfinished": unfinished,
        "total": total,
        "status_file": str(status_file),
    }
    if stopped:
        report["stopped"] = "Slurm step failed; see chunk error log"
    temporary_report = Path(f"{report_file}.tmp")
    temporary_report.write_text(json.dumps(report, indent=2) + "\n")
    os.replace(temporary_report, report_file)
    print(
        "Streaming NPZ batch complete: "
        f"successful={successful}, failed={failed}, stopped={stopped}"
    )
    return int(failed > 0 or stopped)


def run_streaming_chunk(manifest_path, entries_json, result_path):
    """Compute prepared structures locally within a single Slurm step."""
    manifest = json.loads(Path(manifest_path).read_text())
    entries = json.loads(entries_json)
    workers = _batch_capacity(manifest) // manifest["num_nodes"]
    prefix = Path(manifest_path).with_suffix("")

    def run_one(entry):
        index, filename = entry
        label = f"{manifest['base_label']}_idx{index}"
        stdout_path = Path(f"{prefix}_{index:06d}.out")
        stderr_path = Path(f"{prefix}_{index:06d}.err")
        args = _replace_task_structure(
            manifest["task_arguments"], filename, label
        )
        environment = os.environ.copy()
        environment["OMP_NUM_THREADS"] = str(manifest["task_cores"])
        try:
            with (
                stdout_path.open("w") as stdout,
                stderr_path.open("w") as stderr,
            ):
                result = subprocess.run(
                    [
                        sys.executable,
                        "-c",
                        "import sys; from chemsmart.cli.run import run; run(sys.argv[1:])",
                        *args,
                    ],
                    stdout=stdout,
                    stderr=stderr,
                    env=environment,
                    check=False,
                )
            returncode = result.returncode
            details = {}
        except Exception as error:
            returncode = 1
            details = {"error": str(error)}
        finally:
            Path(filename).unlink(missing_ok=True)
        return {
            "index": index,
            "label": label,
            "status": "successful" if returncode == 0 else "failed",
            "returncode": returncode,
            "stdout": str(stdout_path),
            "stderr": str(stderr_path),
            **details,
        }

    with ThreadPoolExecutor(max_workers=workers) as executor:
        records = list(executor.map(run_one, entries))
    output = Path(result_path)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(json.dumps(records))
    os.replace(temporary, output)
    return 0


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
    """Run many manifest tasks inside a bounded number of Slurm steps."""
    manifest_path = Path(manifest_path)
    with manifest_path.open() as stream:
        manifest = json.load(stream)

    tasks = manifest["tasks"]
    if not tasks:
        print("Batch complete: successful=0, failed=0, unfinished=0")
        return 0

    _batch_capacity(manifest)
    prefix = manifest_path.with_suffix("")
    status_file = Path(f"{prefix}.status.jsonl")
    report_file = Path(f"{prefix}.report.json")
    status_lock = threading.Lock()

    # Keep successful work when the same batch is submitted again.  The status
    # journal is more current than the final report after an interrupted run.
    latest = {}
    if status_file.exists():
        with status_file.open() as stream:
            for line in stream:
                try:
                    record = json.loads(line)
                except json.JSONDecodeError:
                    continue
                latest[record["index"]] = record

    def update_status(index, task, status, **details):
        record = {
            "index": index,
            "label": task["label"],
            "status": status,
            **details,
        }
        _append_status(status_file, record, status_lock)

    todo = []
    for index, task in enumerate(tasks, start=1):
        label = _safe_label(task["label"])
        task["stdout"] = str(
            prefix.parent / f"{prefix.name}_{index:05d}_{label}.out"
        )
        task["stderr"] = str(
            prefix.parent / f"{prefix.name}_{index:05d}_{label}.err"
        )
        if latest.get(index, {}).get("status") == "successful":
            continue
        update_status(index, task, "pending")
        todo.append(index)

    # Each step owns one whole node and executes up to 128 molecules before
    # exiting, reducing the number of steps far below the molecule count.
    chunks = [todo[i : i + 128] for i in range(0, len(todo), 128)]

    def run_chunk(number, indices):
        result_path = Path(f"{prefix}.chunk_{number:05d}.json")
        result_path.unlink(missing_ok=True)
        step_error = Path(f"{prefix}.chunk_{number:05d}.err")
        command = [
            "srun",
            "--exclusive",
            "--nodes=1",
            "--ntasks=1",
            f"--cpus-per-task={manifest['cores_per_node']}",
            f"--mem={manifest['memory_per_node_gb']}G",
        ]
        if manifest.get("gpus_per_node", 0):
            command.append(f"--gres=gpu:{manifest['gpus_per_node']}")
        command.extend(
            [
                sys.executable,
                "-c",
                "import sys; from chemsmart.jobs.batch import run_batch_chunk; "
                "raise SystemExit(run_batch_chunk(*sys.argv[1:]))",
                str(manifest_path),
                json.dumps(indices),
                str(result_path),
            ]
        )
        try:
            with step_error.open("w") as stderr:
                result = subprocess.run(command, check=False, stderr=stderr)
            if result.returncode != 0 or not result_path.exists():
                raise RuntimeError(
                    f"Slurm step failed for chunk {number}; see {step_error}"
                )
            return json.loads(result_path.read_text())
        except Exception as error:
            logger.error("Batch execution stopped: %s", error)
            return None

    stopped = False
    with ThreadPoolExecutor(max_workers=manifest["num_nodes"]) as executor:
        pending = {}
        chunk_iterator = iter(enumerate(chunks, start=1))
        for _ in range(min(manifest["num_nodes"], len(chunks))):
            number, indices = next(chunk_iterator)
            pending[executor.submit(run_chunk, number, indices)] = indices
        while pending:
            done, _ = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                indices = pending.pop(future)
                records = future.result()
                if records is None:
                    stopped = True
                    continue
                for record in records:
                    _append_status(status_file, record, status_lock)
                if not stopped:
                    next_chunk = next(chunk_iterator, None)
                    if next_chunk is not None:
                        number, next_indices = next_chunk
                        pending[
                            executor.submit(run_chunk, number, next_indices)
                        ] = next_indices

    latest = {}
    with status_file.open() as stream:
        for line in stream:
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                continue
            latest[record["index"]] = record

    records = [latest[index] for index in range(1, len(tasks) + 1)]
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


def run_batch_chunk(manifest_path, indices_json, result_path):
    """Run a group of molecules locally within one allocated Slurm node."""
    manifest = json.loads(Path(manifest_path).read_text())
    indices = json.loads(indices_json)
    tasks = manifest["tasks"]
    workers = _batch_capacity(manifest) // manifest["num_nodes"]
    prefix = Path(manifest_path).with_suffix("")

    def run_one(index):
        task = tasks[index - 1]
        label = _safe_label(task["label"])
        stdout_path = prefix.parent / f"{prefix.name}_{index:05d}_{label}.out"
        stderr_path = prefix.parent / f"{prefix.name}_{index:05d}_{label}.err"
        details = {}
        try:
            with (
                stdout_path.open("w") as stdout,
                stderr_path.open("w") as stderr,
            ):
                result = subprocess.run(
                    [sys.executable, task["runscript"]],
                    check=False,
                    stdout=stdout,
                    stderr=stderr,
                )
            status = "successful" if result.returncode == 0 else "failed"
            if result.returncode:
                details["error"] = (
                    f"Task exited with status {result.returncode}"
                )
        except Exception as error:
            status = "failed"
            details["error"] = str(error)
        return {
            "index": index,
            "label": task["label"],
            "status": status,
            "stdout": str(stdout_path),
            "stderr": str(stderr_path),
            **details,
        }

    with ThreadPoolExecutor(max_workers=workers) as executor:
        records = list(executor.map(run_one, indices))
    output = Path(result_path)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(json.dumps(records))
    os.replace(temporary, output)
    return 0
