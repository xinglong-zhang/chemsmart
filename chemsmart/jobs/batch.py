import json
import logging
import os
import re
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

logger = logging.getLogger(__name__)


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
