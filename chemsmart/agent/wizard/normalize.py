"""Normalization helpers for wizard scheduler survey data."""

from __future__ import annotations

from chemsmart.agent.wizard.parsers import QueueFacts

SCHEDULER_SUBMIT = {
    "SLURM": "sbatch",
    "PBS": "qsub",
    "LSF": "bsub <",
    "SGE": "qsub",
}
FALLBACK_NUM_HOURS = 24
WALLTIME_SANITY_CEILING_HOURS = 168


def choose_queue(queues: list[QueueFacts]) -> str | None:
    """Choose the default queue name from parsed scheduler facts."""

    for queue in queues:
        if queue.default:
            return queue.name

    candidates = [
        queue
        for queue in queues
        if queue.enabled
        and queue.started
        and (queue.gpus_per_node is None or queue.gpus_per_node <= 0)
    ]
    if not candidates:
        return None
    return min(
        candidates, key=lambda queue: (len(queue.name), queue.name)
    ).name


def normalize_walltime(queue: QueueFacts) -> int:
    """Return a queue-aware walltime in hours for the selected queue."""

    if (
        queue.default_walltime_hours is not None
        and queue.default_walltime_hours > 0
    ):
        return queue.default_walltime_hours
    if queue.max_walltime_hours is not None and queue.max_walltime_hours > 0:
        return min(queue.max_walltime_hours, WALLTIME_SANITY_CEILING_HOURS)
    return FALLBACK_NUM_HOURS


def normalize_resources(queue: QueueFacts) -> dict[str, int | None]:
    """Normalize queue resources to server field names."""

    cores = queue.default_cores
    gpus = queue.gpus_per_node if (queue.gpus_per_node or 0) > 0 else 0
    return {
        "mem_gb": queue.default_mem_gb,
        "cores": cores,
        "gpus": gpus,
        "threads": cores,
    }
