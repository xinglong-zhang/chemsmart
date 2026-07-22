"""Batch manifest helpers for scheduler array submissions.

Writes ``chemsmart_batch_<label>.json`` at submit time describing top-level
``BatchJob`` children and their per-task CLI args. Array execution uses the
rewritten CLI in each ``chemsmart_run_array_<task_id>.py`` runscript, not a
runtime manifest load.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Callable, Mapping, Optional, Sequence

logger = logging.getLogger(__name__)

RewriteCliFn = Callable[
    [Sequence[str], Optional[Mapping[str, Any]]],
    list[str],
]


def batch_manifest_filename(batch_label: str) -> str:
    """Return the manifest filename for *batch_label*."""
    return f"chemsmart_batch_{batch_label}.json"


def get_job_batch_entry(job: Any) -> Optional[dict[str, Any]]:
    """Return ``job.batch_entry`` when it is a mapping, else ``None``."""
    try:
        entry = job.batch_entry
    except AttributeError:
        return None
    if not isinstance(entry, Mapping):
        return None
    return dict(entry)


def set_job_batch_entry(job: Any, entry: Mapping[str, Any]) -> None:
    """Attach an explicit batch-entry mapping on *job*."""
    job.batch_entry = dict(entry)


def build_manifest_children(
    jobs: Sequence[Any],
    shared_cli_args: Sequence[str],
    rewrite_cli: Optional[RewriteCliFn] = None,
) -> list[dict[str, Any]]:
    """Build manifest child records (1-based task ids) from *jobs*.

    When a child has ``batch_entry`` and *rewrite_cli* is provided, store the
    rewritten per-task CLI. Otherwise store the shared CLI list.
    """
    children: list[dict[str, Any]] = []
    for task_id, job in enumerate(jobs, start=1):
        entry = get_job_batch_entry(job)
        child: dict[str, Any] = {
            "task_id": task_id,
            "label": job.label,
        }
        if entry is not None:
            child["batch_entry"] = dict(entry)
            if rewrite_cli is not None:
                child["cli_args"] = rewrite_cli(shared_cli_args, entry)
            else:
                child["cli_args"] = list(shared_cli_args)
        else:
            child["cli_args"] = list(shared_cli_args)
        children.append(child)
    return children


def write_batch_manifest(
    *,
    batch_label: str,
    program: str,
    children: Sequence[Mapping[str, Any]],
    directory: Optional[str | Path] = None,
) -> Path:
    """Write ``chemsmart_batch_<label>.json`` and return its path."""
    payload = {
        "batch_label": batch_label,
        "program": program,
        "children": [dict(child) for child in children],
    }
    path = Path(directory or ".") / batch_manifest_filename(batch_label)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n")
    logger.info("Wrote batch manifest: %s", path)
    return path


def load_batch_manifest(path: str | Path) -> dict[str, Any]:
    """Load a batch manifest JSON file."""
    with open(path, "r") as handle:
        return json.load(handle)


def load_batch_manifest_entry(
    path: str | Path,
    task_id: int,
) -> dict[str, Any]:
    """Return the manifest child record for 1-based *task_id*."""
    payload = load_batch_manifest(path)
    for child in payload.get("children", []):
        if int(child["task_id"]) == int(task_id):
            return dict(child)
    raise KeyError(f"No manifest entry for task_id={task_id} in {path}")


def resolve_array_cli_args(
    jobs: Sequence[Any],
    shared_cli_args: Sequence[str],
    rewrite_cli: Optional[RewriteCliFn] = None,
) -> list[str] | list[list[str]]:
    """Return shared CLI args or a per-task list when batch entries exist.

    Homogeneous batches (no ``batch_entry``) keep a single shared CLI list.
    Heterogeneous batches require *rewrite_cli* and return one CLI list per
    child for array runscripts.
    """
    entries = [get_job_batch_entry(job) for job in jobs]
    if not any(entry is not None for entry in entries):
        return list(shared_cli_args)
    if rewrite_cli is None:
        raise ValueError(
            "Heterogeneous BatchJob children have batch_entry but no "
            "rewrite_cli callback was provided for per-task CLI args."
        )
    return [rewrite_cli(shared_cli_args, entry) for entry in entries]
