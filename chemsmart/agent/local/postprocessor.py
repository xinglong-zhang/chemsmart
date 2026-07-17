"""Inference-time post-processing guardrails for local planner JSON.

The public entry point is :func:`postprocess`, which returns a repaired deep copy of
an already-parsed planner dictionary.  The function is intentionally conservative:
it fixes only deterministic v4 formatting/fidelity defects and does not add/remove
steps or otherwise reinterpret the requested workflow.
"""

from __future__ import annotations

import copy
import os
import re
from difflib import get_close_matches
from typing import Any

PATH_RE = re.compile(
    r"\b(?:examples|inputs)/[^\s,;:'\"`]+?\.xyz\b", re.IGNORECASE
)
BAD_STEP_FIELD_RE = re.compile(r"^\$step\d+\.(molecule|settings|job)$")
ALLOWED_RECOMMEND_REFS = {"$step2.functional", "$step2.basis"}

KIND_WHITELIST: tuple[str, ...] = (
    "gaussian.sp",
    "gaussian.opt",
    "gaussian.ts",
    "gaussian.freq",
    "gaussian.irc",
    "gaussian.scan",
    "gaussian.resp",
    "gaussian.nci",
    "gaussian.dias",
    "gaussian.tddft",
    "gaussian.modred",
    "gaussian.crest",
    "gaussian.traj",
    "gaussian.qmmm",
    "gaussian.qrc",
    "gaussian.wbi",
    "orca.sp",
    "orca.opt",
    "orca.ts",
    "orca.irc",
    "orca.scan",
    "orca.modred",
    "orca.neb",
    "orca.qmmm",
    "orca.qrc",
)

KIND_SHORT: dict[str, str] = {
    "gaussian.opt": "gopt",
    "gaussian.sp": "gsp",
    "gaussian.ts": "gts",
    "gaussian.freq": "gfreq",
    "gaussian.irc": "girc",
    "gaussian.scan": "gscan",
    "gaussian.resp": "gresp",
    "gaussian.nci": "gnci",
    "gaussian.dias": "gdias",
    "gaussian.tddft": "gtddft",
    "gaussian.modred": "gmodred",
    "gaussian.crest": "gcrest",
    "gaussian.traj": "gtraj",
    "gaussian.qmmm": "gqmmm",
    "gaussian.qrc": "gqrc",
    "gaussian.wbi": "gwbi",
    "orca.opt": "oopt",
    "orca.sp": "osp",
    "orca.ts": "ots",
    "orca.irc": "oirc",
    "orca.scan": "oscan",
    "orca.modred": "omodred",
    "orca.neb": "oneb",
    "orca.qmmm": "oqmmm",
    "orca.qrc": "oqrc",
}


def postprocess(plan: dict[str, Any], user_query: str) -> dict[str, Any]:
    """Return a repaired copy of a chemsmart planner output.

    Args:
        plan: Parsed planner JSON dictionary.
        user_query: Original user request, used as the authority for molecule
            file paths.

    Returns:
        A new dictionary with deterministic filepath, step-reference, label, and
        kind-whitelist repairs applied.  The input object is never mutated.
    """

    repaired = copy.deepcopy(plan)
    jobs = repaired.get("jobs")
    if isinstance(jobs, list):
        _repair_compact_job_filepaths(jobs, user_query)

    steps = repaired.get("steps")
    if not isinstance(steps, list):
        return repaired

    _repair_filepaths(steps, user_query)
    _repair_step_refs(repaired)
    _repair_job_kinds_and_labels(steps)
    return repaired


def _repair_compact_job_filepaths(jobs: list[Any], user_query: str) -> None:
    paths = PATH_RE.findall(user_query or "")
    if not paths:
        return

    source_job_by_id: dict[int, dict[str, Any]] = {}
    for index, job in enumerate(jobs):
        if not isinstance(job, dict):
            continue
        job_id = job.get("id")
        if isinstance(job_id, int):
            source_job_by_id[job_id] = job
        if isinstance(job.get("file"), str):
            job["file"] = paths[index] if index < len(paths) else paths[0]

    for job in jobs:
        if not isinstance(job, dict) or "file" in job:
            continue
        geom_from = job.get("geom_from")
        if isinstance(geom_from, int) and geom_from in source_job_by_id:
            continue
        job["file"] = paths[0]


def _repair_filepaths(steps: list[Any], user_query: str) -> None:
    paths = PATH_RE.findall(user_query or "")
    if not paths:
        return

    build_molecule_steps = [
        step
        for step in steps
        if isinstance(step, dict)
        and step.get("tool") == "build_molecule"
        and isinstance(step.get("args"), dict)
    ]
    for index, step in enumerate(build_molecule_steps):
        path = paths[index] if index < len(paths) else paths[0]
        step["args"]["filepath"] = path


def _repair_step_refs(value: Any) -> Any:
    if isinstance(value, dict):
        for key, child in list(value.items()):
            value[key] = _repair_step_refs(child)
        return value
    if isinstance(value, list):
        for index, child in enumerate(value):
            value[index] = _repair_step_refs(child)
        return value
    if (
        isinstance(value, str)
        and BAD_STEP_FIELD_RE.match(value)
        and value not in ALLOWED_RECOMMEND_REFS
    ):
        return value.split(".", 1)[0]
    return value


def _repair_job_kinds_and_labels(steps: list[Any]) -> None:
    filepaths_by_step: dict[int, str] = {}
    job_indices: list[int] = []

    for one_based_index, step in enumerate(steps, start=1):
        if not isinstance(step, dict) or not isinstance(
            step.get("args"), dict
        ):
            continue
        args = step["args"]
        if step.get("tool") == "build_molecule" and isinstance(
            args.get("filepath"), str
        ):
            filepaths_by_step[one_based_index] = args["filepath"]
        elif step.get("tool") == "build_job":
            args["kind"] = _nearest_valid_kind(args.get("kind"))
            job_indices.append(one_based_index)

    for one_based_index in job_indices:
        step = steps[one_based_index - 1]
        args = step["args"]
        kind = args.get("kind", "gaussian.opt")
        filepath = _filepath_for_job(args, filepaths_by_step)
        if filepath:
            stem = os.path.splitext(os.path.basename(filepath))[0].lower()
        else:
            stem = "molecule"
        args["label"] = f"{stem}_{KIND_SHORT.get(kind, 'gopt')}_001"


def _filepath_for_job(
    args: dict[str, Any], filepaths_by_step: dict[int, str]
) -> str | None:
    for key in ("molecule", "reactant", "product"):
        ref = args.get(key)
        if isinstance(ref, str):
            match = re.match(r"^\$step(\d+)$", ref)
            if match and int(match.group(1)) in filepaths_by_step:
                return filepaths_by_step[int(match.group(1))]
    if filepaths_by_step:
        return filepaths_by_step[min(filepaths_by_step)]
    return None


def _nearest_valid_kind(kind: Any) -> str:
    if not isinstance(kind, str):
        return "gaussian.opt"
    normalized = kind.strip().lower()
    if normalized in KIND_WHITELIST:
        return normalized

    prefix_matches = [
        valid
        for valid in KIND_WHITELIST
        if valid.startswith(normalized) or normalized.startswith(valid)
    ]
    if prefix_matches:
        return prefix_matches[0]

    if "." in normalized:
        program, suffix = normalized.split(".", 1)
        suffix_matches = [
            valid
            for valid in KIND_WHITELIST
            if valid.startswith(f"{program}.")
            and valid.split(".", 1)[1].startswith(suffix)
        ]
        if suffix_matches:
            return suffix_matches[0]

    close = get_close_matches(normalized, KIND_WHITELIST, n=1, cutoff=0.82)
    return close[0] if close else "gaussian.opt"
