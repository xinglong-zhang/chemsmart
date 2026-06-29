from __future__ import annotations

import re
from typing import Any

from chemsmart.agent.harness.extractors import extract_gaussian_route
from chemsmart.agent.harness.invariants.gaussian_ts import (
    check_gaussian_ts_route,
)
from chemsmart.agent.harness.models import HarnessResult, InvariantResult

_STEP_REF_RE = re.compile(r"^\$step(?P<index>\d+)(?:\..*)?$")


def evaluate_harness(
    plan: Any,
    dry_run_results: list[dict[str, Any]],
) -> HarnessResult:
    rule_results: list[InvariantResult] = []
    for result_index, (kind, result) in enumerate(
        _iter_dry_run_kinds(plan, dry_run_results)
    ):
        if kind != "gaussian.ts":
            continue
        route = extract_gaussian_route(result.get("content"))
        rule_results.append(
            check_gaussian_ts_route(
                route,
                inputfile=_string_or_none(result.get("inputfile")),
                result_index=result_index,
            )
        )
    return HarnessResult.from_rule_results(rule_results)


def _iter_dry_run_kinds(
    plan: Any,
    dry_run_results: list[dict[str, Any]],
) -> list[tuple[str | None, dict[str, Any]]]:
    if isinstance(plan, dict):
        steps = list(plan.get("steps", []) or [])
    else:
        steps = list(getattr(plan, "steps", []) or [])
    dry_run_steps = [
        step for step in steps if _step_tool(step) == "dry_run_input"
    ]
    pairs: list[tuple[str | None, dict[str, Any]]] = []
    for dry_run_step, result in zip(dry_run_steps, dry_run_results):
        pairs.append((_dry_run_kind(steps, dry_run_step), result))
    return pairs


def _dry_run_kind(steps: list[Any], dry_run_step: Any) -> str | None:
    job_ref = _step_args(dry_run_step).get("job")
    if not isinstance(job_ref, str):
        return None
    match = _STEP_REF_RE.match(job_ref)
    if match is None:
        return None
    job_index = int(match.group("index")) - 1
    if job_index < 0 or job_index >= len(steps):
        return None
    build_job_step = steps[job_index]
    if _step_tool(build_job_step) != "build_job":
        return None
    kind = _step_args(build_job_step).get("kind")
    if not isinstance(kind, str):
        return None
    return kind.strip().lower()


def _step_tool(step: Any) -> str | None:
    if isinstance(step, dict):
        value = step.get("tool")
    else:
        value = getattr(step, "tool", None)
    return value if isinstance(value, str) else None


def _step_args(step: Any) -> dict[str, Any]:
    if isinstance(step, dict):
        value = step.get("args")
    else:
        value = getattr(step, "args", None)
    return value if isinstance(value, dict) else {}


def _string_or_none(value: Any) -> str | None:
    return value if isinstance(value, str) else None
