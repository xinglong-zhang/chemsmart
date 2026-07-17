"""Public facade for checks over safely generated chemistry inputs."""

from __future__ import annotations

import re
from typing import Any

from chemsmart.agent.harness.generated_common import (
    electron_multiplicity_evidence,
    electron_multiplicity_issues,
    expand_range,
    reject,
)
from chemsmart.agent.harness.generated_gaussian import (
    coordinate_in_route as gaussian_coordinate_in_route,
)
from chemsmart.agent.harness.generated_gaussian import (
    coordinate_issues as gaussian_coordinate_issues,
)
from chemsmart.agent.harness.generated_gaussian import (
    dias_issues as gaussian_dias_issues,
)
from chemsmart.agent.harness.generated_gaussian import (
    qmmm_issues as gaussian_qmmm_issues,
)
from chemsmart.agent.harness.generated_gaussian import (
    tddft_issues as gaussian_tddft_issues,
)
from chemsmart.agent.harness.generated_gaussian import (
    ts_extra_issues as gaussian_ts_extra_issues,
)
from chemsmart.agent.harness.generated_orca import (
    coordinate_issues as orca_coordinate_issues,
)
from chemsmart.agent.harness.generated_orca import (
    neb_issues as orca_neb_issues,
)
from chemsmart.agent.harness.generated_orca import (
    qmmm_issues as orca_qmmm_issues,
)
from chemsmart.agent.harness.intent import ObservedIntent
from chemsmart.agent.harness.invariants.gaussian_ts import (
    check_gaussian_ts_route,
)
from chemsmart.agent.harness.invariants.route_checks import (
    check_irc_route,
    check_orca_ts_route,
)
from chemsmart.agent.harness.models import InvariantIssue


def check_generated_input_invariants(
    command: str,
    generated_inputs: list[dict[str, Any]] | tuple[dict[str, Any], ...],
    *,
    cwd: str | None = None,
) -> tuple[InvariantIssue, ...]:
    observed = ObservedIntent.from_command(command, cwd=cwd)
    issues: list[InvariantIssue] = []
    for index, generated in enumerate(generated_inputs):
        issues.extend(
            _generated_input_issues(
                observed,
                generated,
                index=index,
                cwd=cwd,
            )
        )
    return tuple(issues)


def _generated_input_issues(
    observed: ObservedIntent,
    generated: dict[str, Any],
    *,
    index: int,
    cwd: str | None,
) -> list[InvariantIssue]:
    kind = observed.kind or ""
    route = str(generated.get("route") or "")
    content = str(generated.get("content_tail") or "")
    evidence = {
        "kind": kind,
        "path": generated.get("path"),
        "route": route,
        "result_index": index,
    }
    issues = electron_multiplicity_issues(generated, evidence)
    issues.extend(_route_issues(kind, route, observed.chemistry, evidence))
    issues.extend(
        _coordinate_issues(
            kind,
            route,
            content,
            observed,
            evidence,
        )
    )
    issues.extend(
        _workflow_issues(
            kind,
            route,
            content,
            observed,
            generated,
            evidence,
            cwd=cwd,
        )
    )
    return issues


def _route_issues(
    kind: str,
    route: str,
    chemistry: dict[str, Any],
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    issues: list[InvariantIssue] = []
    path = str(evidence.get("path") or "")
    result_index = int(evidence.get("result_index") or 0)
    if kind == "gaussian.ts":
        issues.extend(
            check_gaussian_ts_route(
                route,
                inputfile=path,
                result_index=result_index,
            ).issues
        )
        issues.extend(gaussian_ts_extra_issues(route, chemistry, evidence))
    elif kind == "orca.ts":
        issues.extend(
            check_orca_ts_route(
                route,
                inputfile=path,
                result_index=result_index,
            ).issues
        )
    if kind.endswith(".irc"):
        issues.extend(
            check_irc_route(
                route,
                software=kind.split(".", 1)[0],
                inputfile=path,
                result_index=result_index,
            ).issues
        )
    if kind.endswith(".sp") and re.search(
        r"\b(?:opt|freq|optts|scants)\b",
        route,
        flags=re.IGNORECASE,
    ):
        issues.append(
            reject(
                "input.sp.unrequested_route",
                (
                    "single-point input contains an optimization/frequency "
                    "keyword"
                ),
                evidence,
            )
        )
    return issues


def _coordinate_issues(
    kind: str,
    route: str,
    content: str,
    observed: ObservedIntent,
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    chemistry = observed.chemistry
    if kind in {"gaussian.scan", "gaussian.modred"}:
        issues: list[InvariantIssue] = []
        if gaussian_coordinate_in_route(route):
            issues.append(
                reject(
                    "input.gaussian.coordinate_in_route",
                    (
                        "Gaussian ModRedundant coordinate directives must be "
                        "input-body rows, not opt route options"
                    ),
                    evidence,
                )
            )
        issues.extend(
            gaussian_coordinate_issues(kind, content, chemistry, evidence)
        )
        return issues
    if kind in {"orca.scan", "orca.modred"}:
        return orca_coordinate_issues(
            kind,
            content,
            {
                **chemistry,
                "charge": observed.charge,
                "multiplicity": observed.multiplicity,
            },
            evidence,
        )
    return []


def _workflow_issues(
    kind: str,
    route: str,
    content: str,
    observed: ObservedIntent,
    generated: dict[str, Any],
    evidence: dict[str, Any],
    *,
    cwd: str | None,
) -> list[InvariantIssue]:
    chemistry = observed.chemistry
    issues: list[InvariantIssue] = []
    if kind == "gaussian.tddft":
        issues.extend(gaussian_tddft_issues(route, chemistry, evidence))
    if kind == "orca.neb":
        issues.extend(orca_neb_issues(route, content, chemistry, evidence))
    if kind.endswith(".qmmm"):
        issues.extend(
            _qmmm_issues(
                kind,
                route,
                content,
                observed,
                generated,
                evidence,
                cwd=cwd,
            )
        )
    if kind == "gaussian.dias":
        issues.extend(gaussian_dias_issues(chemistry, evidence))
    return issues


def _qmmm_issues(
    kind: str,
    route: str,
    content: str,
    observed: ObservedIntent,
    generated: dict[str, Any],
    evidence: dict[str, Any],
    *,
    cwd: str | None,
) -> list[InvariantIssue]:
    chemistry = observed.chemistry
    issues: list[InvariantIssue] = []
    lowered = f"{route}\n{content}".lower()
    expected_marker = "oniom" if kind.startswith("gaussian") else "qmmm"
    if expected_marker not in lowered:
        issues.append(
            reject(
                f"input.{kind}.marker",
                (
                    f"generated {kind} input lacks the "
                    f"{expected_marker.upper()} marker"
                ),
                evidence,
            )
        )
    high = chemistry.get("high_level_atoms")
    low = chemistry.get("low_level_atoms")
    if high and low:
        overlap = set(expand_range(str(high))) & set(expand_range(str(low)))
        if overlap:
            issues.append(
                reject(
                    f"input.{kind}.region_overlap",
                    "QM/MM high- and low-level regions overlap",
                    {**evidence, "overlap": sorted(overlap)},
                )
            )
    if kind == "orca.qmmm":
        issues.extend(orca_qmmm_issues(route, content, chemistry, evidence))
    elif kind == "gaussian.qmmm":
        issues.extend(
            gaussian_qmmm_issues(
                route,
                chemistry,
                generated,
                evidence,
                project=observed.project,
                cwd=cwd,
            )
        )
    return issues


__all__ = [
    "check_generated_input_invariants",
    "electron_multiplicity_evidence",
]
