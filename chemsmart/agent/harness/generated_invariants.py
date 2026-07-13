"""Kind-specific checks over safely generated Gaussian and ORCA inputs."""

from __future__ import annotations

import ast
import re
from typing import Any

from chemsmart.agent.harness.intent import ObservedIntent
from chemsmart.agent.harness.invariants.gaussian_ts import check_gaussian_ts_route
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
    kind = observed.kind or ""
    issues: list[InvariantIssue] = []
    for index, generated in enumerate(generated_inputs):
        route = str(generated.get("route") or "")
        content = str(generated.get("content_tail") or "")
        evidence = {
            "kind": kind,
            "path": generated.get("path"),
            "route": route,
            "result_index": index,
        }
        if kind == "gaussian.ts":
            issues.extend(check_gaussian_ts_route(route, inputfile=str(generated.get("path") or ""), result_index=index).issues)
            issues.extend(_gaussian_ts_extra_issues(route, observed.chemistry, evidence))
        elif kind == "orca.ts":
            issues.extend(check_orca_ts_route(route, inputfile=str(generated.get("path") or ""), result_index=index).issues)
        if kind.endswith(".irc"):
            issues.extend(check_irc_route(route, software=kind.split(".", 1)[0], inputfile=str(generated.get("path") or ""), result_index=index).issues)
        if kind.endswith(".sp") and re.search(r"\b(?:opt|freq|optts|scants)\b", route, flags=re.IGNORECASE):
            issues.append(_reject("input.sp.unrequested_route", "single-point input contains an optimization/frequency keyword", evidence))
        if kind == "gaussian.tddft":
            issues.extend(_gaussian_td_issues(route, observed.chemistry, evidence))
        if kind in {"gaussian.scan", "gaussian.modred"}:
            issues.extend(_gaussian_coordinate_issues(kind, content, observed.chemistry, evidence))
        if kind in {"orca.scan", "orca.modred"}:
            issues.extend(_orca_coordinate_issues(kind, content, observed.chemistry, evidence))
        if kind == "orca.neb":
            issues.extend(_orca_neb_issues(route, content, observed.chemistry, evidence))
        if kind.endswith(".qmmm"):
            issues.extend(_qmmm_issues(kind, route, content, observed.chemistry, evidence))
        if kind == "gaussian.dias":
            issues.extend(_dias_issues(observed.chemistry, evidence))
    return tuple(issues)


def _gaussian_td_issues(route: str, chemistry: dict[str, Any], evidence: dict[str, Any]) -> list[InvariantIssue]:
    issues: list[InvariantIssue] = []
    match = re.search(r"\bTD\s*\(([^)]*)\)", route, flags=re.IGNORECASE)
    if not match:
        return [_reject("input.gaussian.tddft.route", "Gaussian TD-DFT input is missing a TD(...) route block", evidence)]
    body = match.group(1).lower()
    for key in ("states", "nstates", "root"):
        expected = chemistry.get(key)
        if expected is None:
            continue
        if key == "states":
            present = str(expected).lower() in body
        else:
            present = re.search(rf"\b{key}\s*=\s*{re.escape(str(expected))}\b", body) is not None
        if not present:
            issues.append(_reject(f"input.gaussian.tddft.{key}", f"TD-DFT route does not preserve requested {key}", {**evidence, "expected": expected}))
    if chemistry.get("eqsolv") is True and "eqsolv" not in body:
        issues.append(_reject("input.gaussian.tddft.eqsolv", "TD-DFT route does not contain requested eqsolv", evidence))
    return issues


def _gaussian_coordinate_issues(kind: str, content: str, chemistry: dict[str, Any], evidence: dict[str, Any]) -> list[InvariantIssue]:
    rows = [line.strip() for line in content.splitlines() if re.match(r"^[BAD]\s+\d", line.strip(), flags=re.IGNORECASE)]
    issues: list[InvariantIssue] = []
    required = "S" if kind.endswith(".scan") else "F"
    if not any(re.search(rf"\s{required}(?:\s|$)", row, flags=re.IGNORECASE) for row in rows):
        issues.append(_reject(f"input.{kind}.directive", f"generated Gaussian {kind.rsplit('.', 1)[1]} input lacks a {required} coordinate directive", {**evidence, "coordinate_rows": rows}))
        return issues
    expected_coordinates = _coordinate_groups(chemistry.get("coordinates"))
    if expected_coordinates:
        directive_rows = [
            row
            for row in rows
            if re.search(rf"\s{required}(?:\s|$)", row, flags=re.IGNORECASE)
        ]
        for group in expected_coordinates:
            if not _gaussian_coordinate_present(directive_rows, group):
                issues.append(
                    _reject(
                        f"input.{kind}.coordinate_atoms",
                        "generated Gaussian coordinate directive does not preserve requested atoms",
                        {**evidence, "expected": group, "coordinate_rows": directive_rows},
                    )
                )
    _check_requested_numbers(issues, rows, chemistry, evidence, prefix=f"input.{kind}")
    constraints = _coordinate_groups(chemistry.get("constrained_coordinates"))
    if constraints:
        frozen_rows = [
            row
            for row in rows
            if re.search(r"\sF(?:\s|$)", row, flags=re.IGNORECASE)
        ]
        for group in constraints:
            if not _gaussian_coordinate_present(frozen_rows, group):
                issues.append(
                    _reject(
                        f"input.{kind}.constraint",
                        "requested frozen coordinate is absent from generated input",
                        {**evidence, "expected": group, "coordinate_rows": frozen_rows},
                    )
                )
    return issues


def _orca_coordinate_issues(kind: str, content: str, chemistry: dict[str, Any], evidence: dict[str, Any]) -> list[InvariantIssue]:
    lowered = content.lower()
    issues: list[InvariantIssue] = []
    if kind.endswith(".scan") and "scan" not in lowered:
        issues.append(_reject("input.orca.scan.block", "generated ORCA scan input lacks a Scan directive", evidence))
    if kind.endswith(".modred") and not any(marker in lowered for marker in ("constraints", "{ b", "{ a", "{ d")):
        issues.append(_reject("input.orca.modred.constraints", "generated ORCA modred input lacks a constraints block", evidence))
    expected_coordinates = _coordinate_groups(chemistry.get("coordinates"))
    if expected_coordinates:
        for group in expected_coordinates:
            if not _orca_coordinate_present(content, group):
                issues.append(
                    _reject(
                        f"input.{kind}.coordinate_atoms",
                        "generated ORCA coordinate block does not preserve requested atoms",
                        {**evidence, "expected": group},
                    )
                )
    _check_requested_numbers(issues, [content], chemistry, evidence, prefix=f"input.{kind}")
    constraints = _coordinate_groups(chemistry.get("constrained_coordinates"))
    if constraints:
        for group in constraints:
            if not _orca_frozen_coordinate_present(content, group):
                issues.append(
                    _reject(
                        f"input.{kind}.constraint",
                        "requested ORCA frozen coordinate is absent from generated input",
                        {**evidence, "expected": group},
                    )
                )
    return issues


def _orca_neb_issues(route: str, content: str, chemistry: dict[str, Any], evidence: dict[str, Any]) -> list[InvariantIssue]:
    issues: list[InvariantIssue] = []
    expected_job = chemistry.get("joboption")
    if expected_job and str(expected_job).lower() not in route.lower():
        issues.append(_reject("input.orca.neb.joboption", "ORCA route does not preserve requested NEB variant", {**evidence, "expected": expected_job}))
    endpoint = chemistry.get("ending_xyzfile")
    if endpoint and str(endpoint) not in content:
        issues.append(_reject("input.orca.neb.endpoint", "ORCA NEB block does not preserve the requested product endpoint", {**evidence, "expected": endpoint}))
    nimages = chemistry.get("nimages")
    if nimages is not None and re.search(rf"\bNImages\s+{re.escape(str(nimages))}\b", content, flags=re.IGNORECASE) is None:
        issues.append(_reject("input.orca.neb.nimages", "ORCA NEB block does not preserve requested image count", {**evidence, "expected": nimages}))
    return issues


def _qmmm_issues(kind: str, route: str, content: str, chemistry: dict[str, Any], evidence: dict[str, Any]) -> list[InvariantIssue]:
    issues: list[InvariantIssue] = []
    lowered = f"{route}\n{content}".lower()
    expected_marker = "oniom" if kind.startswith("gaussian") else "qmmm"
    if expected_marker not in lowered:
        issues.append(_reject(f"input.{kind}.marker", f"generated {kind} input lacks the {expected_marker.upper()} marker", evidence))
    high = chemistry.get("high_level_atoms")
    low = chemistry.get("low_level_atoms")
    if high and low:
        overlap = set(_expand_range(str(high))) & set(_expand_range(str(low)))
        if overlap:
            issues.append(_reject(f"input.{kind}.region_overlap", "QM/MM high- and low-level regions overlap", {**evidence, "overlap": sorted(overlap)}))
    if kind == "orca.qmmm":
        low_method = chemistry.get("low_level_method")
        if low_method and str(low_method).lower() not in lowered:
            issues.append(
                _reject(
                    "input.orca.qmmm.low_level_method",
                    "generated ORCA QM/MM input does not preserve the requested low-level method",
                    {**evidence, "expected": low_method},
                )
            )
    return issues


def _dias_issues(chemistry: dict[str, Any], evidence: dict[str, Any]) -> list[InvariantIssue]:
    fragment = chemistry.get("fragment_indices")
    indices = _expand_range(str(fragment)) if fragment else []
    if not indices:
        return [
            _reject(
                "input.gaussian.dias.fragment_partition",
                "Gaussian DIAS command has no valid fragment-1 atom partition",
                evidence,
            )
        ]
    return []


def _gaussian_ts_extra_issues(route: str, chemistry: dict[str, Any], evidence: dict[str, Any]) -> list[InvariantIssue]:
    issues: list[InvariantIssue] = []
    normalized_route = re.sub(r"\s+", "", route).lower()
    for key in ("route_parameters", "opt_options"):
        value = chemistry.get(key)
        if not isinstance(value, str) or not value.strip():
            continue
        normalized_value = re.sub(r"\s+", "", value).lower()
        if normalized_value not in normalized_route:
            issues.append(
                _reject(
                    f"input.gaussian.ts.{key}",
                    "Gaussian TS route does not preserve a requested extra keyword",
                    {**evidence, "expected": value, "source": key},
                )
            )
    return issues


def _check_requested_numbers(issues: list[InvariantIssue], chunks: list[str], chemistry: dict[str, Any], evidence: dict[str, Any], *, prefix: str) -> None:
    text = "\n".join(chunks)
    numbers = re.findall(r"-?\d+(?:\.\d+)?", text)
    for key in ("num_steps", "step_size", "dist_start", "dist_end"):
        expected = chemistry.get(key)
        expected_numbers = re.findall(r"-?\d+(?:\.\d+)?", str(expected))
        if expected is not None and not all(
            number in numbers for number in expected_numbers
        ):
            issues.append(_reject(f"{prefix}.{key}", f"generated input does not preserve requested {key}", {**evidence, "expected": expected}))


def _expand_range(value: str) -> list[int]:
    values: list[int] = []
    for item in re.split(r"[,\s]+", value.strip()):
        if not item:
            continue
        if "-" in item:
            start, end = item.split("-", 1)
            if start.isdigit() and end.isdigit():
                values.extend(range(int(start), int(end) + 1))
        elif item.isdigit():
            values.append(int(item))
    return values


def _coordinate_groups(value: Any) -> list[tuple[int, ...]]:
    if value is None:
        return []
    try:
        parsed = ast.literal_eval(value) if isinstance(value, str) else value
    except (SyntaxError, ValueError):
        return []
    if not isinstance(parsed, list) or not parsed:
        return []
    groups = parsed if isinstance(parsed[0], list) else [parsed]
    normalized: list[tuple[int, ...]] = []
    for group in groups:
        if isinstance(group, list) and all(
            isinstance(atom, int) and not isinstance(atom, bool) for atom in group
        ):
            normalized.append(tuple(group))
    return normalized


def _coordinate_prefix(group: tuple[int, ...], *, orca: bool = False) -> str:
    tag = {2: "B", 3: "A", 4: "D"}.get(len(group), "")
    atoms = [atom - 1 for atom in group] if orca else list(group)
    return f"{tag} {' '.join(str(atom) for atom in atoms)}".strip()


def _gaussian_coordinate_present(rows: list[str], group: tuple[int, ...]) -> bool:
    prefix = _coordinate_prefix(group)
    return bool(prefix) and any(
        re.match(rf"^{re.escape(prefix)}(?:\s|$)", row, flags=re.IGNORECASE)
        for row in rows
    )


def _orca_coordinate_present(content: str, group: tuple[int, ...]) -> bool:
    prefix = _coordinate_prefix(group, orca=True)
    return bool(prefix) and re.search(
        rf"{re.escape(prefix)}\s*=",
        content,
        flags=re.IGNORECASE,
    ) is not None


def _orca_frozen_coordinate_present(content: str, group: tuple[int, ...]) -> bool:
    prefix = _coordinate_prefix(group, orca=True)
    return bool(prefix) and re.search(
        rf"\{{\s*{re.escape(prefix)}\s+C\s*\}}",
        content,
        flags=re.IGNORECASE,
    ) is not None


def _reject(rule_id: str, message: str, evidence: dict[str, Any]) -> InvariantIssue:
    return InvariantIssue(rule_id=rule_id, severity="reject", message=message, evidence=evidence)


__all__ = ["check_generated_input_invariants"]
