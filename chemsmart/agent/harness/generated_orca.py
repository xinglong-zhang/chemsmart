"""ORCA-specific generated-input invariant checks."""

from __future__ import annotations

import re
from typing import Any

from chemsmart.agent.harness.generated_common import (
    check_requested_numbers,
    coordinate_groups,
    coordinate_prefix,
    numeric_sequence,
    reject,
    sequence_value,
)
from chemsmart.agent.harness.models import InvariantIssue
from chemsmart.agent.harness.scan_values import scan_value_matches


def coordinate_issues(
    kind: str,
    content: str,
    chemistry: dict[str, Any],
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    lowered = content.lower()
    issues: list[InvariantIssue] = []
    if kind.endswith(".scan") and "scan" not in lowered:
        issues.append(
            reject(
                "input.orca.scan.block",
                "generated ORCA scan input lacks a Scan directive",
                evidence,
            )
        )
    if kind.endswith(".modred") and not any(
        marker in lowered for marker in ("constraints", "{ b", "{ a", "{ d")
    ):
        issues.append(
            reject(
                "input.orca.modred.constraints",
                "generated ORCA modred input lacks a constraints block",
                evidence,
            )
        )
    if kind.endswith(".modred"):
        issues.extend(
            _modred_route_issues(str(evidence.get("route") or ""), evidence)
        )

    expected_coordinates = coordinate_groups(chemistry.get("coordinates"))
    if expected_coordinates:
        issues.extend(
            _coordinate_atom_issues(
                kind,
                content,
                expected_coordinates,
                evidence,
            )
        )
        if kind.endswith(".modred"):
            issues.extend(
                _constraint_set_issues(
                    content,
                    expected_coordinates,
                    evidence,
                )
            )
        issues.extend(_xyz_state_issues(kind, content, chemistry, evidence))

    if kind.endswith(".scan"):
        issues.extend(_scan_value_issues(content, chemistry, evidence))
    else:
        check_requested_numbers(
            issues,
            [content],
            chemistry,
            evidence,
            prefix=f"input.{kind}",
        )
    issues.extend(
        _frozen_coordinate_issues(kind, content, chemistry, evidence)
    )
    return issues


def neb_issues(
    route: str,
    content: str,
    chemistry: dict[str, Any],
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    issues: list[InvariantIssue] = []
    expected_job = chemistry.get("joboption")
    if expected_job and str(expected_job).lower() not in route.lower():
        issues.append(
            reject(
                "input.orca.neb.joboption",
                "ORCA route does not preserve requested NEB variant",
                {**evidence, "expected": expected_job},
            )
        )
    endpoint = chemistry.get("ending_xyzfile")
    if endpoint and str(endpoint) not in content:
        issues.append(
            reject(
                "input.orca.neb.endpoint",
                (
                    "ORCA NEB block does not preserve the requested product "
                    "endpoint"
                ),
                {**evidence, "expected": endpoint},
            )
        )
    nimages = chemistry.get("nimages")
    if (
        nimages is not None
        and re.search(
            rf"\bNImages\s+{re.escape(str(nimages))}\b",
            content,
            flags=re.IGNORECASE,
        )
        is None
    ):
        issues.append(
            reject(
                "input.orca.neb.nimages",
                "ORCA NEB block does not preserve requested image count",
                {**evidence, "expected": nimages},
            )
        )
    return issues


def qmmm_issues(
    route: str,
    content: str,
    chemistry: dict[str, Any],
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    low_method = chemistry.get("low_level_method")
    lowered = f"{route}\n{content}".lower()
    if low_method and str(low_method).lower() not in lowered:
        return [
            reject(
                "input.orca.qmmm.low_level_method",
                (
                    "generated ORCA QM/MM input does not preserve the "
                    "requested low-level method"
                ),
                {**evidence, "expected": low_method},
            )
        ]
    return []


def _modred_route_issues(
    route: str,
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    issues: list[InvariantIssue] = []
    if re.search(r"\bopt\b", route, flags=re.IGNORECASE) is None:
        issues.append(
            reject(
                "input.orca.modred.route_opt",
                "generated ORCA modred route is missing the required Opt keyword",
                evidence,
            )
        )
    forbidden = sorted(
        {
            match.group(0).lower()
            for match in re.finditer(
                r"\b(?:scan|optts|scants|irc)\b",
                route,
                flags=re.IGNORECASE,
            )
        }
    )
    if forbidden:
        issues.append(
            reject(
                "input.orca.modred.forbidden_route",
                (
                    "generated ORCA modred route contains a scan, TS, or "
                    "IRC keyword"
                ),
                {**evidence, "forbidden": forbidden},
            )
        )
    return issues


def _coordinate_atom_issues(
    kind: str,
    content: str,
    expected_coordinates: list[tuple[int, ...]],
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    issues: list[InvariantIssue] = []
    for group in expected_coordinates:
        present = (
            _frozen_coordinate_present(content, group)
            if kind.endswith(".modred")
            else _coordinate_present(content, group)
        )
        if not present:
            issues.append(
                reject(
                    f"input.{kind}.coordinate_atoms",
                    (
                        "generated ORCA coordinate block does not preserve "
                        "requested atoms"
                    ),
                    {**evidence, "expected": group},
                )
            )
    return issues


def _constraint_set_issues(
    content: str,
    expected_coordinates: list[tuple[int, ...]],
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    expected_rows = [
        ("BAD"[len(group) - 2], tuple(index - 1 for index in group))
        for group in expected_coordinates
    ]
    observed_rows, unexpected_rows = _constraint_rows(content)
    if observed_rows == expected_rows and not unexpected_rows:
        return []
    return [
        reject(
            "input.orca.modred.constraint_set",
            (
                "generated ORCA constraints do not exactly preserve the "
                "requested internal-coordinate rows and order"
            ),
            {
                **evidence,
                "expected_rows": expected_rows,
                "observed_rows": observed_rows,
                "unexpected_rows": unexpected_rows,
            },
        )
    ]


def _xyz_state_issues(
    kind: str,
    content: str,
    chemistry: dict[str, Any],
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    charge = chemistry.get("charge")
    multiplicity = chemistry.get("multiplicity")
    if charge is None or multiplicity is None:
        return []
    if (
        re.search(
            rf"^\*\s+xyz\s+{re.escape(str(charge))}\s+"
            rf"{re.escape(str(multiplicity))}\s*$",
            content,
            flags=re.IGNORECASE | re.MULTILINE,
        )
        is not None
    ):
        return []
    return [
        reject(
            f"input.{kind}.xyz_state",
            ("generated ORCA input does not preserve charge and multiplicity"),
            {
                **evidence,
                "charge": charge,
                "multiplicity": multiplicity,
            },
        )
    ]


def _frozen_coordinate_issues(
    kind: str,
    content: str,
    chemistry: dict[str, Any],
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    issues: list[InvariantIssue] = []
    constraints = coordinate_groups(chemistry.get("constrained_coordinates"))
    for group in constraints:
        if not _frozen_coordinate_present(content, group):
            issues.append(
                reject(
                    f"input.{kind}.constraint",
                    (
                        "requested ORCA frozen coordinate is absent from "
                        "generated input"
                    ),
                    {**evidence, "expected": group},
                )
            )
    return issues


def _scan_value_issues(
    content: str,
    chemistry: dict[str, Any],
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    """Assert scan values on the exact generated ORCA coordinate row."""

    rows = _scan_rows(content)
    issues: list[InvariantIssue] = []
    groups = coordinate_groups(chemistry.get("coordinates"))
    starts = numeric_sequence(chemistry.get("dist_start"))
    ends = numeric_sequence(chemistry.get("dist_end"))
    steps = numeric_sequence(chemistry.get("num_steps"))
    for index, group in enumerate(groups):
        tag = {2: "B", 3: "A", 4: "D"}.get(len(group))
        expected_atoms = tuple(atom - 1 for atom in group)
        row = next(
            (
                candidate
                for candidate in rows
                if candidate["tag"] == tag
                and candidate["atoms"] == expected_atoms
            ),
            None,
        )
        if row is None:
            continue
        for key, expected_values, observed in (
            ("dist_start", starts, row["dist_start"]),
            ("dist_end", ends, row["dist_end"]),
            ("num_steps", steps, row["num_steps"]),
        ):
            expected = sequence_value(expected_values, index)
            if expected is None:
                continue
            if not scan_value_matches(key, expected, observed):
                issues.append(
                    reject(
                        f"input.orca.scan.{key}",
                        (
                            "generated ORCA scan row does not preserve "
                            f"requested {key}"
                        ),
                        {
                            **evidence,
                            "expected": expected,
                            "observed": observed,
                            "coordinate": group,
                            "scan_rows": rows,
                        },
                    )
                )
    return issues


def _scan_rows(content: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    pattern = re.compile(
        r"^\s*([BAD])\s+([0-9]+(?:\s+[0-9]+){1,3})\s*=\s*"
        r"([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*,\s*(\d+)\b",
        flags=re.IGNORECASE | re.MULTILINE,
    )
    for match in pattern.finditer(content):
        rows.append(
            {
                "tag": match.group(1).upper(),
                "atoms": tuple(int(value) for value in match.group(2).split()),
                "dist_start": float(match.group(3)),
                "dist_end": float(match.group(4)),
                "num_steps": int(match.group(5)),
            }
        )
    return rows


def _constraint_rows(
    content: str,
) -> tuple[list[tuple[str, tuple[int, ...]]], list[str]]:
    """Parse every ORCA brace constraint, including unsupported extras."""

    rows: list[tuple[str, tuple[int, ...]]] = []
    unexpected: list[str] = []
    for match in re.finditer(r"\{([^{}\n]+)\}", content):
        raw = " ".join(match.group(1).split())
        parsed = re.fullmatch(
            r"([BAD])\s+([0-9]+(?:\s+[0-9]+)*)\s+C",
            raw,
            flags=re.IGNORECASE,
        )
        if parsed is None:
            unexpected.append("{" + raw + "}")
            continue
        rows.append(
            (
                parsed.group(1).upper(),
                tuple(int(value) for value in parsed.group(2).split()),
            )
        )
    return rows, unexpected


def _coordinate_present(content: str, group: tuple[int, ...]) -> bool:
    prefix = coordinate_prefix(group, orca=True)
    return (
        bool(prefix)
        and re.search(
            rf"{re.escape(prefix)}\s*=",
            content,
            flags=re.IGNORECASE,
        )
        is not None
    )


def _frozen_coordinate_present(
    content: str,
    group: tuple[int, ...],
) -> bool:
    prefix = coordinate_prefix(group, orca=True)
    return (
        bool(prefix)
        and re.search(
            rf"\{{\s*{re.escape(prefix)}\s+C\s*\}}",
            content,
            flags=re.IGNORECASE,
        )
        is not None
    )


__all__ = ["coordinate_issues", "neb_issues", "qmmm_issues"]
