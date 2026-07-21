"""Shared helpers for generated-input invariant checks."""

from __future__ import annotations

import ast
import re
from typing import Any

from chemsmart.agent.harness.models import InvariantIssue
from chemsmart.utils.periodictable import PeriodicTable

_PERIODIC_TABLE = PeriodicTable()


def reject(
    rule_id: str,
    message: str,
    evidence: dict[str, Any],
) -> InvariantIssue:
    return InvariantIssue(
        rule_id=rule_id,
        severity="reject",
        message=message,
        evidence=evidence,
    )


def electron_multiplicity_issues(
    generated: dict[str, Any],
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    state = electron_multiplicity_evidence(generated)
    if state is None or state["valid"]:
        return []
    return [
        reject(
            "input.state.electron_multiplicity_parity",
            (
                "generated input charge, composition, and multiplicity "
                "have incompatible electron-count parity"
            ),
            {**evidence, **state},
        )
    ]


def electron_multiplicity_evidence(
    generated: dict[str, Any],
) -> dict[str, Any] | None:
    """Return an auditable parity receipt for an explicit geometry."""

    charge = generated.get("charge")
    multiplicity = generated.get("multiplicity")
    counts = generated.get("element_counts")
    if not isinstance(charge, int) or not isinstance(multiplicity, int):
        return None
    if not isinstance(counts, dict) or not counts:
        return None
    try:
        electrons = (
            sum(
                _PERIODIC_TABLE.to_atomic_number(str(symbol)) * int(count)
                for symbol, count in counts.items()
            )
            - charge
        )
    except (TypeError, ValueError):
        return None
    return {
        "charge": charge,
        "multiplicity": multiplicity,
        "element_counts": dict(sorted(counts.items())),
        "electron_count": electrons,
        "valid": (
            electrons >= 0
            and multiplicity >= 1
            and (electrons + multiplicity) % 2 == 1
        ),
    }


def numeric_sequence(value: Any) -> list[float]:
    if value is None:
        return []
    return [
        float(token) for token in re.findall(r"-?\d+(?:\.\d+)?", str(value))
    ]


def sequence_value(values: list[float], index: int) -> float | None:
    if not values:
        return None
    if len(values) == 1:
        return values[0]
    return values[index] if index < len(values) else None


def check_requested_numbers(
    issues: list[InvariantIssue],
    chunks: list[str],
    chemistry: dict[str, Any],
    evidence: dict[str, Any],
    *,
    prefix: str,
) -> None:
    text = "\n".join(chunks)
    numbers = re.findall(r"-?\d+(?:\.\d+)?", text)
    for key in ("num_steps", "step_size", "dist_start", "dist_end"):
        expected = chemistry.get(key)
        expected_numbers = re.findall(
            r"-?\d+(?:\.\d+)?",
            str(expected),
        )
        if expected is not None and not all(
            number in numbers for number in expected_numbers
        ):
            issues.append(
                reject(
                    f"{prefix}.{key}",
                    f"generated input does not preserve requested {key}",
                    {**evidence, "expected": expected},
                )
            )


def expand_range(value: str) -> list[int]:
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


def coordinate_groups(value: Any) -> list[tuple[int, ...]]:
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
            isinstance(atom, int) and not isinstance(atom, bool)
            for atom in group
        ):
            normalized.append(tuple(group))
    return normalized


def coordinate_prefix(
    group: tuple[int, ...],
    *,
    orca: bool = False,
) -> str:
    tag = {2: "B", 3: "A", 4: "D"}.get(len(group), "")
    atoms = [atom - 1 for atom in group] if orca else list(group)
    return f"{tag} {' '.join(str(atom) for atom in atoms)}".strip()


__all__ = [
    "check_requested_numbers",
    "coordinate_groups",
    "coordinate_prefix",
    "electron_multiplicity_evidence",
    "electron_multiplicity_issues",
    "expand_range",
    "numeric_sequence",
    "reject",
    "sequence_value",
]
