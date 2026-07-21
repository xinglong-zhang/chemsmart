"""Gaussian-specific generated-input invariant checks."""

from __future__ import annotations

import re
from typing import Any

import yaml

from chemsmart.agent.harness.generated_common import (
    check_requested_numbers,
    coordinate_groups,
    coordinate_prefix,
    expand_range,
    numeric_sequence,
    reject,
    sequence_value,
)
from chemsmart.agent.harness.models import InvariantIssue
from chemsmart.agent.harness.scan_values import scan_value_matches
from chemsmart.settings.workspace_project import workspace_project_path


def coordinate_in_route(route: str) -> bool:
    opt_blocks = re.findall(
        r"\bopt\s*=\s*\(([^)]*)\)",
        route,
        re.IGNORECASE,
    )
    pattern = re.compile(
        r"(?:^|,)\s*[BAD]\s*,\s*\d+(?:\s*,\s*\d+){1,3}"
        r"\s*,\s*[FS](?:\s*,|\s*$)",
        re.IGNORECASE,
    )
    return any(pattern.search(block) is not None for block in opt_blocks)


def tddft_issues(
    route: str,
    chemistry: dict[str, Any],
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    issues: list[InvariantIssue] = []
    match = re.search(r"\bTD\s*\(([^)]*)\)", route, flags=re.IGNORECASE)
    if not match:
        return [
            reject(
                "input.gaussian.tddft.route",
                "Gaussian TD-DFT input is missing a TD(...) route block",
                evidence,
            )
        ]
    body = match.group(1).lower()
    for key in ("states", "nstates", "root"):
        expected = chemistry.get(key)
        if expected is None:
            continue
        if key == "states":
            present = str(expected).lower() in body
        else:
            present = (
                re.search(
                    rf"\b{key}\s*=\s*{re.escape(str(expected))}\b",
                    body,
                )
                is not None
            )
        if not present:
            issues.append(
                reject(
                    f"input.gaussian.tddft.{key}",
                    f"TD-DFT route does not preserve requested {key}",
                    {**evidence, "expected": expected},
                )
            )
    if chemistry.get("eqsolv") is True and "eqsolv" not in body:
        issues.append(
            reject(
                "input.gaussian.tddft.eqsolv",
                "TD-DFT route does not contain requested eqsolv",
                evidence,
            )
        )
    return issues


def coordinate_issues(
    kind: str,
    content: str,
    chemistry: dict[str, Any],
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    rows = [
        line.strip()
        for line in content.splitlines()
        if re.match(r"^[BAD]\s+\d", line.strip(), flags=re.IGNORECASE)
    ]
    issues: list[InvariantIssue] = []
    required = "S" if kind.endswith(".scan") else "F"
    if not any(
        re.search(rf"\s{required}(?:\s|$)", row, flags=re.IGNORECASE)
        for row in rows
    ):
        issues.append(
            reject(
                f"input.{kind}.directive",
                (
                    f"generated Gaussian {kind.rsplit('.', 1)[1]} input "
                    f"lacks a {required} coordinate directive"
                ),
                {**evidence, "coordinate_rows": rows},
            )
        )
        return issues
    expected_coordinates = coordinate_groups(chemistry.get("coordinates"))
    if expected_coordinates:
        directive_rows = [
            row
            for row in rows
            if re.search(
                rf"\s{required}(?:\s|$)",
                row,
                flags=re.IGNORECASE,
            )
        ]
        for group in expected_coordinates:
            if not _coordinate_present(directive_rows, group):
                issues.append(
                    reject(
                        f"input.{kind}.coordinate_atoms",
                        (
                            "generated Gaussian coordinate directive does "
                            "not preserve requested atoms"
                        ),
                        {
                            **evidence,
                            "expected": group,
                            "coordinate_rows": directive_rows,
                        },
                    )
                )
    if kind.endswith(".scan"):
        issues.extend(_scan_value_issues(rows, chemistry, evidence))
    else:
        check_requested_numbers(
            issues,
            rows,
            chemistry,
            evidence,
            prefix=f"input.{kind}",
        )
    constraints = coordinate_groups(chemistry.get("constrained_coordinates"))
    if constraints:
        frozen_rows = [
            row
            for row in rows
            if re.search(r"\sF(?:\s|$)", row, flags=re.IGNORECASE)
        ]
        for group in constraints:
            if not _coordinate_present(frozen_rows, group):
                issues.append(
                    reject(
                        f"input.{kind}.constraint",
                        (
                            "requested frozen coordinate is absent from "
                            "generated input"
                        ),
                        {
                            **evidence,
                            "expected": group,
                            "coordinate_rows": frozen_rows,
                        },
                    )
                )
    return issues


def qmmm_issues(
    route: str,
    chemistry: dict[str, Any],
    generated: dict[str, Any],
    evidence: dict[str, Any],
    *,
    project: str | None,
    cwd: str | None,
) -> list[InvariantIssue]:
    issues = _qmmm_partition_issues(chemistry, generated, evidence)
    issues.extend(_qmmm_frequency_policy_issues(route, project, cwd, evidence))
    return issues


def dias_issues(
    chemistry: dict[str, Any],
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    fragment = chemistry.get("fragment_indices")
    indices = expand_range(str(fragment)) if fragment else []
    if not indices:
        return [
            reject(
                "input.gaussian.dias.fragment_partition",
                "Gaussian DIAS command has no valid fragment-1 atom partition",
                evidence,
            )
        ]
    return []


def ts_extra_issues(
    route: str,
    chemistry: dict[str, Any],
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    issues: list[InvariantIssue] = []
    normalized_route = re.sub(r"\s+", "", route).lower()
    for key in ("route_parameters", "opt_options"):
        value = chemistry.get(key)
        if not isinstance(value, str) or not value.strip():
            continue
        normalized_value = re.sub(r"\s+", "", value).lower()
        if normalized_value not in normalized_route:
            issues.append(
                reject(
                    f"input.gaussian.ts.{key}",
                    (
                        "Gaussian TS route does not preserve a requested "
                        "extra keyword"
                    ),
                    {**evidence, "expected": value, "source": key},
                )
            )
    return issues


def _scan_value_issues(
    rows: list[str],
    chemistry: dict[str, Any],
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    parsed_rows = _scan_rows(rows)
    issues: list[InvariantIssue] = []
    groups = coordinate_groups(chemistry.get("coordinates"))
    steps = numeric_sequence(chemistry.get("num_steps"))
    increments = numeric_sequence(chemistry.get("step_size"))
    for index, group in enumerate(groups):
        tag = {2: "B", 3: "A", 4: "D"}.get(len(group))
        row = next(
            (
                candidate
                for candidate in parsed_rows
                if candidate["tag"] == tag and candidate["atoms"] == group
            ),
            None,
        )
        if row is None:
            continue
        for key, expected_values, observed in (
            ("num_steps", steps, row["num_steps"]),
            ("step_size", increments, row["step_size"]),
        ):
            expected = sequence_value(expected_values, index)
            if expected is None:
                continue
            if not scan_value_matches(key, expected, observed):
                issues.append(
                    reject(
                        f"input.gaussian.scan.{key}",
                        (
                            "generated Gaussian scan row does not preserve "
                            f"requested {key}"
                        ),
                        {
                            **evidence,
                            "expected": expected,
                            "observed": observed,
                            "coordinate": group,
                            "scan_rows": parsed_rows,
                        },
                    )
                )
    return issues


def _scan_rows(rows: list[str]) -> list[dict[str, Any]]:
    parsed: list[dict[str, Any]] = []
    pattern = re.compile(
        r"^([BAD])\s+([0-9]+(?:\s+[0-9]+){1,3})\s+S\s+"
        r"(\d+)\s+([-+0-9.eE]+)\b",
        flags=re.IGNORECASE,
    )
    for row in rows:
        match = pattern.match(row)
        if match is None:
            continue
        parsed.append(
            {
                "tag": match.group(1).upper(),
                "atoms": tuple(int(value) for value in match.group(2).split()),
                "num_steps": int(match.group(3)),
                "step_size": float(match.group(4)),
            }
        )
    return parsed


def _coordinate_present(
    rows: list[str],
    group: tuple[int, ...],
) -> bool:
    prefix = coordinate_prefix(group)
    return bool(prefix) and any(
        re.match(
            rf"^{re.escape(prefix)}(?:\s|$)",
            row,
            flags=re.IGNORECASE,
        )
        for row in rows
    )


def _qmmm_frequency_policy_issues(
    route: str,
    project: str | None,
    cwd: str | None,
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    """Honor an explicit ONIOM frequency policy in workspace YAML."""

    if not project:
        return []
    path = workspace_project_path(project, "gaussian", cwd=cwd)
    if not path.is_file():
        return []
    try:
        payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError):
        return []
    qmmm = payload.get("qmmm") if isinstance(payload, dict) else None
    policy = qmmm.get("freq") if isinstance(qmmm, dict) else None
    if not isinstance(policy, bool):
        return []
    has_frequency = (
        re.search(r"\bfreq(?:\b|=)", route, flags=re.IGNORECASE) is not None
    )
    if has_frequency == policy:
        return []
    expected = "include" if policy else "omit"
    return [
        reject(
            "input.gaussian.qmmm.frequency_policy",
            (
                "generated Gaussian ONIOM route must "
                f"{expected} frequency analysis per qmmm.freq"
            ),
            {
                **evidence,
                "project": project,
                "project_path": str(path),
                "qmmm_freq": policy,
            },
        )
    ]


def _qmmm_partition_issues(
    chemistry: dict[str, Any],
    generated: dict[str, Any],
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    issues: list[InvariantIssue] = []
    atom_layers = generated.get("atom_layers")
    layer_atoms = generated.get("layer_atoms")
    if not isinstance(atom_layers, list) or not atom_layers:
        return [
            reject(
                "input.gaussian.qmmm.layer_markers",
                (
                    "generated Gaussian ONIOM coordinates lack "
                    "atom-layer markers"
                ),
                evidence,
            )
        ]
    if any(layer not in {"H", "M", "L"} for layer in atom_layers):
        issues.append(
            reject(
                "input.gaussian.qmmm.partition_coverage",
                (
                    "every generated Gaussian ONIOM atom must have one "
                    "H/M/L layer marker"
                ),
                {**evidence, "atom_layers": atom_layers},
            )
        )
    if not isinstance(layer_atoms, dict):
        layer_atoms = {
            layer: [
                index
                for index, observed in enumerate(atom_layers, start=1)
                if observed == layer
            ]
            for layer in ("H", "M", "L")
        }

    expected_high = set(
        expand_range(str(chemistry.get("high_level_atoms") or ""))
    )
    expected_medium = set(
        expand_range(str(chemistry.get("medium_level_atoms") or ""))
    )
    explicit_low = chemistry.get("low_level_atoms")
    expected_low = set(expand_range(str(explicit_low or "")))
    all_atoms = set(range(1, len(atom_layers) + 1))
    if explicit_low in (None, "", []):
        expected_low = all_atoms - expected_high - expected_medium
    for label, marker, expected in (
        ("high", "H", expected_high),
        ("medium", "M", expected_medium),
        ("low", "L", expected_low),
    ):
        observed = set(int(index) for index in layer_atoms.get(marker, []))
        if observed != expected:
            issues.append(
                reject(
                    f"input.gaussian.qmmm.{label}_level_atoms",
                    (
                        f"generated Gaussian ONIOM {label}-level partition "
                        "does not preserve the command"
                    ),
                    {
                        **evidence,
                        "expected": sorted(expected),
                        "observed": sorted(observed),
                    },
                )
            )

    pairs = generated.get("charge_multiplicity_pairs")
    if isinstance(pairs, list) and pairs:
        issues.extend(_qmmm_state_pair_issues(chemistry, pairs, evidence))
    return issues


def _qmmm_state_pair_issues(
    chemistry: dict[str, Any],
    pairs: list[Any],
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    issues: list[InvariantIssue] = []

    def assert_pair(
        label: str,
        pair_index: int,
        charge_key: str,
        multiplicity_key: str,
    ) -> None:
        expected_charge = chemistry.get(charge_key)
        expected_mult = chemistry.get(multiplicity_key)
        if expected_charge is None and expected_mult is None:
            return
        observed = pairs[pair_index] if pair_index < len(pairs) else None
        expected = [
            int(expected_charge) if expected_charge is not None else None,
            int(expected_mult) if expected_mult is not None else None,
        ]
        valid = isinstance(observed, list) and len(observed) == 2
        if valid and expected[0] is not None:
            valid = int(observed[0]) == expected[0]
        if valid and expected[1] is not None:
            valid = int(observed[1]) == expected[1]
        if not valid:
            issues.append(
                reject(
                    f"input.gaussian.qmmm.{label}_state",
                    (
                        f"generated Gaussian ONIOM {label} "
                        "charge/multiplicity does not preserve the command"
                    ),
                    {
                        **evidence,
                        "expected": expected,
                        "observed": observed,
                    },
                )
            )

    assert_pair("total", 0, "charge_total", "mult_total")
    has_medium = bool(
        expand_range(str(chemistry.get("medium_level_atoms") or ""))
    )
    if has_medium:
        assert_pair(
            "intermediate",
            1,
            "charge_intermediate",
            "mult_intermediate",
        )
        assert_pair("high", 3, "charge_high", "mult_high")
    else:
        assert_pair("high", 1, "charge_high", "mult_high")
    return issues


__all__ = [
    "coordinate_in_route",
    "coordinate_issues",
    "dias_issues",
    "qmmm_issues",
    "tddft_issues",
    "ts_extra_issues",
]
