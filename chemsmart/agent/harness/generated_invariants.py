"""Kind-specific checks over safely generated Gaussian and ORCA inputs."""

from __future__ import annotations

import ast
import re
from typing import Any

import yaml

from chemsmart.agent.harness.intent import ObservedIntent
from chemsmart.agent.harness.invariants.gaussian_ts import check_gaussian_ts_route
from chemsmart.agent.harness.invariants.route_checks import (
    check_irc_route,
    check_orca_ts_route,
)
from chemsmart.agent.harness.models import InvariantIssue
from chemsmart.settings.workspace_project import workspace_project_path
from chemsmart.utils.periodictable import PeriodicTable


_PERIODIC_TABLE = PeriodicTable()


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
        issues.extend(_electron_multiplicity_issues(generated, evidence))
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
            issues.extend(
                _orca_coordinate_issues(
                    kind,
                    content,
                    {
                        **observed.chemistry,
                        "charge": observed.charge,
                        "multiplicity": observed.multiplicity,
                    },
                    evidence,
                )
            )
        if kind == "orca.neb":
            issues.extend(_orca_neb_issues(route, content, observed.chemistry, evidence))
        if kind.endswith(".qmmm"):
            issues.extend(
                _qmmm_issues(
                    kind,
                    route,
                    content,
                    observed.chemistry,
                    generated,
                    evidence,
                    project=observed.project,
                    cwd=cwd,
                )
            )
        if kind == "gaussian.dias":
            issues.extend(_dias_issues(observed.chemistry, evidence))
    return tuple(issues)


def _electron_multiplicity_issues(
    generated: dict[str, Any],
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    state = electron_multiplicity_evidence(generated)
    if state is None or state["valid"]:
        return []
    return [
        _reject(
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
    """Return an auditable parity receipt for an explicit generated geometry."""

    charge = generated.get("charge")
    multiplicity = generated.get("multiplicity")
    counts = generated.get("element_counts")
    if not isinstance(charge, int) or not isinstance(multiplicity, int):
        return None
    if not isinstance(counts, dict) or not counts:
        return None
    try:
        electrons = sum(
            _PERIODIC_TABLE.to_atomic_number(str(symbol)) * int(count)
            for symbol, count in counts.items()
        ) - charge
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
    if kind.endswith(".scan"):
        issues.extend(
            _gaussian_scan_value_issues(
                rows,
                chemistry,
                evidence,
            )
        )
    else:
        _check_requested_numbers(
            issues,
            rows,
            chemistry,
            evidence,
            prefix=f"input.{kind}",
        )
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


def _gaussian_scan_value_issues(
    rows: list[str],
    chemistry: dict[str, Any],
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    parsed_rows = _gaussian_scan_rows(rows)
    issues: list[InvariantIssue] = []
    groups = _coordinate_groups(chemistry.get("coordinates"))
    steps = _numeric_sequence(chemistry.get("num_steps"))
    increments = _numeric_sequence(chemistry.get("step_size"))
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
            expected = _sequence_value(expected_values, index)
            if expected is None:
                continue
            matches = (
                int(expected) == int(observed)
                if key == "num_steps"
                else abs(float(expected) - float(observed)) <= 1e-9
            )
            if not matches:
                issues.append(
                    _reject(
                        f"input.gaussian.scan.{key}",
                        f"generated Gaussian scan row does not preserve requested {key}",
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


def _gaussian_scan_rows(rows: list[str]) -> list[dict[str, Any]]:
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


def _orca_coordinate_issues(kind: str, content: str, chemistry: dict[str, Any], evidence: dict[str, Any]) -> list[InvariantIssue]:
    lowered = content.lower()
    issues: list[InvariantIssue] = []
    if kind.endswith(".scan") and "scan" not in lowered:
        issues.append(_reject("input.orca.scan.block", "generated ORCA scan input lacks a Scan directive", evidence))
    if kind.endswith(".modred") and not any(marker in lowered for marker in ("constraints", "{ b", "{ a", "{ d")):
        issues.append(_reject("input.orca.modred.constraints", "generated ORCA modred input lacks a constraints block", evidence))
    if kind.endswith(".modred"):
        route = str(evidence.get("route") or "")
        if re.search(r"\bopt\b", route, flags=re.IGNORECASE) is None:
            issues.append(
                _reject(
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
                _reject(
                    "input.orca.modred.forbidden_route",
                    "generated ORCA modred route contains a scan, TS, or IRC keyword",
                    {**evidence, "forbidden": forbidden},
                )
            )
    expected_coordinates = _coordinate_groups(chemistry.get("coordinates"))
    if expected_coordinates:
        for group in expected_coordinates:
            present = (
                _orca_frozen_coordinate_present(content, group)
                if kind.endswith(".modred")
                else _orca_coordinate_present(content, group)
            )
            if not present:
                issues.append(
                    _reject(
                        f"input.{kind}.coordinate_atoms",
                        "generated ORCA coordinate block does not preserve requested atoms",
                        {**evidence, "expected": group},
                    )
                )
        if kind.endswith(".modred"):
            expected_rows = [
                ("BAD"[len(group) - 2], tuple(index - 1 for index in group))
                for group in expected_coordinates
            ]
            observed_rows, unexpected_rows = _orca_constraint_rows(content)
            if observed_rows != expected_rows or unexpected_rows:
                issues.append(
                    _reject(
                        "input.orca.modred.constraint_set",
                        (
                            "generated ORCA constraints do not exactly preserve "
                            "the requested internal-coordinate rows and order"
                        ),
                        {
                            **evidence,
                            "expected_rows": expected_rows,
                            "observed_rows": observed_rows,
                            "unexpected_rows": unexpected_rows,
                        },
                    )
                )
        charge = chemistry.get("charge")
        multiplicity = chemistry.get("multiplicity")
        if charge is not None and multiplicity is not None and re.search(
            rf"^\*\s+xyz\s+{re.escape(str(charge))}\s+{re.escape(str(multiplicity))}\s*$",
            content,
            flags=re.IGNORECASE | re.MULTILINE,
        ) is None:
            issues.append(
                _reject(
                    f"input.{kind}.xyz_state",
                    "generated ORCA input does not preserve charge and multiplicity",
                    {**evidence, "charge": charge, "multiplicity": multiplicity},
                )
            )
    if kind.endswith(".scan"):
        issues.extend(
            _orca_scan_value_issues(
                content,
                chemistry,
                evidence,
            )
        )
    else:
        _check_requested_numbers(
            issues,
            [content],
            chemistry,
            evidence,
            prefix=f"input.{kind}",
        )
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


def _orca_scan_value_issues(
    content: str,
    chemistry: dict[str, Any],
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    """Assert ORCA scan values on the exact generated coordinate row.

    A global number search can be satisfied accidentally by atom coordinates,
    comments, or memory settings.  Parse the ``%geom Scan`` row so a requested
    13-point scan cannot pass when the generated input silently uses 12.
    """

    rows = _orca_scan_rows(content)
    issues: list[InvariantIssue] = []
    groups = _coordinate_groups(chemistry.get("coordinates"))
    starts = _numeric_sequence(chemistry.get("dist_start"))
    ends = _numeric_sequence(chemistry.get("dist_end"))
    steps = _numeric_sequence(chemistry.get("num_steps"))
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
            expected = _sequence_value(expected_values, index)
            if expected is None:
                continue
            matches = (
                int(expected) == int(observed)
                if key == "num_steps"
                else abs(float(expected) - float(observed)) <= 1e-9
            )
            if not matches:
                issues.append(
                    _reject(
                        f"input.orca.scan.{key}",
                        f"generated ORCA scan row does not preserve requested {key}",
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


def _orca_scan_rows(content: str) -> list[dict[str, Any]]:
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


def _numeric_sequence(value: Any) -> list[float]:
    if value is None:
        return []
    return [float(token) for token in re.findall(r"-?\d+(?:\.\d+)?", str(value))]


def _sequence_value(values: list[float], index: int) -> float | None:
    if not values:
        return None
    if len(values) == 1:
        return values[0]
    return values[index] if index < len(values) else None


def _orca_constraint_rows(
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


def _orca_frozen_coordinate_rows(content: str) -> list[tuple[str, tuple[int, ...]]]:
    """Return the valid internal-coordinate subset for compatibility."""

    return _orca_constraint_rows(content)[0]


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


def _qmmm_issues(
    kind: str,
    route: str,
    content: str,
    chemistry: dict[str, Any],
    generated: dict[str, Any],
    evidence: dict[str, Any],
    *,
    project: str | None,
    cwd: str | None,
) -> list[InvariantIssue]:
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
    if kind == "gaussian.qmmm":
        issues.extend(
            _gaussian_qmmm_partition_issues(chemistry, generated, evidence)
        )
        issues.extend(_gaussian_qmmm_frequency_policy_issues(route, project, cwd, evidence))
    return issues


def _gaussian_qmmm_frequency_policy_issues(
    route: str,
    project: str | None,
    cwd: str | None,
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    """Honor an explicit ONIOM frequency policy in the selected workspace YAML."""

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
    has_frequency = re.search(r"\bfreq(?:\b|=)", route, flags=re.IGNORECASE) is not None
    if has_frequency == policy:
        return []
    expected = "include" if policy else "omit"
    return [
        _reject(
            "input.gaussian.qmmm.frequency_policy",
            f"generated Gaussian ONIOM route must {expected} frequency analysis per qmmm.freq",
            {**evidence, "project": project, "project_path": str(path), "qmmm_freq": policy},
        )
    ]


def _gaussian_qmmm_partition_issues(
    chemistry: dict[str, Any],
    generated: dict[str, Any],
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    issues: list[InvariantIssue] = []
    atom_layers = generated.get("atom_layers")
    layer_atoms = generated.get("layer_atoms")
    if not isinstance(atom_layers, list) or not atom_layers:
        return [
            _reject(
                "input.gaussian.qmmm.layer_markers",
                "generated Gaussian ONIOM coordinates lack atom-layer markers",
                evidence,
            )
        ]
    if any(layer not in {"H", "M", "L"} for layer in atom_layers):
        issues.append(
            _reject(
                "input.gaussian.qmmm.partition_coverage",
                "every generated Gaussian ONIOM atom must have one H/M/L layer marker",
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

    expected_high = set(_expand_range(str(chemistry.get("high_level_atoms") or "")))
    expected_medium = set(
        _expand_range(str(chemistry.get("medium_level_atoms") or ""))
    )
    explicit_low = chemistry.get("low_level_atoms")
    expected_low = set(_expand_range(str(explicit_low or "")))
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
                _reject(
                    f"input.gaussian.qmmm.{label}_level_atoms",
                    f"generated Gaussian ONIOM {label}-level partition does not preserve the command",
                    {
                        **evidence,
                        "expected": sorted(expected),
                        "observed": sorted(observed),
                    },
                )
            )

    pairs = generated.get("charge_multiplicity_pairs")
    if isinstance(pairs, list) and pairs:
        issues.extend(_gaussian_qmmm_state_pair_issues(chemistry, pairs, evidence))
    return issues


def _gaussian_qmmm_state_pair_issues(
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
                _reject(
                    f"input.gaussian.qmmm.{label}_state",
                    f"generated Gaussian ONIOM {label} charge/multiplicity does not preserve the command",
                    {**evidence, "expected": expected, "observed": observed},
                )
            )

    assert_pair("total", 0, "charge_total", "mult_total")
    has_medium = bool(_expand_range(str(chemistry.get("medium_level_atoms") or "")))
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


__all__ = [
    "check_generated_input_invariants",
    "electron_multiplicity_evidence",
]
