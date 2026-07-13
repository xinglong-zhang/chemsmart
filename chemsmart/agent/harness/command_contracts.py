"""Deterministic contracts for synthesized Gaussian and ORCA commands.

These checks cover values that Click accepts as strings but the downstream
chemsmart settings/job builders interpret more narrowly.  Keeping them in a
small, dependency-light module gives command repair stable rule IDs before a
safe fake execution fails with a generic traceback.
"""

from __future__ import annotations

import ast
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Literal

import yaml

from chemsmart.settings.workspace_project import workspace_project_path
from chemsmart.utils.utils import get_list_from_string_range

ContractSeverity = Literal["warn", "reject"]

_MISSING = object()
_COORDINATE_FLAGS = ("-c", "--coordinates")
_SCAN_FLAGS = {
    "gaussian": {
        "step_size": ("-s", "--step-size"),
        "num_steps": ("-n", "--num-steps"),
    },
    "orca": {
        "dist_start": ("-x", "--dist-start"),
        "dist_end": ("-y", "--dist-end"),
        "num_steps": ("-n", "--num-steps"),
    },
}


@dataclass(frozen=True)
class CommandContractIssue:
    rule_id: str
    severity: ContractSeverity
    message: str
    evidence: dict[str, Any] = field(default_factory=dict)
    missing_info: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def check_command_contracts(
    *,
    program: str,
    job: str,
    program_tokens: list[str],
    job_tokens: list[str],
    cwd: str | Path | None = None,
) -> tuple[CommandContractIssue, ...]:
    """Return semantic contract violations for one computational command."""

    normalized_program = str(program).strip().lower()
    normalized_job = str(job).strip().lower()
    issues: list[CommandContractIssue] = []

    if normalized_job == "qmmm":
        issues.append(
            _issue(
                "cmd.contract.qmmm_parent_job",
                (
                    "qmmm is a nested chemsmart subcommand and requires a "
                    "parent calculation such as opt, ts, sp, scan, or modred"
                ),
                {"program": normalized_program, "job": normalized_job},
                ("parent computational job before qmmm",),
            )
        )
        return tuple(issues)

    if normalized_program == "gaussian" and normalized_job == "traj":
        traj_tokens = job_tokens[1:]
        if not _has_option(traj_tokens, ("-j", "--jobtype")):
            issues.append(
                _issue(
                    "cmd.contract.traj_jobtype_required",
                    (
                        "Gaussian traj requires -j/--jobtype so the runtime "
                        "can select project settings for each trajectory frame"
                    ),
                    {"program": normalized_program, "job": normalized_job},
                    ("Gaussian trajectory job type (-j/--jobtype)",),
                )
            )
        has_grouping = _has_option(traj_tokens, ("-g", "--grouping-strategy"))
        has_count = _has_option(
            traj_tokens,
            ("-ns", "--num-structures-to-run"),
        )
        if has_grouping and has_count:
            issues.append(
                _issue(
                    "cmd.contract.traj_selection_conflict",
                    (
                        "Gaussian traj cannot combine conformer grouping with "
                        "an explicit number of structures to run"
                    ),
                    {"program": normalized_program, "job": normalized_job},
                    (
                        "choose either -g/--grouping-strategy or "
                        "-ns/--num-structures-to-run",
                    ),
                )
            )
        count = _option_value(traj_tokens, ("-ns", "--num-structures-to-run"))
        if count is not _MISSING and not _is_positive_integer_literal(count):
            issues.append(
                _issue(
                    "cmd.contract.traj_structure_count",
                    "Gaussian traj structure count must be a positive integer",
                    {"program": normalized_program, "job": normalized_job, "value": count},
                    ("positive trajectory structure count",),
                )
            )
        proportion = _option_value(
            traj_tokens,
            ("-x", "--proportion-structures-to-use"),
        )
        if proportion is not _MISSING and not _is_unit_interval_literal(proportion):
            issues.append(
                _issue(
                    "cmd.contract.traj_proportion_range",
                    "Gaussian traj proportion must satisfy 0 < value <= 1",
                    {
                        "program": normalized_program,
                        "job": normalized_job,
                        "value": proportion,
                    },
                    ("trajectory proportion in the interval (0, 1]",),
                )
            )

    if normalized_program == "gaussian" and normalized_job == "crest":
        crest_tokens = job_tokens[1:]
        has_grouping = _has_option(
            crest_tokens,
            ("-g", "--grouping-strategy"),
        )
        has_count = _has_option(
            crest_tokens,
            ("-nc", "--num-confs-to-run"),
        )
        if has_grouping and has_count:
            issues.append(
                _issue(
                    "cmd.contract.crest_selection_conflict",
                    (
                        "Gaussian crest cannot combine conformer grouping "
                        "with an explicit number of conformers to run"
                    ),
                    {"program": normalized_program, "job": normalized_job},
                    (
                        "choose either -g/--grouping-strategy or "
                        "-nc/--num-confs-to-run",
                    ),
                )
            )
        count = _option_value(crest_tokens, ("-nc", "--num-confs-to-run"))
        if count is not _MISSING and not _is_positive_integer_literal(count):
            issues.append(
                _issue(
                    "cmd.contract.crest_conformer_count",
                    "Gaussian crest conformer count must be a positive integer",
                    {
                        "program": normalized_program,
                        "job": normalized_job,
                        "value": count,
                    },
                    ("positive conformer count",),
                )
            )

    if normalized_job in {"scan", "modred"}:
        issues.extend(
            _coordinate_contract_issues(
                normalized_program,
                normalized_job,
                job_tokens[1:],
                program_tokens=program_tokens,
                cwd=cwd,
            )
        )

    if normalized_program == "gaussian" and normalized_job == "dias":
        issues.extend(
            _dias_contract_issues(
                job_tokens[1:],
                program_tokens=program_tokens,
                cwd=cwd,
            )
        )

    qmmm_index = _token_index(job_tokens[1:], "qmmm")
    if qmmm_index is not None:
        # The index is relative to job_tokens[1:].
        qmmm_tokens = job_tokens[qmmm_index + 2 :]
        issues.extend(
            _qmmm_contract_issues(
                normalized_program,
                normalized_job,
                program_tokens,
                qmmm_tokens,
                cwd=cwd,
            )
        )

    if normalized_program == "gaussian" and normalized_job == "td":
        issue = _gaussian_td_project_issue(program_tokens, cwd=cwd)
        if issue is not None:
            issues.append(issue)

    return tuple(issues)


def _coordinate_contract_issues(
    program: str,
    job: str,
    tokens: list[str],
    *,
    program_tokens: list[str],
    cwd: str | Path | None,
) -> list[CommandContractIssue]:
    issues: list[CommandContractIssue] = []
    coordinate_value = _option_value(tokens, _COORDINATE_FLAGS)

    required: dict[str, tuple[str, ...]] = {"coordinates": _COORDINATE_FLAGS}
    if job == "scan":
        required.update(_SCAN_FLAGS.get(program, {}))
    missing = [
        f"{name} ({'/'.join(aliases)})"
        for name, aliases in required.items()
        if not _has_option(tokens, aliases)
    ]
    if missing:
        issues.append(
            _issue(
                "cmd.contract.scan_required_parameters"
                if job == "scan"
                else "cmd.contract.modred_coordinates_required",
                f"{program} {job} is missing required coordinate parameters",
                {"program": program, "job": job, "missing": missing},
                tuple(missing),
            )
        )
        return issues

    coordinates, coordinate_error = _parse_coordinate_literal(coordinate_value)
    if coordinate_error is not None:
        issues.append(
            _issue(
                "cmd.contract.coordinate_literal",
                coordinate_error,
                {
                    "program": program,
                    "job": job,
                    "coordinates": coordinate_value,
                },
                ("coordinates as [i,j] or [[i,j],[k,l,m]]",),
            )
        )
        return issues

    if coordinates is not None:
        issues.extend(
            _atom_bound_issues(
                coordinates,
                program_tokens,
                cwd=cwd,
                rule_id="cmd.contract.coordinate_atom_bounds",
                label=f"{program} {job} coordinate",
            )
        )

    if job != "scan" or coordinates is None:
        return issues

    parameter_values: dict[str, list[int | float]] = {}
    for name, aliases in _SCAN_FLAGS[program].items():
        value = _option_value(tokens, aliases)
        values, value_error = _parse_numeric_sequence(
            value,
            integers=name == "num_steps",
        )
        if value_error is not None:
            issues.append(
                _issue(
                    "cmd.contract.scan_parameter_literal",
                    f"{name} {value_error}",
                    {
                        "program": program,
                        "job": job,
                        "parameter": name,
                        "value": value,
                    },
                    (f"valid {name} scalar or list",),
                )
            )
            continue
        if values is not None:
            parameter_values[name] = values

    coordinate_count = len(coordinates)
    mismatched = {
        name: len(values)
        for name, values in parameter_values.items()
        if len(values) not in {1, coordinate_count}
    }
    if mismatched:
        issues.append(
            _issue(
                "cmd.contract.scan_parameter_cardinality",
                (
                    "scan parameter lists must contain one value or one value "
                    "per scanned coordinate"
                ),
                {
                    "program": program,
                    "job": job,
                    "coordinate_count": coordinate_count,
                    "parameter_counts": mismatched,
                },
            )
        )
    return issues


def _dias_contract_issues(
    tokens: list[str],
    *,
    program_tokens: list[str],
    cwd: str | Path | None,
) -> list[CommandContractIssue]:
    aliases = ("-i", "--fragment-indices")
    if not _has_option(tokens, aliases):
        return [
            _issue(
                "cmd.contract.dias_fragment_indices_required",
                "Gaussian DIAS requires fragment-1 atom indices",
                {"program": "gaussian", "job": "dias"},
                ("fragment-1 indices (-i/--fragment-indices)",),
            )
        ]

    value = _option_value(tokens, aliases)
    message = _range_literal_error(value)
    if message is not None:
        return [
            _issue(
                "cmd.contract.dias_fragment_indices",
                message,
                {"program": "gaussian", "job": "dias", "value": value},
                ("flat fragment-1 comma/range indices such as 1-6,9",),
            )
        ]
    issues = _atom_bound_issues(
        [get_list_from_string_range(str(value))],
        program_tokens,
        cwd=cwd,
        rule_id="cmd.contract.dias_fragment_atom_bounds",
        label="Gaussian DIAS fragment-1",
    )
    atom_count = _input_atom_count(program_tokens, cwd=cwd)
    if atom_count is not None and len(get_list_from_string_range(str(value))) == atom_count:
        issues.append(
            _issue(
                "cmd.contract.dias_fragment_complement_empty",
                "Gaussian DIAS fragment-1 cannot contain every atom",
                {
                    "atom_count": atom_count,
                    "fragment_indices": value,
                },
                ("a fragment-1 range that leaves at least one atom for fragment 2",),
            )
        )
    return issues


def _qmmm_contract_issues(
    program: str,
    parent_job: str,
    program_tokens: list[str],
    qmmm_tokens: list[str],
    *,
    cwd: str | Path | None,
) -> list[CommandContractIssue]:
    issues: list[CommandContractIssue] = []
    high_atoms_aliases = ("-ha", "--high-level-atoms")
    if not _has_option(qmmm_tokens, high_atoms_aliases):
        issues.append(
            _issue(
                "cmd.contract.qmmm_high_level_atoms_required",
                "QM/MM requires an explicit high-level atom region",
                {"program": program, "parent_job": parent_job},
                ("high-level atom indices (-ha/--high-level-atoms)",),
            )
        )
    else:
        value = _option_value(qmmm_tokens, high_atoms_aliases)
        message = _range_literal_error(value)
        if message is not None:
            issues.append(
                _issue(
                    "cmd.contract.qmmm_high_level_atoms",
                    message,
                    {
                        "program": program,
                        "parent_job": parent_job,
                        "value": value,
                    },
                    ("high-level atoms as comma/range indices",),
                )
            )

    region_values: dict[str, str] = {}
    for name, aliases in (
        ("high", high_atoms_aliases),
        (
            "medium",
            (
                "-ma",
                "--medium-level-atoms",
                "-ia",
                "--intermediate-level-atoms",
            ),
        ),
        ("low", ("-la", "--low-level-atoms")),
    ):
        value = _option_value(qmmm_tokens, aliases)
        if isinstance(value, str) and _range_literal_error(value) is None:
            region_values[name] = value

    for name, value in region_values.items():
        issues.extend(
            _atom_bound_issues(
                [get_list_from_string_range(value)],
                program_tokens,
                cwd=cwd,
                rule_id="cmd.contract.qmmm_atom_bounds",
                label=f"QM/MM {name}-level region",
            )
        )

    expanded_regions = {
        name: set(get_list_from_string_range(value))
        for name, value in region_values.items()
    }
    for left, right in (
        ("high", "medium"),
        ("high", "low"),
        ("medium", "low"),
    ):
        overlap = expanded_regions.get(left, set()) & expanded_regions.get(
            right,
            set(),
        )
        if overlap:
            issues.append(
                _issue(
                    "cmd.contract.qmmm_region_overlap",
                    "QM/MM atom regions must not overlap",
                    {
                        "program": program,
                        "parent_job": parent_job,
                        "regions": [left, right],
                        "overlap": sorted(overlap),
                    },
                    ("non-overlapping QM/MM atom partitions",),
                )
            )

    has_program_charge = _has_option(program_tokens, ("-c", "--charge"))
    has_program_mult = _has_option(program_tokens, ("-m", "--multiplicity"))
    has_total_charge = _has_option(qmmm_tokens, ("-ct", "--charge-total"))
    has_total_mult = _has_option(qmmm_tokens, ("-mt", "--mult-total"))
    missing: list[str] = []
    # Gaussian's qmmm builder does not merge program-level charge/multiplicity
    # into GaussianQMMMJobSettings. ORCA does, so it can use either scope.
    charge_is_explicit = has_total_charge or (
        program == "orca" and has_program_charge
    )
    mult_is_explicit = has_total_mult or (
        program == "orca" and has_program_mult
    )
    if not charge_is_explicit:
        missing.append(
            "total charge (-ct/--charge-total after qmmm)"
            if program == "gaussian"
            else "total charge (-c before parent job or -ct after qmmm)"
        )
    if not mult_is_explicit:
        missing.append(
            "total multiplicity (-mt/--mult-total after qmmm)"
            if program == "gaussian"
            else "total multiplicity (-m before parent job or -mt after qmmm)"
        )
    if missing:
        issues.append(
            _issue(
                "cmd.contract.qmmm_total_state_required",
                "QM/MM requires an unambiguous total charge and multiplicity",
                {
                    "program": program,
                    "parent_job": parent_job,
                    "missing": missing,
                },
                tuple(missing),
            )
        )

    qmmm_jobtype = _option_value(qmmm_tokens, ("-j", "--jobtype"))
    needs_mm_method = not isinstance(qmmm_jobtype, str) or qmmm_jobtype.upper() in {
        "QMMM",
        "QM/MM",
        "QM/QM2/MM",
        "IONIC-CRYSTAL-QMMM",
    }
    has_command_mm_method = _has_option(
        qmmm_tokens,
        ("-lm", "--low-level-method", "--low-level-force-field"),
    )
    has_project_mm_method = _project_has_qmmm_low_level_method(
        program,
        program_tokens,
        cwd=cwd,
    )
    if (
        program == "orca"
        and needs_mm_method
        and not has_command_mm_method
        and not has_project_mm_method
    ):
        issues.append(
            _issue(
                "cmd.contract.qmmm_low_level_method_required",
                (
                    "ORCA QMMM jobs with an MM layer require a low-level "
                    "method/force-field in the command or workspace project YAML"
                ),
                {
                    "program": program,
                    "parent_job": parent_job,
                    "qmmm_jobtype": qmmm_jobtype
                    if isinstance(qmmm_jobtype, str)
                    else "QMMM (default)",
                },
                (
                    "ORCA QM/MM low-level method (-lm/--low-level-method or "
                    "qmmm.low_level_method in project YAML)",
                ),
            )
        )
    return issues


def _gaussian_td_project_issue(
    program_tokens: list[str],
    *,
    cwd: str | Path | None,
) -> CommandContractIssue | None:
    project = _option_value(program_tokens, ("-p", "--project"))
    if not isinstance(project, str) or not project.strip():
        return None
    path = workspace_project_path(project, "gaussian", cwd=cwd)
    if not path.is_file():
        return None
    try:
        payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        return _issue(
            "cmd.contract.project_yaml_invalid",
            f"workspace project YAML could not be loaded: {exc}",
            {"program": "gaussian", "project": project, "path": str(path)},
            ("valid workspace project YAML",),
        )
    if isinstance(payload, dict) and isinstance(payload.get("td"), dict):
        return None
    return _issue(
        "cmd.contract.gaussian_td_project_settings",
        (
            "Gaussian td requires a top-level td settings block in the "
            "selected workspace project YAML"
        ),
        {"program": "gaussian", "project": project, "path": str(path)},
        ("top-level td settings in the Gaussian project YAML",),
    )


def _project_has_qmmm_low_level_method(
    program: str,
    program_tokens: list[str],
    *,
    cwd: str | Path | None,
) -> bool:
    project = _option_value(program_tokens, ("-p", "--project"))
    if not isinstance(project, str) or not project.strip():
        return False
    path = workspace_project_path(project, program, cwd=cwd)
    if not path.is_file():
        return False
    try:
        payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError):
        return False
    if not isinstance(payload, dict) or not isinstance(payload.get("qmmm"), dict):
        return False
    qmmm = payload["qmmm"]
    return bool(
        qmmm.get("low_level_method")
        or qmmm.get("low_level_force_field")
        or qmmm.get("mm_force_field")
    )


def _parse_coordinate_literal(
    value: str | object,
) -> tuple[list[list[int]] | None, str | None]:
    if not isinstance(value, str):
        return None, "coordinates must have a value"
    try:
        parsed = ast.literal_eval(value)
    except (SyntaxError, ValueError):
        return None, "coordinates must be a Python/JSON list literal"
    if not isinstance(parsed, list) or not parsed:
        return None, "coordinates must be a non-empty list"

    if all(isinstance(item, int) and not isinstance(item, bool) for item in parsed):
        groups: list[list[Any]] = [parsed]
    elif all(isinstance(item, list) for item in parsed):
        groups = parsed
    else:
        return None, "coordinates cannot mix atom indices and nested lists"

    normalized: list[list[int]] = []
    for group in groups:
        if not 2 <= len(group) <= 4:
            return None, "each coordinate must contain 2, 3, or 4 atom indices"
        if not all(
            isinstance(item, int) and not isinstance(item, bool) and item > 0
            for item in group
        ):
            return None, "coordinate atom indices must be positive integers"
        normalized.append(list(group))
    return normalized, None


def _parse_numeric_sequence(
    value: str | object,
    *,
    integers: bool,
) -> tuple[list[int | float] | None, str | None]:
    if not isinstance(value, str):
        return None, "must have a value"
    try:
        parsed = ast.literal_eval(value)
    except (SyntaxError, ValueError):
        return None, "must be a numeric scalar or list literal"
    values = parsed if isinstance(parsed, list) else [parsed]
    if not values:
        return None, "cannot be an empty list"
    if any(isinstance(item, bool) for item in values):
        return None, "must contain numeric values"
    if integers:
        if not all(isinstance(item, int) and item > 0 for item in values):
            return None, "must contain positive integers"
        return list(values), None
    if not all(isinstance(item, (int, float)) for item in values):
        return None, "must contain numeric values"
    return list(values), None


def _range_literal_error(value: str | object) -> str | None:
    if not isinstance(value, str) or not value.strip():
        return "atom indices must have a value"
    if "[" in value or "]" in value or any(char.isspace() for char in value):
        return "atom indices must use a flat comma/range form without brackets or spaces"
    try:
        indices = get_list_from_string_range(value)
    except (TypeError, ValueError):
        return "atom indices must use a valid comma/range form such as 1-6,9"
    if not indices or any(index <= 0 for index in indices):
        return "atom indices must be positive and 1-indexed"
    if len(indices) != len(set(indices)):
        return "atom indices must not contain duplicates or overlapping ranges"
    return None


def _atom_bound_issues(
    groups: list[list[int]],
    program_tokens: list[str],
    *,
    cwd: str | Path | None,
    rule_id: str,
    label: str,
) -> list[CommandContractIssue]:
    atom_count = _input_atom_count(program_tokens, cwd=cwd)
    if atom_count is None:
        return []
    invalid = sorted(
        {
            index
            for group in groups
            for index in group
            if index < 1 or index > atom_count
        }
    )
    if not invalid:
        return []
    return [
        _issue(
            rule_id,
            f"{label} contains atom indices outside the input structure",
            {
                "atom_count": atom_count,
                "invalid_indices": invalid,
                "input": _option_value(
                    program_tokens,
                    ("-f", "--filename"),
                ),
            },
            (f"atom indices in the inclusive range 1-{atom_count}",),
        )
    ]


def _input_atom_count(
    program_tokens: list[str],
    *,
    cwd: str | Path | None,
) -> int | None:
    source = _option_value(program_tokens, ("-f", "--filename"))
    if not isinstance(source, str):
        return None
    path = Path(source).expanduser()
    if not path.is_absolute():
        path = Path(cwd or Path.cwd()) / path
    if path.suffix.lower() != ".xyz" or not path.is_file():
        return None
    try:
        first_line = path.read_text(encoding="utf-8").splitlines()[0].strip()
        count = int(first_line)
    except (OSError, UnicodeDecodeError, IndexError, ValueError):
        return None
    return count if count > 0 else None


def _is_positive_integer_literal(value: str | object) -> bool:
    try:
        parsed = ast.literal_eval(value) if isinstance(value, str) else value
    except (SyntaxError, ValueError):
        return False
    return isinstance(parsed, int) and not isinstance(parsed, bool) and parsed > 0


def _is_unit_interval_literal(value: str | object) -> bool:
    try:
        parsed = ast.literal_eval(value) if isinstance(value, str) else value
    except (SyntaxError, ValueError):
        return False
    return (
        isinstance(parsed, (int, float))
        and not isinstance(parsed, bool)
        and 0 < float(parsed) <= 1
    )


def _token_index(tokens: list[str], target: str) -> int | None:
    try:
        return tokens.index(target)
    except ValueError:
        return None


def _has_option(tokens: list[str], aliases: tuple[str, ...]) -> bool:
    return _option_value(tokens, aliases) is not _MISSING


def _option_value(tokens: list[str], aliases: tuple[str, ...]) -> str | object:
    for index, token in enumerate(tokens):
        for alias in aliases:
            if token == alias:
                if index + 1 < len(tokens):
                    return tokens[index + 1]
                return None
            if alias.startswith("--") and token.startswith(f"{alias}="):
                return token.split("=", 1)[1]
    return _MISSING


def _issue(
    rule_id: str,
    message: str,
    evidence: dict[str, Any],
    missing_info: tuple[str, ...] = (),
) -> CommandContractIssue:
    return CommandContractIssue(
        rule_id=rule_id,
        severity="reject",
        message=message,
        evidence=evidence,
        missing_info=missing_info,
    )


__all__ = ["CommandContractIssue", "check_command_contracts"]
