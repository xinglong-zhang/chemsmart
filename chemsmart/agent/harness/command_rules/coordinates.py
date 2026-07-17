"""Scan and ModRedundant command contracts shared by Gaussian and ORCA."""

from __future__ import annotations

from pathlib import Path

from chemsmart.agent.harness.command_rules.models import (
    CommandContractIssue,
    reject,
)
from chemsmart.agent.harness.command_rules.structures import (
    atom_bound_issues,
    parse_coordinate_literal,
    parse_numeric_sequence,
)
from chemsmart.agent.harness.command_rules.tokens import (
    has_option,
    option_value,
)


COORDINATE_FLAGS = ("-c", "--coordinates")
SCAN_FLAGS = {
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


def coordinate_contract_issues(
    program: str,
    job: str,
    tokens: list[str],
    *,
    program_tokens: list[str],
    cwd: str | Path | None,
) -> list[CommandContractIssue]:
    issues: list[CommandContractIssue] = []
    coordinate_value = option_value(tokens, COORDINATE_FLAGS)

    missing = _missing_coordinate_parameters(program, job, tokens)
    if missing:
        issues.append(
            reject(
                "cmd.contract.scan_required_parameters"
                if job == "scan"
                else "cmd.contract.modred_coordinates_required",
                f"{program} {job} is missing required coordinate parameters",
                {"program": program, "job": job, "missing": missing},
                tuple(missing),
            )
        )
        return issues

    coordinates, coordinate_error = parse_coordinate_literal(coordinate_value)
    if coordinate_error is not None:
        issues.append(
            reject(
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
            atom_bound_issues(
                coordinates,
                program_tokens,
                cwd=cwd,
                rule_id="cmd.contract.coordinate_atom_bounds",
                label=f"{program} {job} coordinate",
            )
        )

    if job != "scan" or coordinates is None:
        return issues

    parameter_values, parameter_issues = _parse_scan_parameters(
        program, job, tokens
    )
    issues.extend(parameter_issues)
    coordinate_count = len(coordinates)
    mismatched = {
        name: len(values)
        for name, values in parameter_values.items()
        if len(values) not in {1, coordinate_count}
    }
    if mismatched:
        issues.append(
            reject(
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


def _missing_coordinate_parameters(
    program: str, job: str, tokens: list[str]
) -> list[str]:
    required: dict[str, tuple[str, ...]] = {"coordinates": COORDINATE_FLAGS}
    if job == "scan":
        required.update(SCAN_FLAGS.get(program, {}))
    return [
        f"{name} ({'/'.join(aliases)})"
        for name, aliases in required.items()
        if not has_option(tokens, aliases)
    ]


def _parse_scan_parameters(
    program: str, job: str, tokens: list[str]
) -> tuple[dict[str, list[int | float]], list[CommandContractIssue]]:
    parameter_values: dict[str, list[int | float]] = {}
    issues: list[CommandContractIssue] = []
    for name, aliases in SCAN_FLAGS[program].items():
        value = option_value(tokens, aliases)
        values, value_error = parse_numeric_sequence(
            value,
            integers=name == "num_steps",
        )
        if value_error is not None:
            issues.append(
                reject(
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
        elif values is not None:
            parameter_values[name] = values
    return parameter_values, issues


__all__ = ["coordinate_contract_issues"]
