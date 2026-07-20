"""QM/MM command contracts shared by Gaussian and ORCA."""

from __future__ import annotations

from pathlib import Path

import yaml

from chemsmart.agent.harness.command_rules.models import (
    CommandContractIssue,
    reject,
)
from chemsmart.agent.harness.command_rules.structures import (
    atom_bound_issues,
    range_literal_error,
)
from chemsmart.agent.harness.command_rules.tokens import (
    has_option,
    option_value,
)
from chemsmart.settings.workspace_project import workspace_project_path
from chemsmart.utils.utils import get_list_from_string_range


def qmmm_contract_issues(
    program: str,
    parent_job: str,
    program_tokens: list[str],
    qmmm_tokens: list[str],
    *,
    cwd: str | Path | None,
) -> list[CommandContractIssue]:
    issues: list[CommandContractIssue] = []
    high_atoms_aliases = ("-ha", "--high-level-atoms")
    if not has_option(qmmm_tokens, high_atoms_aliases):
        issues.append(
            reject(
                "cmd.contract.qmmm_high_level_atoms_required",
                "QM/MM requires an explicit high-level atom region",
                {"program": program, "parent_job": parent_job},
                ("high-level atom indices (-ha/--high-level-atoms)",),
            )
        )
    else:
        value = option_value(qmmm_tokens, high_atoms_aliases)
        message = range_literal_error(value)
        if message is not None:
            issues.append(
                reject(
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

    region_values = _region_values(qmmm_tokens, high_atoms_aliases)
    for name, value in region_values.items():
        issues.extend(
            atom_bound_issues(
                [get_list_from_string_range(value)],
                program_tokens,
                cwd=cwd,
                rule_id="cmd.contract.qmmm_atom_bounds",
                label=f"QM/MM {name}-level region",
            )
        )
    issues.extend(_region_overlap_issues(program, parent_job, region_values))
    issues.extend(
        _total_state_issues(
            program,
            parent_job,
            program_tokens,
            qmmm_tokens,
        )
    )
    mm_issue = _low_level_method_issue(
        program,
        parent_job,
        program_tokens,
        qmmm_tokens,
        cwd=cwd,
    )
    if mm_issue is not None:
        issues.append(mm_issue)
    return issues


def _region_values(
    qmmm_tokens: list[str],
    high_atoms_aliases: tuple[str, ...],
) -> dict[str, str]:
    values: dict[str, str] = {}
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
        value = option_value(qmmm_tokens, aliases)
        if isinstance(value, str) and range_literal_error(value) is None:
            values[name] = value
    return values


def _region_overlap_issues(
    program: str,
    parent_job: str,
    region_values: dict[str, str],
) -> list[CommandContractIssue]:
    expanded = {
        name: set(get_list_from_string_range(value))
        for name, value in region_values.items()
    }
    issues: list[CommandContractIssue] = []
    for left, right in (
        ("high", "medium"),
        ("high", "low"),
        ("medium", "low"),
    ):
        overlap = expanded.get(left, set()) & expanded.get(right, set())
        if overlap:
            issues.append(
                reject(
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
    return issues


def _total_state_issues(
    program: str,
    parent_job: str,
    program_tokens: list[str],
    qmmm_tokens: list[str],
) -> list[CommandContractIssue]:
    has_program_charge = has_option(program_tokens, ("-c", "--charge"))
    has_program_mult = has_option(program_tokens, ("-m", "--multiplicity"))
    has_total_charge = has_option(qmmm_tokens, ("-ct", "--charge-total"))
    has_total_mult = has_option(qmmm_tokens, ("-mt", "--mult-total"))
    charge_is_explicit = has_total_charge or (
        program == "orca" and has_program_charge
    )
    mult_is_explicit = has_total_mult or (
        program == "orca" and has_program_mult
    )
    missing: list[str] = []
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
    if not missing:
        return []
    return [
        reject(
            "cmd.contract.qmmm_total_state_required",
            "QM/MM requires an unambiguous total charge and multiplicity",
            {"program": program, "parent_job": parent_job, "missing": missing},
            tuple(missing),
        )
    ]


def _low_level_method_issue(
    program: str,
    parent_job: str,
    program_tokens: list[str],
    qmmm_tokens: list[str],
    *,
    cwd: str | Path | None,
) -> CommandContractIssue | None:
    qmmm_jobtype = option_value(qmmm_tokens, ("-j", "--jobtype"))
    needs_mm_method = not isinstance(
        qmmm_jobtype, str
    ) or qmmm_jobtype.upper() in {
        "QMMM",
        "QM/MM",
        "QM/QM2/MM",
        "IONIC-CRYSTAL-QMMM",
    }
    has_command_mm_method = has_option(
        qmmm_tokens,
        ("-lm", "--low-level-method", "--low-level-force-field"),
    )
    if (
        program != "orca"
        or not needs_mm_method
        or has_command_mm_method
        or _project_has_low_level_method(program, program_tokens, cwd=cwd)
    ):
        return None
    return reject(
        "cmd.contract.qmmm_low_level_method_required",
        (
            "ORCA QMMM jobs with an MM layer require a low-level "
            "method/force-field in the command or workspace project YAML"
        ),
        {
            "program": program,
            "parent_job": parent_job,
            "qmmm_jobtype": (
                qmmm_jobtype
                if isinstance(qmmm_jobtype, str)
                else "QMMM (default)"
            ),
        },
        (
            "ORCA QM/MM low-level method (-lm/--low-level-method or "
            "qmmm.low_level_method in project YAML)",
        ),
    )


def _project_has_low_level_method(
    program: str,
    program_tokens: list[str],
    *,
    cwd: str | Path | None,
) -> bool:
    project = option_value(program_tokens, ("-p", "--project"))
    if not isinstance(project, str) or not project.strip():
        return False
    path = workspace_project_path(project, program, cwd=cwd)
    if not path.is_file():
        return False
    try:
        payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError):
        return False
    if not isinstance(payload, dict) or not isinstance(
        payload.get("qmmm"), dict
    ):
        return False
    qmmm = payload["qmmm"]
    return bool(
        qmmm.get("low_level_method")
        or qmmm.get("low_level_force_field")
        or qmmm.get("mm_force_field")
    )


__all__ = ["qmmm_contract_issues"]
