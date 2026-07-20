"""Gaussian-only command contracts."""

from __future__ import annotations

from pathlib import Path

import yaml

from chemsmart.agent.harness.command_rules.models import (
    CommandContractIssue,
    reject,
)
from chemsmart.agent.harness.command_rules.structures import (
    atom_bound_issues,
    input_atom_count,
    range_literal_error,
)
from chemsmart.agent.harness.command_rules.tokens import (
    MISSING,
    has_option,
    is_positive_integer_literal,
    is_unit_interval_literal,
    option_value,
)
from chemsmart.settings.workspace_project import workspace_project_path
from chemsmart.utils.utils import get_list_from_string_range


def selection_contract_issues(
    job: str,
    job_tokens: list[str],
) -> list[CommandContractIssue]:
    if job == "traj":
        return _trajectory_contract_issues(job_tokens[1:])
    if job == "crest":
        return _crest_contract_issues(job_tokens[1:])
    return []


def _trajectory_contract_issues(
    tokens: list[str],
) -> list[CommandContractIssue]:
    issues: list[CommandContractIssue] = []
    if not has_option(tokens, ("-j", "--jobtype")):
        issues.append(
            reject(
                "cmd.contract.traj_jobtype_required",
                (
                    "Gaussian traj requires -j/--jobtype so the runtime "
                    "can select project settings for each trajectory frame"
                ),
                {"program": "gaussian", "job": "traj"},
                ("Gaussian trajectory job type (-j/--jobtype)",),
            )
        )
    has_grouping = has_option(tokens, ("-g", "--grouping-strategy"))
    has_count = has_option(tokens, ("-ns", "--num-structures-to-run"))
    if has_grouping and has_count:
        issues.append(
            reject(
                "cmd.contract.traj_selection_conflict",
                (
                    "Gaussian traj cannot combine conformer grouping with "
                    "an explicit number of structures to run"
                ),
                {"program": "gaussian", "job": "traj"},
                (
                    "choose either -g/--grouping-strategy or "
                    "-ns/--num-structures-to-run",
                ),
            )
        )
    count = option_value(tokens, ("-ns", "--num-structures-to-run"))
    if count is not MISSING and not is_positive_integer_literal(count):
        issues.append(
            reject(
                "cmd.contract.traj_structure_count",
                "Gaussian traj structure count must be a positive integer",
                {"program": "gaussian", "job": "traj", "value": count},
                ("positive trajectory structure count",),
            )
        )
    proportion = option_value(tokens, ("-x", "--proportion-structures-to-use"))
    if proportion is not MISSING and not is_unit_interval_literal(proportion):
        issues.append(
            reject(
                "cmd.contract.traj_proportion_range",
                "Gaussian traj proportion must satisfy 0 < value <= 1",
                {
                    "program": "gaussian",
                    "job": "traj",
                    "value": proportion,
                },
                ("trajectory proportion in the interval (0, 1]",),
            )
        )
    return issues


def _crest_contract_issues(tokens: list[str]) -> list[CommandContractIssue]:
    issues: list[CommandContractIssue] = []
    has_grouping = has_option(tokens, ("-g", "--grouping-strategy"))
    has_count = has_option(tokens, ("-nc", "--num-confs-to-run"))
    if has_grouping and has_count:
        issues.append(
            reject(
                "cmd.contract.crest_selection_conflict",
                (
                    "Gaussian crest cannot combine conformer grouping "
                    "with an explicit number of conformers to run"
                ),
                {"program": "gaussian", "job": "crest"},
                (
                    "choose either -g/--grouping-strategy or "
                    "-nc/--num-confs-to-run",
                ),
            )
        )
    count = option_value(tokens, ("-nc", "--num-confs-to-run"))
    if count is not MISSING and not is_positive_integer_literal(count):
        issues.append(
            reject(
                "cmd.contract.crest_conformer_count",
                "Gaussian crest conformer count must be a positive integer",
                {"program": "gaussian", "job": "crest", "value": count},
                ("positive conformer count",),
            )
        )
    return issues


def dias_contract_issues(
    tokens: list[str],
    *,
    program_tokens: list[str],
    cwd: str | Path | None,
) -> list[CommandContractIssue]:
    aliases = ("-i", "--fragment-indices")
    if not has_option(tokens, aliases):
        return [
            reject(
                "cmd.contract.dias_fragment_indices_required",
                "Gaussian DIAS requires fragment-1 atom indices",
                {"program": "gaussian", "job": "dias"},
                ("fragment-1 indices (-i/--fragment-indices)",),
            )
        ]

    value = option_value(tokens, aliases)
    message = range_literal_error(value)
    if message is not None:
        return [
            reject(
                "cmd.contract.dias_fragment_indices",
                message,
                {"program": "gaussian", "job": "dias", "value": value},
                ("flat fragment-1 comma/range indices such as 1-6,9",),
            )
        ]
    fragment = get_list_from_string_range(str(value))
    issues = atom_bound_issues(
        [fragment],
        program_tokens,
        cwd=cwd,
        rule_id="cmd.contract.dias_fragment_atom_bounds",
        label="Gaussian DIAS fragment-1",
    )
    atom_count = input_atom_count(program_tokens, cwd=cwd)
    if atom_count is not None and len(fragment) == atom_count:
        issues.append(
            reject(
                "cmd.contract.dias_fragment_complement_empty",
                "Gaussian DIAS fragment-1 cannot contain every atom",
                {"atom_count": atom_count, "fragment_indices": value},
                (
                    "a fragment-1 range that leaves at least one atom for "
                    "fragment 2",
                ),
            )
        )
    return issues


def td_project_issue(
    program_tokens: list[str],
    *,
    cwd: str | Path | None,
) -> CommandContractIssue | None:
    project = option_value(program_tokens, ("-p", "--project"))
    if not isinstance(project, str) or not project.strip():
        return None
    path = workspace_project_path(project, "gaussian", cwd=cwd)
    if not path.is_file():
        return None
    try:
        payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        return reject(
            "cmd.contract.project_yaml_invalid",
            f"workspace project YAML could not be loaded: {exc}",
            {"program": "gaussian", "project": project, "path": str(path)},
            ("valid workspace project YAML",),
        )
    if isinstance(payload, dict) and isinstance(payload.get("td"), dict):
        return None
    return reject(
        "cmd.contract.gaussian_td_project_settings",
        (
            "Gaussian td requires a top-level td settings block in the "
            "selected workspace project YAML"
        ),
        {"program": "gaussian", "project": project, "path": str(path)},
        ("top-level td settings in the Gaussian project YAML",),
    )


__all__ = [
    "dias_contract_issues",
    "selection_contract_issues",
    "td_project_issue",
]
