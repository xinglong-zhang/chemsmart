"""Public facade for deterministic Gaussian and ORCA command contracts.

The program- and job-specific rules live in ``command_rules``. This module
keeps the established import path and dispatch order stable for callers.
"""

from __future__ import annotations

from pathlib import Path

from chemsmart.agent.harness.command_rules.coordinates import (
    coordinate_contract_issues,
)
from chemsmart.agent.harness.command_rules.gaussian import (
    dias_contract_issues,
    selection_contract_issues,
    td_project_issue,
)
from chemsmart.agent.harness.command_rules.models import (
    CommandContractIssue,
    ContractSeverity,
    reject,
)
from chemsmart.agent.harness.command_rules.qmmm import qmmm_contract_issues
from chemsmart.agent.harness.command_rules.tokens import token_index
from chemsmart.agent.harness.command_rules.xtb import solvent_contract_issues


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
        return (
            reject(
                "cmd.contract.qmmm_parent_job",
                (
                    "qmmm is a nested chemsmart subcommand and requires a "
                    "parent calculation such as opt, ts, sp, scan, or modred"
                ),
                {"program": normalized_program, "job": normalized_job},
                ("parent computational job before qmmm",),
            ),
        )

    if normalized_program == "gaussian":
        issues.extend(selection_contract_issues(normalized_job, job_tokens))

    if normalized_job in {"scan", "modred"}:
        issues.extend(
            coordinate_contract_issues(
                normalized_program,
                normalized_job,
                job_tokens[1:],
                program_tokens=program_tokens,
                cwd=cwd,
            )
        )

    if normalized_program == "gaussian" and normalized_job == "dias":
        issues.extend(
            dias_contract_issues(
                job_tokens[1:],
                program_tokens=program_tokens,
                cwd=cwd,
            )
        )

    qmmm_index = token_index(job_tokens[1:], "qmmm")
    if qmmm_index is not None:
        issues.extend(
            qmmm_contract_issues(
                normalized_program,
                normalized_job,
                program_tokens,
                job_tokens[qmmm_index + 2 :],
                cwd=cwd,
            )
        )

    if normalized_program == "gaussian" and normalized_job == "td":
        issue = td_project_issue(program_tokens, cwd=cwd)
        if issue is not None:
            issues.append(issue)

    if normalized_program == "xtb":
        issues.extend(solvent_contract_issues(program_tokens))

    return tuple(issues)


# Preserve the historical public class identity for repr, pickling, and
# introspection while keeping its implementation dependency-light.
CommandContractIssue.__module__ = __name__


__all__ = [
    "CommandContractIssue",
    "ContractSeverity",
    "check_command_contracts",
]
