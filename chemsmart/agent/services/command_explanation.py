"""User-facing rendering for deterministically parsed ChemSmart commands."""

from __future__ import annotations

import os
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from chemsmart.agent.model_command_parser import ParsedModelCommand

_JOB_LABELS = {
    "com": "Gaussian input regeneration",
    "crest": "CREST conformer search",
    "dias": "distortion-interaction analysis",
    "inp": "ORCA input handling",
    "irc": "intrinsic reaction coordinate",
    "link": "Gaussian linked job",
    "modred": "constrained optimization",
    "neb": "nudged elastic band",
    "nci": "non-covalent interaction analysis",
    "opt": "geometry optimization",
    "qmmm": "QM/MM calculation",
    "qrc": "quasi-reaction coordinate scan",
    "resp": "RESP charge calculation",
    "scan": "coordinate scan",
    "sp": "single-point energy",
    "td": "TD-DFT excited-state calculation",
    "traj": "trajectory frame workflow",
    "ts": "transition-state search",
    "userjob": "user-defined Gaussian job",
    "wbi": "Wiberg bond index calculation",
}


def format_model_command_explanation(
    command: str,
    *,
    cwd: str | os.PathLike[str] | None = None,
) -> str:
    from chemsmart.agent.model_command_parser import parse_model_command

    return format_parsed_model_command(parse_model_command(command, cwd=cwd))


def format_parsed_model_command(parsed: ParsedModelCommand) -> str:
    if parsed.parse_error:
        return "\n".join(
            [
                "deterministic command parser:",
                f"- parse status: `error` ({parsed.parse_error})",
                f"- workspace: `{parsed.workspace}`",
                "",
                f"Summary: This command could not be deterministically parsed: {parsed.parse_error}.",
            ]
        )
    action_label = (
        "submit to an HPC/server queue"
        if parsed.action == "sub"
        else "run locally"
    )
    dry_run_text = "yes" if parsed.dry_run else "no"
    program_text = parsed.program or "unknown"
    job_text = parsed.job or "unknown"
    job_label = _JOB_LABELS.get(job_text, job_text)
    project_text = parsed.project or "not specified"
    server_text = parsed.server or (
        "auto/default server" if parsed.action == "sub" else "local/default"
    )
    method_bits = _method_bits(parsed)
    method_text = ", ".join(method_bits) if method_bits else "unresolved"
    lines = _base_lines(
        parsed,
        action_label=action_label,
        dry_run_text=dry_run_text,
        program_text=program_text,
        job_text=job_text,
        job_label=job_label,
        project_text=project_text,
        server_text=server_text,
    )
    _append_optional_lines(lines, parsed, method_text)
    target = parsed.filename or "the selected input"
    if parsed.action == "sub":
        summary_action = (
            f"submit a {program_text} {job_label} job for `{target}` "
            f"to `{server_text}`"
        )
    else:
        summary_action = (
            f"run a {program_text} {job_label} job for `{target}` locally"
        )
    lines.append("")
    lines.append(
        "Summary: This command will "
        f"{summary_action} from `{parsed.workspace}` using project "
        f"`{project_text}` with "
        f"{method_text}; dry-run is {dry_run_text}."
    )
    return "\n".join(lines)


def _method_bits(parsed: ParsedModelCommand) -> list[str]:
    values = (
        ("ab initio", parsed.ab_initio),
        ("functional", parsed.functional),
        ("basis", parsed.basis),
        ("auxiliary basis", parsed.aux_basis),
        ("extrapolation basis", parsed.extrapolation_basis),
    )
    return [f"{label} `{value}`" for label, value in values if value]


def _base_lines(
    parsed: ParsedModelCommand,
    *,
    action_label: str,
    dry_run_text: str,
    program_text: str,
    job_text: str,
    job_label: str,
    project_text: str,
    server_text: str,
) -> list[str]:
    return [
        "deterministic command parser:",
        f"- workspace: `{parsed.workspace}`",
        f"- execution: `{parsed.action}` ({action_label})",
        f"- program: `{program_text}`",
        f"- job: `{job_text}` ({job_label})",
        f"- server: `{server_text}`",
        f"- dry run requested by command: `{dry_run_text}`",
        f"- molecule/input file: `{parsed.filename or 'not specified'}`",
        f"- label: `{parsed.label or 'runtime-derived'}`",
        f"- charge/multiplicity: `{parsed.charge or 'runtime/default'}` / `{parsed.multiplicity or 'runtime/default'}`",
        f"- project: `{project_text}`",
    ]


def _append_optional_lines(
    lines: list[str],
    parsed: ParsedModelCommand,
    method_text: str,
) -> None:
    _append_context_lines(lines, parsed, method_text)
    _append_job_lines(lines, parsed)


def _append_context_lines(
    lines: list[str],
    parsed: ParsedModelCommand,
    method_text: str,
) -> None:
    if parsed.project_p_flag_meaning:
        lines.append(f"- `-p` meaning: {parsed.project_p_flag_meaning}")
    if parsed.top_level_program:
        lines.append(
            "- top-level `-p/--program`: "
            f"`{parsed.top_level_program}` (output-file processing target, not project)"
        )
    db_bits = _database_bits(parsed)
    if db_bits:
        lines.append(f"- database selection: `{', '.join(db_bits)}`")
    lines.append(f"- resolved method: {method_text}")
    if parsed.solvent_model or parsed.solvent_id:
        lines.append(
            "- resolved solvent: "
            f"`{parsed.solvent_model or 'model default'}` / `{parsed.solvent_id or 'id default'}`"
        )
    controls = _route_control_bits(parsed)
    if controls:
        lines.append(f"- resolved route controls: `{', '.join(controls)}`")


def _append_job_lines(
    lines: list[str],
    parsed: ParsedModelCommand,
) -> None:
    if parsed.route_parameters:
        lines.append(f"- route parameters: `{parsed.route_parameters}`")
    if parsed.opt_options:
        lines.append(f"- optimization route options: `{parsed.opt_options}`")
    if parsed.resources:
        text = ", ".join(
            f"{key}={value}" for key, value in sorted(parsed.resources.items())
        )
        lines.append(f"- resources: `{text}`")
    if parsed.structural_options:
        text = ", ".join(
            f"{key}={value}"
            for key, value in sorted(parsed.structural_options.items())
        )
        lines.append(f"- job-specific options: `{text}`")
    if parsed.warnings:
        lines.append("- parser warnings:")
        lines.extend(f"  - {warning}" for warning in parsed.warnings)


def _database_bits(parsed: ParsedModelCommand) -> list[str]:
    values = (
        ("record_index", parsed.record_index),
        ("record_id", parsed.record_id),
        ("structure_index", parsed.structure_index),
        ("structure_id", parsed.structure_id),
        ("molecule_id", parsed.molecule_id),
    )
    return [f"{name}={value}" for name, value in values if value]


def _route_control_bits(parsed: ParsedModelCommand) -> list[str]:
    values = (
        ("defgrid", parsed.defgrid),
        ("scf_tol", parsed.scf_tol),
        ("scf_algorithm", parsed.scf_algorithm),
    )
    return [f"{name}={value}" for name, value in values if value]


__all__ = ["format_model_command_explanation", "format_parsed_model_command"]
