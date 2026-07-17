"""Deterministic explanation for commands produced by the local agent model."""

from __future__ import annotations

import os
import shlex
import contextlib
import io
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from chemsmart.agent.command_models import (
    ParsedModelCommand as ParsedModelCommand,
)
from chemsmart.agent.services.command_explanation import (
    format_parsed_model_command as format_parsed_model_command,
)

_TOP_LEVEL_COMMANDS = {"run", "sub"}
_PROGRAMS = {"gaussian", "orca", "nciplot", "mol"}
_GAUSSIAN_SUBCOMMANDS = {
    "com",
    "crest",
    "dias",
    "irc",
    "link",
    "modred",
    "nci",
    "opt",
    "qmmm",
    "qrc",
    "resp",
    "scan",
    "sp",
    "td",
    "traj",
    "ts",
    "userjob",
    "wbi",
}
_ORCA_SUBCOMMANDS = {
    "inp",
    "irc",
    "modred",
    "neb",
    "opt",
    "qmmm",
    "qrc",
    "scan",
    "sp",
    "ts",
}
_SUBCOMMANDS_BY_PROGRAM = {
    "gaussian": _GAUSSIAN_SUBCOMMANDS,
    "orca": _ORCA_SUBCOMMANDS,
}


@dataclass(frozen=True)
class _OptionSpec:
    name: str
    takes_value: bool = True
    flag_value: str | None = None


_RUNNER_OPTIONS = {
    "-s": _OptionSpec("server"),
    "--server": _OptionSpec("server"),
    "-n": _OptionSpec("num_cores"),
    "--num-cores": _OptionSpec("num_cores"),
    "-g": _OptionSpec("num_gpus"),
    "--num-gpus": _OptionSpec("num_gpus"),
    "-m": _OptionSpec("mem_gb"),
    "--mem-gb": _OptionSpec("mem_gb"),
    "-p": _OptionSpec("top_level_program"),
    "--program": _OptionSpec("top_level_program"),
    "-q": _OptionSpec("queue"),
    "--queue": _OptionSpec("queue"),
    "-t": _OptionSpec("time_hours"),
    "--time-hours": _OptionSpec("time_hours"),
    "--fake": _OptionSpec("fake", takes_value=False, flag_value="true"),
    "--no-fake": _OptionSpec("fake", takes_value=False, flag_value="false"),
    "--scratch": _OptionSpec("scratch", takes_value=False, flag_value="true"),
    "--no-scratch": _OptionSpec(
        "scratch", takes_value=False, flag_value="false"
    ),
    "--delete-scratch": _OptionSpec(
        "delete_scratch", takes_value=False, flag_value="true"
    ),
    "--no-delete-scratch": _OptionSpec(
        "delete_scratch", takes_value=False, flag_value="false"
    ),
    "--test": _OptionSpec("test", takes_value=False, flag_value="true"),
    "--no-test": _OptionSpec("test", takes_value=False, flag_value="false"),
    "--print-command": _OptionSpec(
        "print_command", takes_value=False, flag_value="true"
    ),
    "--no-print-command": _OptionSpec(
        "print_command", takes_value=False, flag_value="false"
    ),
    "-v": _OptionSpec("verbose", takes_value=False, flag_value="true"),
    "--verbose": _OptionSpec("verbose", takes_value=False, flag_value="true"),
    "--no-verbose": _OptionSpec(
        "verbose", takes_value=False, flag_value="false"
    ),
}

_PROGRAM_OPTIONS = {
    "-p": _OptionSpec("project"),
    "--project": _OptionSpec("project"),
    "-f": _OptionSpec("filename"),
    "--filename": _OptionSpec("filename"),
    "-l": _OptionSpec("label"),
    "--label": _OptionSpec("label"),
    "-a": _OptionSpec("append_label"),
    "--append-label": _OptionSpec("append_label"),
    "-i": _OptionSpec("index"),
    "--index": _OptionSpec("index"),
    "-c": _OptionSpec("charge"),
    "--charge": _OptionSpec("charge"),
    "-m": _OptionSpec("multiplicity"),
    "--multiplicity": _OptionSpec("multiplicity"),
    "-x": _OptionSpec("functional"),
    "--functional": _OptionSpec("functional"),
    "-b": _OptionSpec("basis"),
    "--basis": _OptionSpec("basis"),
    "-B": _OptionSpec("aux_basis"),
    "--aux-basis": _OptionSpec("aux_basis"),
    "--extrapolation-basis": _OptionSpec("extrapolation_basis"),
    "--defgrid": _OptionSpec("defgrid"),
    "--scf-tol": _OptionSpec("scf_tol"),
    "--scf-algorithm": _OptionSpec("scf_algorithm"),
    "--ri": _OptionSpec("record_index"),
    "--record-index": _OptionSpec("record_index"),
    "--rid": _OptionSpec("record_id"),
    "--record-id": _OptionSpec("record_id"),
    "--si": _OptionSpec("structure_index"),
    "--structure-index": _OptionSpec("structure_index"),
    "--sid": _OptionSpec("structure_id"),
    "--structure-id": _OptionSpec("structure_id"),
    "--mid": _OptionSpec("molecule_id"),
    "--molecule-id": _OptionSpec("molecule_id"),
    "-r": _OptionSpec("route_parameters"),
    "--additional-route-parameters": _OptionSpec("route_parameters"),
    "-o": _OptionSpec("opt_options"),
    "--additional-opt-options": _OptionSpec("opt_options"),
    "--remove-solvent": _OptionSpec(
        "remove_solvent", takes_value=False, flag_value="true"
    ),
    "--no-remove-solvent": _OptionSpec(
        "remove_solvent", takes_value=False, flag_value="false"
    ),
    "-sm": _OptionSpec("solvent_model"),
    "--solvent-model": _OptionSpec("solvent_model"),
    "-si": _OptionSpec("solvent_id"),
    "--solvent-id": _OptionSpec("solvent_id"),
    "-so": _OptionSpec("solvent_options"),
    "--solvent-options": _OptionSpec("solvent_options"),
    "-A": _OptionSpec("append_additional_info"),
    "--append-additional-info": _OptionSpec("append_additional_info"),
    "--title": _OptionSpec("title"),
    "-t": _OptionSpec("title"),
}

_SUBCOMMAND_OPTIONS = {
    "-f": _OptionSpec("freeze_atoms"),
    "--freeze-atoms": _OptionSpec("freeze_atoms"),
    "-c": _OptionSpec("coordinates"),
    "--coordinates": _OptionSpec("coordinates"),
    "-s": _OptionSpec("step_size_or_solv"),
    "--step-size": _OptionSpec("step_size"),
    "-n": _OptionSpec("num_steps_or_every_n_points"),
    "--num-steps": _OptionSpec("num_steps"),
    "--fragment-indices": _OptionSpec("fragment_indices"),
    "--every-n-points": _OptionSpec("every_n_points"),
    "--states": _OptionSpec("states"),
    "--root": _OptionSpec("root"),
    "--nstates": _OptionSpec("nstates"),
    "--eqsolv": _OptionSpec("eqsolv"),
    "--recalc-hess": _OptionSpec("recalc_hess"),
    "--trust-radius": _OptionSpec("trust_radius"),
    "--tssearch-type": _OptionSpec("tssearch_type"),
    "-e": _OptionSpec("ending_xyzfile"),
    "--ending-xyzfile": _OptionSpec("ending_xyzfile"),
    "--nimages": _OptionSpec("nimages"),
    "--joboption": _OptionSpec("joboption"),
    "-x": _OptionSpec("dist_start"),
    "-y": _OptionSpec("dist_end"),
    "-cc": _OptionSpec("constrained_coordinates"),
    "--constrained-coordinates": _OptionSpec("constrained_coordinates"),
    "--high-level-atoms": _OptionSpec("high_level_atoms"),
    "--low-level-atoms": _OptionSpec("low_level_atoms"),
    "-S": _OptionSpec("skip_completed", takes_value=False, flag_value="true"),
    "-R": _OptionSpec("skip_completed", takes_value=False, flag_value="false"),
    "-j": _OptionSpec("jobtype"),
    "--jobtype": _OptionSpec("jobtype"),
    "-ns": _OptionSpec("num_structures_to_run"),
    "--num-structures-to-run": _OptionSpec("num_structures_to_run"),
    "-g": _OptionSpec("grouping_strategy"),
    "--grouping-strategy": _OptionSpec("grouping_strategy"),
    "--proportion-structures-to-use": _OptionSpec(
        "proportion_structures_to_use"
    ),
    "--skip-completed": _OptionSpec(
        "skip_completed", takes_value=False, flag_value="true"
    ),
    "--no-skip-completed": _OptionSpec(
        "skip_completed", takes_value=False, flag_value="false"
    ),
}

# Gaussian reuses short option letters across nested job commands.  They must
# be read relative to the selected job: ``td -r 2`` is an excited-state root,
# while ``userjob -r pop=(nbo)`` is a custom Gaussian route.  A global map
# silently corrupts one of those meanings.
# ``qmmm``/ONIOM carries per-layer charge and multiplicity as subcommand
# options (``gaussian qmmm`` and ``orca qmmm`` share these letters). They are
# distinct quantities from the program-level ``-c``/``-m`` model-system charge,
# so a request that pins the total or a layer state can only be asserted if the
# parser exposes them instead of dropping them as unrecognized. See the real
# CLI in ``chemsmart/cli/{gaussian,orca}/qmmm.py``.
_QMMM_CHARGE_MULT_OPTIONS: dict[str, _OptionSpec] = {
    "-ct": _OptionSpec("charge_total"),
    "--charge-total": _OptionSpec("charge_total"),
    "-mt": _OptionSpec("mult_total"),
    "--mult-total": _OptionSpec("mult_total"),
    "-ci": _OptionSpec("charge_intermediate"),
    "--charge-intermediate": _OptionSpec("charge_intermediate"),
    "-mi": _OptionSpec("mult_intermediate"),
    "--mult-intermediate": _OptionSpec("mult_intermediate"),
    "-ch": _OptionSpec("charge_high"),
    "--charge-high": _OptionSpec("charge_high"),
    "-mh": _OptionSpec("mult_high"),
    "--mult-high": _OptionSpec("mult_high"),
}

# Per-layer atom-region options (short forms plus the long forms not already in
# ``_SUBCOMMAND_OPTIONS``). Gaussian ONIOM uses high/medium/low; ORCA uses
# high/intermediate/active. Source: ``chemsmart/cli/{gaussian,orca}/qmmm.py``.
_GAUSSIAN_QMMM_OPTIONS: dict[str, _OptionSpec] = {
    **_QMMM_CHARGE_MULT_OPTIONS,
    "-ha": _OptionSpec("high_level_atoms"),
    "-ma": _OptionSpec("medium_level_atoms"),
    "--medium-level-atoms": _OptionSpec("medium_level_atoms"),
    "-la": _OptionSpec("low_level_atoms"),
    "-ba": _OptionSpec("bonded_atoms"),
    "--bonded-atoms": _OptionSpec("bonded_atoms"),
}
_ORCA_QMMM_OPTIONS: dict[str, _OptionSpec] = {
    **_QMMM_CHARGE_MULT_OPTIONS,
    "-ha": _OptionSpec("high_level_atoms"),
    "-ia": _OptionSpec("intermediate_level_atoms"),
    "--intermediate-level-atoms": _OptionSpec("intermediate_level_atoms"),
    "-a": _OptionSpec("active_atoms"),
    "--active-atoms": _OptionSpec("active_atoms"),
}

_JOB_OPTION_OVERRIDES: dict[tuple[str, str], dict[str, _OptionSpec]] = {
    ("gaussian", "td"): {
        "-s": _OptionSpec("states"),
        "-r": _OptionSpec("root"),
        "-n": _OptionSpec("nstates"),
        "-e": _OptionSpec("eqsolv"),
        "--states": _OptionSpec("states"),
        "--root": _OptionSpec("root"),
        "--nstates": _OptionSpec("nstates"),
        "--eqsolv": _OptionSpec("eqsolv"),
    },
    ("gaussian", "userjob"): {
        "-r": _OptionSpec("route_parameters"),
        "--route": _OptionSpec("route_parameters"),
    },
    ("gaussian", "qmmm"): dict(_GAUSSIAN_QMMM_OPTIONS),
    ("orca", "qmmm"): dict(_ORCA_QMMM_OPTIONS),
    # ORCA TS overloads ``-s`` as ``--recalc-hess`` (steps between Hessian
    # recalculations); the generic ``-s`` is a scan step size, which a TS search
    # never takes. Source: ``chemsmart/cli/orca/ts.py``.
    ("orca", "ts"): {
        "-s": _OptionSpec("recalc_hess"),
        "--recalc-hess": _OptionSpec("recalc_hess"),
    },
}


def parse_model_command(
    command: str, *, cwd: str | os.PathLike[str] | None = None
) -> ParsedModelCommand:
    workspace = str(Path(cwd or os.getcwd()).resolve())
    warnings: list[str] = []
    try:
        tokens = shlex.split(command)
    except ValueError as exc:
        return ParsedModelCommand(
            command=command,
            workspace=workspace,
            tokens=[],
            parse_error=f"tokenization failed: {exc}",
        )

    if not tokens:
        return ParsedModelCommand(
            command=command,
            workspace=workspace,
            tokens=tokens,
            parse_error="empty command",
        )

    entrypoint = tokens[0]
    if entrypoint != "chemsmart":
        warnings.append("command does not start with the chemsmart entrypoint")

    action_index = _find_first(tokens, _TOP_LEVEL_COMMANDS, start=1)
    if action_index is None:
        return ParsedModelCommand(
            command=command,
            workspace=workspace,
            tokens=tokens,
            entrypoint=entrypoint,
            parse_error="missing chemsmart run/sub action",
            warnings=warnings,
        )
    action = tokens[action_index]

    program_index = _find_program_index(tokens, action_index + 1)
    if program_index is None:
        return ParsedModelCommand(
            command=command,
            workspace=workspace,
            tokens=tokens,
            entrypoint=entrypoint,
            action=action,
            parse_error="missing computational program after run/sub options",
            warnings=warnings,
        )
    program = tokens[program_index]

    runner_opts, runner_warnings = _parse_options(
        tokens[action_index + 1 : program_index], _RUNNER_OPTIONS
    )
    warnings.extend(runner_warnings)

    subcommand_index = _find_subcommand_index(
        tokens, program, program_index + 1
    )
    if subcommand_index is None:
        subcommand_index = len(tokens)
        job = None
        warnings.append(f"no recognized {program} job subcommand was found")
    else:
        job = tokens[subcommand_index]

    # The real CLI nests ``qmmm`` under a parent job (``opt qmmm``, ``ts qmmm``,
    # ``sp qmmm`` …). ``_find_subcommand_index`` stops at the parent token, so
    # the parent remains the job (and ``qmmm`` is detected as a required token
    # by the intent oracle), but the per-layer charge/mult/region options that
    # follow the nested ``qmmm`` token must be parsed with the qmmm option specs
    # so they land in ``structural_options`` instead of being reported as
    # unrecognized parent-job options. Source:
    # ``chemsmart/cli/{gaussian,orca}/qmmm.py``.
    subcommand_option_index = subcommand_index
    subcommand_specs_job = job
    if job is not None and job != "qmmm":
        nested_qmmm_index = _find_first(
            tokens, {"qmmm"}, start=subcommand_index + 1
        )
        if nested_qmmm_index is not None:
            subcommand_specs_job = "qmmm"
            subcommand_option_index = nested_qmmm_index

    program_opts, program_warnings = _parse_options(
        tokens[program_index + 1 : subcommand_index], _PROGRAM_OPTIONS
    )
    warnings.extend(program_warnings)
    subcommand_specs = dict(_SUBCOMMAND_OPTIONS)
    subcommand_specs.update(
        _JOB_OPTION_OVERRIDES.get((program, subcommand_specs_job or ""), {})
    )
    subcommand_opts, subcommand_warnings = _parse_options(
        tokens[subcommand_option_index + 1 :], subcommand_specs
    )
    warnings.extend(subcommand_warnings)
    if (
        program == "gaussian"
        and job == "traj"
        and "dist_start" in subcommand_opts
    ):
        # ``traj`` owns its post-subcommand ``-x`` as a selection proportion;
        # ORCA scan retains the distinct distance-start meaning.
        subcommand_opts["proportion_structures_to_use"] = subcommand_opts.pop(
            "dist_start"
        )

    project = program_opts.get("project")
    resolved, resolve_warning = _resolve_project_method(
        program=program,
        project=project,
        job=job,
        overrides=program_opts,
    )
    if resolve_warning:
        warnings.append(resolve_warning)

    server = runner_opts.get("server")
    dry_run = (
        runner_opts.get("fake") == "true" or runner_opts.get("test") == "true"
    )
    structural_options = {
        key: value
        for key, value in subcommand_opts.items()
        if value is not None
        and key not in {"skip_completed", "route_parameters"}
    }
    resources = {
        key: value
        for key, value in runner_opts.items()
        if key
        in {
            "num_cores",
            "num_gpus",
            "mem_gb",
            "queue",
            "time_hours",
            "scratch",
            "delete_scratch",
        }
    }

    return ParsedModelCommand(
        command=command,
        workspace=workspace,
        tokens=tokens,
        entrypoint=entrypoint,
        action=action,
        program=program,
        job=job,
        server=server,
        dry_run=dry_run,
        project=project,
        project_p_flag_meaning=(
            "program-level -p/--project for gaussian/orca project settings"
            if project is not None
            else None
        ),
        top_level_program=runner_opts.get("top_level_program"),
        filename=program_opts.get("filename"),
        record_index=program_opts.get("record_index"),
        record_id=program_opts.get("record_id"),
        structure_index=program_opts.get("structure_index"),
        structure_id=program_opts.get("structure_id"),
        molecule_id=program_opts.get("molecule_id"),
        label=program_opts.get("label"),
        charge=program_opts.get("charge"),
        multiplicity=program_opts.get("multiplicity"),
        functional=resolved.get("functional"),
        ab_initio=resolved.get("ab_initio"),
        basis=resolved.get("basis"),
        aux_basis=resolved.get("aux_basis"),
        extrapolation_basis=resolved.get("extrapolation_basis"),
        defgrid=resolved.get("defgrid"),
        scf_tol=resolved.get("scf_tol"),
        scf_algorithm=resolved.get("scf_algorithm"),
        solvent_model=resolved.get("solvent_model"),
        solvent_id=resolved.get("solvent_id"),
        route_parameters=(
            program_opts.get("route_parameters")
            or subcommand_opts.get("route_parameters")
        ),
        opt_options=program_opts.get("opt_options"),
        structural_options=structural_options,
        resources=resources,
        warnings=warnings,
    )


def format_model_command_explanation(
    command: str,
    *,
    cwd: str | os.PathLike[str] | None = None,
) -> str:
    return format_parsed_model_command(parse_model_command(command, cwd=cwd))


def _find_first(
    tokens: list[str], choices: set[str], *, start: int
) -> int | None:
    for index in range(start, len(tokens)):
        if tokens[index] in choices:
            return index
    return None


def _find_program_index(tokens: list[str], start: int) -> int | None:
    index = start
    while index < len(tokens):
        token = tokens[index]
        if token in _PROGRAMS:
            return index
        spec, _value = _option_spec_for_token(token, _RUNNER_OPTIONS)
        if spec is not None:
            index += 2 if spec.takes_value and "=" not in token else 1
        else:
            index += 1
    return None


def _find_subcommand_index(
    tokens: list[str], program: str, start: int
) -> int | None:
    subcommands = _SUBCOMMANDS_BY_PROGRAM.get(program, set())
    index = start
    while index < len(tokens):
        token = tokens[index]
        if token in subcommands:
            return index
        spec, _value = _option_spec_for_token(token, _PROGRAM_OPTIONS)
        if spec is not None:
            index += 2 if spec.takes_value and "=" not in token else 1
        else:
            index += 1
    return None


def _parse_options(
    tokens: list[str], specs: dict[str, _OptionSpec]
) -> tuple[dict[str, str], list[str]]:
    values: dict[str, str] = {}
    warnings: list[str] = []
    index = 0
    while index < len(tokens):
        token = tokens[index]
        spec, inline_value = _option_spec_for_token(token, specs)
        if spec is None:
            if token.startswith("-"):
                warnings.append(f"unrecognized option `{token}`")
            index += 1
            continue
        if spec.takes_value:
            if inline_value is not None:
                values[spec.name] = inline_value
                index += 1
            elif index + 1 < len(tokens):
                values[spec.name] = tokens[index + 1]
                index += 2
            else:
                warnings.append(f"option `{token}` is missing a value")
                index += 1
        else:
            values[spec.name] = spec.flag_value or "true"
            index += 1
    return values, warnings


def _option_spec_for_token(
    token: str, specs: dict[str, _OptionSpec]
) -> tuple[_OptionSpec | None, str | None]:
    if token in specs:
        return specs[token], None
    if token.startswith("--") and "=" in token:
        flag, value = token.split("=", 1)
        if flag in specs:
            return specs[flag], value
    return None, None


def _resolve_project_method(
    *,
    program: str,
    project: str | None,
    job: str | None,
    overrides: dict[str, str],
) -> tuple[dict[str, str | None], str | None]:
    resolved: dict[str, str | None] = {
        "ab_initio": None,
        "functional": None,
        "basis": None,
        "aux_basis": None,
        "extrapolation_basis": None,
        "defgrid": None,
        "scf_tol": None,
        "scf_algorithm": None,
        "solvent_model": None,
        "solvent_id": None,
    }
    if project is None or program not in {"gaussian", "orca"}:
        _apply_method_overrides(resolved, overrides)
        return resolved, None

    try:
        with _quiet_project_settings_resolution():
            settings = _load_job_settings(
                program=program, project=project, job=job
            )
    except Exception as exc:  # pragma: no cover - message content varies
        _apply_method_overrides(resolved, overrides)
        return (
            resolved,
            f"project settings could not be resolved for `{program}:{project}`: {exc}",
        )

    for key in (
        "ab_initio",
        "functional",
        "basis",
        "aux_basis",
        "extrapolation_basis",
        "defgrid",
        "scf_tol",
        "scf_algorithm",
        "solvent_model",
        "solvent_id",
    ):
        value = getattr(settings, key, None)
        if value is not None:
            resolved[key] = str(value)
    _apply_method_overrides(resolved, overrides)
    if overrides.get("remove_solvent") == "true":
        resolved["solvent_model"] = None
        resolved["solvent_id"] = None
    return resolved, None


def _load_job_settings(*, program: str, project: str, job: str | None) -> Any:
    if program == "gaussian":
        from chemsmart.settings.gaussian import GaussianProjectSettings

        project_settings = GaussianProjectSettings.from_project(project)
    else:
        from chemsmart.settings.orca import ORCAProjectSettings

        project_settings = ORCAProjectSettings.from_project(project)

    method_name = _project_method_name(program=program, job=job)
    method = getattr(project_settings, method_name, None)
    if method is None:
        method = getattr(project_settings, "main_settings", None)
    if method is None:
        raise AttributeError(f"project settings has no `{method_name}`")
    return method()


def _project_method_name(*, program: str, job: str | None) -> str:
    if program == "gaussian" and job in {"dias", "resp"}:
        return "sp_settings"
    if program == "gaussian" and job == "td":
        return "td_settings"
    if job in {
        "opt",
        "modred",
        "ts",
        "irc",
        "scan",
        "nci",
        "sp",
        "wbi",
        "qmmm",
        "neb",
    }:
        return f"{job}_settings"
    return "main_settings"


def _apply_method_overrides(
    resolved: dict[str, str | None], overrides: dict[str, str]
) -> None:
    for key in (
        "functional",
        "basis",
        "aux_basis",
        "extrapolation_basis",
        "defgrid",
        "scf_tol",
        "scf_algorithm",
        "solvent_model",
        "solvent_id",
    ):
        if overrides.get(key) is not None:
            resolved[key] = overrides[key]


@contextlib.contextmanager
def _quiet_project_settings_resolution():
    logger_names = (
        "chemsmart.jobs.settings",
        "chemsmart.settings.gaussian",
        "chemsmart.settings.orca",
    )
    loggers = [logging.getLogger(name) for name in logger_names]
    previous_levels = [logger.level for logger in loggers]
    previous_disabled = [logger.disabled for logger in loggers]
    try:
        for logger in loggers:
            logger.setLevel(logging.CRITICAL + 1)
            logger.disabled = True
        with contextlib.redirect_stdout(io.StringIO()):
            with contextlib.redirect_stderr(io.StringIO()):
                yield
    finally:
        for logger, level, disabled in zip(
            loggers, previous_levels, previous_disabled
        ):
            logger.setLevel(level)
            logger.disabled = disabled
