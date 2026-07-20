"""Generated-input rendering grounded in an equivalent ChemSmart CLI command."""

from __future__ import annotations

import json
import os
import shlex
from types import SimpleNamespace
from typing import Any

from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.gaussian.writer import GaussianInputWriter
from chemsmart.jobs.job import Job
from chemsmart.jobs.orca.job import ORCAJob
from chemsmart.jobs.orca.writer import ORCAInputWriter

SUPPORTED_SUBMIT_JOBTYPES = {"opt", "ts", "sp", "irc", "scan"}
_PROJECT_RUNTIME_SETTING_FIELDS = {
    "functional",
    "basis",
    "semiempirical",
    "solvent_model",
    "solvent_id",
    "custom_solvent",
    "heavy_elements",
    "heavy_elements_basis",
    "ab_initio",
    "dispersion",
    "aux_basis",
    "extrapolation_basis",
    "defgrid",
    "scf_tol",
    "scf_algorithm",
    "scf_maxiter",
    "scf_convergence",
}


def dry_run_input(job: Job) -> dict[str, Any]:
    """Render a job input file and return its absolute path and contents."""

    command = reconstruct_run_cli_command(job)
    target_directory = os.path.abspath(job.folder)
    job.set_folder(target_directory)
    if job.jobrunner is None:
        job.jobrunner = SimpleNamespace(num_cores=12, mem_gb=16)
    if isinstance(job, GaussianJob):
        input_writer = GaussianInputWriter(job=job)
    elif isinstance(job, ORCAJob):
        input_writer = ORCAInputWriter(job=job)
    else:
        raise ValueError(
            "dry_run_input only supports GaussianJob and ORCAJob instances"
        )
    input_writer.write(target_directory=target_directory)
    inputfile = os.path.abspath(job.inputfile)
    with open(inputfile) as file:
        content = file.read()
    result = {
        "inputfile": inputfile,
        "content": content,
        "command": command,
        "cli_grounded": command is not None,
    }
    if command is None:
        result["cli_grounding_issue"] = (
            "dry_run_input could not reconstruct an equivalent "
            "chemsmart CLI command for this job"
        )
    return result


def reconstruct_run_cli_command(job: Job) -> str | None:
    """Return an equivalent user-facing ``chemsmart run`` command."""

    try:
        argv = _reconstruct_run_cli_args(job)
    except Exception:
        return None
    if not argv:
        return None
    return " ".join(shlex.quote(str(part)) for part in argv)


def _reconstruct_run_cli_args(job: Job) -> list[str] | None:
    source_filepath = getattr(job, "_agent_source_filepath", None)
    if not isinstance(source_filepath, str) or not source_filepath.strip():
        return None
    settings = getattr(job, "settings", None)
    program_name = getattr(job, "PROGRAM", None)
    if not isinstance(program_name, str):
        return None
    program = program_name.lower()
    if program not in {"gaussian", "orca"}:
        return None
    cli_jobtype = _run_cli_jobtype_name(job)
    argv = ["chemsmart", "run", program]
    project = getattr(job, "_agent_project_name", None)
    project_backed = isinstance(project, str) and bool(project.strip())
    if project_backed:
        argv.extend(["-p", project.strip()])
    argv.extend(
        _program_cli_args(
            settings,
            force_route_freq=_jobtype_name_for_settings(job) == "freq",
            project_backed=project_backed,
        )
    )
    argv.extend(["-f", source_filepath])
    index = _agent_source_index_for_cli(job)
    if index is not None:
        argv.extend(["-i", index])
    label = getattr(job, "label", None)
    if isinstance(label, str) and label.strip():
        argv.extend(["-l", label])
    argv.append(cli_jobtype)
    argv.extend(_jobtype_cli_args(job))
    return argv


def _program_cli_args(
    settings: Any,
    *,
    force_route_freq: bool = False,
    project_backed: bool = False,
) -> list[str]:
    if settings is None:
        return []
    flag_map = [
        ("charge", "-c"),
        ("multiplicity", "-m"),
        ("functional", "-x"),
        ("basis", "-b"),
        ("semiempirical", "--semiempirical"),
        ("additional_route_parameters", "--additional-route-parameters"),
        ("additional_opt_options_in_route", "--additional-opt-options"),
        ("solvent_model", "-sm"),
        ("solvent_id", "-si"),
        ("custom_solvent", "--custom-solvent"),
        ("heavy_elements", "--heavy-elements"),
        ("heavy_elements_basis", "--heavy-elements-basis"),
        ("title", "--title"),
        ("ab_initio", "--ab-initio"),
        ("dispersion", "-D"),
        ("aux_basis", "-B"),
        ("extrapolation_basis", "-e"),
        ("defgrid", "--defgrid"),
        ("scf_tol", "--scf-tol"),
        ("scf_algorithm", "--scf-algorithm"),
        ("scf_maxiter", "--scf-maxiter"),
        ("scf_convergence", "--scf-convergence"),
    ]
    argv: list[str] = []
    for name, flag in flag_map:
        if project_backed and name in _PROJECT_RUNTIME_SETTING_FIELDS:
            continue
        if not hasattr(settings, name):
            continue
        value = getattr(settings, name)
        if value is None or value is False:
            continue
        if isinstance(value, (list, tuple)):
            value = ",".join(str(item) for item in value)
        if name == "additional_route_parameters" and force_route_freq:
            tokens = str(value).replace(",", " ").split()
            if not any(token.lower() == "freq" for token in tokens):
                value = f"{value} freq"
            force_route_freq = False
        argv.extend([flag, str(value)])
    if force_route_freq:
        argv.extend(["--additional-route-parameters", "freq"])
    return argv


def _agent_source_index_for_cli(job: Job) -> str | None:
    index = getattr(job, "_agent_source_index", None)
    if not isinstance(index, str):
        return None
    return None if index == "-1" else index


def _jobtype_cli_args(job: Job) -> list[str]:
    jobtype = _run_cli_jobtype_name(job)
    settings = getattr(job, "settings", None)
    if settings is None:
        return []
    if jobtype == "scan":
        return _dict_to_cli_args(_scan_cli_overrides(settings))
    if jobtype == "modred":
        modred = getattr(settings, "modred", None)
        if isinstance(modred, dict) and modred.get("coords") is not None:
            return ["--coordinates", _json_cli_value(modred["coords"])]
    return []


def _run_cli_jobtype_name(job: Job) -> str:
    jobtype = _jobtype_name_for_settings(job)
    if jobtype == "singlepoint":
        return "sp"
    if jobtype == "freq":
        return "opt"
    if jobtype not in SUPPORTED_SUBMIT_JOBTYPES:
        supported = ", ".join(sorted({*SUPPORTED_SUBMIT_JOBTYPES, "freq"}))
        raise ValueError(
            f"Unsupported run jobtype {jobtype!r}. Supported: {supported}"
        )
    return jobtype


def _jobtype_name_for_settings(job: Job) -> str:
    settings = getattr(job, "settings", None)
    jobtype = getattr(settings, "jobtype", None)
    if not isinstance(jobtype, str):
        raise ValueError(f"Job {job!r} is missing settings.jobtype.")
    return jobtype.lower()


def _dict_to_cli_args(values: dict[str, Any]) -> list[str]:
    argv: list[str] = []
    flag_map = {
        "coordinates": "--coordinates",
        "num_steps": "--num-steps",
        "step_size": "--step-size",
        "constrained_coordinates": "--constrained-coordinates",
        "dist_start": "--dist-start",
        "dist_end": "--dist-end",
    }
    for key, value in values.items():
        flag = flag_map.get(key)
        if flag and value is not None:
            argv.extend([flag, str(value)])
    return argv


def _scan_cli_overrides(settings: Any) -> dict[str, Any]:
    modred = getattr(settings, "modred", None)
    if not isinstance(modred, dict):
        return {}
    overrides: dict[str, Any] = {}
    coords = modred.get("coords")
    if coords is not None:
        overrides["coordinates"] = _json_cli_value(coords)
    if modred.get("num_steps") is not None:
        overrides["num_steps"] = _scalar_or_json_cli_value(modred["num_steps"])
    if modred.get("step_size") is not None:
        overrides["step_size"] = _scalar_or_json_cli_value(modred["step_size"])
    if modred.get("constrained_coordinates") is not None:
        overrides["constrained_coordinates"] = _json_cli_value(
            modred["constrained_coordinates"]
        )
    if modred.get("dist_start") is not None:
        overrides["dist_start"] = _scalar_or_json_cli_value(
            modred["dist_start"]
        )
    if modred.get("dist_end") is not None:
        overrides["dist_end"] = _scalar_or_json_cli_value(modred["dist_end"])
    return overrides


def _json_cli_value(value: Any) -> str:
    return json.dumps(value, separators=(",", ":"))


def _scalar_or_json_cli_value(value: Any) -> str:
    if isinstance(value, (list, tuple)) and len(value) == 1:
        return str(value[0])
    return _json_cli_value(value)


__all__ = [
    "SUPPORTED_SUBMIT_JOBTYPES",
    "dry_run_input",
    "reconstruct_run_cli_command",
]
