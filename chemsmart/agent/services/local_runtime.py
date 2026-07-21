"""Local runtime validation, execution, and result extraction."""

from __future__ import annotations

import contextlib
import os
import re
import shutil
import traceback
from typing import Any

from chemsmart.agent.runtime.result_parsing import inspect_output
from chemsmart.agent.services.server_selection import coerce_server
from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.orca.output import ORCAOutput
from chemsmart.io.xtb.folder import XTBFolder
from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.job import Job
from chemsmart.jobs.orca.job import ORCAJob
from chemsmart.jobs.runner import JobRunner
from chemsmart.jobs.xtb.job import XTBJob
from chemsmart.settings.server import Server
from chemsmart.utils.periodictable import PeriodicTable

_REMOTE_UNKNOWN_SERVER_FIELDS = [
    "server.queue required",
    "server.account required",
    "server.scratch_dir required",
    "server.modules_or_executable_path required",
]
_REMOTE_UNKNOWN_SCRATCH = "scratch_dir writable on HPC"
_REMOTE_UNKNOWN_MODULES = "module load succeeds on HPC"
_REMOTE_UNKNOWN_QUEUE = "queue accepts jobs"
_REMOTE_UNKNOWN_ACCOUNT = "account has remaining hours"
_REMOTE_UNKNOWN_SSH = "ssh login reachable"
_UNRESOLVED_ENVVAR_PATTERN = re.compile(
    r"(\$(?:\{)?[A-Za-z_][A-Za-z0-9_]*(?:\})?)|(%[^%]+%)"
)
_PERIODIC_TABLE = PeriodicTable()


def validate_runtime(job: Job, server=None) -> dict[str, Any]:
    """Check local/runtime prerequisites and remote unknowns for a job."""

    local_issues = _validate_job_fields(job)
    remote_unknown = []
    if server is None:
        remote_unknown.extend(_REMOTE_UNKNOWN_SERVER_FIELDS)
        return _runtime_validation_result(local_issues, remote_unknown)
    if not isinstance(server, Server):
        server = _resolve_server(server, local_issues)
        if server is None:
            return _runtime_validation_result(local_issues, remote_unknown)
    server_config = getattr(server, "config", None)
    if not isinstance(server_config, dict):
        local_issues.append("server.config invalid")
        return _runtime_validation_result(local_issues, remote_unknown)
    queue_name = server_config.get("QUEUE_NAME")
    account_name = server_config.get("ACCOUNT", server_config.get("PROJECT"))
    _validate_server_fields(
        server_config,
        queue_name,
        account_name,
        local_issues,
        remote_unknown,
    )
    runner = _build_runner(job, server, local_issues)
    scratch_dir, executable_path, modules = _runtime_resources(
        runner,
        local_issues,
    )
    has_modules = _validate_resolved_resources(
        scratch_dir,
        executable_path,
        modules,
        local_issues,
    )
    _append_remote_unknowns(
        remote_unknown,
        scratch_dir=scratch_dir,
        has_modules=has_modules,
        queue_name=queue_name,
        account_name=account_name,
    )
    return _runtime_validation_result(local_issues, remote_unknown)


def _resolve_server(server, local_issues: list[str]) -> Server | None:
    try:
        return coerce_server(server)
    except Exception as exc:
        local_issues.append(f"server.config invalid: {exc}")
        return None


def _validate_server_fields(
    config: dict[str, Any],
    queue_name: Any,
    account_name: Any,
    local_issues: list[str],
    remote_unknown: list[str],
) -> None:
    if not _is_non_empty_string(queue_name):
        local_issues.append("server.queue missing")
    if not _is_non_empty_string(account_name):
        remote_unknown.append(_REMOTE_UNKNOWN_ACCOUNT)
    scratch_value = config.get("SCRATCH_DIR")
    if scratch_value is not None and not isinstance(scratch_value, str):
        local_issues.append("server.scratch_dir invalid")


def _validate_resolved_resources(
    scratch_dir: str | None,
    executable_path: Any,
    modules: Any,
    local_issues: list[str],
) -> bool:
    if scratch_dir is not None and _has_unresolved_envvars(scratch_dir):
        local_issues.append("server.scratch_dir unresolved")
    if modules is not None and not isinstance(modules, str):
        local_issues.append("server.modules invalid")
    if executable_path is not None and not isinstance(executable_path, str):
        local_issues.append("server.executable_path invalid")
    has_modules = _is_non_empty_string(modules)
    has_executable = _is_non_empty_string(executable_path)
    if not has_modules and not has_executable:
        local_issues.append("server.modules_or_executable_path missing")
    if has_executable and shutil.which(executable_path) is None:
        local_issues.append("server.executable_path missing")
    return has_modules


def _build_runner(
    job: Job,
    server: Server,
    local_issues: list[str],
) -> JobRunner | None:
    try:
        return JobRunner.from_job(job=job, server=server, scratch=False)
    except Exception:
        local_issues.append("jobrunner unavailable")
        return None


def _runtime_resources(
    runner: JobRunner | None,
    local_issues: list[str],
) -> tuple[str | None, Any, Any]:
    if runner is None:
        return None, None, None
    scratch_dir = None
    try:
        scratch_dir = runner._resolve_scratch_dir_candidate()
        if scratch_dir is not None:
            scratch_dir = runner._normalize_path_string(scratch_dir)
    except Exception:
        local_issues.append("server.scratch_dir invalid")
    try:
        executable = runner.executable
        return scratch_dir, executable.get_executable(), executable.modules
    except Exception:
        return scratch_dir, None, None


def _append_remote_unknowns(
    remote_unknown: list[str],
    *,
    scratch_dir: str | None,
    has_modules: bool,
    queue_name: Any,
    account_name: Any,
) -> None:
    if _is_non_empty_string(scratch_dir):
        remote_unknown.append(_REMOTE_UNKNOWN_SCRATCH)
    if has_modules:
        remote_unknown.append(_REMOTE_UNKNOWN_MODULES)
    if _is_non_empty_string(queue_name):
        remote_unknown.append(_REMOTE_UNKNOWN_QUEUE)
    if _is_non_empty_string(account_name):
        remote_unknown.append(_REMOTE_UNKNOWN_ACCOUNT)
    remote_unknown.append(_REMOTE_UNKNOWN_SSH)


def run_local(job: Job) -> dict[str, Any]:
    """Execute a job locally and summarize the generated output artifacts."""

    folder = os.path.abspath(job.folder)
    os.makedirs(folder, exist_ok=True)
    job.set_folder(folder)
    stdout_path = os.path.join(folder, f"{job.label}.stdout")
    stderr_path = os.path.join(folder, f"{job.label}.stderr")
    job.local = True
    returncode = _execute_job(job, stdout_path, stderr_path)
    summary = _summarize_local_output(job) if returncode == 0 else {}
    return {
        "ok": returncode == 0,
        "returncode": returncode,
        "stdout_path": stdout_path,
        "stderr_path": stderr_path,
        "output_summary": summary,
    }


def _execute_job(job: Job, stdout_path: str, stderr_path: str) -> int:
    with (
        open(stdout_path, "w") as stdout_handle,
        open(stderr_path, "w") as stderr_handle,
    ):
        try:
            with (
                contextlib.redirect_stdout(stdout_handle),
                contextlib.redirect_stderr(stderr_handle),
            ):
                job.run()
        except Exception as exc:
            traceback.print_exc(file=stderr_handle)
            return int(getattr(exc, "returncode", 1) or 1)
    return 0


def extract_optimized_geometry(job: Job) -> Molecule:
    """Extract the final optimized geometry from a completed job log."""

    if isinstance(job, XTBJob):
        return _extract_xtb_geometry(job)
    logfile = _resolve_geometry_logfile(job)
    if not os.path.exists(logfile):
        raise FileNotFoundError(f"Output log not found: {logfile}")
    with open(logfile, encoding="utf-8", errors="ignore") as file:
        lines = file.readlines()
    if isinstance(job, GaussianJob):
        symbols, positions = _parse_gaussian_geometry(lines)
    elif isinstance(job, ORCAJob):
        symbols, positions = _parse_orca_geometry(lines)
    else:
        raise ValueError(
            "extract_optimized_geometry only supports GaussianJob, ORCAJob, "
            "and XTBJob instances"
        )
    return Molecule(
        symbols=symbols,
        positions=positions,
        charge=getattr(job.settings, "charge", None),
        multiplicity=getattr(job.settings, "multiplicity", None),
    )


def _extract_xtb_geometry(job: Job) -> Molecule:
    # xTB does not print a final-geometry table in its .out file; the
    # optimizer writes the relaxed structure to xtbopt.* in the job folder
    # (same format as the input, .xyz for chemsmart-driven jobs).
    folder = XTBFolder(folder=str(job.folder))
    geometry_path = folder._xtbopt_geometry()
    if geometry_path is None:
        raise FileNotFoundError(
            "No optimized geometry (xtbopt.*) found in "
            f"{job.folder}; only completed xtb opt jobs write one"
        )
    molecule = Molecule.from_filepath(
        filepath=geometry_path,
        index="-1",
        return_list=False,
    )
    molecule.charge = getattr(job.settings, "charge", None)
    molecule.multiplicity = getattr(job.settings, "multiplicity", None)
    # xtbopt.* is a real on-disk file, so it can ground a follow-up
    # `chemsmart run/sub ... -f` command directly.
    setattr(molecule, "_agent_source_filepath", str(geometry_path))
    setattr(molecule, "_agent_source_index", "-1")
    return molecule


def _validate_job_fields(job: Job) -> list[str]:
    issues = []
    if getattr(job, "molecule", None) is None:
        issues.append("job.molecule missing")
    if getattr(job, "settings", None) is None:
        issues.append("job.settings missing")
    if not _is_non_empty_string(getattr(job, "label", None)):
        issues.append("job.label missing")
    return issues


def _runtime_validation_result(
    local_issues: list[str],
    remote_unknown: list[str],
) -> dict[str, Any]:
    status = "fail" if local_issues else "partial" if remote_unknown else "ok"
    return {
        "ok": status,
        "local_ok": not local_issues,
        "local_issues": local_issues,
        "remote_unknown": remote_unknown,
    }


def _is_non_empty_string(value: Any) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _resolve_geometry_logfile(job: Job) -> str:
    outputfile = getattr(job, "outputfile", None)
    if isinstance(outputfile, str) and os.path.exists(outputfile):
        return outputfile
    fallbacks = [
        os.path.join(job.folder, f"{job.label}.log"),
        os.path.join(job.folder, f"{job.label}.out"),
    ]
    for path in fallbacks:
        if os.path.exists(path):
            return path
    return outputfile if isinstance(outputfile, str) else fallbacks[0]


def _parse_gaussian_geometry(
    lines: list[str],
) -> tuple[list[str], list[list[float]]]:
    blocks: list[tuple[list[str], list[list[float]]]] = []
    for index, line in enumerate(lines):
        if "Standard orientation:" not in line:
            continue
        symbols = []
        positions = []
        for coord_line in lines[index + 5 :]:
            stripped = coord_line.strip()
            if not stripped or set(stripped) == {"-"}:
                break
            parts = stripped.split()
            if len(parts) < 6:
                break
            symbols.append(_PERIODIC_TABLE.to_symbol(int(parts[1])))
            positions.append([float(value) for value in parts[3:6]])
        if symbols:
            blocks.append((symbols, positions))
    if not blocks:
        raise ValueError("No optimized Gaussian geometry found in output log")
    return blocks[-1]


def _parse_orca_geometry(
    lines: list[str],
) -> tuple[list[str], list[list[float]]]:
    blocks: list[tuple[list[str], list[list[float]]]] = []
    for index, line in enumerate(lines):
        if "CARTESIAN COORDINATES (ANGSTROEM)" not in line:
            continue
        symbols = []
        positions = []
        for coord_line in lines[index + 2 :]:
            stripped = coord_line.strip()
            if (
                not stripped
                or set(stripped) == {"-"}
                or "CARTESIAN COORDINATES" in stripped
            ):
                break
            parts = stripped.split()
            if len(parts) < 4:
                break
            symbols.append(parts[0])
            positions.append([float(value) for value in parts[1:4]])
        if symbols:
            blocks.append((symbols, positions))
    if not blocks:
        raise ValueError("No optimized ORCA geometry found in output log")
    return blocks[-1]


def _has_unresolved_envvars(path: str) -> bool:
    return bool(_UNRESOLVED_ENVVAR_PATTERN.search(path))


def _summarize_local_output(job: Job) -> dict[str, Any]:
    if isinstance(job, XTBJob):
        return _summarize_xtb_output(job)
    if isinstance(job, GaussianJob):
        parser_cls = Gaussian16Output
        output_path = job.outputfile
    elif isinstance(job, ORCAJob):
        parser_cls = ORCAOutput
        output_path = job.outputfile
    else:
        return {}
    if output_path is None or not os.path.exists(output_path):
        return {}
    try:
        parser = parser_cls(output_path)
        energies = list(getattr(parser, "energies", []) or [])
        energy = energies[-1] if energies else None
        converged = bool(getattr(parser, "normal_termination", False))
        frequencies = list(
            getattr(parser, "vibrational_frequencies", []) or []
        )
        imag_freqs = [freq for freq in frequencies if freq < 0]
        structures = list(getattr(parser, "all_structures", []) or [])
        geometry_count = len(structures)
        if (
            energy is None
            and not converged
            and not imag_freqs
            and not geometry_count
        ):
            return {}
        return {
            "energy": energy,
            "converged": converged,
            "imag_freqs": imag_freqs,
            "optimized_geometry_count": geometry_count,
        }
    except Exception:
        return {}


def _summarize_xtb_output(job: Job) -> dict[str, Any]:
    # xTB has no Gaussian/ORCA output parser class; reuse the deterministic
    # xtb-aware parser that reads `TOTAL ENERGY ... Eh`, `* finished run`,
    # and the eigval frequency rows. Without this the model gets an empty
    # summary from a genuine run and cannot report the energy that the
    # pre-optimization chain depends on.
    output_path = getattr(job, "outputfile", None)
    if output_path is None or not os.path.exists(output_path):
        return {}
    inspected = inspect_output(output_path, program="xtb")
    energy = inspected.get("energy")
    converged = bool(inspected.get("normal_termination"))
    imag_freqs = list(inspected.get("imag_freqs") or [])
    if (
        energy is None
        and not converged
        and not imag_freqs
        and not inspected.get("frequency_count")
    ):
        return {}
    return {
        "energy": energy,
        "converged": converged,
        "imag_freqs": imag_freqs,
        "frequency_count": inspected.get("frequency_count"),
        "optimized_geometry_count": (
            1 if inspected.get("optimization_converged") else 0
        ),
    }


__all__ = ["extract_optimized_geometry", "run_local", "validate_runtime"]
