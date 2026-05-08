from __future__ import annotations

import contextlib
import importlib
import os
import re
import shutil
import traceback
from types import SimpleNamespace
from typing import Any, Literal

import yaml

from chemsmart.agent.transport import LocalDryRunTransport, SubmitTransport
from chemsmart.cli.sub import sub as sub_cli
from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.orca.output import ORCAOutput
from chemsmart.jobs.gaussian.irc import GaussianIRCJob
from chemsmart.jobs.gaussian.job import GaussianGeneralJob, GaussianJob
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.scan import GaussianScanJob
from chemsmart.jobs.gaussian.settings import GaussianJobSettings
from chemsmart.jobs.gaussian.singlepoint import GaussianSinglePointJob
from chemsmart.jobs.gaussian.ts import GaussianTSJob
from chemsmart.jobs.gaussian.writer import GaussianInputWriter
from chemsmart.jobs.job import Job
from chemsmart.jobs.orca.irc import ORCAIRCJob
from chemsmart.jobs.orca.job import ORCAGeneralJob, ORCAJob
from chemsmart.jobs.orca.opt import ORCAOptJob
from chemsmart.jobs.orca.scan import ORCAScanJob
from chemsmart.jobs.orca.settings import ORCAJobSettings
from chemsmart.jobs.orca.singlepoint import ORCASinglePointJob
from chemsmart.jobs.orca.ts import ORCATSJob
from chemsmart.jobs.orca.writer import ORCAInputWriter
from chemsmart.jobs.runner import JobRunner
from chemsmart.settings.server import Server
from chemsmart.settings.user import ChemsmartUserSettings
from chemsmart.utils.cli import CtxObjArguments
from chemsmart.utils.periodictable import PeriodicTable

_TASK_PROJECT_MAP = {
    "opt": ["dft_default", "organics"],
    "opt+freq": ["dft_default", "organics"],
    "ts": ["tspaths"],
}
_TASK_JOBTYPE_MAP = {
    "opt": "opt",
    "opt+freq": "opt",
    "ts": "ts",
}
_PERIODIC_TABLE = PeriodicTable()


def _resolve_optional_job_class(
    module_path: str,
    class_name: str,
    fallback,
):
    try:
        module = importlib.import_module(module_path)
        return getattr(module, class_name)
    except (ImportError, AttributeError):
        return fallback


_GAUSSIAN_FREQ_JOB_CLASS = _resolve_optional_job_class(
    "chemsmart.jobs.gaussian.freq",
    "GaussianFreqJob",
    GaussianGeneralJob,
)
_ORCA_FREQ_JOB_CLASS = _resolve_optional_job_class(
    "chemsmart.jobs.orca.freq",
    "ORCAFreqJob",
    ORCAGeneralJob,
)
JobKind = Literal[
    "gaussian.opt",
    "gaussian.ts",
    "gaussian.freq",
    "gaussian.sp",
    "gaussian.irc",
    "gaussian.scan",
    "orca.opt",
    "orca.ts",
    "orca.freq",
    "orca.sp",
    "orca.irc",
    "orca.scan",
]
_CANONICAL_JOB_KINDS: tuple[JobKind, ...] = (
    "gaussian.opt",
    "gaussian.ts",
    "gaussian.freq",
    "gaussian.sp",
    "gaussian.irc",
    "gaussian.scan",
    "orca.opt",
    "orca.ts",
    "orca.freq",
    "orca.sp",
    "orca.irc",
    "orca.scan",
)
_JOBTYPE_BY_KIND_SUFFIX = {
    "opt": "opt",
    "ts": "ts",
    "freq": "freq",
    "sp": "sp",
    "irc": "irc",
    "scan": "scan",
}
_JOB_CLASS_BY_KIND = {
    "gaussian.opt": GaussianOptJob,
    "gaussian.ts": GaussianTSJob,
    "gaussian.freq": _GAUSSIAN_FREQ_JOB_CLASS,
    "gaussian.sp": GaussianSinglePointJob,
    "gaussian.irc": GaussianIRCJob,
    "gaussian.scan": GaussianScanJob,
    "orca.opt": ORCAOptJob,
    "orca.ts": ORCATSJob,
    "orca.freq": _ORCA_FREQ_JOB_CLASS,
    "orca.sp": ORCASinglePointJob,
    "orca.irc": ORCAIRCJob,
    "orca.scan": ORCAScanJob,
}
_JOB_KIND_ALIASES = {
    "gaussian.singlepoint": "gaussian.sp",
    "orca.singlepoint": "orca.sp",
}
_SUPPORTED_JOB_KINDS = _CANONICAL_JOB_KINDS
_SUPPORTED_SUBMIT_JOBTYPES = {"opt", "ts", "sp", "irc", "scan"}
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
_ORCA_AB_INITIO_KEYWORDS = (
    "MP2",
    "MP3",
    "MP4",
    "CCSD",
    "DLPNO",
    "CASSCF",
    "NEVPT2",
    "MRCI",
)
_ORCA_AB_INITIO_EXACT_METHODS = {"HF", "RHF", "UHF", "ROHF"}


def build_molecule(filepath: str, index: str = "-1") -> Molecule:
    """Load one molecule from a structure file using chemsmart parsing."""
    return Molecule.from_filepath(
        filepath=filepath,
        index=index,
        return_list=False,
    )


def build_gaussian_settings(
    functional,
    basis,
    charge=0,
    multiplicity=1,
    solvent_model=None,
    solvent_id=None,
    heavy_elements=None,
    heavy_elements_basis=None,
    title=None,
    freq=False,
    numfreq=False,
    additional_opt_options_in_route=None,
    additional_route_parameters=None,
    scan_definition: str | None = None,
    **extras,
) -> GaussianJobSettings:
    """Build Gaussian settings.

    `scan_definition` is the ModRedundant coordinate scan specification for
    Gaussian scan jobs, for example ``D 1 2 3 4 S 10 36.0`` for a 10-step
    dihedral scan in 36° increments. It is required for ``gaussian.scan``
    workflows and is translated into ``settings.modred``.
    """
    modred = extras.pop("modred", None)
    if scan_definition is not None:
        if modred is not None:
            raise ValueError(
                "Provide either scan_definition or modred, not both."
            )
        modred = _parse_gaussian_scan_definition(scan_definition)

    return GaussianJobSettings(
        functional=functional,
        basis=basis,
        charge=charge,
        multiplicity=multiplicity,
        solvent_model=solvent_model,
        solvent_id=solvent_id,
        heavy_elements=heavy_elements,
        heavy_elements_basis=heavy_elements_basis,
        title=title,
        freq=freq,
        numfreq=numfreq,
        additional_opt_options_in_route=additional_opt_options_in_route,
        additional_route_parameters=additional_route_parameters,
        modred=modred,
        **extras,
    )


def build_orca_settings(
    functional,
    basis,
    charge=0,
    multiplicity=1,
    solvent_model=None,
    solvent_id=None,
    heavy_elements=None,
    heavy_elements_basis=None,
    ab_initio=None,
    dispersion=None,
    aux_basis=None,
    extrapolation_basis=None,
    defgrid=None,
    scf_tol=None,
    scf_algorithm=None,
    scf_maxiter=None,
    scf_convergence=None,
    gbw=True,
    freq=False,
    numfreq=False,
    dipole=False,
    quadrupole=False,
    mdci_cutoff=None,
    mdci_density=None,
    jobtype=None,
    title=None,
    additional_solvent_options=None,
    solventfilename=None,
    additional_route_parameters=None,
    route_to_be_written=None,
    modred=None,
    gen_genecp_file=None,
    light_elements_basis=None,
    custom_solvent=None,
    forces=False,
    input_string=None,
    invert_constraints=False,
    **extras,
) -> ORCAJobSettings:
    """Build validated ORCA job settings from planner-supplied fields."""
    if (
        ab_initio is None
        and isinstance(functional, str)
        and _looks_like_orca_ab_initio_method(functional)
    ):
        ab_initio = functional
        functional = None

    return ORCAJobSettings(
        ab_initio=ab_initio,
        functional=functional,
        dispersion=dispersion,
        basis=basis,
        aux_basis=aux_basis,
        extrapolation_basis=extrapolation_basis,
        defgrid=defgrid,
        scf_tol=scf_tol,
        scf_algorithm=scf_algorithm,
        scf_maxiter=scf_maxiter,
        scf_convergence=scf_convergence,
        charge=charge,
        multiplicity=multiplicity,
        gbw=gbw,
        freq=freq,
        numfreq=numfreq,
        dipole=dipole,
        quadrupole=quadrupole,
        mdci_cutoff=mdci_cutoff,
        mdci_density=mdci_density,
        jobtype=jobtype,
        title=title,
        solvent_model=solvent_model,
        solvent_id=solvent_id,
        additional_solvent_options=additional_solvent_options,
        solventfilename=solventfilename,
        additional_route_parameters=additional_route_parameters,
        route_to_be_written=route_to_be_written,
        modred=modred,
        gen_genecp_file=gen_genecp_file,
        heavy_elements=heavy_elements,
        heavy_elements_basis=heavy_elements_basis,
        light_elements_basis=light_elements_basis,
        custom_solvent=custom_solvent,
        forces=forces,
        input_string=input_string,
        invert_constraints=invert_constraints,
        **extras,
    )


def _looks_like_orca_ab_initio_method(method: Any) -> bool:
    if not isinstance(method, str):
        return False

    normalized = method.strip().upper()
    if not normalized:
        return False
    if normalized in _ORCA_AB_INITIO_EXACT_METHODS:
        return True
    return any(keyword in normalized for keyword in _ORCA_AB_INITIO_KEYWORDS)


def _parse_gaussian_scan_definition(scan_definition: str) -> dict[str, Any]:
    lines = [
        line.strip()
        for chunk in scan_definition.splitlines()
        for line in chunk.split(";")
        if line.strip()
    ]
    if not lines:
        raise ValueError("scan_definition cannot be empty.")

    atom_count_by_coordinate = {"B": 2, "A": 3, "D": 4}
    coords: list[list[int]] = []
    num_steps: list[int] = []
    step_sizes: list[float] = []
    constrained_coordinates: list[list[int]] = []

    for line in lines:
        tokens = line.split()
        coordinate_type = tokens[0].upper()
        expected_atoms = atom_count_by_coordinate.get(coordinate_type)
        if expected_atoms is None:
            raise ValueError(
                "scan_definition must start with B, A, or D followed by "
                "1-based atom indices."
            )
        if len(tokens) < expected_atoms + 2:
            raise ValueError(
                f"scan_definition {line!r} is incomplete. Expected a "
                "coordinate, indices, and a scan directive."
            )

        try:
            coordinate = [
                int(token) for token in tokens[1 : 1 + expected_atoms]
            ]
        except ValueError as exc:
            raise ValueError(
                f"scan_definition {line!r} has non-integer atom indices."
            ) from exc

        directive = tokens[1 + expected_atoms].upper()
        if directive == "F":
            constrained_coordinates.append(coordinate)
            continue

        if directive != "S" or len(tokens) != expected_atoms + 4:
            raise ValueError(
                "scan_definition entries must look like "
                "'D 1 2 3 4 S 10 36.0' or 'B 1 2 F'."
            )

        try:
            steps = int(tokens[2 + expected_atoms])
            step_size = float(tokens[3 + expected_atoms])
        except ValueError as exc:
            raise ValueError(
                f"scan_definition {line!r} has an invalid step count "
                "or step size."
            ) from exc

        coords.append(coordinate)
        num_steps.append(steps)
        step_sizes.append(step_size)

    if not coords:
        raise ValueError(
            "scan_definition must include at least one scan entry with "
            "an 'S <steps> <step_size>' directive."
        )

    modred: dict[str, Any] = {
        "coords": coords,
        "num_steps": num_steps,
        "step_size": step_sizes,
    }
    if constrained_coordinates:
        modred["constrained_coordinates"] = constrained_coordinates
    return modred


def build_job(
    kind: JobKind,
    molecule: Molecule,
    settings,
    label: str | None = None,
    jobrunner=None,
) -> Job:
    """Instantiate a chemsmart job object for a canonical agent job kind."""
    normalized_kind = (kind or "").strip().lower()
    normalized_kind = _JOB_KIND_ALIASES.get(normalized_kind, normalized_kind)
    job_class = _JOB_CLASS_BY_KIND.get(normalized_kind)
    if job_class is None:
        supported_kinds = ", ".join(_SUPPORTED_JOB_KINDS)
        raise ValueError(
            f"Unknown job kind {kind!r}. Supported kinds: {supported_kinds}"
        )

    job_settings = settings.copy()
    expected_settings_class = job_class.settings_class()
    if not isinstance(job_settings, expected_settings_class):
        job_settings = expected_settings_class(**job_settings.__dict__)
    kind_suffix = normalized_kind.split(".", maxsplit=1)[1]
    job_settings.jobtype = _JOBTYPE_BY_KIND_SUFFIX[kind_suffix]
    if (
        kind_suffix == "freq"
        and not job_settings.freq
        and not job_settings.numfreq
    ):
        job_settings.freq = True
    if normalized_kind == "gaussian.scan" and not job_settings.modred:
        raise ValueError(
            "gaussian.scan requires scan_definition in "
            "build_gaussian_settings, for example "
            "'D 1 2 3 4 S 10 36.0' or 'B 1 2 S 10 0.05'."
        )

    return job_class(
        molecule=molecule,
        settings=job_settings,
        label=label,
        jobrunner=jobrunner,
    )


def dry_run_input(job: Job) -> dict[str, str]:
    """Render a job input file and return its absolute path and contents."""
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
    return {"inputfile": inputfile, "content": content}


def validate_runtime(
    job: Job,
    server=None,
) -> dict[str, Any]:
    """Check local/runtime prerequisites and remote unknowns for a job."""
    local_issues = _validate_job_fields(job)
    remote_unknown = []

    if server is None:
        remote_unknown.extend(_REMOTE_UNKNOWN_SERVER_FIELDS)
        return _runtime_validation_result(local_issues, remote_unknown)

    server_config = getattr(server, "config", None)
    if not isinstance(server_config, dict):
        local_issues.append("server.config invalid")
        return _runtime_validation_result(local_issues, remote_unknown)

    queue_name = _get_server_queue_name(server_config)
    account_name = _get_server_account_name(server_config)
    scratch_value = server_config.get("SCRATCH_DIR")

    if not _is_non_empty_string(queue_name):
        local_issues.append("server.queue missing")
    if not _is_non_empty_string(account_name):
        local_issues.append("server.account missing")
    if scratch_value is None:
        local_issues.append("server.scratch_dir missing")
    elif not isinstance(scratch_value, str):
        local_issues.append("server.scratch_dir invalid")

    runner = None
    try:
        runner = JobRunner.from_job(job=job, server=server, scratch=False)
    except Exception:
        local_issues.append("jobrunner unavailable")

    executable_path = None
    modules = None
    resolved_scratch_dir = None

    if runner is not None:
        try:
            resolved_scratch_dir = runner._resolve_scratch_dir_candidate()
            if resolved_scratch_dir is not None:
                resolved_scratch_dir = runner._normalize_path_string(
                    resolved_scratch_dir
                )
        except Exception:
            local_issues.append("server.scratch_dir invalid")

        try:
            executable = runner.executable
            executable_path = executable.get_executable()
            modules = executable.modules
        except Exception:
            executable = None
    else:
        executable = None

    if resolved_scratch_dir is not None and _has_unresolved_envvars(
        resolved_scratch_dir
    ):
        local_issues.append("server.scratch_dir unresolved")

    if modules is not None and not isinstance(modules, str):
        local_issues.append("server.modules invalid")
    if executable_path is not None and not isinstance(executable_path, str):
        local_issues.append("server.executable_path invalid")

    has_modules = _is_non_empty_string(modules)
    has_executable_path = _is_non_empty_string(executable_path)
    if not has_modules and not has_executable_path:
        local_issues.append("server.modules_or_executable_path missing")

    if has_executable_path and shutil.which(executable_path) is None:
        local_issues.append("server.executable_path missing")

    if _is_non_empty_string(resolved_scratch_dir):
        remote_unknown.append(_REMOTE_UNKNOWN_SCRATCH)
    if has_modules:
        remote_unknown.append(_REMOTE_UNKNOWN_MODULES)
    if _is_non_empty_string(queue_name):
        remote_unknown.append(_REMOTE_UNKNOWN_QUEUE)
    if _is_non_empty_string(account_name):
        remote_unknown.append(_REMOTE_UNKNOWN_ACCOUNT)
    remote_unknown.append(_REMOTE_UNKNOWN_SSH)

    return _runtime_validation_result(local_issues, remote_unknown)


def run_local(job: Job) -> dict[str, Any]:
    """Execute a job locally and summarize the generated output artifacts."""
    job_folder = os.path.abspath(job.folder)
    os.makedirs(job_folder, exist_ok=True)
    job.set_folder(job_folder)

    stdout_path = os.path.join(job_folder, f"{job.label}.stdout")
    stderr_path = os.path.join(job_folder, f"{job.label}.stderr")

    job.local = True
    returncode = 0

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
            returncode = int(getattr(exc, "returncode", 1) or 1)

    output_summary = {}
    if returncode == 0:
        output_summary = _summarize_local_output(job)

    return {
        "ok": returncode == 0,
        "returncode": returncode,
        "stdout_path": stdout_path,
        "stderr_path": stderr_path,
        "output_summary": output_summary,
    }


def extract_optimized_geometry(job: Job) -> Molecule:
    """Extract the final optimized geometry from a completed job log."""
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
            "extract_optimized_geometry only supports GaussianJob and "
            "ORCAJob instances"
        )

    return Molecule(
        symbols=symbols,
        positions=positions,
        charge=getattr(job.settings, "charge", None),
        multiplicity=getattr(job.settings, "multiplicity", None),
    )


def submit_hpc(
    job: Job,
    server=None,
    transport: SubmitTransport | None = None,
    execute: bool = False,
) -> dict[str, Any]:
    """Generate and optionally submit an HPC script for a prepared job."""
    selected_transport = _select_submit_transport(
        transport=transport,
        execute=execute,
    )
    duplicate_check = _check_duplicate_submission(job)
    if duplicate_check["duplicate"]:
        return {
            "transport": selected_transport.__class__.__name__,
            "script_path": None,
            "script_bytes": None,
            "command_executed": None,
            "job_id": None,
            "duplicate_check": duplicate_check,
        }

    server_obj = _coerce_server(server)
    job_folder = os.path.abspath(job.folder)
    os.makedirs(job_folder, exist_ok=True)
    job.set_folder(job_folder)

    dry_run_input(job)
    cli_args = _reconstruct_submit_cli_args(job, server_obj)

    with _pushd(job_folder):
        server_obj.submit(job=job, test=True, cli_args=cli_args)

    submitter = server_obj.get_submitter(job)
    script_path = os.path.abspath(
        os.path.join(job_folder, submitter.submit_script)
    )
    with open(script_path, "rb") as file:
        script_bytes = file.read()

    submit_result = selected_transport.submit(
        script_path=script_path,
        working_dir=job_folder,
        server=server_obj,
    )
    return {
        "transport": selected_transport.__class__.__name__,
        "script_path": script_path,
        "script_bytes": script_bytes,
        "command_executed": submit_result["command_executed"],
        "job_id": submit_result["job_id"],
        "duplicate_check": duplicate_check,
    }


def recommend_method(
    task: str,
    charge: int = 0,
    multiplicity: int = 1,
    atomic_numbers: list[int] | None = None,
    project_hint: str | None = None,
) -> dict[str, Any]:
    """Return a conservative project-based method recommendation or no-match."""
    user_settings = ChemsmartUserSettings()
    available_project_paths = _get_available_project_paths(user_settings)
    available_projects = sorted(available_project_paths)
    normalized_task = (task or "").strip().lower()
    heavy_symbols = _get_heavy_symbols(atomic_numbers)

    if project_hint in available_project_paths:
        return _matched_recommendation(
            project_name=project_hint,
            project_path=available_project_paths[project_hint],
            task=normalized_task,
            rationale=("matched project_hint rule: " f"{project_hint}.yaml"),
            available_projects=available_projects,
        )

    if charge != 0:
        return _no_match_response(
            available_projects=available_projects,
            reason=f"charge={charge}",
        )

    if multiplicity > 1:
        return _no_match_response(
            available_projects=available_projects,
            reason=f"multiplicity={multiplicity}",
        )

    heavy_element_projects = _get_heavy_element_projects(
        available_project_paths=available_project_paths,
        task=normalized_task,
        heavy_symbols=heavy_symbols,
    )
    if heavy_symbols and not heavy_element_projects:
        return _no_match_response(
            available_projects=available_projects,
            reason=(
                "no project lists heavy_elements for "
                f"{', '.join(heavy_symbols)}"
            ),
        )

    task_candidates = [
        project_name
        for project_name in _TASK_PROJECT_MAP.get(normalized_task, [])
        if project_name in available_project_paths
    ]
    if heavy_symbols:
        task_candidates = [
            project_name
            for project_name in task_candidates
            if project_name in heavy_element_projects
        ]

    if task_candidates:
        project_name = task_candidates[0]
        return _matched_recommendation(
            project_name=project_name,
            project_path=available_project_paths[project_name],
            task=normalized_task,
            rationale=(
                f"matched task rule for '{normalized_task}': "
                f"{project_name}.yaml"
            ),
            available_projects=available_projects,
        )

    if heavy_element_projects:
        project_name = heavy_element_projects[0]
        return _matched_recommendation(
            project_name=project_name,
            project_path=available_project_paths[project_name],
            task=normalized_task,
            rationale=(
                "matched heavy_elements rule for "
                f"{', '.join(heavy_symbols)}: {project_name}.yaml"
            ),
            available_projects=available_projects,
        )

    return _no_match_response(
        available_projects=available_projects,
        reason=f"no task rule matched for '{normalized_task}'",
    )


def _get_available_project_paths(
    user_settings: ChemsmartUserSettings,
) -> dict[str, str]:
    project_paths = {}
    for filepath in user_settings.gaussian_project_yaml_files:
        project_name = os.path.basename(filepath).removesuffix(".yaml")
        if project_name == "defaults":
            continue
        project_paths[project_name] = filepath
    return project_paths


def _get_heavy_symbols(
    atomic_numbers: list[int] | None,
) -> list[str]:
    if not atomic_numbers:
        return []

    heavy_symbols = {
        _PERIODIC_TABLE.to_symbol(atomic_number)
        for atomic_number in atomic_numbers
        if atomic_number >= 54
    }
    return _PERIODIC_TABLE.sorted_periodic_table_list(list(heavy_symbols))


def _get_heavy_element_projects(
    available_project_paths: dict[str, str],
    task: str,
    heavy_symbols: list[str],
) -> list[str]:
    if not heavy_symbols:
        return []

    matching_projects = []
    for project_name, project_path in available_project_paths.items():
        settings = _get_project_settings(project_path, task)
        project_heavy_elements = set(
            _normalize_heavy_elements(settings.get("heavy_elements")) or []
        )
        if set(heavy_symbols).issubset(project_heavy_elements):
            matching_projects.append(project_name)

    return sorted(matching_projects)


def _matched_recommendation(
    project_name: str,
    project_path: str,
    task: str,
    rationale: str,
    available_projects: list[str],
) -> dict[str, Any]:
    settings = _get_project_settings(project_path, task)
    return {
        "match": project_name,
        "functional": settings.get("functional"),
        "basis": settings.get("basis"),
        "solvent_model": settings.get("solvent_model"),
        "solvent_id": settings.get("solvent_id"),
        "heavy_elements": _normalize_heavy_elements(
            settings.get("heavy_elements")
        ),
        "heavy_elements_basis": settings.get("heavy_elements_basis"),
        "rationale": rationale,
        "available_projects": available_projects,
    }


def _no_match_response(
    available_projects: list[str],
    reason: str,
) -> dict[str, Any]:
    return {
        "match": None,
        "functional": None,
        "basis": None,
        "solvent_model": None,
        "solvent_id": None,
        "heavy_elements": None,
        "heavy_elements_basis": None,
        "rationale": (
            f"no rule matched: {reason}; pick from available_projects"
        ),
        "available_projects": available_projects,
    }


def _get_project_settings(project_path: str, task: str) -> dict[str, Any]:
    project_settings = _load_yaml(project_path)
    gas_settings = project_settings.get("gas")
    solv_settings = project_settings.get("solv")
    jobtype = _TASK_JOBTYPE_MAP.get(task, "opt")

    if jobtype == "sp":
        settings = solv_settings or gas_settings or project_settings
    else:
        settings = gas_settings or solv_settings or project_settings

    if not isinstance(settings, dict):
        return {}
    return settings


def _normalize_heavy_elements(
    heavy_elements: list[str] | str | None,
) -> list[str] | None:
    if heavy_elements is None:
        return None
    if isinstance(heavy_elements, str):
        heavy_elements = heavy_elements.replace(",", " ").split()
    return [_PERIODIC_TABLE.to_element(element) for element in heavy_elements]


def _load_yaml(filepath: str) -> dict[str, Any]:
    with open(filepath) as file:
        return yaml.safe_load(file) or {}


def _validate_job_fields(job: Job) -> list[str]:
    issues = []
    if getattr(job, "molecule", None) is None:
        issues.append("job.molecule missing")
    if getattr(job, "settings", None) is None:
        issues.append("job.settings missing")

    label = getattr(job, "label", None)
    if not _is_non_empty_string(label):
        issues.append("job.label missing")
    return issues


def _runtime_validation_result(
    local_issues: list[str],
    remote_unknown: list[str],
) -> dict[str, Any]:
    if local_issues:
        status = "fail"
    elif remote_unknown:
        status = "partial"
    else:
        status = "ok"

    return {
        "ok": status,
        "local_ok": not local_issues,
        "local_issues": local_issues,
        "remote_unknown": remote_unknown,
    }


def _get_server_queue_name(server_config: dict[str, Any]) -> Any:
    return server_config.get("QUEUE_NAME")


def _get_server_account_name(server_config: dict[str, Any]) -> Any:
    return server_config.get("ACCOUNT", server_config.get("PROJECT"))


def _is_non_empty_string(value: Any) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _resolve_geometry_logfile(job: Job) -> str:
    outputfile = getattr(job, "outputfile", None)
    if isinstance(outputfile, str) and os.path.exists(outputfile):
        return outputfile

    fallback_paths = [
        os.path.join(job.folder, f"{job.label}.log"),
        os.path.join(job.folder, f"{job.label}.out"),
    ]
    for path in fallback_paths:
        if os.path.exists(path):
            return path

    if isinstance(outputfile, str):
        return outputfile
    return fallback_paths[0]


def _parse_gaussian_geometry(
    lines: list[str],
) -> tuple[list[str], list[list[float]]]:
    blocks: list[tuple[list[str], list[list[float]]]] = []

    for index, line in enumerate(lines):
        if "Standard orientation:" not in line:
            continue

        symbols = []
        positions = []
        coord_start = index + 5
        for coord_line in lines[coord_start:]:
            stripped = coord_line.strip()
            if not stripped or set(stripped) == {"-"}:
                break
            parts = stripped.split()
            if len(parts) < 6:
                break
            atomic_number = int(parts[1])
            symbols.append(_PERIODIC_TABLE.to_symbol(atomic_number))
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
        coord_start = index + 2
        for coord_line in lines[coord_start:]:
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
    parser_cls = None
    output_path = None

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
        all_structures = list(getattr(parser, "all_structures", []) or [])
        optimized_geometry_count = len(all_structures)

        if (
            energy is None
            and not converged
            and not imag_freqs
            and optimized_geometry_count == 0
        ):
            return {}

        return {
            "energy": energy,
            "converged": converged,
            "imag_freqs": imag_freqs,
            "optimized_geometry_count": optimized_geometry_count,
        }
    except Exception:
        return {}


def _select_submit_transport(
    transport: SubmitTransport | None,
    execute: bool,
) -> SubmitTransport:
    if not execute:
        return LocalDryRunTransport()
    if transport is None:
        return LocalDryRunTransport()
    return transport


def _check_duplicate_submission(job: Job) -> dict[str, Any]:
    try:
        Server._check_running_jobs(job)
    except SystemExit as exc:
        return {
            "duplicate": True,
            "message": str(exc),
        }

    return {
        "duplicate": False,
        "message": None,
    }


def _coerce_server(server) -> Server:
    if server is None:
        server = _default_submit_server_name()
    if isinstance(server, Server):
        return server
    if isinstance(server, str):
        return Server.from_servername(server)
    raise TypeError("server must be a chemsmart.settings.server.Server or str")


def _default_submit_server_name() -> str:
    user_settings = ChemsmartUserSettings()
    available_servers = sorted(user_settings.all_available_servers)
    if len(available_servers) == 1:
        return available_servers[0]
    if not available_servers:
        raise ValueError(
            "submit_hpc requires server when no configured servers are available"
        )
    raise ValueError(
        "submit_hpc requires server when multiple configured servers are "
        f"available: {', '.join(available_servers)}"
    )


def _reconstruct_submit_cli_args(job: Job, server: Server) -> list[str]:
    program_command = _program_submit_command(job)
    jobtype_command = _jobtype_submit_command(job, program_command)
    subcommands = [
        _click_command_to_ctx_obj(
            sub_cli,
            {
                "server": _server_cli_name(server),
                "num_cores": (
                    getattr(job.jobrunner, "num_cores", None)
                    if getattr(job, "jobrunner", None) is not None
                    else None
                ),
                "num_gpus": (
                    getattr(job.jobrunner, "num_gpus", None)
                    if getattr(job, "jobrunner", None) is not None
                    else None
                ),
                "mem_gb": (
                    getattr(job.jobrunner, "mem_gb", None)
                    if getattr(job, "jobrunner", None) is not None
                    else None
                ),
                "time_hours": None,
                "queue": None,
                "verbose": None,
                "test": None,
                "print_command": None,
            },
            parent=None,
        ),
        _click_command_to_ctx_obj(
            program_command,
            {
                "filename": os.path.abspath(job.inputfile),
                "label": job.label,
            },
            parent="sub",
        ),
        _click_command_to_ctx_obj(
            jobtype_command,
            {},
            parent=program_command.name,
        ),
    ]
    return CtxObjArguments(
        subcommands,
        entry_point="sub",
    ).reconstruct_command_line()[1:]


def _program_submit_command(job: Job):
    program_name = getattr(job, "PROGRAM", None)
    if not isinstance(program_name, str):
        raise ValueError(f"Job {job!r} does not define PROGRAM.")
    try:
        return sub_cli.commands[program_name.lower()]
    except KeyError as exc:
        raise ValueError(
            f"Unsupported submit program for {job!r}: {program_name!r}"
        ) from exc


def _jobtype_submit_command(job: Job, program_command):
    jobtype = _jobtype_name(job)
    try:
        return program_command.commands[jobtype]
    except KeyError as exc:
        raise ValueError(
            f"Unsupported submit jobtype {jobtype!r} for {job!r}"
        ) from exc


def _jobtype_name(job: Job) -> str:
    job_settings = getattr(job, "settings", None)
    jobtype = getattr(job_settings, "jobtype", None)
    if not isinstance(jobtype, str):
        raise ValueError(f"Job {job!r} is missing settings.jobtype.")
    jobtype = jobtype.lower()
    if jobtype == "singlepoint":
        return "sp"
    if jobtype not in _SUPPORTED_SUBMIT_JOBTYPES:
        supported = ", ".join(sorted(_SUPPORTED_SUBMIT_JOBTYPES))
        raise ValueError(
            f"Unsupported submit jobtype {jobtype!r}. Supported: {supported}"
        )
    return jobtype


def _click_command_to_ctx_obj(command, overrides: dict[str, Any], parent):
    kwargs = {}
    for param in command.params:
        kwargs[param.name] = {
            "value": overrides.get(param.name, param.default),
            "nargs": param.nargs,
            "is_multiple": param.multiple,
            "type": param.type,
            "is_flag": param.is_flag,
            "secondary_opts": param.secondary_opts,
        }
    return {
        "name": command.name,
        "kwargs": kwargs,
        "parent": parent,
    }


def _server_cli_name(server: Server) -> str | None:
    server_name = getattr(server, "name", None)
    if not isinstance(server_name, str):
        return None
    basename = os.path.basename(server_name)
    if basename.endswith(".yaml"):
        return basename.removesuffix(".yaml")
    return basename


@contextlib.contextmanager
def _pushd(path: str):
    cwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(cwd)
