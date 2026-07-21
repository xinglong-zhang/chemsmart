from __future__ import annotations

import contextlib
import importlib
import os
from typing import Any, Literal

from chemsmart.agent.harness.workflow_state import current_workflow_state
from chemsmart.agent.services.job_cli import (
    SUBMIT_JOBTYPES_BY_PROGRAM,
)
from chemsmart.agent.services.job_cli import (
    SUPPORTED_SUBMIT_JOBTYPES as _SUPPORTED_SUBMIT_JOBTYPES,
)
from chemsmart.agent.services.job_cli import dry_run_input as _dry_run_input
from chemsmart.agent.services.local_runtime import (
    extract_optimized_geometry as _extract_optimized_geometry,
)
from chemsmart.agent.services.local_runtime import run_local as _run_local
from chemsmart.agent.services.local_runtime import (
    validate_runtime as _validate_runtime,
)
from chemsmart.agent.services.method_recommendation import (
    recommend_method as _recommend_method,
)
from chemsmart.agent.services.scan_directives import (
    parse_gaussian_scan_definition as _parse_gaussian_scan_definition,
)
from chemsmart.agent.services.scan_directives import (
    partition_modred_route_directives as _partition_modred_route_directives,
)
from chemsmart.agent.services.server_selection import (
    coerce_server as _coerce_server,
)
from chemsmart.agent.transport import LocalDryRunTransport, SubmitTransport
from chemsmart.cli.sub import sub as sub_cli
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.irc import GaussianIRCJob
from chemsmart.jobs.gaussian.job import GaussianGeneralJob
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.scan import GaussianScanJob
from chemsmart.jobs.gaussian.settings import GaussianJobSettings
from chemsmart.jobs.gaussian.singlepoint import GaussianSinglePointJob
from chemsmart.jobs.gaussian.ts import GaussianTSJob
from chemsmart.jobs.job import Job
from chemsmart.jobs.orca.irc import ORCAIRCJob
from chemsmart.jobs.orca.job import ORCAGeneralJob
from chemsmart.jobs.orca.opt import ORCAOptJob
from chemsmart.jobs.orca.scan import ORCAScanJob
from chemsmart.jobs.orca.settings import ORCAJobSettings
from chemsmart.jobs.orca.singlepoint import ORCASinglePointJob
from chemsmart.jobs.orca.ts import ORCATSJob
from chemsmart.jobs.runner import JobRunner
from chemsmart.jobs.xtb.hess import XTBHessJob
from chemsmart.jobs.xtb.opt import XTBOptJob
from chemsmart.jobs.xtb.settings import XTBJobSettings
from chemsmart.jobs.xtb.singlepoint import XTBSinglePointJob
from chemsmart.settings.server import Server
from chemsmart.utils.cli import CtxObjArguments


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
    "xtb.opt",
    "xtb.sp",
    "xtb.hess",
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
    "xtb.opt",
    "xtb.sp",
    "xtb.hess",
)
_JOBTYPE_BY_KIND_SUFFIX = {
    "opt": "opt",
    "ts": "ts",
    "freq": "freq",
    "sp": "sp",
    "irc": "irc",
    "scan": "scan",
    "hess": "hess",
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
    "xtb.opt": XTBOptJob,
    "xtb.sp": XTBSinglePointJob,
    "xtb.hess": XTBHessJob,
}
_JOB_KIND_ALIASES = {
    "gaussian.singlepoint": "gaussian.sp",
    "orca.singlepoint": "orca.sp",
    "xtb.singlepoint": "xtb.sp",
    "xtb.freq": "xtb.hess",
}
_SUPPORTED_JOB_KINDS = _CANONICAL_JOB_KINDS
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


def dry_run_input(job: Job) -> dict[str, Any]:
    """Render a job input file and return its absolute path and contents."""

    return _dry_run_input(job)


def recommend_method(
    task: str,
    charge: int = 0,
    multiplicity: int = 1,
    atomic_numbers: list[int] | None = None,
    project_hint: str | None = None,
) -> dict[str, Any]:
    """Return a conservative project-based method recommendation or no-match."""

    return _recommend_method(
        task=task,
        charge=charge,
        multiplicity=multiplicity,
        atomic_numbers=atomic_numbers,
        project_hint=project_hint,
    )


def validate_runtime(job: Job, server=None) -> dict[str, Any]:
    """Check local/runtime prerequisites and remote unknowns for a job."""

    return _validate_runtime(job, server)


def run_local(job: Job) -> dict[str, Any]:
    """Execute a job locally and summarize the generated output artifacts."""

    return _run_local(job)


def extract_optimized_geometry(job: Job) -> Molecule:
    """Extract the final optimized geometry from a completed job log."""

    return _extract_optimized_geometry(job)


def build_molecule(filepath: str, index: str = "-1") -> Molecule:
    """Load one molecule from a structure file using chemsmart parsing."""
    molecule = Molecule.from_filepath(
        filepath=filepath,
        index=index,
        return_list=False,
    )
    # Molecule objects do not preserve their source path. The agent harness
    # needs it to keep every generated input grounded in an equivalent
    # `chemsmart run/sub ...` command.
    setattr(molecule, "_agent_source_filepath", filepath)
    setattr(molecule, "_agent_source_index", index)
    return molecule


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
    if scan_definition is not None and additional_opt_options_in_route:
        coordinate_directives, route_options = (
            _partition_modred_route_directives(additional_opt_options_in_route)
        )
        if coordinate_directives:
            scan_definition = "\n".join(
                (scan_definition, *coordinate_directives)
            )
            additional_opt_options_in_route = (
                ",".join(route_options) if route_options else None
            )
    if scan_definition is not None:
        if modred is not None:
            raise ValueError(
                "Provide either scan_definition or modred, not both."
            )
        modred = _parse_gaussian_scan_definition(scan_definition)

    # `ts`, `calcfc`, `noeigentest` are the transition-state route options the
    # runtime auto-derives for a TS job (see jobs/gaussian/settings.py:627-631,
    # which always renders `opt=(ts,calcfc,noeigentest,<extras>)`). They are
    # never valid *manual* extras — supplying one duplicates the opt keyword,
    # e.g. a caller passing `ts` yields `opt=(ts,calcfc,noeigentest,ts)`, which
    # the runtime harness rejects (`gaussian.ts.route`). Strip them defensively
    # so a redundant token from any caller does not corrupt the route. Genuine
    # extras such as `maxstep=8` or `calcall` are preserved.
    if additional_opt_options_in_route:
        _runtime_ts_opts = {"ts", "calcfc", "noeigentest"}
        _raw = additional_opt_options_in_route
        if isinstance(_raw, (list, tuple)):
            _raw = " ".join(str(x) for x in _raw)
        _kept = [
            tok
            for tok in str(_raw).replace(",", " ").split()
            if tok and tok.lower() not in _runtime_ts_opts
        ]
        additional_opt_options_in_route = ",".join(_kept) if _kept else None

    settings = GaussianJobSettings(
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
    _attach_selected_project(settings, "gaussian")
    return settings


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

    settings = ORCAJobSettings(
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
    _attach_selected_project(settings, "orca")
    return settings


def build_xtb_settings(
    gfn_version="gfn2",
    charge=0,
    multiplicity=1,
    optimization_level=None,
    solvent_model=None,
    solvent_id=None,
    grad=False,
    jobtype=None,
    title=None,
    **extras,
) -> XTBJobSettings:
    """Build validated xTB job settings from planner-supplied fields.

    xTB has no functional or basis set: the GFN Hamiltonian is the method,
    and implicit solvation needs both a model (alpb, gbsa) and a solvent
    identifier. Each field maps one-to-one onto a ``chemsmart run xtb``
    flag so the grounded command stays equivalent to the built job.
    """

    if solvent_model is not None and solvent_id is None:
        raise ValueError(
            "xTB solvation requires solvent_id alongside solvent_model; "
            "a model alone runs in the gas phase."
        )
    if solvent_id is not None and solvent_model is None:
        raise ValueError(
            "xTB solvation requires solvent_model alongside solvent_id, "
            "for example solvent_model='alpb'."
        )
    settings = XTBJobSettings(
        gfn_version=gfn_version,
        charge=charge,
        multiplicity=multiplicity,
        jobtype=jobtype,
        title=title,
        grad=grad,
        solvent_model=solvent_model,
        solvent_id=solvent_id,
        **(
            {}
            if optimization_level is None
            else {"optimization_level": optimization_level}
        ),
        **extras,
    )
    _attach_selected_project(settings, "xtb")
    return settings


def _looks_like_orca_ab_initio_method(method: Any) -> bool:
    if not isinstance(method, str):
        return False

    normalized = method.strip().upper()
    if not normalized:
        return False
    if normalized in _ORCA_AB_INITIO_EXACT_METHODS:
        return True
    return any(keyword in normalized for keyword in _ORCA_AB_INITIO_KEYWORDS)


def _attach_selected_project(settings: Any, program: str) -> None:
    selected = current_workflow_state().project
    if selected is None or selected.program != program:
        return
    setattr(settings, "_agent_project_name", selected.name)
    setattr(settings, "_agent_project_path", selected.path)
    setattr(settings, "_agent_project_sha256", selected.sha256)


def build_job(
    kind: JobKind,
    molecule: Molecule,
    settings,
    label: str | None = None,
    jobrunner=None,
) -> Job:
    """Instantiate a chemsmart job object for a canonical agent job kind."""
    # jobrunner must be a proper runner object; reject strings/dicts passed by mistake

    if jobrunner is not None and not isinstance(jobrunner, JobRunner):
        jobrunner = None

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

    job = job_class(
        molecule=molecule,
        settings=job_settings,
        label=label,
        jobrunner=jobrunner,
    )
    source_filepath = getattr(molecule, "_agent_source_filepath", None)
    source_index = getattr(molecule, "_agent_source_index", None)
    if source_filepath is not None:
        setattr(job, "_agent_source_filepath", source_filepath)
    if source_index is not None:
        setattr(job, "_agent_source_index", source_index)
    setattr(job, "_agent_kind", normalized_kind)
    for attribute in (
        "_agent_project_name",
        "_agent_project_path",
        "_agent_project_sha256",
    ):
        value = getattr(settings, attribute, None)
        if value:
            setattr(job, attribute, value)
    return job


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
    program = str(getattr(job, "PROGRAM", "") or "").lower()
    supported = SUBMIT_JOBTYPES_BY_PROGRAM.get(
        program, _SUPPORTED_SUBMIT_JOBTYPES
    )
    if jobtype not in supported:
        names = ", ".join(sorted(supported))
        raise ValueError(
            f"Unsupported submit jobtype {jobtype!r}. Supported: {names}"
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
