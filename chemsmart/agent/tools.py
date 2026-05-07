from __future__ import annotations

import importlib
import os
from types import SimpleNamespace
from typing import Any

import yaml

from chemsmart.io.molecules.structure import Molecule
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
from chemsmart.settings.user import ChemsmartUserSettings
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
_JOBTYPE_BY_KIND_SUFFIX = {
    "opt": "opt",
    "ts": "ts",
    "freq": "freq",
    "sp": "sp",
    "singlepoint": "sp",
    "irc": "irc",
    "scan": "scan",
}
_JOB_CLASS_BY_KIND = {
    "gaussian.opt": GaussianOptJob,
    "gaussian.ts": GaussianTSJob,
    "gaussian.freq": _GAUSSIAN_FREQ_JOB_CLASS,
    "gaussian.sp": GaussianSinglePointJob,
    "gaussian.singlepoint": GaussianSinglePointJob,
    "gaussian.irc": GaussianIRCJob,
    "gaussian.scan": GaussianScanJob,
    "orca.opt": ORCAOptJob,
    "orca.ts": ORCATSJob,
    "orca.freq": _ORCA_FREQ_JOB_CLASS,
    "orca.sp": ORCASinglePointJob,
    "orca.singlepoint": ORCASinglePointJob,
    "orca.irc": ORCAIRCJob,
    "orca.scan": ORCAScanJob,
}
_SUPPORTED_JOB_KINDS = tuple(_JOB_CLASS_BY_KIND)


def build_molecule(filepath: str, index: str = "-1") -> Molecule:
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
    **extras,
) -> GaussianJobSettings:
    return GaussianJobSettings(
        functional=functional,
        basis=basis,
        charge=charge,
        multiplicity=multiplicity,
        solvent_model=solvent_model,
        solvent_id=solvent_id,
        heavy_elements=heavy_elements,
        heavy_elements_basis=heavy_elements_basis,
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


def build_job(
    kind: str,
    molecule: Molecule,
    settings,
    label: str | None = None,
    jobrunner=None,
) -> Job:
    normalized_kind = (kind or "").strip().lower()
    job_class = _JOB_CLASS_BY_KIND.get(normalized_kind)
    if job_class is None:
        supported_kinds = ", ".join(_SUPPORTED_JOB_KINDS)
        raise ValueError(
            f"Unknown job kind {kind!r}. Supported kinds: {supported_kinds}"
        )

    job_settings = settings.copy()
    kind_suffix = normalized_kind.split(".", maxsplit=1)[1]
    job_settings.jobtype = _JOBTYPE_BY_KIND_SUFFIX[kind_suffix]
    if (
        kind_suffix == "freq"
        and not job_settings.freq
        and not job_settings.numfreq
    ):
        job_settings.freq = True

    return job_class(
        molecule=molecule,
        settings=job_settings,
        label=label,
        jobrunner=jobrunner,
    )


def dry_run_input(job: Job) -> dict[str, str]:
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


def recommend_method(
    task: str,
    charge: int = 0,
    multiplicity: int = 1,
    atomic_numbers: list[int] | None = None,
    project_hint: str | None = None,
) -> dict[str, Any]:
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
