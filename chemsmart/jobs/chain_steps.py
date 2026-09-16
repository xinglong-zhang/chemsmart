"""Chain-step registry and ChainJob builder."""

import shlex
from dataclasses import dataclass

from chemsmart.jobs.chain import ChainJob, JobPhase
from chemsmart.jobs.crest.conformers import CRESTConformerSearchJob
from chemsmart.jobs.crest.job import CRESTJob
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.singlepoint import GaussianSinglePointJob
from chemsmart.jobs.gaussian.ts import GaussianTSJob
from chemsmart.jobs.orca.opt import ORCAOptJob
from chemsmart.jobs.orca.singlepoint import ORCASinglePointJob
from chemsmart.jobs.orca.ts import ORCATSJob
from chemsmart.jobs.xtb.hess import XTBHessJob
from chemsmart.jobs.xtb.opt import XTBOptJob
from chemsmart.jobs.xtb.singlepoint import XTBSinglePointJob
from chemsmart.settings.chain import ChainProjectSettings
from chemsmart.settings.crest import CRESTProjectSettings
from chemsmart.settings.gaussian import GaussianProjectSettings
from chemsmart.settings.orca import ORCAProjectSettings
from chemsmart.settings.xtb import XTBProjectSettings

__all__ = [
    "ChainStep",
    "ChainStepSpec",
    "CHAIN_PROGRAM_SETTINGS",
    "CHAIN_STEP_SPECS",
    "CHAIN_NESTED_ONLY_JOBS",
    "get_chain_step_spec",
    "parse_chain_step",
    "molecule_from_completed_job",
    "build_chain_job",
]

# Workflow jobs that are not ``-s`` pipeline steps.
CHAIN_NESTED_ONLY_JOBS = frozenset({"pka", "fukui", "redox", "reaction"})


@dataclass(frozen=True)
class ChainStep:
    """One pipeline step: program name, job type, and optional CLI extras."""

    program: str
    job: str
    extra_option_tokens: tuple = ()
    extra_settings: tuple = ()
    extra_keywords: tuple = ()


@dataclass(frozen=True)
class ChainStepSpec:
    """Resolved chain step: job class and project-settings class."""

    job_class: type
    project_settings_class: type


_CHAIN_PROGRAMS = {
    "crest": (
        CRESTProjectSettings,
        {"conformers": CRESTConformerSearchJob},
    ),
    "xtb": (
        XTBProjectSettings,
        {
            "opt": XTBOptJob,
            "sp": XTBSinglePointJob,
            "hess": XTBHessJob,
        },
    ),
    "gaussian": (
        GaussianProjectSettings,
        {
            "opt": GaussianOptJob,
            "ts": GaussianTSJob,
            "sp": GaussianSinglePointJob,
        },
    ),
    "orca": (
        ORCAProjectSettings,
        {
            "opt": ORCAOptJob,
            "ts": ORCATSJob,
            "sp": ORCASinglePointJob,
        },
    ),
}

CHAIN_PROGRAM_SETTINGS = {
    program: settings_class
    for program, (settings_class, _jobs) in _CHAIN_PROGRAMS.items()
}
CHAIN_STEP_SPECS = {
    program: jobs
    for program, (_settings_class, jobs) in _CHAIN_PROGRAMS.items()
}

_JOB_SETTINGS = {
    "conformers": lambda settings: settings.conformer_settings(),
    "opt": lambda settings: settings.opt_settings(),
    "ts": lambda settings: settings.ts_settings(),
    "sp": lambda settings: settings.sp_settings(),
    "hess": lambda settings: settings.hess_settings(),
}


def parse_chain_step(token):
    """Parse a CLI step token into a :class:`ChainStep`.

    Accepts ``PROGRAM:JOB`` or ``PROGRAM: [OPTIONS] JOB``, e.g.
    ``gaussian:opt`` and ``gaussian: -o maxstep=8,maxsize=12 opt``.

    Args:
        token (str): One ``-s/--steps`` value.

    Returns:
        ChainStep: Parsed program, job, and extra option tokens.

    Raises:
        ValueError: If the token is missing ``:``, has empty parts, or is
            otherwise malformed.
    """
    token = token.strip()
    if ":" not in token:
        raise ValueError(
            f"Invalid chain step {token!r}; expected PROGRAM:JOB "
            "or PROGRAM: [OPTIONS] JOB."
        )
    program, rest = token.split(":", 1)
    program = program.strip().lower()
    rest = rest.strip()
    if not program or not rest:
        raise ValueError(
            f"Invalid chain step {token!r}; expected PROGRAM:JOB "
            "or PROGRAM: [OPTIONS] JOB."
        )
    try:
        tokens = shlex.split(rest)
    except ValueError as exc:
        raise ValueError(f"Invalid chain step {token!r}: {exc}") from exc
    if not tokens:
        raise ValueError(
            f"Invalid chain step {token!r}; expected PROGRAM:JOB "
            "or PROGRAM: [OPTIONS] JOB."
        )

    known_jobs = set(CHAIN_STEP_SPECS.get(program, ()))
    job_index = None
    job = None
    for index, part in enumerate(tokens):
        name = part.lower()
        if name in known_jobs or name in CHAIN_NESTED_ONLY_JOBS:
            job_index = index
            job = name
            break
    if job is None:
        if len(tokens) == 1 and not tokens[0].startswith("-"):
            job = tokens[0].lower()
            extra_option_tokens = ()
        else:
            raise ValueError(
                f"Invalid chain step {token!r}; expected PROGRAM:JOB "
                "or PROGRAM: [OPTIONS] JOB."
            )
    else:
        extra_option_tokens = tuple(tokens[:job_index])
        suffix = tokens[job_index + 1 :]
        if suffix:
            raise ValueError(
                f"Unsupported options after {job!r} in chain step {token!r}; "
                f"put program options before the job name, e.g. "
                '-s "gaussian: -o maxstep=8,maxsize=12 opt".'
            )
    return ChainStep(
        program=program,
        job=job,
        extra_option_tokens=extra_option_tokens,
    )


def _settings_for_chain_step(program, job, project_settings):
    settings_for = _JOB_SETTINGS.get(job)
    if settings_for is None:
        raise ValueError(
            f"No project settings method for chain step {program} {job!r}."
        )
    return settings_for(project_settings)


def get_chain_step_spec(program, job):
    """Return the registry entry for a ``(program, job)`` pair.

    Args:
        program (str): Chain program name (``crest``, ``xtb``,
            ``gaussian``, or ``orca``).
        job (str): Job name (e.g. ``opt``, ``sp``).

    Returns:
        ChainStepSpec: Job class and project-settings class for the pair.

    Raises:
        ValueError: If the pair is not a pipeline step. Workflow jobs
            (pKa, Fukui, redox, reaction) must be run as chain workflow
            subcommands or on the program CLI with that program's own
            ``-p`` project YAML.
    """
    jobs = CHAIN_STEP_SPECS.get(program)
    job_class = None if jobs is None else jobs.get(job)
    if job_class is not None:
        return ChainStepSpec(
            job_class=job_class,
            project_settings_class=CHAIN_PROGRAM_SETTINGS[program],
        )
    if jobs is not None and job in CHAIN_NESTED_ONLY_JOBS:
        raise ValueError(
            f"Chain pipeline steps do not support {program} {job!r}; "
            f"run it as chemsmart sub {program} -p <project> {job} ..."
        )
    allowed = ", ".join(jobs or ())
    if allowed:
        raise ValueError(
            f"Unsupported chain step job {job!r} for program "
            f"{program!r}. Supported jobs: {allowed}."
        )
    allowed_programs = ", ".join(ChainProjectSettings.PROGRAMS)
    raise ValueError(
        f"Unsupported chain step program {program!r}. "
        f"Allowed programs: {allowed_programs}."
    )


def molecule_from_completed_job(job, fallback=None):
    """Geometry to pass from a completed child to the next chain step.

    CREST uses the best conformer. Other jobs use
    ``optimized_structure()``, falling back to the output molecule when
    the child did not optimize (e.g. single-point).

    If the child is incomplete or no geometry can be read, return
    ``fallback``.
    """
    if not job.is_complete():
        return fallback
    if isinstance(job, CRESTJob):
        output = job._output()
        mol = None if output is None else output.best_conformer
        return fallback if mol is None else mol
    mol = job.optimized_structure()
    if mol is not None:
        return mol
    output = job._output()
    mol = None if output is None else output.molecule
    return fallback if mol is None else mol


def _first_set(*values):
    for value in values:
        if value is not None:
            return value
    return None


def _charge_multiplicity_overrides(molecule, charge, multiplicity):
    mol_charge = None if molecule is None else molecule.charge
    mol_mult = None if molecule is None else molecule.multiplicity
    overrides = {}
    for key, value in (
        ("charge", _first_set(charge, mol_charge)),
        ("multiplicity", _first_set(multiplicity, mol_mult)),
    ):
        if value is not None:
            overrides[key] = value
    return overrides


def _molecule_for_step(index, phases, initial_molecule):
    if index == 0:
        return initial_molecule
    previous_jobs = phases[index - 1].resolve_jobs()
    if not previous_jobs or not all(
        job.is_complete() for job in previous_jobs
    ):
        return None
    return molecule_from_completed_job(previous_jobs[0])


def build_chain_job(
    chain_settings,
    steps,
    molecule,
    label,
    charge=None,
    multiplicity=None,
    jobrunner=None,
    skip_completed=True,
    **kwargs,
):
    """Build a ``ChainJob`` from explicit pipeline steps.

    Each step becomes a ``JobPhase`` whose factory reads geometry from the
    previous completed child. Charge and multiplicity from ``charge`` /
    ``multiplicity`` (and otherwise from ``molecule``) are merged into
    each child with keywords ``("charge", "multiplicity")``.

    Args:
        chain_settings: Loaded chain project settings with program aliases.
        steps: Sequence of :class:`ChainStep` values.
        molecule: Input structure for the first step.
        label (str): Chain label; children are
            ``{label}_{index:02d}_{program}_{job}``.
        charge (int, optional): Chain-level charge override.
        multiplicity (int, optional): Chain-level multiplicity override.
        jobrunner: Runner assigned to the parent ``ChainJob``.
        skip_completed (bool): Forwarded to child jobs.

    Returns:
        ChainJob: Phased pipeline ready for ``run`` / ``sub``.

    Raises:
        ValueError: If there are no steps, ``molecule`` is missing, or a
            step's ``(program, job)`` pair is not a pipeline step.
    """
    if not steps:
        raise ValueError("Chain pipeline requires at least one step.")
    if molecule is None:
        raise ValueError("Chain pipeline requires a molecule.")

    prepared = []
    project_settings_by_program = {}
    for step in steps:
        spec = get_chain_step_spec(step.program, step.job)
        if step.program not in project_settings_by_program:
            project_name = chain_settings.project_for(step.program)
            project_settings_by_program[step.program] = (
                spec.project_settings_class.from_project(project_name)
            )
        prepared.append((step, spec))

    overrides = _charge_multiplicity_overrides(molecule, charge, multiplicity)
    phases = []

    def make_factory(index):
        def factory():
            mol = _molecule_for_step(index, phases, molecule)
            if mol is None:
                return []
            step, spec = prepared[index]
            settings = _settings_for_chain_step(
                step.program,
                step.job,
                project_settings_by_program[step.program],
            )
            merge_values = dict(overrides)
            merge_keywords = list(overrides)
            if step.extra_settings:
                merge_values.update(dict(step.extra_settings))
                merge_keywords.extend(step.extra_keywords)
            if merge_values:
                settings = settings.merge(
                    merge_values, keywords=tuple(merge_keywords)
                )
            child_label = f"{label}_{index:02d}_{step.program}_{step.job}"
            return [
                spec.job_class(
                    molecule=mol,
                    settings=settings,
                    label=child_label,
                    jobrunner=None,
                    skip_completed=skip_completed,
                )
            ]

        return factory

    for index, (step, _spec) in enumerate(prepared):
        phase_name = f"{index:02d}_{step.program}_{step.job}"
        phases.append(
            JobPhase(
                name=phase_name,
                jobs_factory=make_factory(index),
                require_complete=True,
                stop_on_incomplete=True,
                stop_message=(
                    f"{phase_name} incomplete, halting serial execution."
                ),
            )
        )

    return ChainJob(
        molecule=molecule,
        label=label,
        jobrunner=jobrunner,
        skip_completed=skip_completed,
        phases=phases,
        programs=tuple(dict.fromkeys(step.program for step, _ in prepared)),
        **kwargs,
    )
