"""Chain step parsing and ChainJob builder."""

import logging
import os
import shlex
from dataclasses import dataclass, replace
from typing import Callable, Sequence

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.chain import ChainJob, JobPhase
from chemsmart.jobs.crest.job import CRESTJob
from chemsmart.jobs.gaussian.crest import GaussianCrestJob
from chemsmart.jobs.gaussian.irc import GaussianIRCJob
from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.gaussian.link import GaussianLinkJob
from chemsmart.jobs.gaussian.qrc import GaussianQRCJob
from chemsmart.jobs.grouper.job import GrouperJob
from chemsmart.jobs.job import Job
from chemsmart.jobs.orca.job import ORCAJob
from chemsmart.jobs.orca.qrc import ORCAQRCJob
from chemsmart.jobs.runner import JobRunner
from chemsmart.jobs.thermochemistry.job import ThermochemistryJob
from chemsmart.jobs.xtb.job import XTBJob
from chemsmart.settings.chain import ChainProjectSettings
from chemsmart.settings.crest import CRESTProjectSettings
from chemsmart.settings.gaussian import GaussianProjectSettings
from chemsmart.settings.orca import ORCAProjectSettings
from chemsmart.settings.xtb import XTBProjectSettings

logger = logging.getLogger(__name__)

__all__ = [
    "ChainStep",
    "ChainBuildContext",
    "HandoffResult",
    "ChainStepParser",
    "ChainHandoffResolver",
    "ProgramAdapter",
    "PROGRAM_ADAPTERS",
    "ChainPipelineBuilder",
    "CHAIN_PROGRAM_SETTINGS",
    "CHAIN_NESTED_ONLY_JOBS",
    "validate_chain_step",
    "parse_chain_step",
    "molecule_from_completed_job",
    "build_chain_job",
]

# Public models

# Workflow jobs that are not ``-s`` pipeline steps.
CHAIN_NESTED_ONLY_JOBS = frozenset({"pka", "fukui", "redox", "reaction"})

# Program options the chain sets itself. A later ``-f`` is that command's
# own option, such as frozen atoms on ``gaussian opt``.
_CHAIN_OWNED_OPTIONS = frozenset(
    {
        "-f",
        "--filename",
        "-p",
        "--project",
        "-i",
        "--index",
        "--si",
        "--structure-index",
    }
)

# Subcommands that need a prepared file or a trajectory, not one structure.
_CHAIN_STEPS_WITHOUT_STRUCTURE = frozenset(
    {
        ("gaussian", "com"),
        ("gaussian", "traj"),
        ("gaussian", "dias"),
        ("orca", "inp"),
    }
)


@dataclass(frozen=True)
class ChainStep:
    """One pipeline step: program name, job type, and optional CLI extras."""

    program: str
    job: str
    extra_option_tokens: tuple = ()
    suffix_option_tokens: tuple = ()
    num_conformers: int | None = None


@dataclass(frozen=True)
class ChainBuildContext:
    """Inputs shared by every phase of one chain build."""

    chain_settings: ChainProjectSettings
    molecules: tuple[Molecule, ...]
    label: str
    charge: int | None = None
    multiplicity: int | None = None
    jobrunner: JobRunner | None = None
    skip_completed: bool = True


# Program registry


@dataclass(frozen=True)
class ProgramAdapter:
    """Project settings class and Click command for one chain program."""

    settings_class: type
    command_loader: Callable


def _load_crest_command():
    import chemsmart.cli.crest as crest_cli

    return crest_cli.crest


def _load_xtb_command():
    import chemsmart.cli.xtb as xtb_cli

    return xtb_cli.xtb


def _load_gaussian_command():
    import chemsmart.cli.gaussian as gaussian_cli

    return gaussian_cli.gaussian


def _load_orca_command():
    import chemsmart.cli.orca as orca_cli

    return orca_cli.orca


PROGRAM_ADAPTERS = {
    "crest": ProgramAdapter(CRESTProjectSettings, _load_crest_command),
    "xtb": ProgramAdapter(XTBProjectSettings, _load_xtb_command),
    "gaussian": ProgramAdapter(
        GaussianProjectSettings, _load_gaussian_command
    ),
    "orca": ProgramAdapter(ORCAProjectSettings, _load_orca_command),
}


# Step parsing and validation


class ChainStepParser:
    def parse(self, token: str) -> ChainStep:
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

        # Grouper reads the CREST ensemble in one command and has no
        # chain project.
        if program == "grouper":
            known_jobs = self._grouper_command_names()
        else:
            known_jobs = self.pipeline_job_names(program)
        blocked = (
            self._thermochemistry_command_names()
            if program == "thermochemistry"
            else frozenset()
        )
        job_index = None
        job = None
        for index, part in enumerate(tokens):
            name = part.lower()
            if name in blocked:
                raise ValueError(
                    f"Chain step {token!r} cannot use {name!r}. "
                    "A thermochemistry step runs the default analysis "
                    "and takes options only."
                )
            if name in known_jobs or name in CHAIN_NESTED_ONLY_JOBS:
                job_index = index
                job = name
                break
        suffix = ()
        num_conformers = None
        if job is None and program == "thermochemistry":
            job = "thermochemistry"
            extra_option_tokens = tuple(tokens)
        elif job is None:
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
            suffix = tuple(tokens[job_index + 1 :])
            if program == "crest" and job == "conformers" and suffix:
                num_conformers, suffix = self._split_crest_num_conformers(
                    suffix, token
                )
        return ChainStep(
            program=program,
            job=job,
            extra_option_tokens=extra_option_tokens,
            suffix_option_tokens=suffix,
            num_conformers=num_conformers,
        )

    def validate(self, step: ChainStep) -> None:
        program = step.program
        job = step.job
        for token in step.extra_option_tokens:
            name = token.split("=", 1)[0]
            if name in _CHAIN_OWNED_OPTIONS:
                raise ValueError(
                    f"Chain step {program}:{job} cannot set {name}; "
                    "the chain supplies the file, project, and "
                    "structure index."
                )
        if program == "grouper" and job in self._grouper_command_names():
            return
        if program == "thermochemistry":
            if job != "thermochemistry":
                raise ValueError(
                    "A thermochemistry chain step runs the default "
                    "analysis and takes options only."
                )
            return
        if program in PROGRAM_ADAPTERS and job in CHAIN_NESTED_ONLY_JOBS:
            raise ValueError(
                f"Chain pipeline steps do not support {program} {job!r}; "
                f"run it as chemsmart sub {program} -p <project> {job} ..."
            )
        if (program, job) in _CHAIN_STEPS_WITHOUT_STRUCTURE:
            raise ValueError(
                f"Chain pipeline steps do not support {program} {job!r}; "
                "the chain hands off one structure."
            )
        names = self.pipeline_job_names(program)
        if job in names:
            return
        if names:
            allowed = ", ".join(sorted(names))
            raise ValueError(
                f"Unsupported chain step job {job!r} for program "
                f"{program!r}. Supported jobs: {allowed}."
            )
        if program == "grouper":
            allowed = ", ".join(sorted(self._grouper_command_names()))
            raise ValueError(
                f"Unsupported chain step job {job!r} for program "
                f"{program!r}. Supported jobs: {allowed}."
            )
        allowed_programs = ", ".join(ChainProjectSettings.PROGRAMS)
        raise ValueError(
            f"Unsupported chain step program {program!r}. "
            f"Allowed programs: {allowed_programs}."
        )

    def pipeline_job_names(self, program: str) -> frozenset[str]:
        """Subcommands of ``program`` that can be a chain step.

        Workflow commands (pKa, Fukui, redox, reaction) stay on their own
        chain subcommands.
        """
        adapter = PROGRAM_ADAPTERS.get(program)
        if adapter is None:
            return frozenset()
        names = adapter.command_loader().commands
        return frozenset(
            name
            for name in names
            if name not in CHAIN_NESTED_ONLY_JOBS
            and (program, name) not in _CHAIN_STEPS_WITHOUT_STRUCTURE
        )

    def _grouper_command_names(self):
        """Strategy names registered on the grouper command."""
        import chemsmart.cli.grouper as grouper_cli

        return frozenset(grouper_cli.grouper.commands)

    def _thermochemistry_command_names(self):
        """Subcommands of the thermochemistry command."""
        import chemsmart.cli.thermochemistry as thermochemistry_cli

        return frozenset(thermochemistry_cli.thermochemistry.commands)

    def _split_crest_num_conformers(self, tokens, token):
        """Pull ``-nc N`` out of a ``crest:conformers`` step.

        ``N`` chooses the CREST step's handoff file and is not passed to
        CREST. ``1`` selects ``crest_best.xyz``. Any other count, and an
        omitted ``-nc``, select ``crest_conformers.xyz``. Other tokens are
        left for the CREST command.
        """
        names = {"-nc", "--num-confs-to-run"}
        kept = []
        count = None
        index = 0
        while index < len(tokens):
            name = tokens[index]
            if name not in names:
                kept.append(name)
                index += 1
                continue
            if count is not None:
                raise ValueError(f"Repeated -nc in chain step {token!r}.")
            if index + 1 >= len(tokens):
                raise ValueError(
                    "Unsupported options after 'conformers' in chain step "
                    f'{token!r}; use -nc N, e.g. -s "crest:conformers -nc 10".'
                )
            raw = tokens[index + 1]
            try:
                count = int(raw)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid -nc value {raw!r} in chain step {token!r}; "
                    "expected an integer >= 1."
                ) from exc
            if count < 1:
                raise ValueError(
                    f"Invalid -nc value {count} in chain step {token!r}; "
                    "expected an integer >= 1."
                )
            index += 2
        return count, tuple(kept)


# Geometry handoff


@dataclass(frozen=True)
class HandoffResult:
    """Structures for the next phase, and the fan-out label kind."""

    molecules: tuple[Molecule, ...]
    label_kind: str | None = None
    sources: tuple[tuple[str, str] | None, ...] = ()


class _ChainGroupPassthrough(Job):
    """Completed stand-in when CREST produced fewer than two structures."""

    TYPE = "grouper"
    PROGRAM = "grouper"

    def __init__(self, molecule, label, skip_completed=True):
        super().__init__(
            molecule=molecule,
            label=label,
            jobrunner=None,
            skip_completed=skip_completed,
        )
        self.groups = [[molecule]]

    def is_complete(self):
        return True

    def _run(self, **kwargs):
        return None


class ChainHandoffResolver:
    def molecule_from_job(self, job, fallback=None):
        if not job.is_complete():
            return fallback
        if isinstance(job, ThermochemistryJob):
            mol = job.molecule
            return fallback if mol is None else mol
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

    def resolve(self, index, phases, initial_molecule, num_conformers=None):
        """Structures for one step.

        Returns ``None`` when the previous phase is not ready.
        ``label_kind`` is ``g`` for group representatives, ``idx`` for
        structures selected by ``-i``, and ``c`` for other multi-structure
        handoffs. A completed Gaussian crest job
        contributes only its lowest-energy optimized conformer. A
        completed CREST job hands off ``crest_best.xyz`` when ``-nc`` is
        1, and ``crest_conformers.xyz`` otherwise. A
        completed Gaussian IRC job contributes its directional endpoint
        jobs, one per log. Gaussian, ORCA, and xTB jobs also pass their
        output file, so the next step reads the last structure from it.
        A thermochemistry step passes that same log or out file on.
        """
        if index == 0:
            if isinstance(initial_molecule, Molecule):
                molecules = (initial_molecule,)
            else:
                molecules = tuple(initial_molecule)
            kind = "idx" if len(molecules) > 1 else None
            return HandoffResult(molecules=molecules, label_kind=kind)
        jobs = self._completed_previous_jobs(index, phases)
        if jobs is None:
            return None
        if len(jobs) == 1:
            groups = self._groups_from_job(jobs[0])
            if groups is not None:
                representatives = tuple(
                    _with_source_charge(group[0], jobs[0].molecule)
                    for group in groups
                    if group
                )
                return HandoffResult(molecules=representatives, label_kind="g")
            lowest = self._lowest_gaussian_crest_job(jobs[0])
            if lowest is not None:
                jobs = lowest
            elif isinstance(jobs[0], CRESTJob):
                crest_job = jobs[0]
                if num_conformers == 1:
                    selected = self._best_conformer(crest_job)
                    if selected is None:
                        return None
                    return HandoffResult(
                        molecules=tuple(selected),
                        sources=((crest_job.best_file, None),),
                    )
                if not os.path.isfile(crest_job.conformers_file):
                    logger.info(
                        "CREST conformers file is not available yet; "
                        "waiting before the next step."
                    )
                    return None
                frames = (
                    Molecule.from_filepath(
                        filepath=crest_job.conformers_file,
                        index=":",
                        return_list=True,
                    )
                    or []
                )
                molecules = [
                    _with_source_charge(molecule, crest_job.molecule)
                    for molecule in frames
                ]
                if not molecules:
                    return None
                sources = tuple(
                    (crest_job.conformers_file, str(frame_index))
                    for frame_index in range(1, len(molecules) + 1)
                )
                kind = "c" if len(molecules) > 1 else None
                return HandoffResult(
                    molecules=tuple(molecules),
                    label_kind=kind,
                    sources=sources,
                )
        molecules = []
        sources = []
        file_indices = []
        for job in self._expand_handoff_jobs(jobs):
            molecule = self.molecule_from_job(job)
            if molecule is None:
                return None
            molecules.append(molecule)
            file_indices.append(
                job.molecule.structure_index_in_file
                if job.molecule is not None
                else None
            )
            if isinstance(job, (GaussianJob, ORCAJob, XTBJob)):
                sources.append((job.outputfile, "-1"))
            elif isinstance(job, ThermochemistryJob) and job.filename:
                sources.append((job.filename, "-1"))
            else:
                sources.append(None)
        use_file_index = (
            len(molecules) > 1
            and None not in file_indices
            and len(set(file_indices)) == len(file_indices)
        )
        if use_file_index:
            for molecule, file_index in zip(molecules, file_indices):
                molecule.structure_index_in_file = file_index
            kind = "idx"
        else:
            kind = "c" if len(molecules) > 1 else None
        return HandoffResult(
            molecules=tuple(molecules),
            label_kind=kind,
            sources=tuple(sources),
        )

    def _expand_handoff_jobs(self, jobs):
        """Replace a two-direction job with the jobs that write its logs."""
        expanded = []
        for job in jobs:
            if isinstance(job, GaussianIRCJob):
                expanded.extend(job.endpoint_jobs())
            elif isinstance(job, (GaussianQRCJob, ORCAQRCJob)):
                expanded.extend(job.both_qrc_jobs)
            elif isinstance(job, GaussianLinkJob):
                # handle the case where the job is a Gaussian link IRC job
                endpoints = job._get_irc_jobs()
                expanded.extend(endpoints if endpoints else [job])
            else:
                expanded.append(job)
        return expanded

    def _completed_previous_jobs(self, index, phases):
        if index == 0:
            return []
        jobs = phases[index - 1].resolve_jobs()
        if not jobs or not all(job.is_complete() for job in jobs):
            return None
        return jobs

    def _lowest_gaussian_crest_job(self, job):
        """Lowest-energy conformer job of a finished Gaussian crest job.

        Returns ``None`` when ``job`` is not a Gaussian crest job, ``[]``
        when a selected conformer geometry or its energy is not available
        yet, and a one-item list otherwise.
        """
        if not isinstance(job, GaussianCrestJob):
            return None
        children = job.all_conformers_jobs[: job.num_confs_to_opt]
        energies = []
        for child in children:
            molecule = self.molecule_from_job(child)
            if molecule is None or molecule.energy is None:
                logger.info(
                    "Gaussian crest conformer energies are not available "
                    "yet; waiting before the next step."
                )
                return []
            energies.append((molecule.energy, child))
        if not energies:
            return []
        return [min(energies, key=lambda pair: pair[0])[1]]

    def _best_conformer(self, job):
        """Single structure from ``crest_best.xyz``.

        Returns ``None`` when that file is not available yet.
        """
        if not os.path.isfile(job.best_file):
            logger.info(
                "CREST best conformer file is not available yet; "
                "waiting before the next step."
            )
            return None
        frames = (
            Molecule.from_filepath(
                filepath=job.best_file, index=":", return_list=True
            )
            or []
        )
        if not frames:
            return None
        return [_with_source_charge(frames[0], job.molecule)]

    def _groups_from_job(self, job):
        if isinstance(job, _ChainGroupPassthrough):
            return job.groups
        if not isinstance(job, GrouperJob):
            return None
        if job._groups:
            return [list(group) for group in job._groups]
        return self._groups_from_output_files(job)

    def _groups_from_output_files(self, job):
        groups = []
        for path in self._group_xyz_paths(job):
            molecules = (
                Molecule.from_filepath(
                    filepath=path, index=":", return_list=True
                )
                or []
            )
            if molecules:
                groups.append(
                    [
                        _with_source_charge(molecule, job.molecule)
                        for molecule in molecules
                    ]
                )
        return groups

    def _group_xyz_paths(self, job):
        folder = job.output_dir
        if not os.path.isdir(folder):
            return []
        prefix = f"{job.label}_group_"
        numbered = []
        for name in os.listdir(folder):
            if not name.startswith(prefix) or not name.endswith(".xyz"):
                continue
            stem = name[len(prefix) : -4]
            if stem.isdigit():
                numbered.append((int(stem), os.path.join(folder, name)))
        numbered.sort()
        return [path for _number, path in numbered]


def _with_source_charge(molecule, source):
    copied = molecule.copy()
    if source is None:
        return copied
    if source.charge is not None:
        copied.charge = source.charge
    if source.multiplicity is not None:
        copied.multiplicity = source.multiplicity
    return copied


# Command adapters


def _argv_has_option(tokens, *names):
    return any(token in names for token in tokens)


def _with_fanout_label(tokens, suffix):
    """Append ``suffix`` to a step's ``-l``/``--label`` value."""
    updated = []
    index = 0
    while index < len(tokens):
        token = tokens[index]
        if token in ("-l", "--label") and index + 1 < len(tokens):
            updated.extend((token, f"{tokens[index + 1]}{suffix}"))
            index += 2
            continue
        if token.startswith("--label="):
            updated.append(f"{token}{suffix}")
            index += 1
            continue
        updated.append(token)
        index += 1
    return tuple(updated)


def _invoke_cli_command(command, argv, jobrunner):
    import click

    try:
        context = command.make_context(
            command.name, list(argv), obj={"jobrunner": jobrunner}
        )
        result = command.invoke(context)
    except click.ClickException as exc:
        raise ValueError(str(exc)) from exc
    if result is None:
        raise ValueError(f"{command.name} did not return a job.")
    return result


def _program_argv(
    step,
    filename,
    structure_index,
    project_name,
    child_label,
    charge,
    multiplicity,
    skip_completed,
):
    argv = ["-p", project_name, "-f", filename]
    if structure_index is not None:
        argv.extend(["-i", structure_index])
    before = step.extra_option_tokens
    after = step.suffix_option_tokens
    if charge is not None and not _argv_has_option(before, "-c", "--charge"):
        argv.extend(["-c", str(charge)])
    if multiplicity is not None and not _argv_has_option(
        before, "-m", "--multiplicity"
    ):
        argv.extend(["-m", str(multiplicity)])
    if not _argv_has_option(before, "-l", "--label", "-a", "--append-label"):
        argv.extend(["-l", child_label])
    argv.extend(before)
    argv.append(step.job)
    argv.extend(after)
    flags = before + after
    if not skip_completed and not _argv_has_option(
        flags, "-S", "-R", "--skip-completed", "--no-skip-completed"
    ):
        argv.append("--no-skip-completed")
    return argv


def _thermochemistry_invocation(step):
    """Parse a thermochemistry step into settings and Click parameters.

    ``-f``, ``-d``, and ``-p`` are rejected. The chain supplies the
    previous output file, and ``-T`` is required.
    """
    import click

    from chemsmart.cli.thermochemistry.thermochemistry import (
        resolve_entropy_cutoff,
        thermochemistry,
    )
    from chemsmart.jobs.thermochemistry.settings import (
        ThermochemistryJobSettings,
    )

    tokens = list(step.extra_option_tokens) + list(step.suffix_option_tokens)
    try:
        context = thermochemistry.make_context(
            "thermochemistry", tokens, obj={}
        )
    except click.ClickException as exc:
        raise ValueError(str(exc)) from exc
    params = context.params
    if params.get("temperature") is None:
        raise ValueError(
            "A thermochemistry chain step requires -T/--temperature."
        )
    if (
        params.get("directory")
        or params.get("filenames")
        or params.get("program")
        or params.get("filetype")
    ):
        raise ValueError(
            "A thermochemistry chain step cannot set -d/--directory, "
            "-t/--filetype, -f/--filenames, or -p/--program. "
            "The chain supplies the previous output file."
        )
    cutoff_entropy, entropy_method = resolve_entropy_cutoff(
        params["cutoff_entropy_grimme"], params["cutoff_entropy_truhlar"]
    )
    settings = ThermochemistryJobSettings(
        temperature=params["temperature"],
        concentration=params["concentration"],
        pressure=params["pressure"],
        use_weighted_mass=params["weighted"],
        alpha=params["alpha"],
        s_freq_cutoff=cutoff_entropy,
        entropy_method=entropy_method,
        h_freq_cutoff=params["cutoff_enthalpy"],
        energy_units=params["energy_units"],
        outputfile=params["outputfile"],
        overwrite=params["overwrite"],
        check_imaginary_frequencies=params["check_imaginary_frequencies"],
    )
    return settings, params


def _grouper_argv(step, conformers_file, label):
    argv = ["-f", conformers_file]
    before = step.extra_option_tokens
    if not _argv_has_option(before, "-l", "--label", "-a", "--append-label"):
        argv.extend(["-l", label])
    argv.extend(before)
    argv.append(step.job)
    argv.extend(step.suffix_option_tokens)
    return argv


def _jobs_from_cli_result(result):
    """Keep the charge the command set.

    Charge and multiplicity are carried on from the previous step.
    """
    jobs = result if isinstance(result, list) else [result]
    for job in jobs:
        if job.molecule.charge is None:
            job.molecule.charge = job.settings.charge
        if job.molecule.multiplicity is None:
            job.molecule.multiplicity = job.settings.multiplicity
    return jobs


def _step_slug(step):
    """Phase and child-label token for one step."""
    if step.program == "thermochemistry":
        return step.program
    return f"{step.program}_{step.job}"


def _child_label(chain_label, index, step, ordinal, count, kind):
    base = f"{chain_label}_{index:02d}_{_step_slug(step)}"
    if count > 1 and kind:
        return f"{base}_{kind}{ordinal}"
    return base


def _handoff_child_label(chain_label, index, step, molecule, ordinal, handoff):
    """Child label for one structure in a fan-out."""
    label_ordinal = ordinal
    if (
        handoff.label_kind == "idx"
        and molecule.structure_index_in_file is not None
    ):
        label_ordinal = molecule.structure_index_in_file
    child_label = _child_label(
        chain_label,
        index,
        step,
        label_ordinal,
        len(handoff.molecules),
        handoff.label_kind,
    )
    return label_ordinal, child_label


# Pipeline builder


def _first_set(*values):
    for value in values:
        if value is not None:
            return value
    return None


class ChainPipelineBuilder:
    def __init__(self, context: ChainBuildContext):
        self.context = context
        self.phases = []
        self.handoff = ChainHandoffResolver()
        self._thermochemistry = {}

    def build(self, steps: Sequence[ChainStep], **kwargs) -> ChainJob:
        if not steps:
            raise ValueError("Chain pipeline requires at least one step.")
        if not self.context.molecules:
            raise ValueError("Chain pipeline requires a molecule.")

        prepared, project_names = self._validate_pipeline(steps)
        self._prepared = prepared
        self._project_names = project_names
        parent = (
            self.context.molecules[-1]
            if len(self.context.molecules) == 1
            else None
        )
        self.phases = []

        def make_factory(index):
            def factory():
                return self._jobs_for_phase(index)

            return factory

        for index, step in enumerate(prepared):
            phase_name = f"{index:02d}_{_step_slug(step)}"
            self.phases.append(
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

        # Grouper runs inside the chain. ``programs`` lists the external
        # programs whose server environment the submission script loads.
        external_programs = tuple(
            dict.fromkeys(
                step.program
                for step in prepared
                if step.program in PROGRAM_ADAPTERS
            )
        )
        return ChainJob(
            molecule=parent,
            label=self.context.label,
            jobrunner=self.context.jobrunner,
            skip_completed=self.context.skip_completed,
            phases=self.phases,
            programs=external_programs,
            **kwargs,
        )

    def _validate_pipeline(self, steps):
        parser = ChainStepParser()
        prepared = []
        project_names = {}
        settings = self.context.chain_settings
        for step in steps:
            parser.validate(step)
            if step.program in PROGRAM_ADAPTERS:
                if step.program not in project_names:
                    project_names[step.program] = settings.project_for(
                        step.program
                    )
            prepared.append(step)

        for index, step in enumerate(prepared):
            if step.program == "thermochemistry":
                previous = None if index == 0 else prepared[index - 1]
                if previous is None or previous.program not in {
                    "gaussian",
                    "orca",
                    "xtb",
                    "thermochemistry",
                }:
                    raise ValueError(
                        "thermochemistry must follow a gaussian, orca, "
                        "or xtb step that writes a .log or .out file."
                    )
                self._thermochemistry[index] = _thermochemistry_invocation(
                    step
                )
            if not self._needs_conformers_file(step):
                continue
            previous = None if index == 0 else prepared[index - 1]
            if (
                previous is None
                or previous.program != "crest"
                or previous.job != "conformers"
            ):
                raise ValueError(
                    f"{step.program}:{step.job} must follow crest:conformers."
                )
        return prepared, project_names

    def _jobs_for_phase(self, index):
        step = self._prepared[index]
        previous_step = None if index == 0 else self._prepared[index - 1]
        num_conformers = (
            None if previous_step is None else previous_step.num_conformers
        )
        initial = self.context.molecules
        handoff = self.handoff.resolve(
            index, self.phases, initial, num_conformers
        )
        if handoff is None or not handoff.molecules:
            return []
        if self._needs_conformers_file(step):
            source = handoff.sources[0] if handoff.sources else None
            if source is None:
                return []
            filename = source[0]
            if step.program == "grouper":
                return self._grouper_jobs(
                    step, filename, handoff.molecules, index
                )
            return self._gaussian_crest_jobs(
                step, filename, handoff.molecules[0]
            )
        if step.program == "thermochemistry":
            return self._thermochemistry_jobs(step, handoff, index)
        return self._program_jobs(step, handoff, index)

    def _program_jobs(self, step, handoff, index):
        """Build one program command per structure.

        The chain passes the previous job's output file, or writes the
        structure to an xyz file when there is none. The first step gets
        the chain charge and multiplicity, or each structure's own when
        the chain omits them. A later log or out supplies its own. An xyz
        supplies the charge carried on its molecule.
        """
        command = PROGRAM_ADAPTERS[step.program].command_loader()
        molecules = handoff.molecules
        count = len(molecules)
        sources = handoff.sources or (None,) * count
        jobs = []
        for ordinal, (molecule, source) in enumerate(
            zip(molecules, sources), start=1
        ):
            label_ordinal, child_label = _handoff_child_label(
                self.context.label, index, step, molecule, ordinal, handoff
            )
            if source is None:
                filename = os.path.abspath(f"{child_label}.xyz")
                molecule.write_xyz(filename, mode="w")
                structure_index = None
            else:
                filename, structure_index = source
            argv_step = step
            if (
                count > 1
                and handoff.label_kind
                and (
                    _argv_has_option(step.extra_option_tokens, "-l", "--label")
                    or any(
                        token.startswith("--label=")
                        for token in step.extra_option_tokens
                    )
                )
            ):
                argv_step = replace(
                    step,
                    extra_option_tokens=_with_fanout_label(
                        step.extra_option_tokens,
                        f"_{handoff.label_kind}{label_ordinal}",
                    ),
                )
            if index == 0:
                charge = _first_set(self.context.charge, molecule.charge)
                multiplicity = _first_set(
                    self.context.multiplicity, molecule.multiplicity
                )
            elif filename.endswith(".xyz"):
                charge, multiplicity = molecule.charge, molecule.multiplicity
            else:
                charge, multiplicity = None, None
            argv = _program_argv(
                argv_step,
                filename,
                structure_index,
                self._project_names[step.program],
                child_label,
                charge,
                multiplicity,
                self.context.skip_completed,
            )
            result = _invoke_cli_command(command, argv, self.context.jobrunner)
            children = _jobs_from_cli_result(result)
            if handoff.label_kind == "idx":
                for child in children:
                    child.molecule.structure_index_in_file = (
                        molecule.structure_index_in_file
                    )
            jobs.extend(children)
        return jobs

    def _thermochemistry_jobs(self, step, handoff, index):
        """One thermochemistry job per previous .log or .out file."""
        settings, params = self._thermochemistry[index]
        molecules = handoff.molecules
        count = len(molecules)
        sources = handoff.sources or (None,) * count
        if len(sources) != count:
            sources = (None,) * count
        flags = step.extra_option_tokens + step.suffix_option_tokens
        if _argv_has_option(
            flags, "-S", "-R", "--skip-completed", "--no-skip-completed"
        ):
            skip_completed = params["skip_completed"]
        else:
            skip_completed = self.context.skip_completed
        user_label = params.get("label")
        append_label = params.get("append_label")
        jobs = []
        for ordinal, (molecule, source) in enumerate(
            zip(molecules, sources), start=1
        ):
            filename = None if source is None else source[0]
            if filename is None or not str(filename).endswith(
                (".log", ".out")
            ):
                raise ValueError(
                    "thermochemistry needs the previous step's "
                    ".log or .out file."
                )
            label_ordinal, child_label = _handoff_child_label(
                self.context.label, index, step, molecule, ordinal, handoff
            )
            if user_label:
                child_label = user_label
                if count > 1 and handoff.label_kind:
                    child_label = (
                        f"{user_label}_{handoff.label_kind}{label_ordinal}"
                    )
            elif append_label:
                child_label = f"{child_label}_{append_label}"
            job_settings = settings.copy()
            if count > 1 and job_settings.outputfile and handoff.label_kind:
                root, ext = os.path.splitext(job_settings.outputfile)
                job_settings.outputfile = (
                    f"{root}_{handoff.label_kind}{label_ordinal}{ext}"
                )
            job = ThermochemistryJob(
                filename=filename,
                molecule=molecule,
                settings=job_settings,
                label=child_label,
                jobrunner=self.context.jobrunner,
                skip_completed=skip_completed,
            )
            if job.settings.outputfile is None:
                job.settings.outputfile = job.outputfile
            jobs.append(job)
        return jobs

    def _grouper_jobs(self, step, conformers_file, molecules, index):
        if len(molecules) < 2:
            logger.info(
                "Fewer than 2 CREST conformers; skipping grouper and "
                "passing the available structure on."
            )
            if not molecules:
                return []
            label = f"{self.context.label}_{index:02d}_grouper_{step.job}"
            return [
                _ChainGroupPassthrough(
                    molecule=molecules[0],
                    label=label,
                    skip_completed=self.context.skip_completed,
                )
            ]
        from chemsmart.cli.grouper.grouper import grouper

        job = _invoke_cli_command(
            grouper,
            _grouper_argv(step, conformers_file, self.context.label),
            self.context.jobrunner,
        )
        job.skip_completed = self.context.skip_completed
        return [job]

    def _gaussian_crest_jobs(self, step, conformers_file, molecule):
        argv = _program_argv(
            step,
            conformers_file,
            None,
            self._project_names[step.program],
            self.context.label,
            molecule.charge,
            molecule.multiplicity,
            self.context.skip_completed,
        )
        command = PROGRAM_ADAPTERS["gaussian"].command_loader()
        return _jobs_from_cli_result(
            _invoke_cli_command(command, argv, self.context.jobrunner)
        )

    @staticmethod
    def _needs_conformers_file(step):
        return step.program == "grouper" or (
            step.program == "gaussian" and step.job == "crest"
        )


# Compatibility API

CHAIN_PROGRAM_SETTINGS = {
    name: adapter.settings_class for name, adapter in PROGRAM_ADAPTERS.items()
}

_DEFAULT_PARSER = ChainStepParser()
_DEFAULT_HANDOFF = ChainHandoffResolver()


def parse_chain_step(token):
    """Parse a CLI step token into a :class:`ChainStep`.

    Accepts ``PROGRAM:JOB`` or ``PROGRAM: [OPTIONS] JOB``, e.g.
    ``gaussian:opt`` and ``gaussian: -o maxstep=8,maxsize=12 opt``.
    Options before the job name are program options. Options after it
    are that command's options. ``crest:conformers -nc 1`` hands off
    ``crest_best.xyz``. Omitting ``-nc``, or any other ``-nc N``, hands
    off ``crest_conformers.xyz``. ``-nc`` is not passed to CREST.
    ``thermochemistry: [OPTIONS]`` has no job name. It runs the default
    thermochemistry analysis, for example
    ``thermochemistry: -T 353.15 -c 1.0 -csg 100 -ch 100``. ``-T`` is
    required. ``-c`` on that step is concentration in mol/L.

    Args:
        token (str): One ``-s/--steps`` value.

    Returns:
        ChainStep: Parsed program, job, and extra option tokens.

    Raises:
        ValueError: If the token is missing ``:``, has empty parts, or is
            otherwise malformed.
    """
    return _DEFAULT_PARSER.parse(token)


def validate_chain_step(program, job):
    """Require ``(program, job)`` to be a pipeline step.

    Args:
        program (str): Chain program name (``crest``, ``xtb``,
            ``gaussian``, ``orca``, ``grouper``, or ``thermochemistry``).
        job (str): Job name (e.g. ``opt``, ``sp``, ``irc``). A grouper
            job is a strategy registered on the grouper command. A
            program job is a subcommand of that program, other than
            pKa, Fukui, redox, reaction, ``gaussian:com``,
            ``gaussian:traj``, ``gaussian:dias``, and ``orca:inp``.
            A thermochemistry job is ``thermochemistry``: the default
            analysis, with options and no separate job token.

    Raises:
        ValueError: If the pair is not a pipeline step. Workflow jobs
            (pKa, Fukui, redox, reaction) must be run as chain workflow
            subcommands or on the program CLI with that program's own
            ``-p`` project YAML.
    """
    return _DEFAULT_PARSER.validate(ChainStep(program, job))


def molecule_from_completed_job(job, fallback=None):
    """Geometry to pass from a completed child to the next chain step.

    CREST uses the best conformer. Other jobs use
    ``optimized_structure()``, falling back to the output molecule when
    the child did not optimize (e.g. single-point).

    If the child is incomplete or no geometry can be read, return
    ``fallback``.
    """
    return _DEFAULT_HANDOFF.molecule_from_job(job, fallback)


def build_chain_job(
    chain_settings,
    steps,
    molecules,
    label,
    charge=None,
    multiplicity=None,
    jobrunner=None,
    skip_completed=True,
    **kwargs,
):
    """Build a ``ChainJob`` from explicit pipeline steps.

    Each step becomes a ``JobPhase`` whose factory reads geometry from the
    previous completed child. ``crest:conformers -nc 1`` passes
    ``crest_best.xyz``. Omitting ``-nc``, or any other ``-nc N``, passes
    ``crest_conformers.xyz``. ``gaussian:crest`` and ``grouper`` read that
    same file. The step after Gaussian crest receives that job's
    lowest-energy optimized conformer. The step after Gaussian IRC
    receives one structure from each directional log. A thermochemistry
    step reads the previous Gaussian, ORCA, or xTB ``.log`` or ``.out``
    and does not use a chain project alias. ``-T`` is required.
    ``charge`` and ``multiplicity`` (otherwise taken from each starting
    structure) are passed only to the first program command. One starting
    structure is stored on the chain job. More than one starting structure
    generates one chain job per structure; each child stores its own
    and is labeled with that structure's 1-based
    ``structure_index_in_file`` (``_idxN``). Later steps keep the index.
    A later Gaussian, ORCA, or xTB step reads charge and multiplicity from
    the previous output file.

    Args:
        chain_settings: Loaded chain project settings with program aliases.
        steps: Sequence of :class:`ChainStep` values.
        molecules: Starting structures. A single :class:`Molecule` is one
            structure. With one structure, the chain job stores
            ``molecules[-1]``. With more than one, the chain job stores
            none; each child stores its own. Each structure selected by
            ``-i`` needs ``structure_index_in_file`` set to its 1-based
            index in the file.
        label (str): Chain label; children are
            ``{label}_{index:02d}_{program}_{job}``.
        charge (int, optional): Chain-level charge override.
        multiplicity (int, optional): Chain-level multiplicity override.
        jobrunner: Runner assigned to the parent ``ChainJob``.
        skip_completed (bool): Forwarded to child jobs.

    Returns:
        ChainJob: Phased pipeline ready for ``run`` / ``sub``.

    Raises:
        ValueError: If there are no steps, ``molecules`` is missing, a
            step's ``(program, job)`` pair is not a pipeline step, or
            ``gaussian:crest`` / ``grouper`` does not follow
            ``crest:conformers``.
    """
    if isinstance(molecules, Molecule):
        molecules = (molecules,)
    else:
        molecules = () if molecules is None else tuple(molecules)
    if not molecules:
        raise ValueError("Chain pipeline requires a molecule.")
    context = ChainBuildContext(
        chain_settings=chain_settings,
        molecules=molecules,
        label=label,
        charge=charge,
        multiplicity=multiplicity,
        jobrunner=jobrunner,
        skip_completed=skip_completed,
    )
    return ChainPipelineBuilder(context).build(steps, **kwargs)
