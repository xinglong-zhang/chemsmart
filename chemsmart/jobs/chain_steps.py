"""Chain step parsing and ChainJob builder."""

import logging
import os
import shlex
from dataclasses import dataclass
from typing import Callable, Sequence

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.chain import ChainJob, JobPhase
from chemsmart.jobs.crest.job import CRESTJob
from chemsmart.jobs.gaussian.crest import GaussianCrestJob
from chemsmart.jobs.grouper.job import GrouperJob
from chemsmart.jobs.job import Job
from chemsmart.jobs.runner import JobRunner
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
    molecule: Molecule
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
        job_index = None
        job = None
        for index, part in enumerate(tokens):
            name = part.lower()
            if name in known_jobs or name in CHAIN_NESTED_ONLY_JOBS:
                job_index = index
                job = name
                break
        suffix = ()
        num_conformers = None
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
        if program == "grouper" and job in self._grouper_command_names():
            return
        if program in PROGRAM_ADAPTERS and job in CHAIN_NESTED_ONLY_JOBS:
            raise ValueError(
                f"Chain pipeline steps do not support {program} {job!r}; "
                f"run it as chemsmart sub {program} -p <project> {job} ..."
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
            name for name in names if name not in CHAIN_NESTED_ONLY_JOBS
        )

    def _grouper_command_names(self):
        """Strategy names registered on the grouper command."""
        import chemsmart.cli.grouper as grouper_cli

        return frozenset(grouper_cli.grouper.commands)

    def _split_crest_num_conformers(self, tokens, token):
        """Pull ``-nc N`` out of a ``crest:conformers`` step.

        ``N`` selects conformers for the next step and is not passed to
        CREST. Other tokens are left for the CREST command.
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
        ``label_kind`` is ``g`` for group representatives and ``c`` for
        other multi-structure handoffs. A completed Gaussian crest job
        contributes only its lowest-energy optimized conformer. A
        completed CREST job contributes its lowest ``num_conformers``
        frames, or the lowest frame when that value is omitted.
        """
        if index == 0:
            return HandoffResult(molecules=(initial_molecule,))
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
            ensemble = self._gaussian_crest_molecules(jobs[0])
            if ensemble is not None:
                return HandoffResult(molecules=tuple(ensemble))
            if isinstance(jobs[0], CRESTJob):
                selected = self._selected_crest_conformers(
                    jobs[0], num_conformers
                )
                if selected is None:
                    return None
                kind = "c" if len(selected) > 1 else None
                return HandoffResult(
                    molecules=tuple(selected), label_kind=kind
                )
        molecules = []
        for job in jobs:
            molecule = self.molecule_from_job(job)
            if molecule is None:
                return None
            molecules.append(molecule)
        kind = "c" if len(molecules) > 1 else None
        return HandoffResult(molecules=tuple(molecules), label_kind=kind)

    def conformers_file(self, crest_job, num_conformers):
        """CREST xyz path for a following gaussian-crest or grouper step.

        With ``-nc N``, write the selected conformers to a new file.
        Without it, return ``crest_conformers.xyz``.
        """
        if crest_job is None or not os.path.isfile(crest_job.conformers_file):
            return None
        if num_conformers is None:
            return crest_job.conformers_file
        selected = self._selected_crest_conformers(crest_job, num_conformers)
        if not selected:
            return None
        folder = os.path.dirname(crest_job.conformers_file)
        path = os.path.join(folder, f"crest_lowest_{len(selected)}.xyz")
        for index, molecule in enumerate(selected):
            molecule.write_xyz(path, mode="w" if index == 0 else "a")
        return path

    def _completed_previous_jobs(self, index, phases):
        if index == 0:
            return []
        jobs = phases[index - 1].resolve_jobs()
        if not jobs or not all(job.is_complete() for job in jobs):
            return None
        return jobs

    def _gaussian_crest_molecules(self, job):
        """Lowest-energy optimized conformer from a finished Gaussian crest job.

        Returns ``None`` when ``job`` is not a Gaussian crest job, ``[]``
        when a selected conformer geometry or its energy is not available
        yet, and a one-item list otherwise.
        """
        if not isinstance(job, GaussianCrestJob):
            return None
        children = job.all_conformers_jobs[: job.num_confs_to_opt]
        molecules = []
        for child in children:
            molecule = self.molecule_from_job(child)
            if molecule is None or molecule.energy is None:
                logger.info(
                    "Gaussian crest conformer energies are not available "
                    "yet; waiting before the next step."
                )
                return []
            molecules.append(molecule)
        if not molecules:
            return []
        return [min(molecules, key=lambda molecule: molecule.energy)]

    def _selected_crest_conformers(self, job, num_conformers):
        """Lowest CREST conformers to pass to the next step.

        ``num_conformers`` defaults to 1. The result is that many
        lowest-energy frames of ``crest_conformers.xyz``. Returns ``None``
        when that file is not available yet.
        """
        count = 1 if num_conformers is None else num_conformers
        if not os.path.isfile(job.conformers_file):
            logger.info(
                "CREST conformers file is not available yet; "
                "waiting before selecting conformers."
            )
            return None
        frames = (
            Molecule.from_filepath(
                filepath=job.conformers_file, index=":", return_list=True
            )
            or []
        )
        molecules = [
            _with_source_charge(molecule, job.molecule) for molecule in frames
        ]
        if not molecules:
            return None
        if all(molecule.energy is not None for molecule in molecules):
            molecules = sorted(molecules, key=lambda molecule: molecule.energy)
        if len(molecules) < count:
            logger.info(
                "CREST produced %s conformer(s); requested %s. "
                "Using the available conformers.",
                len(molecules),
                count,
            )
        return molecules[:count]

    def crest_source(self, jobs):
        for job in jobs:
            if isinstance(job, CRESTJob) and os.path.isfile(
                job.conformers_file
            ):
                return job
        return None

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
    project_name,
    child_label,
    charge,
    multiplicity,
    skip_completed,
):
    argv = ["-p", project_name, "-f", filename]
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


def _grouper_argv(step, conformers_file, label):
    argv = ["-f", conformers_file]
    before = step.extra_option_tokens
    if not _argv_has_option(before, "-l", "--label", "-a", "--append-label"):
        argv.extend(["-l", label])
    argv.extend(before)
    argv.append(step.job)
    argv.extend(step.suffix_option_tokens)
    return argv


def _jobs_from_cli_result(result, molecule):
    jobs = result if isinstance(result, list) else [result]
    for job in jobs:
        if molecule.charge is not None:
            job.molecule.charge = molecule.charge
        if molecule.multiplicity is not None:
            job.molecule.multiplicity = molecule.multiplicity
    return jobs


def _child_label(chain_label, index, step, ordinal, count, kind):
    base = f"{chain_label}_{index:02d}_{step.program}_{step.job}"
    if count > 1 and kind:
        return f"{base}_{kind}{ordinal}"
    return base


# Pipeline builder


def _first_set(*values):
    """Return the first value that is not ``None``.

    Args:
        *values: Candidates in priority order.

    Returns:
        The first set value, or ``None`` when every candidate is ``None``.
    """
    for value in values:
        if value is not None:
            return value
    return None


class ChainPipelineBuilder:
    def __init__(self, context: ChainBuildContext):
        self.context = context
        self.phases = []
        self.handoff = ChainHandoffResolver()

    def build(self, steps: Sequence[ChainStep], **kwargs) -> ChainJob:
        if not steps:
            raise ValueError("Chain pipeline requires at least one step.")
        if self.context.molecule is None:
            raise ValueError("Chain pipeline requires a molecule.")

        prepared, project_names = self._validate_pipeline(steps)
        self._prepared = prepared
        self._project_names = project_names
        molecule = self.context.molecule
        self._charge = _first_set(self.context.charge, molecule.charge)
        self._multiplicity = _first_set(
            self.context.multiplicity, molecule.multiplicity
        )
        self.phases = []

        def make_factory(index):
            def factory():
                return self._jobs_for_phase(index)

            return factory

        for index, step in enumerate(prepared):
            phase_name = f"{index:02d}_{step.program}_{step.job}"
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

        return ChainJob(
            molecule=molecule,
            label=self.context.label,
            jobrunner=self.context.jobrunner,
            skip_completed=self.context.skip_completed,
            phases=self.phases,
            programs=tuple(dict.fromkeys(step.program for step in prepared)),
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
        if self._needs_conformers_file(step):
            previous_jobs = self.handoff._completed_previous_jobs(
                index, self.phases
            )
            if previous_jobs is None:
                return []
            crest_job = self.handoff.crest_source(previous_jobs)
            if step.program == "grouper":
                return self._grouper_jobs(
                    step, crest_job, index, num_conformers
                )
            return self._gaussian_crest_jobs(step, crest_job, num_conformers)

        handoff = self.handoff.resolve(
            index, self.phases, self.context.molecule, num_conformers
        )
        if handoff is None or not handoff.molecules:
            return []
        return self._program_jobs(
            step, handoff.molecules, index, handoff.label_kind
        )

    def _program_jobs(self, step, molecules, index, kind):
        """Build one program command per structure.

        The chain writes the structure and supplies project, charge, and
        multiplicity. The program command builds the job.
        """
        command = PROGRAM_ADAPTERS[step.program].command_loader()
        count = len(molecules)
        jobs = []
        for ordinal, molecule in enumerate(molecules, start=1):
            child_label = _child_label(
                self.context.label, index, step, ordinal, count, kind
            )
            filename = os.path.abspath(f"{child_label}.xyz")
            molecule.write_xyz(filename, mode="w")
            argv = _program_argv(
                step,
                filename,
                self._project_names[step.program],
                child_label,
                self._charge,
                self._multiplicity,
                self.context.skip_completed,
            )
            result = _invoke_cli_command(command, argv, self.context.jobrunner)
            jobs.extend(_jobs_from_cli_result(result, molecule))
        return jobs

    def _grouper_jobs(self, step, crest_job, index, num_conformers):
        conformers_file = self.handoff.conformers_file(
            crest_job, num_conformers
        )
        if conformers_file is None:
            logger.info(
                "CREST conformers file is not available yet; "
                "waiting before grouper."
            )
            return []
        frames = (
            Molecule.from_filepath(
                filepath=conformers_file, index=":", return_list=True
            )
            or []
        )
        molecules = [
            _with_source_charge(molecule, crest_job.molecule)
            for molecule in frames
        ]
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

    def _gaussian_crest_jobs(self, step, crest_job, num_conformers):
        conformers_file = self.handoff.conformers_file(
            crest_job, num_conformers
        )
        if conformers_file is None:
            logger.info(
                "CREST conformers file is not available yet; "
                "waiting before gaussian crest."
            )
            return []
        argv = _program_argv(
            step,
            conformers_file,
            self._project_names[step.program],
            self.context.label,
            self._charge,
            self._multiplicity,
            self.context.skip_completed,
        )
        command = PROGRAM_ADAPTERS["gaussian"].command_loader()
        return [_invoke_cli_command(command, argv, self.context.jobrunner)]

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
    are that command's options. ``crest:conformers`` accepts ``-nc N``
    after the job name; that selects the lowest N conformers for the
    next step and is not passed to CREST. Omitting it selects the
    lowest conformer.

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
            ``gaussian``, ``orca``, or ``grouper``).
        job (str): Job name (e.g. ``opt``, ``sp``, ``irc``). A grouper
            job is a strategy registered on the grouper command. A
            program job is a subcommand of that program, other than
            pKa, Fukui, redox, and reaction.

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
    previous completed child. ``crest:conformers -nc N`` passes the lowest
    N conformers to the next step. Omitting ``-nc`` passes the lowest
    conformer. ``gaussian:crest`` and ``grouper`` read that same
    selection from the preceding CREST step, or every frame of
    ``crest_conformers.xyz`` when ``-nc`` is omitted. The step after
    Gaussian crest receives that job's lowest-energy optimized conformer.
    Charge and multiplicity from ``charge`` / ``multiplicity`` (and
    otherwise from ``molecule``) are passed to each program command. Each
    program step is that program's own subcommand. The chain supplies the
    geometry.

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
        ValueError: If there are no steps, ``molecule`` is missing, a
            step's ``(program, job)`` pair is not a pipeline step, or
            ``gaussian:crest`` / ``grouper`` does not follow
            ``crest:conformers``.
    """
    context = ChainBuildContext(
        chain_settings=chain_settings,
        molecule=molecule,
        label=label,
        charge=charge,
        multiplicity=multiplicity,
        jobrunner=jobrunner,
        skip_completed=skip_completed,
    )
    return ChainPipelineBuilder(context).build(steps, **kwargs)
