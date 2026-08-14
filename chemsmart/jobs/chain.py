"""Ordered multi-job workflows.

``ChainJob`` runs an ordered list of ``JobPhase`` objects. Each phase contains
one or more child jobs that execute before the next phase starts. Later phases
may be built from earlier phase results via ``jobs_factory``. How those results
are interpreted (optimized geometry, wavefunction, energies, etc.) is left to
the concrete workflow. pKa and Fukui jobs use ``ChainJob``.
"""

import logging
from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.job import Job
from chemsmart.jobs.runner import decide_phase_transition, run_phase_jobs

logger = logging.getLogger(__name__)


@dataclass
class JobPhase:
    """One step in a ``ChainJob`` workflow.

    Attributes:
        name (str): Label used in logs and phase-transition messages.
        jobs (Sequence, optional): Child jobs to run. Ignored when
            ``jobs_factory`` is set.
        jobs_factory (Callable, optional): Returns the jobs for this phase.
            Use this when jobs depend on earlier phase outputs.
        require_complete (bool): If True, the chain stops when this phase
            is incomplete after running.
        stop_on_incomplete (bool): If True, stop the phase after the first
            incomplete child job.
        skip_if (Callable, optional): When this returns True, the phase is
            omitted from both execution and completeness checks.
        before_run (Callable, optional): Called once before the phase runs.
        stop_message (str, optional): Message logged when the chain stops
            because this phase is incomplete.
    """

    name: str
    jobs: Optional[Sequence] = None
    jobs_factory: Optional[Callable[[], Optional[Sequence]]] = None
    require_complete: bool = True
    stop_on_incomplete: bool = True
    skip_if: Optional[Callable[[], bool]] = None
    before_run: Optional[Callable[[], None]] = None
    stop_message: Optional[str] = None
    _resolved_jobs: Optional[list] = field(
        default=None, init=False, repr=False, compare=False
    )

    def should_skip(self) -> bool:
        if self.skip_if is None:
            return False
        return bool(self.skip_if())

    def resolve_jobs(self):
        if self._resolved_jobs is not None:
            return self._resolved_jobs
        if self.jobs_factory is not None:
            jobs = self.jobs_factory()
        else:
            jobs = self.jobs
        self._resolved_jobs = list(jobs) if jobs else []
        return self._resolved_jobs

    def is_complete(self) -> bool:
        jobs = self.resolve_jobs()
        if not jobs:
            return False
        return all(job.is_complete() for job in jobs)

    def run(self, parent_runner, logger_obj=None) -> None:
        def _before_run():
            if self.before_run is not None:
                self.before_run()
            self._resolved_jobs = None

        run_phase_jobs(
            parent_runner=parent_runner,
            jobs=None,
            jobs_factory=self.resolve_jobs,
            stop_on_incomplete=self.stop_on_incomplete,
            before_run=_before_run,
            logger_obj=logger_obj,
            phase_label=self.name,
        )


class ChainJob(Job):
    """Job that runs child jobs in ordered phases.

    Subclasses prepare child jobs and assign ``self.phases``. Execution,
    phase-to-phase gating, and overall completeness are handled here.
    """

    TYPE = "chain"

    def __init__(
        self,
        molecule,
        label,
        jobrunner=None,
        settings=None,
        phases=None,
        skip_completed=True,
        **kwargs,
    ):
        if molecule is not None and not isinstance(molecule, Molecule):
            raise ValueError(
                f"Molecule must be instance of Molecule for {self}, "
                f"but is {molecule} instead!"
            )

        super().__init__(
            molecule=molecule,
            label=label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            **kwargs,
        )

        self.molecule = molecule.copy() if molecule is not None else None
        self.settings = settings.copy() if settings is not None else None

        if label is None and self.molecule is not None:
            self.label = self.molecule.get_chemical_formula(empirical=True)

        self._phases = list(phases) if phases is not None else []

    @property
    def phases(self):
        return self._phases

    @phases.setter
    def phases(self, phases):
        self._phases = list(phases) if phases is not None else []

    def add_phase(self, phase: JobPhase) -> None:
        self._phases.append(phase)

    def phase_by_name(self, name: str) -> Optional[JobPhase]:
        for phase in self._phases:
            if phase.name == name:
                return phase
        return None

    def run_phase(self, name: str) -> None:
        """Run a single named phase without advancing the rest of the chain."""
        phase = self.phase_by_name(name)
        if phase is None:
            raise ValueError(f"No chain phase named {name!r}.")
        if phase.should_skip():
            return
        phase.run(parent_runner=self.jobrunner, logger_obj=logger)

    def _run(self, **kwargs):
        for phase in self._phases:
            if phase.should_skip():
                continue
            phase.run(parent_runner=self.jobrunner, logger_obj=logger)
            transition = decide_phase_transition(
                phase_name=phase.name,
                require_complete=phase.require_complete,
                is_complete=phase.is_complete(),
                stop_message=phase.stop_message,
            )
            if not transition.proceed:
                if transition.should_raise:
                    raise RuntimeError(transition.message)
                logger.info(transition.message)
                return

    def is_complete(self):
        active_phases = [
            phase for phase in self._phases if not phase.should_skip()
        ]
        if not active_phases:
            return False
        return all(phase.is_complete() for phase in active_phases)
