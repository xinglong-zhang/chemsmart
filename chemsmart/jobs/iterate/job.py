import logging
from typing import TYPE_CHECKING, Type

from chemsmart.jobs.iterate.settings import IterateJobSettings
from chemsmart.jobs.job import Job

if TYPE_CHECKING:
    from chemsmart.jobs.iterate.runner import IterateRunSummary

logger = logging.getLogger(__name__)


class IterateJob(Job):
    """
    Job for generating molecular structures
    by attaching substituents to skeletons.

    This job holds configuration from
    settings (skeleton_list, substituent_list)
    and delegates all execution logic to IterateJobRunner.
    """

    PROGRAM = "Iterate"
    TYPE = "iterate"

    def __init__(
        self,
        settings=None,
        jobrunner=None,
        nprocs=1,
        timeout=120,
        outputfile="iterate_out",
        separate_outputs=False,
        output_directory=None,
        command_line=None,
        **kwargs,
    ):
        """
        Initialize an IterateJob instance.

        Parameters
        ----------
        settings : IterateJobSettings, optional
            Job settings containing skeleton_list, substituent_list, algorithm
        jobrunner : IterateJobRunner, optional
            Job runner instance for executing combinations
        nprocs : int, optional
            Number of processes for parallel execution. Default 1.
        timeout : int, optional
            Timeout in seconds for each worker
            process. Default 120 (2 minutes).
        outputfile : str, optional
            Output filename (without .xyz extension). Default is 'iterate_out'.
        separate_outputs : bool, optional
            If True, save each structure as a separate file. Default False.
        output_directory : str, optional
            Directory for output files if separate_outputs is True.
        command_line : str, optional
            The original CLI command, recorded verbatim in the run report.
        **kwargs
            Additional keyword arguments
        """
        logger.debug("Initializing IterateJob")
        # IterateJob doesn't use a single
        # molecule, pass None and a default label
        super().__init__(
            molecule=None,
            label="iterate_job",
            jobrunner=jobrunner,
            skip_completed=False,  # IterateJob doesn't support skip_completed
            **kwargs,
        )

        self.settings = (
            settings if settings is not None else IterateJobSettings()
        )
        self.nprocs = nprocs
        self.timeout = timeout
        self._outputfile = outputfile
        self.separate_outputs = separate_outputs
        self.output_directory = output_directory
        self.command_line = command_line

    @classmethod
    def settings_class(cls) -> Type[IterateJobSettings]:
        """Return the settings class for iterate jobs."""
        return IterateJobSettings

    @property
    def outputfile(self) -> str:
        """
        Get the path to the output xyz file.

        Returns
        -------
        str
            Path to the output file with .xyz extension
        """
        # Add .xyz extension if not already present
        if self._outputfile.endswith(".xyz"):
            return self._outputfile
        return f"{self._outputfile}.xyz"

    def run(self, progress_callback=None) -> "IterateRunSummary":
        """
        Run the iterate job by delegating to IterateJobRunner.

        The runner handles all execution logic including:
        - Loading molecules from settings
        - Generating combinations
        - Running combinations with multiprocessing
        - Writing results to output xyz file

        Parameters
        ----------
        progress_callback : callable, optional
            ``callback(completed, total)`` forwarded to the runner for
            display-only progress reporting. When ``None`` (the default) the
            run behaves exactly as before.

        Returns
        -------
        IterateRunSummary
            Success/failure/timeout counts and the paths actually written.
        """
        logger.debug(f"Running IterateJob with {self.nprocs} process(es)")

        # Use jobrunner if available, otherwise create a temporary one
        if self.jobrunner is not None:
            runner = self.jobrunner
        else:
            from chemsmart.jobs.iterate.runner import IterateJobRunner

            runner = IterateJobRunner()

        # Delegate all execution to runner and surface its run summary.
        return runner.run(self, progress_callback=progress_callback)
