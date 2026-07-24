import logging
import os
import re
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Type

from chemsmart.jobs.iterate.settings import IterateJobSettings
from chemsmart.jobs.job import Job
from chemsmart.utils.io import clean_label

if TYPE_CHECKING:
    from chemsmart.jobs.iterate.runner import IterateRunSummary

logger = logging.getLogger(__name__)


class IterateExecutionError(RuntimeError):
    """Raised when an Iterate run completes with a non-zero result."""

    def __init__(self, summary: "IterateRunSummary"):
        self.summary = summary
        if summary.summary_write_error:
            detail = f"report write failed: {summary.summary_write_error}"
        elif summary.error_codes:
            detail = f"error codes: {', '.join(summary.error_codes)}"
        else:
            detail = "see the Iterate run report for details"
        super().__init__(
            f"Iterate terminated with exit code {summary.exit_code} ({detail})."
        )


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
        outputfile=None,
        separate_outputs=False,
        output_directory=None,
        command_line=None,
        show_worker_logs=False,
        skip_completed=True,
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
            Output filename (without .xyz extension). When omitted, use
            ``<configuration_stem>_iterate.xyz``.
        separate_outputs : bool, optional
            If True, save each structure as a separate file. Default False.
        output_directory : str, optional
            Directory for output files if separate_outputs is True.
        command_line : str, optional
            The original CLI command, recorded verbatim in the run report.
        show_worker_logs : bool, optional
            Preserve worker and RDKit logs during debug-mode local runs.
        skip_completed : bool, optional
            Skip execution when the normal report and all expected outputs
            already exist. Default True.
        **kwargs
            Additional keyword arguments
        """
        logger.debug("Initializing IterateJob")
        self.settings = (
            settings if settings is not None else IterateJobSettings()
        )
        label = self._build_label(settings=self.settings)

        super().__init__(
            molecule=None,
            label=label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            **kwargs,
        )

        self.nprocs = nprocs
        self.timeout = timeout
        self._outputfile = outputfile if outputfile is not None else label
        self.separate_outputs = separate_outputs
        self.output_directory = output_directory
        self.command_line = command_line
        self.show_worker_logs = show_worker_logs

    @staticmethod
    def _build_label(settings) -> str:
        """Build ``<config_stem>_iterate`` following other CLI jobs."""
        config_file = getattr(settings, "config_file", None)
        stem = Path(config_file).stem if config_file else "iterate"
        return clean_label(f"{stem}_iterate")

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

    @property
    def reportfile(self) -> str:
        """Return the absolute path of the Iterate run report."""
        filename = f"{self.label}.out"

        if self.separate_outputs:
            report_directory = self.output_directory or self.folder
        else:
            output_directory = os.path.dirname(self.outputfile)
            report_directory = output_directory or self.folder

        if not os.path.isabs(report_directory):
            report_directory = os.path.join(self.folder, report_directory)
        return os.path.join(report_directory, filename)

    def _resolve_output_path(self, path: str) -> str:
        if os.path.isabs(path):
            return path
        return os.path.join(self.folder, path)

    @staticmethod
    def _report_value(lines, key):
        prefix = f"{key}:"
        for line in lines:
            stripped = line.strip()
            if stripped.startswith(prefix):
                return stripped.split(":", 1)[1].strip()
        return None

    @staticmethod
    def _separate_output_paths(lines) -> list[str]:
        in_output_section = False
        in_mapping = False
        output_paths = []
        pattern = re.compile(r"^\s{3}.+?\s{2,}(.+\.xyz)\s*$")

        for line in lines:
            stripped = line.strip()
            if stripped == "OUTPUT STRUCTURES":
                in_output_section = True
                continue
            if in_output_section and stripped == "FINAL SUMMARY":
                break
            if in_output_section and stripped == "Separate output mapping:":
                in_mapping = True
                continue
            if not in_mapping:
                continue
            match = pattern.match(line)
            if match is not None:
                output_paths.append(match.group(1).strip())
        return output_paths

    def is_complete(self) -> bool:
        """Return whether a matching normal report and outputs exist."""
        try:
            if not os.path.isfile(self.reportfile):
                return False
            with open(self.reportfile, "r", errors="replace") as handle:
                lines = handle.readlines()

            nonempty_lines = [line.strip() for line in lines if line.strip()]
            if not nonempty_lines or not nonempty_lines[-1].startswith(
                "Normal termination of CHEMSMART Iterate"
            ):
                return False

            structures_written = self._report_value(
                lines, "Structures written"
            )
            if structures_written is None:
                return False
            structures_written = int(structures_written)
            if structures_written == 0:
                return True

            if not self.separate_outputs:
                return os.path.isfile(
                    self._resolve_output_path(self.outputfile)
                )

            output_paths = self._separate_output_paths(lines)
            if len(output_paths) != structures_written:
                return False
            return all(
                os.path.isfile(self._resolve_output_path(path))
                for path in output_paths
            )
        except (OSError, TypeError, ValueError) as error:
            logger.debug(f"Could not inspect Iterate completion: {error}")
            return False

    def run(
        self, progress_callback=None
    ) -> Optional["IterateRunSummary"]:
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
        IterateRunSummary or None
            Success/failure/timeout counts and the paths actually written, or
            ``None`` when a completed job is skipped.

        Raises
        ------
        IterateExecutionError
            If the completed run reports a non-zero application exit code.
        """
        if self.is_complete() and self.skip_completed:
            logger.info(f"{self} is already complete, not running.")
            return None

        logger.debug(f"Running IterateJob with {self.nprocs} process(es)")

        # Use jobrunner if available, otherwise create a temporary one
        if self.jobrunner is not None:
            runner = self.jobrunner
        else:
            from chemsmart.jobs.iterate.runner import IterateJobRunner

            runner = IterateJobRunner()

        if hasattr(runner, "show_worker_logs"):
            runner.show_worker_logs = self.show_worker_logs

        # The runner writes the report before returning. Convert its application
        # status into a Python exception so CLI/scheduler exit codes reflect an
        # Iterate failure instead of reporting a successfully exited wrapper.
        summary = runner.run(self, progress_callback=progress_callback)
        if summary.exit_code != 0:
            raise IterateExecutionError(summary)
        return summary
