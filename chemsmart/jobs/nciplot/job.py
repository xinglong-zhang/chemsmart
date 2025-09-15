"""
NCIPLOT job implementation.

This module contains the NCIPLOTJob class for running non-covalent
interaction job for creating dens.cube and grad.cube files
using the NCIPLOT program.
"""

import logging
import os
from typing import Type

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.job import Job
from chemsmart.jobs.nciplot.settings import NCIPLOTJobSettings

logger = logging.getLogger(__name__)


class NCIPLOTJob(Job):
    """
    NCIPLOT job for non-covalent interaction job creating dens.cube 
    and grad.cube files.

    This class handles the setup and execution of NCIPLOT calculations for
    creat non-covalent interaction cube files in molecular systems.

    Attributes:
        PROGRAM (str): Program identifier ('NCIPLOT').
        TYPE (str): Job type identifier ('nciplot').
        filenames (list[str] | None): Input filenames (.xyz/.wfn/.wfx).
        molecule (Molecule | None): Molecular structure for analysis.
        settings (NCIPLOTJobSettings): Job configuration options.
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend for running the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    PROGRAM = "NCIPLOT"
    TYPE = "nciplot"

    def __init__(
        self,
        filenames=None,  # accepts multiple files
        molecule=None,
        settings=None,
        label=None,
        jobrunner=None,
        **kwargs,
    ):
        """
        Initialize NCIPLOTJob.

        Args:
            filenames: List of input filenames (.xyz, .wfn, .wfx files)
            molecule: Molecule object for structure-based analysis
            settings: NCIPLOTJobSettings instance
            label: Job label for identification
            jobrunner: Job runner instance
            **kwargs: Additional keyword arguments
        """
        super().__init__(
            molecule=molecule, label=label, jobrunner=jobrunner, **kwargs
        )

        # Validate settings type
        if settings is not None and not isinstance(
            settings, NCIPLOTJobSettings
        ):
            raise ValueError(
                f"Settings must be instance of {NCIPLOTJobSettings} for "
                f"{self}, but is {settings} instead!"
            )

        # Either filenames or molecule must be provided, but not both
        logger.debug(f"filenames: {filenames}, molecule: {molecule}")
        if filenames is not None and molecule is not None:
            raise ValueError(
                "Either filenames or molecule must be provided, but not both!"
            )
        elif filenames is None and molecule is None:
            raise ValueError(
                "Either filenames or molecule must be provided, but both "
                "are None!"
            )
        elif filenames is not None and molecule is None:
            if len(filenames) == 0:
                raise ValueError(
                    "No filenames provided for NCIPLOT job. Please provide "
                    "at least one file."
                )
            elif len(filenames) == 1:
                label = filenames[0].split(".")[0] if label is None else label
            else:
                # Add filenames together
                label = (
                    "_and_".join(
                        [filename.split(".")[0] for filename in filenames]
                    )
                    if label is None
                    else label
                )

        elif molecule is not None and filenames is None:
            if not isinstance(molecule, Molecule):
                raise ValueError(
                    f"Molecule must be instance of Molecule for {self}, but "
                    f"is {molecule} instead!"
                )

        self.molecule = molecule.copy() if molecule is not None else None
        self.settings = settings.copy()
        self.filenames = filenames
        self.label = label

    @classmethod
    def settings_class(cls) -> Type[NCIPLOTJobSettings]:
        """
        Return the settings class for this job type.

        Returns:
            NCIPLOTJobSettings class
        """
        return NCIPLOTJobSettings

    @property
    def inputfile(self):
        """
        Get the input file path for the NCIPLOT job.

        Returns:
            str: Absolute path to the input file
        """
        inputfile = self.label + ".nci"
        return os.path.join(self.folder, inputfile)

    @property
    def outputfile(self):
        """
        Get the output file path for the NCIPLOT job.

        Returns:
            str: Absolute path to the output file
        """
        outputfile = self.label + ".nciout"
        return os.path.join(self.folder, outputfile)

    @property
    def errfile(self):
        """
        Get the error file path for the NCIPLOT job.

        Returns:
            str: Absolute path to the error file
        """
        errfile = self.label + ".ncierr"
        return os.path.join(self.folder, errfile)

    def _backup_files(self, **kwargs):
        """
        Create backup of output files.

        Args:
            **kwargs: Additional arguments for backup operation
        """
        folder = self._create_backup_folder_name()
        logger.debug(f"Backing up files to folder: {folder}")
        self.backup_file(self.outputfile, folder=folder, **kwargs)

    def _output(self):
        """
        Get the output file path if it exists.

        Returns:
            str or None: Absolute path to output file if exists, None otherwise
        """
        if not os.path.exists(self.outputfile):
            logger.debug(f"Output file not found: {self.outputfile}")
            return None
        return os.path.abspath(self.outputfile)

    def _job_is_complete(self):
        """
        Check if the job is complete by checking the output file.

        Returns:
            bool: True if job completed successfully, False otherwise
        """
        file_exists = os.path.exists(self.outputfile)
        end_in_last_line_of_file = False
        
        if file_exists:
            with open(self.outputfile) as f:
                lines = f.readlines()
                if lines:
                    end_in_last_line_of_file = lines[-1].startswith("End")
        return file_exists and end_in_last_line_of_file

    def _run(self, **kwargs):
        """
        Run the NCIPLOT job.

        Args:
            **kwargs: Additional arguments for job execution
        """
        logger.info(
            f"Running NCIPLOTJob {self} with jobrunner {self.jobrunner}"
        )
        self.jobrunner.run(self, **kwargs)
