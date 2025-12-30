"""
xTB job implementation.

This module contains the main xTB job classes for running quantum chemistry
calculations using the xTB program package.
"""

import logging
import os
from typing import Type

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.job import Job
from chemsmart.jobs.xtb.settings import XTBJobSettings
from chemsmart.utils.utils import string2index_1based

logger = logging.getLogger(__name__)


class XTBJob(Job):
    """
    Base class for xTB jobs.

    This class handles the setup, execution, and management of xTB calculations.

    Attributes:
        PROGRAM (str): The program name ('xtb').
        molecule (Molecule): The molecule object to be calculated.
        settings (XTBJobSettings): The settings for the xTB calculation.
        label (str): The label for the job, used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
    """

    PROGRAM = "xtb"

    def __init__(
        self, molecule, settings=None, label=None, jobrunner=None, **kwargs
    ):
        """
        Initialize XTBJob.

        Args:
            molecule: Molecule object for the calculation
            settings: XTBJobSettings instance
            label: Job label for identification
            jobrunner: Job runner instance
            **kwargs: Additional keyword arguments
        """
        if label is None:
            label = molecule.get_chemical_formula(empirical=True)
        self.label = label

        super().__init__(
            molecule=molecule, label=label, jobrunner=jobrunner, **kwargs
        )
        if not isinstance(settings, XTBJobSettings):
            raise ValueError(
                f"Settings must be instance of {XTBJobSettings} for {self}, but is {settings} instead!"
            )
        if not isinstance(molecule, Molecule):
            raise ValueError(
                f"Molecule must be instance of Molecule for {self}, but is {molecule} instead!"
            )
        molecule = molecule.copy()
        settings = settings.copy()

        self.settings = settings
        self.molecule = molecule
        self.label = label

    @classmethod
    def settings_class(cls) -> Type[XTBJobSettings]:
        """
        Return the settings class for this job type.

        Returns:
            XTBJobSettings class
        """
        return XTBJobSettings

    @property
    def xyzfile(self):
        """
        Get the xyz file path for the xTB job.

        Returns:
            str: Absolute path to the xyz file (.xyz).
        """
        xyzfile = self.label + ".xyz"
        return os.path.join(self.folder, xyzfile)

    @property
    def inputfile(self):
        """
        Get the input file path for the xTB job.

        Returns:
            str: Absolute path to the input file (.inp).
        """
        inputfile = self.label + ".inp"
        return os.path.join(self.folder, inputfile)

    @property
    def outputfile(self):
        """
        Get the output file path for the xTB job.

        Returns:
            str: Absolute path to the output file
        """
        outputfile = self.label + ".out"
        return os.path.join(self.folder, outputfile)

    @property
    def errfile(self):
        """
        Get the error file path for the xTB job.

        Returns:
            str: Absolute path to the error file
        """
        errfile = self.label + ".err"
        return os.path.join(self.folder, errfile)

    @property
    def restartfile(self):
        """
        Get the restart file path for the xTB job.

        Returns:
            str: Absolute path to the restart file (xtbrestart).
        """
        return os.path.join(self.folder, "xtbrestart")

    def _backup_files(self, backup_restart=False, **kwargs):
        """
        Backup input and output files to a backup folder.

        Args:
            backup_restart (bool): Whether to backup restart files.
            **kwargs: Additional arguments passed to backup_file.
        """
        folder = self._create_backup_folder_name()
        self.backup_file(self.xyzfile, folder=folder, **kwargs)
        self.backup_file(self.inputfile, folder=folder, **kwargs)
        self.backup_file(self.outputfile, folder=folder, **kwargs)
        if backup_restart:
            self.backup_file(self.restartfile, folder=folder, **kwargs)

    def _output(self):
        """
        Get the output object for the job.

        Returns:
            XTBOutput: The output object if the output file exists, else None.
        """
        if not os.path.exists(self.outputfile):
            return None

        try:
            from chemsmart.io.xtb.output import XTBOutput

            return XTBOutput(folder=self.folder)
        except Exception:
            logger.error(f"Error reading output file: {self.outputfile}")
            return None

    def _run(self, **kwargs):
        """
        Run the job using the assigned jobrunner.

        Args:
            **kwargs: Additional arguments passed to the jobrunner.
        """
        self.jobrunner.run(self, **kwargs)

    @classmethod
    def from_filename(
        cls,
        filename,
        settings=None,
        index="-1",
        label=None,
        keywords=("charge", "multiplicity"),
        **kwargs,
    ):
        """
        Create xTB Job from a file containing molecular structures.

        Args:
            filename: Path to molecular structure file
            settings: xTB job settings (optional)
            index: Molecule index to use (default: "-1" for last)
            label: Job label (optional)
            keywords: Settings keywords to use
            **kwargs: Additional arguments

        Returns:
            XTBJob: Configured xTB job instance
        """
        # get all molecule in a file and give the result as a list
        logger.info(f"Reading images from file: {filename}.")
        molecules = Molecule.from_filepath(
            filepath=filename, index=":", return_list=True
        )
        logger.info(f"Num of images read: {len(molecules)}.")

        # Select specified molecule by index
        molecules = molecules[string2index_1based(index)]

        # only supply last molecule in some jobs; but require all molecule in others e.g., dias job
        return cls(
            molecule=molecules,
            settings=settings,
            label=label,
            **kwargs,
        )

    @classmethod
    def from_pubchem(cls, identifier, settings=None, label=None, **kwargs):
        """
        Create xTB job from PubChem molecular database.

        Args:
            identifier: PubChem compound identifier (name, CID, etc.)
            settings: xTB job settings (optional)
            label: Job label (optional)
            **kwargs: Additional arguments

        Returns:
            XTBJob: Configured xTB job instance
        """
        molecules = Molecule.from_pubchem(identifier=identifier)
        return cls(
            molecule=molecules,
            settings=settings,
            label=label,
            **kwargs,
        )

    @classmethod
    def from_jobtype(
        cls, jobtype, molecule, settings=None, label=None, **kwargs
    ):
        """
        Create specific xTB job type from jobtype string.

        Args:
            jobtype: Type of xTB job ('opt').
            molecule: Molecule object for the calculation
            settings: xTB job settings (optional)
            label: Job label (optional)
            **kwargs: Additional arguments

        Returns:
            XTBJob: Specific xTB job subclass instance

        Raises:
            ValueError: If the job type is invalid.
        """
        if jobtype.lower() == "opt":
            from chemsmart.jobs.xtb.opt import XTBOptJob

            logger.debug(f"Creating XTBOptJob from jobtype: {jobtype}")

            return XTBOptJob(
                molecule=molecule,
                settings=settings,
                label=label,
                **kwargs,
            )
        else:
            raise ValueError(f"Invalid job type: {jobtype}")

    def _determine_folder(self):
        """
        Determine the folder based on the current working directory
        where the job is submitted.
        """
        # Get the current working directory at runtime
        cwd = os.getcwd()
        folder = os.path.abspath(cwd)

        # XTB jobs must be run in their own directory
        if self.label:
            folder = os.path.join(folder, self.label)
            if not os.path.exists(folder):
                os.makedirs(folder)

        return folder

    def base_folder(self):
        """
        Base folder for the job, where input and output files are located.
        This is typically the folder where the job was submitted from.
        """
        return self._determine_folder()
