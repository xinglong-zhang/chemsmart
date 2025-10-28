"""
Gaussian computational chemistry job classes and base functionality.

This module provides the core job classes for Gaussian calculations,
including the base GaussianJob class and specialized variants for
different input types. These classes handle job creation, execution,
file management, and output processing for Gaussian computational
chemistry calculations.

Key classes:
- GaussianJob: Base class for all Gaussian calculations
- GaussianComJob: For running existing .com input files
- GaussianGeneralJob: General-purpose Gaussian job runner
"""

import logging
import os
from typing import Type

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.settings import GaussianJobSettings
from chemsmart.jobs.job import Job
from chemsmart.jobs.runner import JobRunner

logger = logging.getLogger(__name__)


class GaussianJob(Job):
    """
    Base class for all Gaussian computational chemistry jobs.

    Provides common functionality for Gaussian calculations including
    file management, job execution, output processing, and molecular
    structure handling. All specialized Gaussian job types inherit
    from this base class.

    Attributes:
        PROGRAM (str): The computational program name ('Gaussian').
        molecule (Molecule): The molecular structure for calculation.
        settings (GaussianJobSettings): Job configuration parameters.
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend (inherited from Job).
        folder (str): Working directory where files are written (inherited).
    """

    PROGRAM = "Gaussian"

    def __init__(
        self, molecule, settings=None, label=None, jobrunner=None, **kwargs
    ):
        """
        Initialize a Gaussian job with molecule and calculation settings.

        Validates input parameters and sets up the job configuration
        for Gaussian calculations. Handles molecular structure copying
        and settings validation.

        Args:
            molecule (Molecule): Molecular structure for calculation.
            settings (GaussianJobSettings): Job configuration (required).
            label (str, optional): Job identifier for file naming.
            jobrunner (JobRunner, optional): Job execution handler.
            **kwargs: Additional keyword arguments for parent class.

        Raises:
            ValueError: If settings or molecule types are invalid.
        """
        super().__init__(
            molecule=molecule, label=label, jobrunner=jobrunner, **kwargs
        )

        if not isinstance(settings, GaussianJobSettings):
            raise ValueError(
                f"Settings must be instance of {GaussianJobSettings} for {self}, but is {settings} instead!"
            )

        if not isinstance(molecule, Molecule):
            raise ValueError(
                f"Molecule must be instance of Molecule for {self}, but is {molecule} instead!"
            )

        self.molecule = molecule.copy() if molecule is not None else None
        self.settings = settings.copy()

        if label is None:
            label = molecule.get_chemical_formula(empirical=True)
        self.label = label

    @classmethod
    def settings_class(cls) -> Type[GaussianJobSettings]:
        """
        Get the settings class used by this job type.

        Returns the appropriate settings class for configuring
        Gaussian job parameters and calculation options.

        Returns:
            Type[GaussianJobSettings]: Settings class for this job type.
        """
        return GaussianJobSettings

    @property
    def inputfile(self):
        """
        Get the path to the Gaussian input file (.com).

        Constructs the full path to the input file using the job
        label and folder location.

        Returns:
            str: Full path to the Gaussian input file.
        """
        inputfile = self.label + ".com"
        return os.path.join(self.folder, inputfile)

    @property
    def outputfile(self):
        """
        Get the path to the Gaussian output file (.log).

        Constructs the full path to the output file using the job
        label and folder location.

        Returns:
            str: Full path to the Gaussian output file.
        """
        outputfile = self.label + ".log"
        return os.path.join(self.folder, outputfile)

    @property
    def chkfile(self):
        """
        Get the path to the Gaussian checkpoint file (.chk).

        Constructs the full path to the checkpoint file using the job
        label and folder location.

        Returns:
            str: Full path to the Gaussian checkpoint file.
        """
        chkfile = self.label + ".chk"
        return os.path.join(self.folder, chkfile)

    @property
    def errfile(self):
        """
        Get the path to the Gaussian error file (.err).

        Constructs the full path to the error file using the job
        label and folder location.

        Returns:
            str: Full path to the Gaussian error file.
        """
        errfile = self.label + ".err"
        return os.path.join(self.folder, errfile)

    def _backup_files(self, backup_chk=False, **kwargs):
        """
        Create backup copies of important job files.

        Backs up input and output files, and optionally the checkpoint
        file, to a timestamped backup directory for data preservation.

        Args:
            backup_chk (bool): Whether to backup checkpoint files.
            **kwargs: Additional arguments passed to backup_file method.
        """
        folder = self._create_backup_folder_name()
        self.backup_file(self.inputfile, folder=folder, **kwargs)
        self.backup_file(self.outputfile, folder=folder, **kwargs)
        if backup_chk:
            self.backup_file(self.chkfile, folder=folder, **kwargs)

    def _output(self):
        """
        Create and return a Gaussian output parser object.

        Attempts to parse the output file using appropriate Gaussian
        output readers, with fallback for periodic boundary conditions.

        Returns:
            Gaussian16Output or Gaussian16OutputWithPBC or None:
                Parsed output object or None if file doesn't exist
                or parsing fails.
        """
        if not os.path.exists(self.outputfile):
            return None

        try:
            from chemsmart.io.gaussian.output import Gaussian16Output

            return Gaussian16Output(filename=self.outputfile)
        except AttributeError:
            return None
        except ValueError:
            from chemsmart.io.gaussian.output import Gaussian16OutputWithPBC

            return Gaussian16OutputWithPBC(filename=self.outputfile)

    def _run(self, **kwargs):
        """
        Execute the Gaussian job using the assigned jobrunner.

        Logs the job execution details and delegates the actual
        running to the configured jobrunner instance.

        Args:
            **kwargs: Additional keyword arguments passed to jobrunner.
        """
        logger.info(
            f"Running GaussianJob {self} with jobrunner {self.jobrunner}"
        )
        self.jobrunner.run(self, **kwargs)

    @classmethod
    def from_filename(
        cls,
        filename,
        settings=None,
        index="-1",
        label=None,
        jobrunner=None,
        keywords=("charge", "multiplicity"),
        **kwargs,
    ):
        """
        Create a GaussianJob from a file containing molecular data.

        Reads molecular structures from various file formats and creates
        a Gaussian job with the specified settings and configuration.
        Supports multiple molecules and index-based selection.

        Args:
            filename (str): Path to file containing molecular data.
            settings (GaussianJobSettings, optional): Job configuration.
            index (str): Molecule index selection (default "-1").
            label (str, optional): Job identifier for file naming.
            jobrunner (JobRunner, optional): Job execution handler.
            keywords (tuple): Keywords to extract from molecule data.
            **kwargs: Additional keyword arguments.

        Returns:
            GaussianJob: Configured Gaussian job instance.
        """
        logger.info(f"Reading molecules from file: {filename}.")
        molecules = Molecule.from_filepath(
            filepath=filename, index=":", return_list=True
        )
        logger.info(f"Num of molecules read: {len(molecules)}.")

        # use 1-based indexing for the index
        from chemsmart.utils.utils import return_objects_from_string_index

        molecules = return_objects_from_string_index(molecules, index)

        # Create jobrunner if not provided
        if jobrunner is None:
            jobrunner = JobRunner.from_job(
                cls(
                    molecule=molecules,
                    settings=settings,
                    label=label,
                    jobrunner=None,
                    **kwargs,
                ),
                server=kwargs.get("server"),
                scratch=kwargs.get("scratch"),
                fake=kwargs.get("fake", False),
                **kwargs,
            )

        return cls(
            molecule=molecules,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )

    @classmethod
    def from_pubchem(
        cls, identifier, settings=None, label=None, jobrunner=None, **kwargs
    ):
        """
        Create a GaussianJob from a PubChem molecular identifier.

        Downloads molecular structure data from PubChem database
        and creates a Gaussian job with the specified configuration.

        Args:
            identifier (str): PubChem compound identifier (CID, name, etc.).
            settings (GaussianJobSettings, optional): Job configuration.
            label (str, optional): Job identifier for file naming.
            jobrunner (JobRunner, optional): Job execution handler.
            **kwargs: Additional keyword arguments.

        Returns:
            GaussianJob: Configured Gaussian job instance.
        """
        molecules = Molecule.from_pubchem(identifier=identifier)

        # Create jobrunner if not provided
        if jobrunner is None:
            jobrunner = JobRunner.from_job(
                cls(
                    molecule=molecules,
                    settings=settings,
                    label=label,
                    jobrunner=None,
                    **kwargs,
                ),
                server=kwargs.get("server"),
                scratch=kwargs.get("scratch"),
                fake=kwargs.get("fake", False),
                **kwargs,
            )

        return cls(
            molecule=molecules,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )

    @classmethod
    def from_jobtype(
        cls,
        jobtype,
        molecule,
        settings=None,
        label=None,
        jobrunner=None,
        **kwargs,
    ):
        """
        Create a specific GaussianJob subclass based on job type.

        Factory method that creates appropriate job instances based
        on the specified calculation type (e.g., 'opt', 'com', 'g16').

        Args:
            jobtype (str): Type of calculation ('opt', 'com', 'g16').
            molecule (Molecule): Molecular structure for calculation.
            settings (GaussianJobSettings, optional): Job configuration.
            label (str, optional): Job identifier for file naming.
            jobrunner (JobRunner, optional): Job execution handler.
            **kwargs: Additional keyword arguments.

        Returns:
            GaussianJob: Appropriate job subclass instance.

        Raises:
            ValueError: If jobtype is not recognized.
        """
        if jobtype.lower() == "opt":
            from chemsmart.jobs.gaussian.opt import GaussianOptJob

            logger.debug(f"Creating GaussianOptJob from jobtype: {jobtype}")

            # Create jobrunner if not provided
            if jobrunner is None:
                jobrunner = JobRunner.from_job(
                    GaussianOptJob(
                        molecule=molecule,
                        settings=settings,
                        label=label,
                        jobrunner=None,
                        **kwargs,
                    ),
                    server=kwargs.get("server"),
                    scratch=kwargs.get("scratch"),
                    fake=kwargs.get("fake", False),
                    **kwargs,
                )

            return GaussianOptJob(
                molecule=molecule,
                settings=settings,
                label=label,
                jobrunner=jobrunner,
                **kwargs,
            )

        elif jobtype.lower() == "com":
            # Create jobrunner if not provided
            if jobrunner is None:
                jobrunner = JobRunner.from_job(
                    GaussianComJob(
                        molecule=molecule,
                        settings=settings,
                        label=label,
                        jobrunner=None,
                        **kwargs,
                    ),
                    server=kwargs.get("server"),
                    scratch=kwargs.get("scratch"),
                    fake=kwargs.get("fake", False),
                    **kwargs,
                )

            return GaussianComJob(
                molecule=molecule,
                settings=settings,
                label=label,
                jobrunner=jobrunner,
                **kwargs,
            )
        elif jobtype.lower() == "g16":
            # Create jobrunner if not provided
            if jobrunner is None:
                jobrunner = JobRunner.from_job(
                    GaussianGeneralJob(
                        molecule=molecule,
                        settings=settings,
                        label=label,
                        jobrunner=None,
                        **kwargs,
                    ),
                    server=kwargs.get("server"),
                    scratch=kwargs.get("scratch"),
                    fake=kwargs.get("fake", False),
                    **kwargs,
                )

                return GaussianGeneralJob(
                    molecule=molecule,
                    settings=settings,
                    label=label,
                    jobrunner=jobrunner,
                    **kwargs,
                )
            else:
                raise ValueError(f"Invalid job type: {jobtype}")


class GaussianComJob(GaussianJob):
    """
    Gaussian job for running existing .com input files as-is.

    Specialized job class that takes pre-written Gaussian input files
    and runs them without modification. Useful for running existing
    input files or custom calculations with specific parameters.

    Attributes:
        TYPE (str): Job type identifier ('g16com').
    """

    TYPE = "g16com"

    def __init__(
        self, molecule, settings=None, label=None, jobrunner=None, **kwargs
    ):
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )

    @classmethod
    def from_filename(
        cls, filename, settings=None, label=None, jobrunner=None, **kwargs
    ):
        """
        Create a GaussianComJob from an existing .com input file.

        Reads a Gaussian input file and creates a job that will run
        the file exactly as written, preserving all input parameters
        and molecular coordinates.

        Args:
            filename (str): Path to the .com input file.
            settings (GaussianJobSettings, optional): Additional settings.
            label (str, optional): Job identifier (defaults to filename).
            jobrunner (JobRunner, optional): Job execution handler.
            **kwargs: Additional keyword arguments.

        Returns:
            GaussianComJob: Job configured to run the input file.
        """
        # job.label as the filename (without extension) used
        if label is None:
            label = os.path.splitext(os.path.basename(filename))[0]

        # get input lines
        from chemsmart.io.gaussian.input import Gaussian16Input

        g16com_file = Gaussian16Input(filename=filename)
        input_lines = g16com_file.content_lines_string

        # get settings from file
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings

        # store file lines in settings
        settings = GaussianJobSettings.from_filepath(filename)
        settings.input_string = input_lines

        logger.debug(
            f"Supplied file {filename} settings are: \n{settings.__dict__}"
        )
        logger.debug(f"Writing input lines: \n{input_lines}")

        # Create jobrunner if not provided
        if jobrunner is None:
            jobrunner = JobRunner.from_job(
                cls(
                    molecule=None,
                    settings=settings,
                    label=label,
                    jobrunner=None,
                    **kwargs,
                ),
                server=kwargs.get("server"),
                scratch=kwargs.get("scratch"),
                fake=kwargs.get("fake", False),
                **kwargs,
            )

        return cls(
            molecule=None,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )


class GaussianGeneralJob(GaussianJob):
    """
    General-purpose Gaussian job class for standard calculations.

    Subclasses GaussianJob to prevent recursive loops when other
    job types need to run general Gaussian calculations. For example,
    prevents infinite recursion in specialized jobs like GaussianCrestOptJob
    that need to run standard Gaussian calculations internally.

    Attributes:
        TYPE (str): Job type identifier ('g16').
    """

    TYPE = "g16"

    def __init__(
        self, molecule, settings=None, label=None, jobrunner=None, **kwargs
    ):
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
