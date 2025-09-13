"""
ORCA job implementation.

This module contains the main ORCA job classes for running quantum chemistry
calculations using the ORCA program package.
"""

import logging
import os
import shutil
from contextlib import suppress
from typing import Type

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.job import Job
from chemsmart.jobs.orca.settings import ORCAJobSettings
from chemsmart.jobs.runner import JobRunner
from chemsmart.utils.utils import string2index_1based

logger = logging.getLogger(__name__)


class ORCAJob(Job):
    """
    Base ORCA job class.

    This class provides the foundation for all ORCA quantum chemistry
    calculations including setup, execution, and output handling.

    Attributes:
        PROGRAM (str): Program identifier ('ORCA').
        molecule (Molecule): Molecular structure used for the calculation.
        settings (ORCAJobSettings): Configuration options for the job.
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    PROGRAM = "ORCA"

    def __init__(
        self, molecule, settings=None, label=None, jobrunner=None, **kwargs
    ):
        """
        Initialize ORCAJob.

        Args:
            molecule: Molecule object for the calculation
            settings: ORCAJobSettings instance
            label: Job label for identification
            jobrunner: Job runner instance
            **kwargs: Additional keyword arguments
        """
        super().__init__(
            molecule=molecule, label=label, jobrunner=jobrunner, **kwargs
        )

        # Validate settings type
        if not isinstance(settings, ORCAJobSettings):
            raise ValueError(
                f"Settings must be instance of {ORCAJobSettings} for {self}, "
                f"but is {settings} instead!"
            )
        
        # Validate molecule type
        if not isinstance(molecule, Molecule):
            raise ValueError(
                f"Molecule must be instance of Molecule for {self}, but is "
                f"{molecule} instead!"
            )

        # Store validated parameters
        self.molecule = molecule.copy() if molecule is not None else None
        self.settings = settings.copy()

        # Set default label if not provided
        if label is None:
            label = molecule.get_chemical_formula(empirical=True)
        self.label = label

    @classmethod
    def settings_class(cls) -> Type[ORCAJobSettings]:
        """
        Return the settings class for this job type.

        Returns:
            ORCAJobSettings class
        """
        return ORCAJobSettings

    @property
    def inputfile(self):
        """
        Get the input file path for the ORCA job.

        Returns:
            str: Absolute path to the input file
        """
        inputfile = self.label + ".inp"
        return os.path.join(self.folder, inputfile)

    @property
    def outputfile(self):
        """
        Get the output file path for the ORCA job.

        Returns:
            str: Absolute path to the output file
        """
        outputfile = self.label + ".out"
        return os.path.join(self.folder, outputfile)

    @property
    def gbwfile(self):
        """
        Get the GBW file path for the ORCA job.

        Returns:
            str: Absolute path to the GBW file
        """
        orca_gbwfile = self.label + ".gbw"
        return os.path.join(self.folder, orca_gbwfile)

    @property
    def errfile(self):
        """
        Get the error file path for the ORCA job.

        Returns:
            str: Absolute path to the error file
        """
        errfile = self.label + ".err"
        return os.path.join(self.folder, errfile)

    def _backup_files(self, backup_gbw=False, **kwargs):
        """
        Create backup of important files.

        Args:
            backup_gbw: Whether to backup the GBW file
            **kwargs: Additional arguments for backup operation
        """
        folder = self._create_backup_folder_name()
        self.backup_file(self.inputfile, folder=folder, **kwargs)
        self.backup_file(self.outputfile, folder=folder, **kwargs)
        
        if backup_gbw:
            self.backup_file(self.gbwfile, folder=folder, **kwargs)

    def _output(self):
        """
        Get the output object if the output file exists.

        Returns:
            ORCAOutput or None: Output object if file exists, None otherwise
        """
        if not os.path.exists(self.outputfile):
            logger.debug(f"Output file not found: {self.outputfile}")
            return None

        try:
            from chemsmart.io.orca.output import ORCAOutput

            return ORCAOutput(self.outputfile)
        except AttributeError as e:
            logger.error(f"Error creating ORCAOutput object: {e}")
            return None

    def _run(self, **kwargs):
        """
        Run the job using the assigned jobrunner.

        Args:
            **kwargs: Additional arguments for job execution
        """
        logger.info(f"Running ORCAJob {self} with jobrunner {self.jobrunner}")
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
        Create ORCA job from molecular structure file.

        This factory method reads molecular structures from various file
        formats and creates an ORCA job with appropriate settings.

        Args:
            filename: Path to molecular structure file
            settings: ORCA job settings (optional)
            index: Molecule index to use (default: "-1" for last)
            label: Job label (optional)
            jobrunner: Job runner instance (optional)
            keywords: Settings keywords to use
            **kwargs: Additional arguments

        Returns:
            ORCAJob: Configured ORCA job instance
        """

        # Read all molecules from file
        molecules = Molecule.from_filepath(
            filepath=filename, index=":", return_list=True
        )
        logger.info(f"Number of molecules read: {len(molecules)}")
        
        # Select specified molecule by index
        molecules = molecules[string2index_1based(index)]

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

        # Create and return job instance
        # Note: Some jobs use last molecule only, others require all molecules
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
        Create ORCA job from PubChem molecular database.

        This factory method retrieves molecular structures from PubChem
        database and creates an ORCA job for quantum chemistry calculations.

        Args:
            identifier: PubChem compound identifier (name, CID, etc.)
            settings: ORCA job settings (optional)
            label: Job label (optional)
            jobrunner: Job runner instance (optional)
            **kwargs: Additional arguments

        Returns:
            ORCAJob: Configured ORCA job instance
        """

        # Retrieve molecule from PubChem database
        molecule = Molecule.from_pubchem(identifier=identifier)

        # Create jobrunner if not provided
        if jobrunner is None:
            jobrunner = JobRunner.from_job(
                cls(
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

        return cls(
            molecule=molecule,
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
        Create specific ORCA job type from jobtype string.

        This factory method creates specialized ORCA job instances based
        on the specified job type (optimization, input file, general).

        Args:
            jobtype: Type of ORCA job ('opt', 'inp', 'orca')
            molecule: Molecule object for the calculation
            settings: ORCA job settings (optional)
            label: Job label (optional)
            jobrunner: Job runner instance (optional)
            **kwargs: Additional arguments

        Returns:
            ORCAJob: Specific ORCA job subclass instance

        Raises:
            ValueError: If invalid job type is specified
        """
        
        if jobtype.lower() == "opt":
            # Import optimization job class
            from chemsmart.jobs.orca.opt import ORCAOptJob

            logger.debug(f"Creating GaussianOptJob from jobtype: {jobtype}")

            # Create jobrunner if not provided
            if jobrunner is None:
                jobrunner = JobRunner.from_job(
                    ORCAOptJob(
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

            return ORCAOptJob(
                molecule=molecule,
                settings=settings,
                label=label,
                jobrunner=jobrunner,
                **kwargs,
            )
            
        elif jobtype.lower() == "inp":
            # Create jobrunner if not provided
            if jobrunner is None:
                jobrunner = JobRunner.from_job(
                    ORCAInpJob(
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

            return ORCAInpJob(
                molecule=molecule,
                settings=settings,
                label=label,
                jobrunner=jobrunner,
                **kwargs,
            )
            
        elif jobtype.lower() == "orca":
            # Create jobrunner if not provided
            if jobrunner is None:
                jobrunner = JobRunner.from_job(
                    ORCAGeneralJob(
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

            return ORCAGeneralJob(
                molecule=molecule,
                settings=settings,
                label=label,
                jobrunner=jobrunner,
                **kwargs,
            )
        else:
            raise ValueError(f"Invalid job type: {jobtype}")


class ORCAInpJob(ORCAJob):
    """
    ORCA input file job runner.

    This class runs ORCA calculations using existing .inp input files
    without modification, preserving all original input specifications.

    Attributes:
        TYPE (str): Job type identifier ('orcainp').
        molecule (Molecule): Molecular structure associated with the input.
        settings (ORCAJobSettings): Job settings (read from file when using
            from_filename).
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "orcainp"

    def __init__(
        self, molecule, settings=None, label=None, jobrunner=None, **kwargs
    ):
        """
        Initialize ORCA input file job.

        Args:
            molecule: Molecule object for the calculation
            settings: ORCA job settings (optional)
            label: Job label (optional)
            jobrunner: Job runner instance (optional)
            **kwargs: Additional keyword arguments
        """
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
        Create ORCA input job from existing .inp file.

        This factory method reads an existing ORCA input file and creates
        a job that will run the input file exactly as provided.

        Args:
            filename: Path to ORCA .inp input file
            settings: ORCA job settings (optional, will be read from file)
            label: Job label (optional, derived from filename)
            jobrunner: Job runner instance (optional)
            **kwargs: Additional arguments

        Returns:
            ORCAInpJob: Job instance configured from input file

        Raises:
            AssertionError: If file is not .inp format
        """

        logger.debug(f"Checking if {filename} is an ORCA input file.")

        assert filename.endswith(".inp"), "Input file must be .inp file."

        # Read molecule from input file
        molecule = Molecule.from_filepath(filepath=filename)

        # job.label as the filename (without extension) used
        if label is None:
            label = os.path.splitext(os.path.basename(filename))[0]

        # set file
        from chemsmart.io.orca.input import ORCAInput

        orca_file = ORCAInput(filename=filename)
        input_lines = orca_file.content_lines_string

        # Get settings from input file
        from chemsmart.jobs.orca.settings import ORCAJobSettings

        # Store file lines in settings for direct writing
        settings = ORCAJobSettings.from_filepath(filename)
        settings.input_string = input_lines

        logger.debug(
            f"Supplied file {filename} settings are: \n{settings.__dict__}"
        )

        logger.debug(f"Writing input lines: \n{input_lines}")

        # Create jobrunner if not provided
        if jobrunner is None:
            jobrunner = JobRunner.from_job(
                cls(
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

        # Create and return job instance
        return cls(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )

    def _run(self, **kwargs):
        """
        Run the job using the original input file.

        This method overrides the parent _run method to execute the job
        using the supplied .inp file without regenerating input content.

        Args:
            **kwargs: Additional arguments for job execution
        """

        self._copy_input()
        
        # Execute the job using the jobrunner
        self.jobrunner.run(self, **kwargs)

    def _copy_input(self):
        """
        Copy the supplied ORCA .inp file to execution directory.

        This method handles copying the input file to the appropriate
        location (scratch directory or job folder) for execution.
        """

        if (
            self.jobrunner.scratch
            and self.jobrunner.scratch_dir is not None
            and os.path.exists(self.jobrunner.scratch_dir)
        ):
            # Running job in scratch directory
            job_scratch_dir = os.path.join(
                self.jobrunner.scratch_dir, self.label
            )
            
            # Create scratch directory if needed
            with suppress(FileExistsError):
                os.mkdir(job_scratch_dir)
                logger.info(f"Folder in scratch {job_scratch_dir} is made.")
            shutil.copy(self.inputfile, job_scratch_dir)
            scratch_inputfile = os.path.join(
                job_scratch_dir, f"{self.label}.inp"
            )
            assert os.path.exists(
                scratch_inputfile
            ), f"inputfile {scratch_inputfile} is not found"
        elif self.jobrunner.scratch and self.jobrunner.scratch_dir is not None:
            # Scratch directory specified but doesn't exist
            logger.warning(
                f"Scratch directory {self.jobrunner.scratch_dir} does not "
                f"exist! Running job in {self.folder}"
            )
        else:
            # Running job in regular job folder
            logger.info(f"Running job in regular folder: {self.folder}")


class ORCAGeneralJob(ORCAJob):
    """
    General ORCA job implementation.

    This class provides a general ORCA job implementation that subclasses
    ORCAJob. It prevents recursive loops in specialized ORCA job classes
    that need to create and run general ORCA jobs internally.

    Attributes:
        TYPE (str): Job type identifier ('orcajob').
        molecule (Molecule): Molecular structure used for the calculation.
        settings (ORCAJobSettings): Configuration options for the job.
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "orcajob"

    def __init__(self, molecule, settings=None, label=None, **kwargs):
        """
        Initialize general ORCA job.

        Args:
            molecule: Molecule object for the calculation
            settings: ORCA job settings (optional)
            label: Job label (optional)
            **kwargs: Additional keyword arguments
        """
        super().__init__(
            molecule=molecule, settings=settings, label=label, **kwargs
        )
