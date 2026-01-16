"""
Core thermochemistry job implementation.

This module provides the main ThermochemistryJob class for computing
thermochemical properties from quantum chemistry calculations, including
thermal corrections, entropies, and Gibbs free energies.
"""

import logging
import os
from typing import Type

from chemsmart.analysis.thermochemistry import Thermochemistry
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.job import Job
from chemsmart.jobs.runner import JobRunner
from chemsmart.jobs.thermochemistry.settings import ThermochemistryJobSettings

logger = logging.getLogger(__name__)


class ThermochemistryJob(Job):
    """
    Job for computing thermochemical properties from quantum calculations.

    This class handles the calculation of thermodynamic properties such as
    enthalpy, entropy, and Gibbs free energy from frequency calculations
    in Gaussian or ORCA output files.

    Attributes:
        PROGRAM (str): Program identifier ('Thermochemistry').
        TYPE (str): Job type identifier ('thermochemistry').
        filename (str | None): Path to the QC output file (.log/.out).
        molecule (Molecule | None): Molecular structure used for context.
        settings (ThermochemistryJobSettings): Thermochemistry configuration.
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    PROGRAM = "Thermochemistry"
    TYPE = "thermochemistry"

    def __init__(
        self,
        filename=None,
        folder=None,
        molecule=None,
        settings=None,
        label=None,
        jobrunner=None,
        **kwargs,
    ):
        """
        Initialize a thermochemistry job instance.

        Args:
            filename (str, optional): Path to quantum chemistry output file
                (for Gaussian/ORCA)
            folder (str, optional): Path to calculation folder (for xTB)
            molecule (Molecule, optional): Molecule object for the calculation
            settings (ThermochemistryJobSettings, optional): Job configuration
            label (str, optional): Custom label for the job
            jobrunner (JobRunner, optional): Job execution manager
            **kwargs: Additional keyword arguments for parent class

        Raises:
            ValueError: If settings or molecule types are invalid
        """
        logger.debug("Initializing ThermochemistryJob")
        super().__init__(
            molecule=molecule, label=label, jobrunner=jobrunner, **kwargs
        )

        # Validate that either filename or folder is provided, but not both
        if filename is None and folder is None:
            raise ValueError("Either 'filename' or 'folder' must be provided.")
        if filename is not None and folder is not None:
            raise ValueError(
                "Cannot specify both 'filename' and 'folder'. Choose one."
            )

        # Validate file extension for filename
        if filename is not None and not filename.endswith((".log", ".out")):
            raise ValueError(
                f"Unsupported file extension for '{filename}'. "
                f"Only .log or .out files are accepted."
            )

        # Validate settings type
        if settings is not None and not isinstance(
            settings, ThermochemistryJobSettings
        ):
            raise ValueError(
                f"Settings must be instance of {ThermochemistryJobSettings} for {self}, but is {settings} instead!"
            )

        # Validate molecule type
        if molecule is not None and not isinstance(molecule, Molecule):
            logger.error(f"Invalid molecule type: {type(molecule)}")
            raise ValueError(
                f"Molecule must be instance of Molecule for {self}, "
                f"but is {molecule} instead!"
            )

        # Set instance attributes
        self.molecule = molecule.copy() if molecule is not None else None
        self.settings = (
            settings.copy()
            if settings is not None
            else ThermochemistryJobSettings()
        )
        self.filename = filename
        self.folder = folder  # For xTB calculations

        # Generate label if not provided
        if label is None:
            if filename is not None:
                label = os.path.splitext(os.path.basename(filename))[0]
            elif folder is not None:
                label = os.path.basename(os.path.normpath(folder))
            elif molecule is not None:
                label = molecule.get_chemical_formula(empirical=True)
            else:
                label = "thermochemistry_job"
        self.label = label

    @classmethod
    def settings_class(cls) -> Type[ThermochemistryJobSettings]:
        """
        Return the settings class for thermochemistry jobs.

        Returns:
            Type[ThermochemistryJobSettings]: Settings class for
                                            thermochemistry jobs
        """
        return ThermochemistryJobSettings

    @property
    def inputfile(self):
        """
        Get the absolute path to the input file.

        Returns:
            str or None: Absolute path to input file if filename exists,
                        None otherwise
        """
        if self.filename:
            return os.path.abspath(self.filename)
        return None

    @property
    def outputfile(self):
        """
        Get the path to the thermochemistry output file.

        Returns:
            str: Absolute path to the thermochemistry output file
        """
        outputfile = self.label + ".dat"
        return os.path.join(self.folder, outputfile)

    @property
    def errfile(self):
        """
        Get the path to the error file for the thermochemistry job.

        Returns:
            str: Absolute path to the thermochemistry error file
        """
        errfile = self.label + ".err"
        return os.path.join(self.folder, errfile)

    def _backup_files(self, **kwargs):
        """
        Create backup copies of output files before overwriting.

        Args:
            **kwargs: Additional keyword arguments for backup operations
        """
        folder = self._create_backup_folder_name()
        self.backup_file(self.outputfile, folder=folder, **kwargs)

    def _output(self):
        """
        Get the path to the main output file if it exists.

        Returns:
            str or None: Absolute path to output file if it exists,
                        None otherwise
        """
        if not os.path.exists(self.settings.outputfile):
            return None
        return os.path.abspath(self.settings.outputfile)

    def _job_is_complete(self):
        """
        Check if the thermochemistry job has completed successfully.

        Returns:
            bool: True if output file exists, False otherwise
        """
        return os.path.exists(self.settings.outputfile)

    def _run(self, **kwargs):
        """
        Execute the thermochemistry analysis job.

        Args:
            **kwargs: Additional keyword arguments for job execution
        """
        logger.info(
            f"Running ThermochemistryJob {self} with jobrunner "
            f"{self.jobrunner}"
        )
        self.jobrunner.run(self, **kwargs)

    @classmethod
    def from_filename(
        cls,
        filename,
        settings=None,
        label=None,
        jobrunner=None,
        **kwargs,
    ):
        """
        Create a thermochemistry job from a quantum chemistry output file.

        Factory method that creates a ThermochemistryJob instance from
        a Gaussian or ORCA output file containing frequency calculations.

        Args:
            filename (str): Path to quantum chemistry output file
            settings (ThermochemistryJobSettings, optional): Job configuration
            label (str, optional): Custom label for the job
            jobrunner (JobRunner, optional): Job execution manager
            **kwargs: Additional keyword arguments for job configuration

        Returns:
            ThermochemistryJob: Configured thermochemistry job instance

        Raises:
            ValueError: If output file type is not supported
        """

        logger.info(f"Reading molecule from file: {filename}")
        molecule = Molecule.from_filepath(filename)

        # Use default settings if not provided
        if settings is None:
            settings = ThermochemistryJobSettings()

        # Generate label from filename if not provided
        if label is None:
            label = os.path.splitext(os.path.basename(filename))[0]

        # Create jobrunner if not provided
        if jobrunner is None:
            jobrunner = JobRunner.from_job(
                cls(
                    molecule=molecule,
                    settings=settings,
                    label=label,
                    filename=filename,
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
            filename=filename,
            jobrunner=jobrunner,
            **kwargs,
        )

    @classmethod
    def from_folder(
        cls,
        folder,
        settings=None,
        label=None,
        jobrunner=None,
        **kwargs,
    ):
        """
        Create a thermochemistry job from an xTB calculation folder.

        Factory method that creates a ThermochemistryJob instance from
        a folder containing xTB output files with frequency calculations.

        Args:
            folder (str): Path to xTB calculation folder
            settings (ThermochemistryJobSettings, optional): Job configuration
            label (str, optional): Custom label for the job
            jobrunner (JobRunner, optional): Job execution manager
            **kwargs: Additional keyword arguments for job configuration

        Returns:
            ThermochemistryJob: Configured thermochemistry job instance
        """
        logger.info(f"Reading molecule from folder: {folder}")
        molecule = Molecule.from_directorypath(folder, program="xtb")

        # Use default settings if not provided
        if settings is None:
            settings = ThermochemistryJobSettings()

        # Generate label from folder name if not provided
        if label is None:
            label = os.path.basename(os.path.normpath(folder))

        # Create jobrunner if not provided
        if jobrunner is None:
            jobrunner = JobRunner.from_job(
                cls(
                    molecule=molecule,
                    settings=settings,
                    label=label,
                    folder=folder,
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
            folder=folder,
            jobrunner=jobrunner,
            **kwargs,
        )

    def compute_thermochemistry(self):
        """
        Perform thermochemistry calculation and save results.

        Computes thermochemical properties including electronic energy,
        zero-point energy, thermal corrections, enthalpy, entropy, and
        Gibbs free energy from frequency calculation data.

        Raises:
            ValueError: If no input file or folder is provided
            Exception: If calculation fails during processing
        """
        if not self.filename and not self.folder:
            raise ValueError(
                "No input file or folder provided for thermochemistry calculation."
            )

        # Set default output file if not specified
        if self.settings.outputfile is None:
            self.settings.outputfile = self.outputfile

        try:
            thermochemistry = Thermochemistry(
                filename=self.filename,
                folder=self.folder,
                temperature=self.settings.temperature,
                concentration=self.settings.concentration,
                pressure=self.settings.pressure,
                use_weighted_mass=self.settings.use_weighted_mass,
                alpha=self.settings.alpha,
                s_freq_cutoff=self.settings.s_freq_cutoff,
                entropy_method=self.settings.entropy_method,
                h_freq_cutoff=self.settings.h_freq_cutoff,
                energy_units=self.settings.energy_units,
                outputfile=self.settings.outputfile,
                overwrite=self.settings.overwrite,
                check_imaginary_frequencies=(
                    self.settings.check_imaginary_frequencies
                ),
            )
            (
                structure,
                electronic_energy,
                zero_point_energy,
                enthalpy,
                qrrho_enthalpy,
                entropy_times_temperature,
                qrrho_entropy_times_temperature,
                gibbs_free_energy,
                qrrho_gibbs_free_energy,
            ) = thermochemistry.compute_thermochemistry()
            thermochemistry.log_results_to_file(
                structure,
                electronic_energy,
                zero_point_energy,
                enthalpy,
                qrrho_enthalpy,
                entropy_times_temperature,
                qrrho_entropy_times_temperature,
                gibbs_free_energy,
                qrrho_gibbs_free_energy,
                outputfile=self.settings.outputfile,
                overwrite=self.settings.overwrite,
                write_header=self.settings.write_header,
            )

        except Exception as e:
            logger.error(f"Error processing {self.filename}: {e}")
            raise

    def show_results(self):
        """
        Display the thermochemistry calculation results.

        Reads and prints the contents of the output file containing
        the computed thermochemical properties to the console.
        """
        if self.settings.outputfile and os.path.exists(
            self.settings.outputfile
        ):
            with open(self.settings.outputfile, "r") as out:
                print()
                results = out.read()
                print(results)
