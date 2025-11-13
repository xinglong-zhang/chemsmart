"""
Boltzmann-averaged thermochemistry job implementation.

This module provides functionality for computing Boltzmann-weighted
thermochemical properties from multiple conformers or structures,
enabling accurate statistical mechanical calculations.
"""

import logging
import os
from typing import Type

from chemsmart.jobs.runner import JobRunner
from chemsmart.jobs.thermochemistry.job import ThermochemistryJob
from chemsmart.jobs.thermochemistry.settings import ThermochemistryJobSettings

logger = logging.getLogger(__name__)


class BoltzmannAverageThermochemistryJob(ThermochemistryJob):
    """
    Job for computing Boltzmann-weighted thermochemical properties.

    This class handles thermochemical calculations involving multiple
    conformers or structures, computing population-weighted averages
    of thermodynamic properties based on Boltzmann statistics.

    Attributes:
        PROGRAM (str): Program identifier ('Thermochemistry').
        TYPE (str): Job type identifier ('boltzmann').
        files (list[str] | None): List of output files to analyze (.log/.out).
        energy_type (str): Energy used for weighting ('gibbs', 'enthalpy', 'electronic').
        settings (ThermochemistryJobSettings): Thermochemistry configuration.
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    PROGRAM = "Thermochemistry"
    TYPE = "boltzmann"

    def __init__(
        self,
        files=None,
        energy_type="gibbs",
        **kwargs,
    ):
        """
        Initialize Boltzmann-averaged thermochemistry job.

        Args:
            files (list): List of output files for thermochemical analysis
            energy_type (str): Type of energy for Boltzmann weighting
                             ('gibbs', 'enthalpy', 'electronic')
            **kwargs: Additional keyword arguments passed to parent class
        """
        super().__init__(**kwargs)
        self.files = files
        self.energy_type = energy_type

        # Generate default label from common filename prefix if not provided
        if self.label is None:
            logger.debug("Generating label from file list")
            self.label = (
                os.path.commonprefix(
                    [
                        os.path.splitext(os.path.basename(f))[0]
                        for f in self.files
                    ]
                )
                + f"_boltzmann_avg_by_{self.energy_type}"
            )

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
        Get the path to the primary input file.

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
        Get the path to the Boltzmann average output file.

        Returns:
            str: Absolute path to the Boltzmann output file
        """
        outputfile = self.label + "_boltzmann.dat"
        return os.path.join(self.folder, outputfile)

    @property
    def errfile(self):
        """
        Get the path to the error file for the Boltzmann job.

        Returns:
            str: Absolute path to the Boltzmann error file
        """
        errfile = self.label + "_boltzmann.err"
        return os.path.join(self.folder, errfile)

    @classmethod
    def from_files(
        cls,
        files,
        jobrunner=None,
        **kwargs,
    ):
        """
        Create a Boltzmann thermochemistry job from output files.

        Creates a job instance from multiple Gaussian or ORCA output files
        for Boltzmann-weighted thermochemical analysis.

        Args:
            files (list): List of paths to quantum chemistry output files
            jobrunner (JobRunner, optional): Job runner instance for
                                           execution management
            **kwargs: Additional keyword arguments for job configuration

        Returns:
            BoltzmannAverageThermochemistryJob: Configured job instance
        """

        # Create jobrunner if not provided
        if jobrunner is None:
            jobrunner = JobRunner.from_job(
                cls(
                    files=files,
                    **kwargs,
                ),
                server=kwargs.get("server"),
                scratch=kwargs.get("scratch"),
                fake=kwargs.get("fake", False),
                **kwargs,
            )

        return cls(
            files=files,
            jobrunner=jobrunner,
            **kwargs,
        )

    def compute_boltzmann_averages(self):
        """
        Perform Boltzmann-weighted thermochemistry calculation.

        Computes population-weighted thermochemical properties from
        multiple conformers using Boltzmann statistics and saves
        results to the specified output file.

        Raises:
            ValueError: If no input files are provided
            Exception: If calculation fails during processing
        """
        if not self.files:
            raise ValueError(
                "No input file provided for thermochemistry calculation."
            )

        # Set default output file if not specified
        if self.settings.outputfile is None:
            self.settings.outputfile = self.outputfile

        try:
            from chemsmart.analysis.thermochemistry import (
                BoltzmannAverageThermochemistry,
            )

            # Initialize thermochemistry analyzer with job settings
            thermochemistry = BoltzmannAverageThermochemistry(
                files=self.files,
                energy_type=self.energy_type,
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
            ) = thermochemistry.compute_boltzmann_averages()
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
            )

        except Exception as e:
            logger.error(f"Error processing {self.filename}: {e}")
            raise
