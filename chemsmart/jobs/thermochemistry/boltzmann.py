import logging
import os
from typing import Type

from chemsmart.jobs.runner import JobRunner
from chemsmart.jobs.thermochemistry.job import ThermochemistryJob
from chemsmart.jobs.thermochemistry.settings import ThermochemistryJobSettings

logger = logging.getLogger(__name__)


class BoltzmannAverageThermochemistryJob(ThermochemistryJob):
    PROGRAM = "Thermochemistry"
    TYPE = "boltzmann"

    def __init__(
        self,
        files=None,
        energy_type="gibbs",
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.files = files
        self.energy_type = energy_type

        if self.label is None:
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
        return ThermochemistryJobSettings

    @property
    def inputfile(self):
        if self.filename:
            return os.path.abspath(self.filename)
        return None

    @property
    def outputfile(self):
        outputfile = self.label + "_boltzmann.dat"
        return os.path.join(self.folder, outputfile)

    @property
    def errfile(self):
        errfile = self.label + "_boltzmann.err"
        return os.path.join(self.folder, errfile)

    @classmethod
    def from_files(
        cls,
        files,
        jobrunner=None,
        **kwargs,
    ):
        """Create a ThermochemistryJob from a Gaussian or ORCA output file."""
        for file in files:
            if not file.endswith((".log", ".out")):
                raise ValueError(
                    f"Unsupported file extension for '{file}'. Only .log or .out files are accepted."
                )

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
        """Perform the thermochemistry calculation and save results."""
        if not self.files:
            raise ValueError(
                "No input file provided for thermochemistry calculation."
            )

        if self.settings.outputfile is None:
            self.settings.outputfile = self.outputfile

        try:
            from chemsmart.analysis.thermochemistry import (
                BoltzmannAverageThermochemistry,
            )

            thermochemistry = BoltzmannAverageThermochemistry(
                files=self.files,
                energy_type=self.energy_type,
                temperature=self.settings.temperature,
                concentration=self.settings.concentration,
                pressure=self.settings.pressure,
                use_weighted_mass=self.settings.use_weighted_mass,
                alpha=self.settings.alpha,
                s_freq_cutoff=self.settings.s_freq_cutoff,
                h_freq_cutoff=self.settings.h_freq_cutoff,
                energy_units=self.settings.energy_units,
                outputfile=self.settings.outputfile,
                overwrite=self.settings.overwrite,
                check_imaginary_frequencies=self.settings.check_imaginary_frequencies,
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
