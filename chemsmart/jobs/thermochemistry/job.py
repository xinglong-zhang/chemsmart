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
    PROGRAM = "Thermochemistry"
    TYPE = "thermochemistry"

    def __init__(
        self,
        filename=None,
        molecule=None,
        settings=None,
        label=None,
        jobrunner=None,
        **kwargs,
    ):
        super().__init__(
            molecule=molecule, label=label, jobrunner=jobrunner, **kwargs
        )

        if settings is not None and not isinstance(
            settings, ThermochemistryJobSettings
        ):
            raise ValueError(
                f"Settings must be instance of {ThermochemistryJobSettings} for {self}, but is {settings} instead!"
            )

        if molecule is not None and not isinstance(molecule, Molecule):
            raise ValueError(
                f"Molecule must be instance of Molecule for {self}, but is {molecule} instead!"
            )

        self.molecule = molecule.copy() if molecule is not None else None
        self.settings = (
            settings.copy()
            if settings is not None
            else ThermochemistryJobSettings()
        )
        self.filename = filename

        if label is None:
            if filename is not None:
                label = os.path.splitext(os.path.basename(filename))[0]
            elif molecule is not None:
                label = molecule.get_chemical_formula(empirical=True)
            else:
                label = "thermochemistry_job"
        self.label = label

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
        outputfile = self.label + ".dat"
        return os.path.join(self.folder, outputfile)

    @property
    def errfile(self):
        errfile = self.label + ".err"
        return os.path.join(self.folder, errfile)

    def _backup_files(self, **kwargs):
        folder = self._create_backup_folder_name()
        self.backup_file(self.outputfile, folder=folder, **kwargs)

    def _output(self):
        if not os.path.exists(self.settings.outputfile):
            return None
        return os.path.abspath(self.settings.outputfile)

    def _job_is_complete(self):
        return os.path.exists(self.settings.outputfile)

    def _run(self, **kwargs):
        """Run the thermochemistry analysis job."""
        logger.info(
            f"Running ThermochemistryJob {self} with jobrunner {self.jobrunner}"
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
        """Create a ThermochemistryJob from a Gaussian or ORCA output file."""
        if not filename.endswith((".log", ".out")):
            raise ValueError(
                f"Unsupported file extension for '{filename}'. Only .log or .out files are accepted."
            )

        logger.info(f"Reading molecule from file: {filename}.")
        molecule = Molecule.from_filepath(filename)

        if settings is None:
            settings = ThermochemistryJobSettings()

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

    def compute_thermochemistry(self):
        """Perform the thermochemistry calculation and save results."""
        if not self.filename:
            raise ValueError(
                "No input file provided for thermochemistry calculation."
            )

        if self.settings.outputfile is None:
            self.settings.outputfile = self.outputfile

        try:
            thermochemistry = Thermochemistry(
                filename=self.filename,
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
            )

        except Exception as e:
            logger.error(f"Error processing {self.filename}: {e}")
            raise

    def show_results(self):
        with open(self.settings.outputfile, "r") as out:
            print()
            results = out.read()
            print(results)
