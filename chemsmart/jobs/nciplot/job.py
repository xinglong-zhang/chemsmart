import logging
import os
from typing import Type

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.job import Job
from chemsmart.jobs.nciplot.settings import NCIPLOTJobSettings

logger = logging.getLogger(__name__)


class NCIPLOTJob(Job):
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
        super().__init__(
            molecule=molecule, label=label, jobrunner=jobrunner, **kwargs
        )

        if settings is not None and not isinstance(
            settings, NCIPLOTJobSettings
        ):
            raise ValueError(
                f"Settings must be instance of {NCIPLOTJobSettings} for {self}, but is {settings} instead!"
            )

        # either filenames or molecule must be provided, but not both
        print(filenames, molecule)
        if filenames is not None and molecule is not None:
            raise ValueError(
                "Either filenames or molecule must be provided, but not both!"
            )
        elif filenames is None and molecule is None:
            raise ValueError(
                "Either filenames or molecule must be provided, but both are None!"
            )

        if molecule is not None and not isinstance(molecule, Molecule):
            raise ValueError(
                f"Molecule must be instance of Molecule for {self}, but is {molecule} instead!"
            )

        self.molecule = molecule.copy() if molecule is not None else None
        self.settings = settings.copy()
        self.filenames = filenames

        if label is None:
            if self.settings.label is not None:
                label = self.settings.label
            else:
                if molecule is not None:
                    label = molecule.get_chemical_formula(empirical=True)
                else:
                    label = "nciplot_job"
        self.label = label

    @classmethod
    def settings_class(cls) -> Type[NCIPLOTJobSettings]:
        return NCIPLOTJobSettings

    @property
    def inputfile(self):
        inputfile = self.label + ".nci"
        return os.path.join(self.folder, inputfile)

    @property
    def outputfile(self):
        outputfile = self.label + ".nciout"
        return os.path.join(self.folder, outputfile)

    @property
    def errfile(self):
        errfile = self.label + ".ncierr"
        return os.path.join(self.folder, errfile)

    def _backup_files(self, **kwargs):
        folder = self._create_backup_folder_name()
        self.backup_file(self.outputfile, folder=folder, **kwargs)

    def _output(self):
        if not os.path.exists(self.outputfile):
            return None
        return os.path.abspath(self.outputfile)

    def _job_is_complete(self):
        return os.path.exists(self.outputfile)

    def _run(self, **kwargs):
        """Run the thermochemistry analysis job."""
        logger.info(
            f"Running NCIPLOTJob {self} with jobrunner {self.jobrunner}"
        )
        self.jobrunner.run(self, **kwargs)
