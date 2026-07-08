import logging
import os
from typing import Type

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.job import Job
from chemsmart.jobs.xtb.settings import XTBJobSettings
from chemsmart.utils.utils import string2index_1based

logger = logging.getLogger(__name__)


class XTBJob(Job):
    """Base class for xTB jobs."""

    PROGRAM = "xtb"

    def __init__(
        self, molecule, settings=None, label=None, jobrunner=None, **kwargs
    ):
        if label is None:
            label = molecule.get_chemical_formula(empirical=True)
        self.label = label

        super().__init__(
            molecule=molecule, label=label, jobrunner=jobrunner, **kwargs
        )

        if not isinstance(settings, XTBJobSettings):
            raise ValueError(
                f"Settings must be instance of {XTBJobSettings} for {self}, "
                f"but is {settings} instead!"
            )
        if not isinstance(molecule, Molecule):
            raise ValueError(
                f"Molecule must be instance of Molecule for {self}, but is "
                f"{molecule} instead!"
            )

        self.molecule = molecule.copy()
        self.settings = settings.copy()
        self.label = label

    @classmethod
    def settings_class(cls) -> Type[XTBJobSettings]:
        return XTBJobSettings

    @property
    def xyzfile(self):
        return os.path.join(self.folder, f"{self.label}.xyz")

    @property
    def inputfile(self):
        return self.xyzfile

    @property
    def outputfile(self):
        return os.path.join(self.folder, f"{self.label}.out")

    @property
    def errfile(self):
        return os.path.join(self.folder, f"{self.label}.err")

    def _determine_folder(self):
        folder = os.path.abspath(os.getcwd())
        if self.label:
            folder = os.path.join(folder, self.label)
            os.makedirs(folder, exist_ok=True)
        return folder

    def base_folder(self):
        return self._determine_folder()

    def _backup_files(self, **kwargs):
        folder = self._create_backup_folder_name()
        self.backup_file(self.xyzfile, folder=folder, **kwargs)
        self.backup_file(self.outputfile, folder=folder, **kwargs)

    def _output(self):
        if not os.path.exists(self.outputfile):
            logger.debug(f"xTB output file not found: {self.outputfile}")
            return None
        try:
            from chemsmart.io.xtb.output import XTBOutput

            return XTBOutput(folder=self.folder)
        except Exception as exc:
            logger.error(
                f"Error reading xTB output folder {self.folder}: {exc}"
            )
            return None

    def _run(self, **kwargs):
        logger.info(f"Running XTBJob {self} with jobrunner {self.jobrunner}")
        self.jobrunner.run(self, **kwargs)

    @classmethod
    def from_filename(
        cls,
        filename,
        settings=None,
        index="-1",
        label=None,
        jobrunner=None,
        **kwargs,
    ):
        molecules = Molecule.from_filepath(
            filepath=filename, index=":", return_list=True
        )
        molecule = molecules[string2index_1based(index)]
        return cls(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
