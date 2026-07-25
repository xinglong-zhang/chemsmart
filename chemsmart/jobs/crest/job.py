import logging
import os
from typing import Type

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.crest.settings import CRESTJobSettings
from chemsmart.jobs.job import Job
from chemsmart.utils.utils import string2index_1based

logger = logging.getLogger(__name__)


class CRESTJob(Job):
    """Base class for CREST jobs."""

    PROGRAM = "crest"

    def __init__(
        self, molecule, settings=None, label=None, jobrunner=None, **kwargs
    ):
        if label is None:
            label = molecule.get_chemical_formula(empirical=True)
        self.label = label

        super().__init__(
            molecule=molecule, label=label, jobrunner=jobrunner, **kwargs
        )

        if not isinstance(settings, CRESTJobSettings):
            raise ValueError(
                f"Settings must be instance of {CRESTJobSettings} for {self}, "
                f"but is {settings} instead!"
            )
        if not isinstance(molecule, Molecule):
            raise ValueError(
                f"Molecule must be instance of Molecule for {self}, but is "
                f"{molecule} instead!"
            )

        self.molecule = molecule.copy()
        self.settings = settings.copy()

    @classmethod
    def settings_class(cls) -> Type[CRESTJobSettings]:
        return CRESTJobSettings

    @property
    def xyzfile(self):
        return os.path.join(self.folder, f"{self.label}.xyz")

    @property
    def outputfile(self):
        return os.path.join(self.folder, f"{self.label}.out")

    @property
    def errfile(self):
        return os.path.join(self.folder, f"{self.label}.err")

    @property
    def constraints_file(self):
        return os.path.join(self.folder, "constraints.inp")

    @property
    def toml_file(self):
        return os.path.join(self.folder, f"{self.label}.toml")

    @property
    def conformers_file(self):
        return os.path.join(self.folder, "crest_conformers.xyz")

    @property
    def best_file(self):
        return os.path.join(self.folder, "crest_best.xyz")

    @property
    def rotamers_file(self):
        return os.path.join(self.folder, "crest_rotamers.xyz")

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
            logger.debug(f"CREST output file not found: {self.outputfile}")
            return None
        try:
            from chemsmart.io.crest.output import CRESTOutput

            return CRESTOutput(folder=self.folder)
        except Exception as exc:
            logger.error(
                f"Error reading CREST output folder {self.folder}: {exc}"
            )
            return None

    def _run(self, **kwargs):
        logger.info(f"Running CRESTJob {self} with jobrunner {self.jobrunner}")
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
