import logging
import os
from typing import Type

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.job import Job
from chemsmart.jobs.xtb.settings import XTBJobSettings
from chemsmart.utils.utils import string2index_1based

logger = logging.getLogger(__name__)


class XTBJob(Job):
    PROGRAM = "xtb"

    def __init__(self, molecule, settings=None, label=None, **kwargs):
        super().__init__(molecule=molecule, label=label, **kwargs)
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

        if label is None:
            label = molecule.get_chemical_formula(empirical=True)

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
    def inputfile(self):
        inputfile = self.label + ".xyz"
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

    def _backup_files(self, backup_chk=False, **kwargs):
        folder = self._create_backup_folder_name()
        self.backup_file(self.inputfile, folder=folder, **kwargs)
        self.backup_file(self.outputfile, folder=folder, **kwargs)

    def _output(self):
        if not os.path.exists(self.outputfile):
            return None

        try:
            from chemsmart.io.xtb.output import XTBOutput

            return XTBOutput(filename=self.outputfile)
        except Exception:
            logger.error(f"Error reading output file: {self.outputfile}")
            return None

    def _run(self, jobrunner, **kwargs):
        jobrunner.run(self, **kwargs)

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
        # get all molecule in a file and give the result as a list
        logger.info(f"Reading images from file: {filename}.")
        molecules = Molecule.from_filepath(
            filepath=filename, index=":", return_list=True
        )
        logger.info(f"Num of images read: {len(molecules)}.")
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
