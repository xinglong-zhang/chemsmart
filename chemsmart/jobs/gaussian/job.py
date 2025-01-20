import logging
import os
from typing import Type
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.job import Job
from chemsmart.jobs.gaussian.settings import GaussianJobSettings
from chemsmart.utils.utils import string2index_1based

logger = logging.getLogger(__name__)


class GaussianJob(Job):
    PROGRAM = "Gaussian"

    def __init__(self, molecule, settings=None, label=None, **kwargs):
        super().__init__(molecule=molecule, label=label, **kwargs)
        if not isinstance(settings, GaussianJobSettings):
            raise ValueError(
                f"Settings must be instance of {GaussianJobSettings} for {self}, but is {settings} instead!"
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
    def settings_class(cls) -> Type[GaussianJobSettings]:
        return GaussianJobSettings

    @property
    def inputfile(self):
        inputfile = self.label + ".com"
        return os.path.join(self.folder, inputfile)

    @property
    def outputfile(self):
        outputfile = self.label + ".log"
        return os.path.join(self.folder, outputfile)

    @property
    def chkfile(self):
        chkfile = self.label + ".chk"
        return os.path.join(self.folder, chkfile)

    @property
    def errfile(self):
        errfile = self.label + ".err"
        return os.path.join(self.folder, errfile)

    def _backup_files(self, backup_chk=False, **kwargs):
        folder = self._create_backup_folder_name()
        self.backup_file(self.inputfile, folder=folder, **kwargs)
        self.backup_file(self.outputfile, folder=folder, **kwargs)
        if backup_chk:
            self.backup_file(self.chkfile, folder=folder, **kwargs)

    def _output(self):
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

    def _run(self, jobrunner):
        jobrunner.run(self)

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
            from chemsmart.jobs.gaussian.opt import GaussianOptJob

            logger.debug(f"Creating GaussianOptJob from jobtype: {jobtype}")

            return GaussianOptJob(
                molecule=molecule,
                settings=settings,
                label=label,
                **kwargs,
            )
        elif jobtype.lower() == "com":
            return GaussianComJob(
                molecule=molecule,
                settings=settings,
                label=label,
                **kwargs,
            )
        elif jobtype.lower() == "g16":
            return GaussianGeneralJob(
                molecule=molecule,
                settings=settings,
                label=label,
                **kwargs,
            )
        else:
            raise ValueError(f"Invalid job type: {jobtype}")


class GaussianComJob(GaussianJob):
    """Runs any given .com Gaussian input file as is."""

    TYPE = "g16com"

    def __init__(self, molecule, settings=None, label=None, **kwargs):
        super().__init__(
            molecule=molecule, settings=settings, label=label, **kwargs
        )

    @classmethod
    def from_filename(cls, filename, settings=None, label=None, **kwargs):
        # job.label as the filename (without extension) used
        if label is None:
            label = os.path.splitext(os.path.basename(filename))[0]

        # set file
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

        return cls(molecule=None, settings=settings, label=label, **kwargs)


class GaussianGeneralJob(GaussianJob):
    """GaussianGeneralJob subclasses GaussianJob, this is needed to prevent recursive loop.

    For example, recursive loop occurs in class GaussianCrestOptJob(GaussianJob) that
    subclasses GaussianJob and calls and runs GaussianGeneralJob.
    """

    TYPE = "g16"

    def __init__(self, molecule, settings=None, label=None, **kwargs):
        super().__init__(
            molecule=molecule, settings=settings, label=label, **kwargs
        )
