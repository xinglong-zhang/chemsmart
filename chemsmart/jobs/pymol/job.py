import logging
import os

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.job import Job
from chemsmart.utils.utils import string2index_1based

logger = logging.getLogger(__name__)


class PyMOLJob(Job):
    PROGRAM = "PyMOL"

    def __init__(self, molecule=None, label=None, **kwargs):
        super().__init__(molecule=molecule, label=label, **kwargs)
        molecule = molecule.copy()

        self.molecule = molecule
        self.label = label

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
