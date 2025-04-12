import logging
import os

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.job import Job
from chemsmart.utils.utils import string2index_1based

logger = logging.getLogger(__name__)


class PyMOLJob(Job):
    PROGRAM = "PyMOL"

    def __init__(
        self,
        molecule=None,
        label=None,
        pymol_script=None,
        render=None,
        trace=None,
        vdw=None,
        quiet_mode=True,
        command_line_only=True,
        **kwargs,
    ):
        super().__init__(molecule=molecule, label=label, **kwargs)
        self.pymol_script = pymol_script
        self.render = render
        self.trace = trace
        self.vdw = vdw
        self.quiet_mode = quiet_mode
        self.command_line_only = command_line_only
        molecule = molecule.copy()

        self.molecule = molecule
        self.label = label

    @property
    def inputfile(self):
        inputfile = self.label + ".xyz"
        return os.path.join(self.folder, inputfile)

    @property
    def logfile(self):
        logfile = "log." + self.label
        return os.path.join(self.folder, logfile)

    @property
    def outputfile(self):
        outputfile = self.label + ".pse"
        return os.path.join(self.folder, outputfile)

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
        return os.path.abspath(self.outputfile)

    def _job_is_complete(self):
        return os.path.exists(self.outputfile)

    def _run(self, jobrunner):
        jobrunner.run(self)

    @classmethod
    def from_filename(
        cls,
        filename,
        pymol_script=None,
        index="-1",
        label=None,
        render_style=None,
        vdw=None,
        quiet_mode=True,
        command_line_only=True,
        **kwargs,
    ):
        # get all molecule in a file and give the result as a list
        logger.info(f"Reading images from file: {filename}.")
        molecules = Molecule.from_filepath(
            filepath=filename, index=":", return_list=True
        )
        if label is None:
            # by default, if no label is given and the job is read in
            # from a file, the label is set to the file basename
            label = os.path.basename(filename).split(".")[0]

        logger.info(f"Num of images read: {len(molecules)}.")
        molecules = molecules[string2index_1based(index)]

        return cls(
            molecule=molecules,
            label=label,
            pymol_script=pymol_script,
            render=render_style,
            vdw=vdw,
            quiet_mode=quiet_mode,
            command_line_only=command_line_only,
            **kwargs,
        )

    @classmethod
    def from_pubchem(cls, identifier, label=None, **kwargs):
        molecules = Molecule.from_pubchem(identifier=identifier)
        return cls(
            molecule=molecules,
            label=label,
            **kwargs,
        )
