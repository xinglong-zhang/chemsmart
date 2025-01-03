import ase
import logging
import os

from chemsmart.jobs.job import Job

from pyatoms.io.ase.atoms import AtomsWrapper
from pyatoms.jobs.gaussian.settings import GaussianJobSettings
from pyatoms.utils.utils import string2index

logger = logging.getLogger(__name__)


class GaussianJob(Job):
    PROGRAM_TYPE = 'Gaussian'
    _SETTINGS_CLS = GaussianJobSettings

    def __init__(self, folder, atoms, settings=None, label=None, **kwargs):
        super().__init__(folder=folder, write_jobtype_file=False, num_tasks_per_node=1, label=label, **kwargs)
        if not isinstance(settings, self._SETTINGS_CLS):
            raise ValueError(
                f'Settings must be instance of {self._SETTINGS_CLS} for {self}, but is {settings} instead!'
            )

        atoms = atoms.copy()
        settings = settings.copy()

        if label is None:
            label = atoms.get_chemical_formula(empirical=True)

        self.settings = settings
        self.atoms = atoms
        self.label = label

    @property
    def inputfile(self):
        inputfile = self.label + '.com'
        return os.path.join(self.folder, inputfile)

    @property
    def outfile(self):
        outputfile = self.label + '.log'
        return os.path.join(self.folder, outputfile)

    @property
    def chkfile(self):
        chkfile = self.label + '.chk'
        return os.path.join(self.folder, chkfile)

    @property
    def errfile(self):
        errfile = self.label + '.err'
        return os.path.join(self.folder, errfile)

    @property
    def all_intermediate_optimization_points(self):
        intermediate_optimization_points_path = os.path.join(self.label + '_intermediate_opt_points.xyz')
        all_points = self._intermediate_optimization_points()
        ase.io.write(intermediate_optimization_points_path, all_points)
        return all_points

    def optimized_atoms(self):
        output = self._output()
        if output is not None and output.normal_termination:
            return output.get_atoms(index='-1')
        return None

    def _intermediate_optimization_points(self):
        output = self._output()
        if output is None:
            return []
        return output.get_atoms(index=':')

    def _compute_resources_required(self):
        return self._compute_units(num_units=1)

    def backup_files(self, backup_chk=False):
        folder = self._create_backup_folder_name()
        self.backup_file(self.inputfile, folder=folder)
        self.backup_file(self.outfile, folder=folder)
        if backup_chk:
            self.backup_file(self.chkfile, folder=folder)

    def _output(self):
        if not os.path.exists(self.outfile):
            return None

        try:
            from pyatoms.io.gaussian.outputs import Gaussian16Output

            return Gaussian16Output(self.outfile)
        except AttributeError:
            return None
        except ValueError:
            from pyatoms.io.gaussian.outputs import Gaussian16OutputWithPBC

            return Gaussian16OutputWithPBC(self.outfile)
        else:
            return None

    def _job_is_complete(self):
        # private method to check if the job is complete
        if self._output() is None:
            return False
        return self._output().normal_termination

    def is_complete(self):
        return self._job_is_complete()

    def _run(self, jobrunner, queue_manager=None, **kwargs):
        self.settings.apply_on(self, jobrunner, **kwargs)
        jobrunner.run(self)

    @classmethod
    def from_filename(
        cls, folder, filename, settings=None, index='-1', label=None, keywords=('charge', 'multiplicity'), **kwargs
    ):
        # get all atoms in a file and give the result as a list
        logger.info(f'Reading images from file: {filename}.')
        atoms = AtomsWrapper.from_filepath(filepath=filename, index=':', return_list=True)
        logger.info(f'Num of images read: {len(atoms)}.')
        atoms = atoms[string2index(index)]  # python 0-indexed

        # only supply last atoms in some jobs; but require all atoms in others e.g., dias job
        return cls(folder=folder, atoms=atoms, settings=settings, label=label, **kwargs)

    @classmethod
    def from_pubchem(cls, folder, identifier, settings=None, label=None, **kwargs):
        atoms = AtomsWrapper.from_pubchem(identifier=identifier)
        return cls(folder=folder, atoms=atoms, settings=settings, label=label, **kwargs)

    @classmethod
    def from_label(cls, folder, label, atoms=None, settings=None, **kwargs):
        if atoms is None:
            comfile = os.path.join(folder, label + '.com')
            if not os.path.exists(comfile):
                pass
            atoms = AtomsWrapper.from_filepath(comfile)

        if settings is None:
            settings = GaussianJobSettings.default()

        return cls(folder=folder, label=label, atoms=atoms, settings=settings, **kwargs)


class GaussianComJob(GaussianJob):
    """Runs any given .com Gaussian input file as is."""

    TYPE = 'g16com'

    def __init__(self, folder, atoms, settings=None, **kwargs):
        super().__init__(folder=folder, atoms=atoms, settings=settings, **kwargs)

    @classmethod
    def from_filename(cls, folder, filename, settings=None, label=None, **kwargs):
        # job.label as the filename (without extension) used
        if label is None:
            label = os.path.splitext(os.path.basename(filename))[0]

        # set file
        from pyatoms.io.gaussian.inputs import Gaussian16Input

        g16com_file = Gaussian16Input(comfile=filename)

        # get atoms object from supplied .com file
        atoms = AtomsWrapper.from_filepath(filepath=filename)

        # get settings from file
        from pyatoms.jobs.gaussian.settings import GaussianJobSettings

        settings = GaussianJobSettings.from_filepath(filename)
        logger.info(f'Supplied file {filename} settings are: \n{settings.__dict__}')

        # update settings to write the same file - turn off everything that would be appended at the end from route info
        settings.route_to_be_written = g16com_file.route
        settings.modred = None
        settings.gen_genecp = None
        settings.custom_solvent = None

        # append original information as is to the end of input file
        settings.append_additional_info = g16com_file.append_info_after_coords_as_string
        logger.info(f'Rewritten settings are: \n{settings.__dict__}')

        return cls(folder=folder, atoms=atoms, settings=settings, label=label, **kwargs)


class GaussianGeneralJob(GaussianJob):
    """GaussianGeneralJob subclasses GaussianJob, this is needed to prevent recursive loop.

    For example, recursive loop occurs in class GaussianCrestOptJob(GaussianJob) that subclasses GaussianJob
    and calls and runs GaussianGeneralJob.
    """

    TYPE = 'g16'

    def __init__(self, folder, atoms, settings, **kwargs):
        super().__init__(folder=folder, atoms=atoms, settings=settings, **kwargs)
