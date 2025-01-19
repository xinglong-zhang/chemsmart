import ase
import logging
import os
import shutil
from contextlib import suppress
from typing import Type
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.job import Job
from chemsmart.jobs.orca.settings import ORCAJobSettings
from chemsmart.utils.utils import string2index_1based

logger = logging.getLogger(__name__)


class ORCAJob(Job):
    PROGRAM_TYPE = "ORCA"

    def __init__(self, molecule, settings=None, label=None, **kwargs):
        super().__init__(molecule=molecule, label=label, **kwargs)

        if not isinstance(settings, ORCAJobSettings):
            raise ValueError(
                f"Settings must be instance of {ORCAJobSettings} for {self}, but is {settings} instead!"
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
        self.atoms = molecule
        self.label = label

    @classmethod
    def settings_class(cls) -> Type[ORCAJobSettings]:
        return ORCAJobSettings

    @property
    def inputfile(self):
        inputfile = self.label + ".inp"
        return os.path.join(self.folder, inputfile)

    @property
    def outfile(self):
        outputfile = self.label + ".out"
        return os.path.join(self.folder, outputfile)

    @property
    def gbwfile(self):
        orca_gbwfile = self.label + ".gbw"
        return os.path.join(self.folder, orca_gbwfile)

    @property
    def errfile(self):
        errfile = self.label + ".err"
        return os.path.join(self.folder, errfile)

    @property
    def all_intermediate_optimization_points(self):
        intermediate_optimization_points_path = os.path.join(
            self.label + "_intermediate_opt_points.xyz"
        )
        all_points = self._intermediate_optimization_points()
        ase.io.write(intermediate_optimization_points_path, all_points)
        return all_points

    def optimized_structure(self):
        output = self._output()
        if output is not None and output.normal_termination:
            return output.get_atoms(index="-1")
        return None

    def _intermediate_optimization_points(self):
        output = self._output()
        if output is None:
            return []
        return output.all_structures

    def _compute_resources_required(self):
        return self._compute_units(num_units=1)

    def backup_files(self, backup_gbw=False):
        folder = self._create_backup_folder_name()
        self.backup_file(self.inputfile, folder=folder)
        self.backup_file(self.outfile, folder=folder)
        if backup_gbw:
            self.backup_file(self.gbwfile, folder=folder)

    def _output(self):
        if not os.path.exists(self.outfile):
            return None

        try:
            from pyatoms.io.orca.outputs import ORCAOutput

            return ORCAOutput(self.outfile)
        except AttributeError:
            return None
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
        cls,
        folder,
        filename,
        settings=None,
        index="-1",
        label=None,
        keywords=("charge", "multiplicity"),
        **kwargs,
    ):
        # get all atoms in a file and give the result as a list
        logger.info(f"Reading images from file: {filename}.")
        atoms = AtomsWrapper.from_filepath(
            filepath=filename, index=":", return_list=True
        )
        logger.info(f"Num of images read: {len(atoms)}.")
        atoms = atoms[string2index_1based(index)]  # python 0-indexed

        # only supply last atoms in some jobs; but require all atoms in others e.g., dias job
        return cls(
            folder=folder,
            atoms=atoms,
            settings=settings,
            label=label,
            **kwargs,
        )

    @classmethod
    def from_pubchem(
        cls, folder, identifier, settings=None, label=None, **kwargs
    ):
        atoms = AtomsWrapper.from_pubchem(identifier=identifier)
        return cls(
            folder=folder,
            atoms=atoms,
            settings=settings,
            label=label,
            **kwargs,
        )

    @classmethod
    def from_label(cls, folder, label, atoms=None, settings=None, **kwargs):
        if atoms is None:
            inpfile = os.path.join(folder, label + ".inp")
            if not os.path.exists(inpfile):
                pass
            atoms = AtomsWrapper.from_filepath(inpfile)

        if settings is None:
            settings = ORCAJobSettings.default()

        return cls(
            folder=folder,
            label=label,
            atoms=atoms,
            settings=settings,
            **kwargs,
        )


class ORCAInpJob(ORCAJob):
    """Runs any given .inp ORCA input file as is."""

    TYPE = "orcainp"

    def __init__(self, folder, atoms, settings=None, **kwargs):
        super().__init__(
            folder=folder, atoms=atoms, settings=settings, **kwargs
        )

    def _run(self, jobrunner, queue_manager=None, **kwargs):
        """Override the _run method of parent to run the job as is."""
        # instead of applying setting and write the input file, create the input file from the supplied .inp file
        # and run it as is
        self._copy_input(jobrunner=jobrunner)
        jobrunner.run(self)

    def _copy_input(self, jobrunner):
        """Coy the supplied orca .inp file."""
        if (
            jobrunner.scratch
            and jobrunner.scratch_dir is not None
            and os.path.exists(jobrunner.scratch_dir)
        ):
            # if running job in scratch, then copy the input file over to scratch
            job_scratch_dir = os.path.join(jobrunner.scratch_dir, self.label)
            with suppress(FileExistsError):
                os.mkdir(job_scratch_dir)
                logger.info(f"Folder in scratch {job_scratch_dir} is made.")
            shutil.copy(self.inputfile, job_scratch_dir)
            scratch_inputfile = os.path.join(
                job_scratch_dir, self.label + ".inp"
            )
            assert os.path.exists(
                scratch_inputfile
            ), f"inputfile {scratch_inputfile} is not found"
        elif jobrunner.scratch and not os.path.exists(jobrunner.scratch_dir):
            # if running job in scratch, but scratch dir is not found, then run the job in the run folder
            logger.info(
                f"{jobrunner.scratch_dir} does not exist! Running job in {self.folder}."
            )
        else:
            logger.info(f"Running job in {self.folder}.")

    @classmethod
    def from_filename(
        cls, folder, filename, settings=None, label=None, **kwargs
    ):
        # job.label as the filename (without extension) used
        if label is None:
            label = os.path.splitext(os.path.basename(filename))[0]

        # get atoms object from supplied file
        atoms = AtomsWrapper.from_filepath(filepath=filename)

        # get settings from file
        from pyatoms.jobs.orca.settings import ORCAJobSettings

        settings = ORCAJobSettings.from_inpfile(inp_path=filename)
        logger.info(
            f"Supplied file {filename} settings are: \n{settings.__dict__}"
        )
        return cls(
            folder=folder,
            atoms=atoms,
            settings=settings,
            label=label,
            **kwargs,
        )
