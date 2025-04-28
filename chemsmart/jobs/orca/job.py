import logging
import os
import shutil
from contextlib import suppress
from typing import Type

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.job import Job
from chemsmart.jobs.orca.settings import ORCAJobSettings
from chemsmart.jobs.runner import JobRunner
from chemsmart.utils.utils import string2index_1based

logger = logging.getLogger(__name__)


class ORCAJob(Job):
    PROGRAM = "ORCA"

    def __init__(
        self, molecule, settings=None, label=None, jobrunner=None, **kwargs
    ):
        super().__init__(
            molecule=molecule, label=label, jobrunner=jobrunner, **kwargs
        )

        if not isinstance(settings, ORCAJobSettings):
            raise ValueError(
                f"Settings must be instance of {ORCAJobSettings} for {self}, but is {settings} instead!"
            )
        if not isinstance(molecule, Molecule):
            raise ValueError(
                f"Molecule must be instance of Molecule for {self}, but is {molecule} instead!"
            )

        self.molecule = molecule.copy() if molecule is not None else None
        self.settings = settings.copy()

        if label is None:
            label = molecule.get_chemical_formula(empirical=True)
        self.label = label

    @classmethod
    def settings_class(cls) -> Type[ORCAJobSettings]:
        return ORCAJobSettings

    @property
    def inputfile(self):
        inputfile = self.label + ".inp"
        return os.path.join(self.folder, inputfile)

    @property
    def outputfile(self):
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

    def _backup_files(self, backup_gbw=False, **kwargs):
        folder = self._create_backup_folder_name()
        self.backup_file(self.inputfile, folder=folder, **kwargs)
        self.backup_file(self.outputfile, folder=folder, **kwargs)
        if backup_gbw:
            self.backup_file(self.gbwfile, folder=folder, **kwargs)

    def _output(self):
        if not os.path.exists(self.outputfile):
            return None

        try:
            from chemsmart.io.orca.output import ORCAOutput

            return ORCAOutput(self.outputfile)
        except AttributeError:
            return None

    def _run(self, **kwargs):
        """Run the job using the assigned jobrunner."""
        logger.info(f"Running ORCAJob {self} with jobrunner {self.jobrunner}")
        self.jobrunner.run(self, **kwargs)

    @classmethod
    def from_filename(
        cls,
        filename,
        settings=None,
        index="-1",
        label=None,
        jobrunner=None,
        keywords=("charge", "multiplicity"),
        **kwargs,
    ):
        molecules = Molecule.from_filepath(
            filepath=filename, index=":", return_list=True
        )
        logger.info(f"Num of molecules read: {len(molecules)}.")
        molecules = molecules[string2index_1based(index)]

        # Create jobrunner if not provided
        if jobrunner is None:
            jobrunner = JobRunner.from_job(
                cls(
                    molecule=molecules,
                    settings=settings,
                    label=label,
                    jobrunner=None,
                    **kwargs,
                ),
                server=kwargs.get("server"),
                scratch=kwargs.get("scratch"),
                fake=kwargs.get("fake", False),
                **kwargs,
            )

        # only supply last molecule in some jobs; but require all molecule in others e.g., dias job
        return cls(
            molecule=molecules,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )

    @classmethod
    def from_pubchem(
        cls, identifier, settings=None, label=None, jobrunner=None, **kwargs
    ):
        molecule = Molecule.from_pubchem(identifier=identifier)

        # Create jobrunner if not provided
        if jobrunner is None:
            jobrunner = JobRunner.from_job(
                cls(
                    molecule=molecule,
                    settings=settings,
                    label=label,
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
            jobrunner=jobrunner,
            **kwargs,
        )

    @classmethod
    def from_jobtype(
        cls,
        jobtype,
        molecule,
        settings=None,
        label=None,
        jobrunner=None,
        **kwargs,
    ):
        if jobtype.lower() == "opt":
            from chemsmart.jobs.orca.opt import ORCAOptJob

            logger.debug(f"Creating GaussianOptJob from jobtype: {jobtype}")

            # Create jobrunner if not provided
            if jobrunner is None:
                jobrunner = JobRunner.from_job(
                    ORCAOptJob(
                        molecule=molecule,
                        settings=settings,
                        label=label,
                        jobrunner=None,
                        **kwargs,
                    ),
                    server=kwargs.get("server"),
                    scratch=kwargs.get("scratch"),
                    fake=kwargs.get("fake", False),
                    **kwargs,
                )

            return ORCAOptJob(
                molecule=molecule,
                settings=settings,
                label=label,
                jobrunner=jobrunner,
                **kwargs,
            )
        elif jobtype.lower() == "inp":
            # Create jobrunner if not provided
            if jobrunner is None:
                jobrunner = JobRunner.from_job(
                    ORCAInpJob(
                        molecule=molecule,
                        settings=settings,
                        label=label,
                        jobrunner=None,
                        **kwargs,
                    ),
                    server=kwargs.get("server"),
                    scratch=kwargs.get("scratch"),
                    fake=kwargs.get("fake", False),
                    **kwargs,
                )

            return ORCAInpJob(
                molecule=molecule,
                settings=settings,
                label=label,
                jobrunner=jobrunner,
                **kwargs,
            )
        elif jobtype.lower() == "orca":
            # Create jobrunner if not provided
            if jobrunner is None:
                jobrunner = JobRunner.from_job(
                    ORCAGeneralJob(
                        molecule=molecule,
                        settings=settings,
                        label=label,
                        jobrunner=None,
                        **kwargs,
                    ),
                    server=kwargs.get("server"),
                    scratch=kwargs.get("scratch"),
                    fake=kwargs.get("fake", False),
                    **kwargs,
                )

            return ORCAGeneralJob(
                molecule=molecule,
                settings=settings,
                label=label,
                jobrunner=jobrunner,
                **kwargs,
            )
        else:
            raise ValueError(f"Invalid job type: {jobtype}")


class ORCAInpJob(ORCAJob):
    """Runs any given .inp ORCA input file as is."""

    TYPE = "orcainp"

    def __init__(
        self, molecule, settings=None, label=None, jobrunner=None, **kwargs
    ):
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )

    @classmethod
    def from_filename(
        cls, filename, settings=None, label=None, jobrunner=None, **kwargs
    ):

        # first check that the file is orca format with .inp extension
        logger.debug(f"Checking if {filename} is an ORCA input file.")
        assert filename.endswith(".inp"), "Input file must be .inp file."

        molecule = Molecule.from_filepath(filepath=filename)

        # job.label as the filename (without extension) used
        if label is None:
            label = os.path.splitext(os.path.basename(filename))[0]

        # set file
        from chemsmart.io.orca.input import ORCAInput

        orca_file = ORCAInput(filename=filename)
        input_lines = orca_file.content_lines_string

        # get settings from file
        from chemsmart.jobs.orca.settings import ORCAJobSettings

        # store file lines in settings for direct writing of settings.input_string
        settings = ORCAJobSettings.from_filepath(filename)
        settings.input_string = input_lines

        logger.debug(
            f"Supplied file {filename} settings are: \n{settings.__dict__}"
        )

        logger.debug(f"Writing input lines: \n{input_lines}")

        # Create jobrunner if not provided
        if jobrunner is None:
            jobrunner = JobRunner.from_job(
                cls(
                    molecule=molecule,
                    settings=settings,
                    label=label,
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
            jobrunner=jobrunner,
            **kwargs,
        )

    def _run(self, **kwargs):
        """Override the _run method of parent to run the job as is."""
        # instead of applying setting and write the input file, create the input file
        # from the supplied .inp file and run it as is
        self._copy_input()
        self.jobrunner.run(self, **kwargs)

    def _copy_input(self):
        """Copy the supplied orca .inp file."""
        if (
            self.jobrunner.scratch
            and self.jobrunner.scratch_dir is not None
            and os.path.exists(self.jobrunner.scratch_dir)
        ):
            # if running job in scratch, then copy the input file over to scratch
            job_scratch_dir = os.path.join(
                self.jobrunner.scratch_dir, self.label
            )
            with suppress(FileExistsError):
                os.mkdir(job_scratch_dir)
                logger.info(f"Folder in scratch {job_scratch_dir} is made.")
            shutil.copy(self.inputfile, job_scratch_dir)
            scratch_inputfile = os.path.join(
                job_scratch_dir, f"{self.label}.inp"
            )
            assert os.path.exists(
                scratch_inputfile
            ), f"inputfile {scratch_inputfile} is not found"
        elif self.jobrunner.scratch and self.jobrunner.scratch_dir is not None:
            # if running job in scratch, but scratch dir is not found, then run the job in the run folder
            logger.warning(
                f"{self.jobrunner.scratch_dir} does not exist! Running job in {self.folder}."
            )
        else:
            logger.info(f"Running job in {self.folder}.")


class ORCAGeneralJob(ORCAJob):
    """ORCAGeneralJob subclasses ORCAJob, this is needed to prevent recursive loop.

    For example, recursive loop occurs in class ORCACrestJob(ORCAJob) that
    subclasses ORCAJob and calls and runs ORCAGeneralJob.
    """

    TYPE = "orcajob"

    def __init__(self, molecule, settings=None, label=None, **kwargs):
        super().__init__(
            molecule=molecule, settings=settings, label=label, **kwargs
        )
