import logging
import os
from typing import Type

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.settings import GaussianJobSettings
from chemsmart.jobs.job import Job
from chemsmart.jobs.runner import JobRunner

logger = logging.getLogger(__name__)


class GaussianJob(Job):
    PROGRAM = "Gaussian"

    def __init__(
        self, molecule, settings=None, label=None, jobrunner=None, **kwargs
    ):
        super().__init__(
            molecule=molecule, label=label, jobrunner=jobrunner, **kwargs
        )

        if not isinstance(settings, GaussianJobSettings):
            raise ValueError(
                f"Settings must be instance of {GaussianJobSettings} for {self}, but is {settings} instead!"
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

    def _run(self, **kwargs):
        """Run the job using the assigned jobrunner."""
        logger.info(
            f"Running GaussianJob {self} with jobrunner {self.jobrunner}"
        )
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
        """Create a GaussianJob from a file containing molecule data."""
        logger.info(f"Reading molecules from file: {filename}.")
        molecules = Molecule.from_filepath(
            filepath=filename, index=":", return_list=True
        )
        logger.info(f"Num of molecules read: {len(molecules)}.")

        # use 1-based indexing for the index
        from chemsmart.utils.utils import return_objects_from_string_index

        molecules = return_objects_from_string_index(molecules, index)

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
        """Create a GaussianJob from a PubChem identifier."""
        molecules = Molecule.from_pubchem(identifier=identifier)

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

        return cls(
            molecule=molecules,
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
        """Create a GaussianJob based on the specified job type."""
        if jobtype.lower() == "opt":
            from chemsmart.jobs.gaussian.opt import GaussianOptJob

            logger.debug(f"Creating GaussianOptJob from jobtype: {jobtype}")

            # Create jobrunner if not provided
            if jobrunner is None:
                jobrunner = JobRunner.from_job(
                    GaussianOptJob(
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

            return GaussianOptJob(
                molecule=molecule,
                settings=settings,
                label=label,
                jobrunner=jobrunner,
                **kwargs,
            )

        elif jobtype.lower() == "com":
            # Create jobrunner if not provided
            if jobrunner is None:
                jobrunner = JobRunner.from_job(
                    GaussianComJob(
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

            return GaussianComJob(
                molecule=molecule,
                settings=settings,
                label=label,
                jobrunner=jobrunner,
                **kwargs,
            )
        elif jobtype.lower() == "g16":
            # Create jobrunner if not provided
            if jobrunner is None:
                jobrunner = JobRunner.from_job(
                    GaussianGeneralJob(
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

                return GaussianGeneralJob(
                    molecule=molecule,
                    settings=settings,
                    label=label,
                    jobrunner=jobrunner,
                    **kwargs,
                )
            else:
                raise ValueError(f"Invalid job type: {jobtype}")


class GaussianComJob(GaussianJob):
    """Runs any given .com Gaussian input file as is."""

    TYPE = "g16com"

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
        """Create a GaussianComJob from a .com file."""
        # job.label as the filename (without extension) used
        if label is None:
            label = os.path.splitext(os.path.basename(filename))[0]

        # get input lines
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

        # Create jobrunner if not provided
        if jobrunner is None:
            jobrunner = JobRunner.from_job(
                cls(
                    molecule=None,
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
            molecule=None,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )


class GaussianGeneralJob(GaussianJob):
    """GaussianGeneralJob subclasses GaussianJob, this is needed
    to prevent recursive loop.
    For example, recursive loop occurs in class GaussianCrestOptJob(GaussianJob)
    that subclasses GaussianJob and calls and runs GaussianGeneralJob.
    """

    TYPE = "g16"

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
