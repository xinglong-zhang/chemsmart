import datetime
import glob
import logging
import os
import shutil
import time
from abc import abstractmethod
from contextlib import suppress
from typing import Optional

from chemsmart.jobs.runner import JobRunner
from chemsmart.utils.mixins import RegistryMixin

logger = logging.getLogger(__name__)


class Job(RegistryMixin):
    """Class that encapsulates and runs a task.
    Args:
        molecule: The molecule object associated with the job.
        label (str): A label for the job.
        jobrunner (JobRunner): The JobRunner instance to execute the job.
        local (bool): If True, run the job locally. Defaults to False.
        skip_completed (bool): If True, completed jobs will not be rerun. Defaults to True.
        **kwargs: Additional keyword arguments.
    """

    TYPE: Optional[str] = None
    PROGRAM: Optional[str] = None

    def __init__(
        self,
        molecule,
        label,
        jobrunner,
        local=False,
        skip_completed=True,
        **kwargs,
    ):
        if not isinstance(jobrunner, JobRunner):
            raise ValueError(
                f"jobrunner must be an instance of JobRunner. Instead was: {jobrunner}!"
            )
        self._folder = self._determine_folder()  # Private backing attribute
        self.molecule = molecule
        self.label = label
        self.jobrunner = jobrunner
        self.local = local
        self.skip_completed = skip_completed
        self.kwargs = kwargs

    @property
    def folder(self):
        return self._folder

    @folder.setter
    def folder(self, folder):
        self._folder = folder

    def _determine_folder(self):
        """
        Determine the folder based on the current working directory
        where the job is submitted.
        """
        # Get the current working directory at runtime
        cwd = os.getcwd()
        folder = os.path.abspath(cwd)
        return folder

    def base_folder(self):
        if self.folder is None:
            return None
        return os.path.basename(os.path.abspath(self.folder))

    def __repr__(self):
        return (
            f"{self.__class__.__qualname__}<folder={self.folder}, "
            f"label={self.label}, jobrunner={self.jobrunner}>"
        )

    @property
    def joblog(self):
        return os.path.join(self.folder, f"{self.label}.joblog")

    def run(self, **kwargs):
        if self.is_complete() and self.skip_completed:
            logger.info(f"{self} is already complete, not running.")
            return
        self._run(**kwargs)

    @abstractmethod
    def _run(self, **kwargs):
        """
        Run the job using the assigned jobrunner.
        Subclasses can override this method for custom behavior.
        """
        logger.info(f"Running job {self} with jobrunner {self.jobrunner}")
        self.jobrunner.run(self, **kwargs)

    def set_folder(self, folder):
        self.folder = folder

    def _previous_backup_folders(self):
        return sorted(glob.glob(f"{self.folder}/bk*"))

    def _create_backup_folder_name(self, prefix="bk"):
        ts = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        folder = f"{prefix}.{ts}"
        return os.path.join(self.folder, folder)

    def backup(self, **kwargs):
        """Backup the current folder to a new folder.
        Return the new folder path.
        """
        self._backup_files(**kwargs)

    def _backup_files(self):
        raise NotImplementedError

    def backup_file(self, file, folder=None, remove=False):
        if not os.path.exists(file):
            return

        if folder is None:
            folder = self._create_backup_folder_name()

        with suppress(FileExistsError):
            os.mkdir(folder)

        newfilepath = os.path.join(
            folder, os.path.basename(os.path.abspath(file))
        )

        if remove:
            shutil.move(src=file, dst=newfilepath)
        else:
            if os.path.isdir(file):
                if os.path.exists(newfilepath):
                    os.rename(newfilepath, newfilepath + f"_{time.time()}")
                shutil.copytree(src=file, dst=newfilepath)
            else:
                shutil.copy(src=file, dst=newfilepath)

    @property
    def all_intermediate_optimization_points(self):
        import ase

        intermediate_optimization_points_path = os.path.join(
            self.label + "_intermediate_opt_points.xyz"
        )
        all_points = self._intermediate_optimization_points()
        ase.io.write(intermediate_optimization_points_path, all_points)
        return all_points

    def optimized_structure(self):
        output = self._output()
        if output is not None and output.normal_termination:
            return output.optimized_structure
        return None

    def _intermediate_optimization_points(self):
        output = self._output()
        if output is None:
            return []
        return output.all_structures

    def _job_is_complete(self):
        # private method to check if the job is complete
        if self._output() is None:
            return False
        return self._output().normal_termination

    def is_complete(self):
        return self._job_is_complete()

    @classmethod
    def from_molecule(
        cls, molecule, label, server=None, scratch=None, fake=False, **kwargs
    ):
        """
        Create a Job instance with a JobRunner initialized from the job's TYPE.

        Args:
            molecule: The molecule object associated with the job.
            label (str): A label for the job.
            server: The server to run the job on (optional, defaults to current server).
            scratch (bool): Whether to use a scratch directory (optional).
            fake (bool): Whether to use a fake job runner (optional).
            **kwargs: Additional keyword arguments for Job or JobRunner.

        Returns:
            Job: A new Job instance with an appropriate JobRunner.
        """
        jobrunner = JobRunner.from_job(
            cls(molecule=molecule, label=label, jobrunner=None, **kwargs),
            server=server,
            scratch=scratch,
            fake=fake,
            **kwargs,
        )
        return cls(
            molecule=molecule, label=label, jobrunner=jobrunner, **kwargs
        )
