import datetime
import glob
import logging
import os
import shutil
import time
from abc import abstractmethod
from contextlib import suppress
from typing import Optional

from chemsmart.utils.mixins import RegistryMixin

logger = logging.getLogger(__name__)


class Job(RegistryMixin):
    """Class that encapsulates and runs a task.
    Args:
        molecule: The molecule object associated with the job.
        label (str): A label for the job.
        jobrunner (JobRunner): The JobRunner instance to execute the job.
        local (bool): If True, run the job locally. Defaults to False.
        skip_completed (bool): If True, completed
        jobs will not be rerun. Defaults to True.
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

    @staticmethod
    def _propagate_runner(
        runner, job, *, num_cores=None, mem_gb=None, num_gpus=None
    ):
        """Copy *runner* onto *job*, optionally patching resource fields.

        This is the **single, centralised** mechanism for propagating a
        runner to a child job.  Every call site that needs to hand a runner
        to a child should use this helper rather than assigning directly, so
        that:

        * The child always receives its own ``.copy()`` — no shared-reference
          mutations between siblings.
        * Resource overrides (e.g. splitting cores across parallel workers)
          are applied in one place, making future additions (``num_gpus``,
          ``mem_per_cpu``, …) a single-line change.

        Args:
            runner: The parent runner to propagate, or ``None``.
            job (Job): The child job that will receive the runner copy.
            num_cores (int, optional): Override ``runner.num_cores`` on the copy.
            mem_gb (int, optional): Override ``runner.mem_gb`` on the copy.
            num_gpus (int, optional): Override ``runner.num_gpus`` on the copy.

        Returns:
            The runner copy that was assigned, or ``None`` if *runner* was falsy.
        """
        if not runner:
            return None
        child_runner = runner.copy()
        if num_cores is not None:
            child_runner.num_cores = num_cores
        if mem_gb is not None:
            child_runner.mem_gb = mem_gb
        if num_gpus is not None:
            child_runner.num_gpus = num_gpus
        job.jobrunner = child_runner
        return child_runner

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
    def from_molecule(cls, molecule, label, jobrunner=None, **kwargs):
        """
        Create a Job instance with a JobRunner initialized from the job's TYPE.

        Args:
            molecule: The molecule object associated with the job.
            label (str): A label for the job.
            jobrunner: The JobRunner instance to execute the job.
            If None, it will be created based on the job's TYPE.
            **kwargs: Additional keyword arguments for Job or JobRunner.

        Returns:
            Job: A new Job instance with an appropriate JobRunner.
        """
        return cls(
            molecule=molecule, label=label, jobrunner=jobrunner, **kwargs
        )
