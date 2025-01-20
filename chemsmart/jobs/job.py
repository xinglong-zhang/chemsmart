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
        skip_completed (bool): If True, completed jobs will not be rerun. Defaults to True.
    """

    TYPE: Optional[str] = None
    PROGRAM: Optional[str] = None

    def __init__(self, molecule, label, local=False, skip_completed=True):
        self.molecule = molecule
        self.label = label
        self.local = local
        self.skip_completed = skip_completed

    @property
    def folder(self):
        return self._determine_folder()

    def _determine_folder(self):
        """
        Determine the folder based on the current working directory where the job is submitted.
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
        return f"{self.__class__.__qualname__}<folder={self.folder}, label={self.label}>"

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
        raise NotImplementedError

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
        else:  # noqa: PLR5501
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
