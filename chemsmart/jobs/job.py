import datetime
import glob
import logging
import os
import time
import shutil
from contextlib import suppress
from abc import abstractmethod
from chemsmart.utils.mixins import RegistryMixin

logger = logging.getLogger(__name__)


class Job(RegistryMixin):
    """Class that encapsulates and runs a task.

    Args:
        skip_completed (bool): If True, completed jobs will not be rerun. Defaults to True.
    """

    TYPE = NotImplemented
    PROGRAM = NotImplemented
    # ALIASES = []
    # RUNNERS = []

    def __init__(self, molecule, label, local=False, skip_completed=True):
        self.molecule = molecule
        self.label = label
        self.local = local
        self.skip_completed = skip_completed
        self.folder = self._determine_folder()

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

    @abstractmethod
    def is_complete(self):
        raise NotImplementedError

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

    def backup(self, dest_folder=None):
        """Backup the current folder to a new folder.
        Return the new folder path.
        """
        if self.folder is None:
            return None

        shutil.copytree(src=self.folder, dst=dest_folder)
        return dest_folder

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
