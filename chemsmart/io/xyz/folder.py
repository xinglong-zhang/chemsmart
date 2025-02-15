import logging

from chemsmart.utils.mixins import FolderMixin
from chemsmart.utils.periodictable import PeriodicTable

p = PeriodicTable()
logger = logging.getLogger(__name__)


class XYZFolder(FolderMixin):
    def __init__(self, folder):
        self.folder = folder

    @property
    def all_xyzfiles(self):
        return self.get_all_files_in_current_folder_and_subfolders(
            filetype="xyz"
        )

    @property
    def all_xyzfiles_in_current_folder(self):
        return self.get_all_files_in_current_folder(filetype="xyz")
