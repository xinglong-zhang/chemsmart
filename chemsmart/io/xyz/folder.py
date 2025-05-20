import logging

from chemsmart.utils.mixins import BaseFolder
from chemsmart.utils.periodictable import PeriodicTable

p = PeriodicTable()
logger = logging.getLogger(__name__)


class XYZFolder(BaseFolder):
    def __init__(self, folder):
        self.folder = folder

    @property
    def all_xyzfiles(self):
        return self.get_all_files_in_current_folder_and_subfolders_by_suffix(
            filetype="xyz"
        )

    @property
    def all_xyzfiles_in_current_folder(self):
        return self.get_all_files_in_current_folder_by_suffix(filetype="xyz")
