import logging

from chemsmart.utils.mixins import BaseFolder
from chemsmart.utils.periodictable import PeriodicTable

p = PeriodicTable()
logger = logging.getLogger(__name__)


class XYZFolder(BaseFolder):
    """
    Folder handler for XYZ coordinate files.
    
    This class provides utilities for managing collections of XYZ files
    within a directory structure, including recursive search capabilities
    for batch processing of molecular coordinate files.
    """
    
    def __init__(self, folder):
        """
        Initialize XYZ folder handler.
        
        Args:
            folder (str): Path to folder containing XYZ files
        """
        self.folder = folder

    @property
    def all_xyzfiles(self):
        """
        Get all XYZ files in the folder and subfolders.
        
        Returns:
            list: Paths to all .xyz files found recursively
        """
        return self.get_all_files_in_current_folder_and_subfolders_by_suffix(
            filetype="xyz"
        )

    @property
    def all_xyzfiles_in_current_folder(self):
        """
        Get all XYZ files in the current folder only.
        
        Returns:
            list: Paths to all .xyz files in current directory
        """
        return self.get_all_files_in_current_folder_by_suffix(filetype="xyz")
