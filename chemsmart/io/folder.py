from chemsmart.utils.mixins import FolderMixin


class BaseFolder(FolderMixin):
    """
    Base class for folder management with path operations.

    Provides basic folder initialization and validation functionality.
    Combines folder path management with file discovery capabilities
    from the FolderMixin class.
    """

    def __init__(self, folder):
        """
        Initialize the BaseFolder with a specified path.

        Args:
            folder (str): Path to the folder to be managed.
        """
        self.folder = folder
