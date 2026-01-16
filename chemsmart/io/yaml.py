import logging

from chemsmart.utils.mixins import YAMLFileMixin

logger = logging.getLogger(__name__)


class YAMLFile(YAMLFileMixin):
    """
    A class for handling YAML files with file I/O operations.

    This class provides functionality to read and write YAML files,
    inheriting common YAML operations from YAMLFileMixin.
    """

    def __init__(self, filename):
        """
        Initialize YAMLFile with a filename.

        Args:
            filename (str): Path to the YAML file to be processed.
        """
        self.filename = filename
