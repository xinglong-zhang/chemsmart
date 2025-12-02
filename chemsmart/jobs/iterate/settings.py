import copy
import logging

logger = logging.getLogger(__name__)


class IterateJobSettings:

    def __init__(
        self,
        config_file=None,
        algorithm='lagrange_multipliers',
    ):
        """
        Initialize iterate job settings.
        
        Parameters
        ----------
        config_file : str, optional
            Path to the YAML configuration file.
        algorithm : str, optional
            Algorithm to use for position optimization. Default is 'lagrange_multipliers'.
        """
        logger.debug("Initialized iterate job settings.")
        self.config_file = config_file
        self.skeleton_list: list[dict] = []
        self.substituent_list: list[dict] = []
        self.algorithm = algorithm
    
    def copy(self):
        """
        Create a deep copy of the current settings.
        """
        return copy.deepcopy(self)