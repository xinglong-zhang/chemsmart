import copy
import logging

logger = logging.getLogger(__name__)


class IterateJobSettings:

    def __init__(
        self,
        config_file=None,
        algorithm="lagrange_multipliers",
        sphere_direction_samples_num=96,
        axial_rotations_sample_num=6,
    ):
        """
        Initialize iterate job settings.

        Parameters
        ----------
        config_file : str, optional
            Path to the YAML configuration file.
        algorithm : str, optional
            Algorithm to use for position optimization. Default is 'lagrange_multipliers'.
        sphere_direction_samples_num : int, optional
            Number of points to sample on the unit sphere. Default is 96.
        axial_rotations_sample_num : int, optional
            Number of axial rotations per sphere point. Default is 6.
        """
        logger.debug("Initialized iterate job settings.")
        self.config_file = config_file
        self.skeleton_list: list[dict] = []
        self.substituent_list: list[dict] = []
        self.algorithm = algorithm
        self.sphere_direction_samples_num = sphere_direction_samples_num
        self.axial_rotations_sample_num = axial_rotations_sample_num

    def copy(self):
        """
        Create a deep copy of the current settings.
        """
        return copy.deepcopy(self)
