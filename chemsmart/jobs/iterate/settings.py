import copy
import logging

logger = logging.getLogger(__name__)


class IterateJobSettings:

    def __init__(
        self,
        config_file=None,
        method="lagrange_multipliers",
        sphere_direction_samples_num=96,
        axial_rotations_sample_num=6,
        combination_mode="independent",
    ):
        """
        Initialize iterate job settings.

        Parameters
        ----------
        config_file : str, optional
            Path to the TOML configuration file.
        method : str, optional
            Mathematical method to use for position
            optimization. Default is 'lagrange_multipliers'.
        sphere_direction_samples_num : int, optional
            Number of points to sample on the unit sphere. Default is 96.
        axial_rotations_sample_num : int, optional
            Number of axial rotations per sphere point. Default is 6.
        combination_mode : str, optional
            Combination strategy for skeleton slots.
            'independent' (default): each group generates combinations
            independently, results are merged.
            'global': all groups merge into one pool of position
            options, then a single Cartesian product is taken.
        """
        logger.debug("Initialized iterate job settings.")
        self.config_file = config_file
        self.skeleton_list: list[dict] = []
        self.substituent_list: list[dict] = []
        self.method = method
        self.sphere_direction_samples_num = sphere_direction_samples_num
        self.axial_rotations_sample_num = axial_rotations_sample_num
        self.combination_mode = combination_mode

    def copy(self):
        """
        Create a deep copy of the current settings.
        """
        return copy.deepcopy(self)
