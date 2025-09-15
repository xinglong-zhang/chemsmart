"""
NCIPLOT job settings module.

This module contains the NCIPLOTJobSettings class for configuring
NCIPLOT calculations with various analysis parameters.
"""

import logging

logger = logging.getLogger(__name__)


class NCIPLOTJobSettings:
    """
    Settings for NCIPLOT job configuration.

    This class handles all configuration parameters for NCIPLOT calculations
    including thresholds, cutoffs, grid settings, and analysis options.

    Attributes:
        rthres (float): Density threshold for analysis.
        ligand_file_number (int): File number for ligand identification.
        ligand_radius (float): Radius parameter for ligand analysis.
        radius_positions (list | tuple): Spatial positions for radius analysis.
        radius_r (float): Radius value for analysis.
        intercut1 (float): First intermolecular cutoff parameter.
        intercut2 (float): Second intermolecular cutoff parameter.
        increments (list | tuple): Grid increment values.
        fragments (dict): Fragment definitions mapping names/ids to atoms.
        cutoff_density_dat (float): Density cutoff for .dat outputs.
        cutoff_rdg_dat (float): RDG cutoff for .dat outputs.
        cutoff_density_cube (float): Density cutoff for .cube outputs.
        cutoff_rdg_cube (float): RDG cutoff for .cube outputs.
        dgrid (bool): Whether to enable DGRID option.
        integrate (bool): Whether to perform integration analysis.
        ranges (list[list] | tuple): Ranges for integration/analysis.
        grid_quality (str | int): Grid quality setting.
    """

    def __init__(
        self,
        rthres,
        ligand_file_number,
        ligand_radius,
        radius_positions,
        radius_r,
        intercut1,
        intercut2,
        increments,
        fragments,  # dictionary
        cutoff_density_dat,
        cutoff_rdg_dat,
        cutoff_density_cube,
        cutoff_rdg_cube,
        dgrid,
        integrate,
        ranges,  # list of lists
        grid_quality,
    ):
        """
        Initialize NCIPLOT job settings.

        Args:
            rthres: Density threshold for analysis
            ligand_file_number: File number for ligand identification
            ligand_radius: Radius parameter for ligand analysis
            radius_positions: Spatial positions for radius analysis
            radius_r: Radius value for analysis
            intercut1: First intermolecular cutoff parameter
            intercut2: Second intermolecular cutoff parameter
            increments: Grid increment values
            fragments: Dictionary defining molecular fragments
            cutoff_density_dat: Density cutoff for .dat files
            cutoff_rdg_dat: RDG cutoff for .dat files
            cutoff_density_cube: Density cutoff for .cube files
            cutoff_rdg_cube: RDG cutoff for .cube files
            dgrid: Enable DGRID option
            integrate: Enable integration analysis
            ranges: List of analysis ranges
            grid_quality: Quality setting for grid generation
        """
        # Store all configuration parameters
        self.rthres = rthres
        self.ligand_file_number = ligand_file_number
        self.ligand_radius = ligand_radius
        self.radius_positions = radius_positions
        self.radius_r = radius_r
        self.intercut1 = intercut1
        self.intercut2 = intercut2
        self.increments = increments
        self.fragments = fragments
        self.cutoff_density_dat = cutoff_density_dat
        self.cutoff_rdg_dat = cutoff_rdg_dat
        self.cutoff_density_cube = cutoff_density_cube
        self.cutoff_rdg_cube = cutoff_rdg_cube
        self.dgrid = dgrid
        self.integrate = integrate
        self.ranges = ranges
        self.grid_quality = grid_quality

    def copy(self):
        """
        Create a deep copy of the settings.

        Returns:
            NCIPLOTJobSettings: New instance with copied parameters
        """
        return NCIPLOTJobSettings(
            rthres=self.rthres,
            ligand_file_number=self.ligand_file_number,
            ligand_radius=self.ligand_radius,
            radius_positions=self.radius_positions,
            radius_r=self.radius_r,
            intercut1=self.intercut1,
            intercut2=self.intercut2,
            increments=self.increments,
            fragments=self.fragments,
            cutoff_density_dat=self.cutoff_density_dat,
            cutoff_rdg_dat=self.cutoff_rdg_dat,
            cutoff_density_cube=self.cutoff_density_cube,
            cutoff_rdg_cube=self.cutoff_rdg_cube,
            dgrid=self.dgrid,
            integrate=self.integrate,
            ranges=self.ranges,
            grid_quality=self.grid_quality,
        )

    @classmethod
    def from_dict(cls, settings_dict):
        """
        Create settings instance from dictionary.

        Args:
            settings_dict: Dictionary containing setting parameters

        Returns:
            NCIPLOTJobSettings: New instance created from dictionary
        """
        return cls(**settings_dict)
