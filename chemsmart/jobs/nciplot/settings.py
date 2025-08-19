import logging

logger = logging.getLogger(__name__)


class NCIPLOTJobSettings:
    """Settings for NCIPLOTJob."""

    def __init__(
        self,
        filenames,
        label,
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
        self.filenames = filenames
        self.label = label
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
        return NCIPLOTJobSettings(
            filenames=self.filenames,
            label=self.label,
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
        return cls(**settings_dict)
