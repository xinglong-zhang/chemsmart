"""
Gaussian coordinate scan job implementation.

This module provides the GaussianScanJob class for performing
potential energy surface scans using Gaussian. These calculations
systematically vary one or more internal coordinates while optimizing
the remaining degrees of freedom.
"""

from chemsmart.jobs.gaussian.job import GaussianJob


class GaussianScanJob(GaussianJob):
    """
    Gaussian job class for potential energy surface scan calculations.

    Performs systematic scans of potential energy surfaces by varying
    specified internal coordinates (bonds, angles, dihedrals) while
    optimizing the remaining molecular geometry. Useful for studying
    reaction paths and conformational preferences.

    Attributes:
        TYPE (str): Job type identifier ('g16scan').
        molecule (Molecule): Molecular structure used for the scan.
        settings (GaussianJobSettings): Calculation configuration, including
            scan coordinate specifications (e.g., modredundant settings).
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "g16scan"

    def __init__(self, molecule, settings, label, **kwargs):
        """
        Initialize a Gaussian coordinate scan job.

        Sets up a potential energy surface (PES) scan using the supplied
        molecular structure and scan parameters. Scan coordinates are
        typically defined via Gaussian's modredundant syntax in
        `settings.modred` and the route (e.g., `opt=modredundant`).

        Common patterns for `settings.modred`:
        - Constraint optimization (list): `[[i,j], [k,l,m], ...]`
        - Coordinate scan (dict): `{"num_steps": N, "step_size": dx, "coords": [[...], ...]}`

        Args:
            molecule (Molecule): Molecular structure for scanning.
            settings (GaussianJobSettings): Calculation configuration,
                including scan coordinate specifications.
            label (str): Job identifier for file naming.
            **kwargs: Additional keyword arguments for the parent class
                (e.g., `jobrunner`, `skip_completed`).

        Raises:
            ValueError: If `molecule` is not a Molecule or `settings` is not
                a GaussianJobSettings instance (validated by the base class).
        """
        super().__init__(
            molecule=molecule, settings=settings, label=label, **kwargs
        )
