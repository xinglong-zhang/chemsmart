"""
ORCA quick reaction coordinate (QRC) job implementation.

Provides the ORCAQRCJob for running ORCA QRC optimizations
from TS structure file.
"""

import logging

import numpy as np

from chemsmart.jobs.orca.job import ORCAGeneralJob, ORCAJob

logger = logging.getLogger(__name__)


class ORCAQRCJob(ORCAJob):
    """
    ORCA job class for quick reaction coordinate (QRC) calculations.

    QRC Reference: Tetrahedron Letters, 2003, 44, 45, 8233-8236.
    """

    TYPE = "orcaqrc"

    def __init__(
        self,
        molecule,
        settings=None,
        label=None,
        jobrunner=None,
        mode_idx=1,
        amp=0.5,
        nframes=None,
        phase=np.pi / 2,
        normalize=False,
        return_xyz=False,
        skip_completed=True,
        **kwargs,
    ):
        """
        Initialize an ORCA QRC calculation.

        Args:
            molecule: Molecule objects representing the structure for QRC run.
            settings (ORCAJobSettings): Calculation configuration settings.
            label (str, optional): Base label for conformer jobs.
            jobrunner (JobRunner, optional): Job execution handler.
            mode_idx (int): vibrational mode number to be used.
            amp (float): amplitude of displacement along the chosen mode.
            nframes (int): number of image frames to create.
            phase (float): phase of displacement.
            normalize (bool): If true, will normalize displacement such that
                max displacement is 1.
            return_xyz (bool): If True and `nframes` is set, return
                a multi-frame XYZ string.
            **kwargs: Additional keyword arguments for parent class.

        Raises:
            ValueError: If molecule does not have vibrational modes.
        """

        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            **kwargs,
        )
        if not molecule.has_vibrations:
            raise ValueError(f"Molecule {molecule} has no vibrational modes!")
        self.mode_idx = mode_idx
        self.amp = amp
        self.nframes = nframes
        self.phase = phase
        self.normalize = normalize
        self.return_xyz = return_xyz

        self.jobtype = self.settings.__dict__["jobtype"]

    @property
    def both_qrc_jobs(self):
        """
        Both QRC forward and reverse jobs.

        Returns:
            list: List of ORCAGeneralJob objects for QRC.
        """
        return self._prepare_both_qrc_jobs()

    @property
    def qrcf_molecule(self):
        """
        QRC forward molecule with vibrational displacement at given amplitude.
        """
        return self.molecule.vibrationally_displaced(
            mode_idx=self.mode_idx,
            amp=self.amp,
            nframes=self.nframes,
            phase=self.phase,
            normalize=self.normalize,
            return_xyz=self.return_xyz,
        )

    @property
    def qrcr_molecule(self):
        """
        QRC reverse molecule with vibrational displacement at -ve given amplitude.
        """
        return self.molecule.vibrationally_displaced(
            mode_idx=self.mode_idx,
            amp=-self.amp,
            nframes=self.nframes,
            phase=self.phase,
            normalize=self.normalize,
            return_xyz=self.return_xyz,
        )

    def _prepare_both_qrc_jobs(self):
        """
        Create both QRCr and QRCf jobs for the given TS/OPT file.
        Generates two GaussianGeneralJob objects, one for each QRC job.
        """
        jobs = []
        for direction in ["f", "r"]:
            label = f"{self.label}{direction}"
            if self.jobtype is not None:
                label += f"_{self.jobtype}"
            mol = self.molecule
            if direction == "f":
                mol = self.qrcf_molecule
            elif direction == "r":
                mol = self.qrcr_molecule
            logger.debug(f"Molecule created: {mol} with label: {label}")
            jobs.append(
                ORCAGeneralJob(
                    molecule=mol,
                    settings=self.settings,
                    label=label,
                    jobrunner=self.jobrunner,
                    skip_completed=self.skip_completed,
                )
            )
        logger.debug(f"QRC jobs created: {jobs}")
        return jobs

    def _run_both_jobs(self):
        """
        Execute both QRC jobs (forward and reverse) sequentially.
        Runs the QRC forward and reverse jobs for the current molecule.
        """
        # Check if jobs should be run in serial based on jobrunner flag
        if self.jobrunner and self.jobrunner.run_in_serial:
            logger.info("Running QRC jobs in serial mode (one after another)")
        else:
            logger.info("Running QRC jobs using default behavior")

        for job in self.both_qrc_jobs:
            job.run()

    def _run(self):
        """
        Execute the QRC jobs calculation.
        """
        self._run_both_jobs()

    def is_complete(self):
        """
        Check if both QRC jobs are complete.
        """
        return self._both_qrc_jobs_are_complete()

    def _both_qrc_jobs_are_complete(self):
        """
        Verify completion status of both QRC jobs.
        """
        return all(job.is_complete() for job in self.both_qrc_jobs)
