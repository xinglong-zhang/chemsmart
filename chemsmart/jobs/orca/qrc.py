"""
ORCA quick reaction coordinate (QRC) job implementation.

Provides the ORCAQRCJob for running ORCA QRC optimizations
from TS structure file.
"""

import logging

import numpy as np

from chemsmart.jobs.batch import (
    NestableJobMixin,
    run_child_jobs_as_batch,
    run_nestable_job,
)
from chemsmart.jobs.orca.batch import ORCABatchJob
from chemsmart.jobs.orca.job import ORCAGeneralJob, ORCAJob

logger = logging.getLogger(__name__)


class ORCAQRCJob(NestableJobMixin, ORCAJob):
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
        child_index=None,
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
            child_index (int, optional): 1-based nestable child index for
                single-child array tasks.
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
        self.child_index = child_index

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
    def num_array_children(self) -> int:
        """Forward and reverse QRC children."""
        return 2

    def get_array_child_job(self, index: int):
        """Build only the forward (0) or reverse (1) QRC child."""
        self.validate_array_child_index(index)
        direction = "f" if index == 0 else "r"
        return self._prepare_qrc_job(direction)

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
        QRC reverse molecule with vibrational
        displacement at -ve given amplitude.
        """
        return self.molecule.vibrationally_displaced(
            mode_idx=self.mode_idx,
            amp=-self.amp,
            nframes=self.nframes,
            phase=self.phase,
            normalize=self.normalize,
            return_xyz=self.return_xyz,
        )

    def _prepare_qrc_job(self, direction: str):
        """Create one QRC child job for *direction* ``f`` or ``r``."""
        label = f"{self.label}{direction}"
        if self.jobtype is not None:
            label += f"_{self.jobtype}"
        if direction == "f":
            mol = self.qrcf_molecule
        elif direction == "r":
            mol = self.qrcr_molecule
        else:
            raise ValueError(
                f"Invalid QRC direction {direction!r}; expected 'f' or 'r'."
            )
        logger.debug(f"Molecule created: {mol} with label: {label}")
        return ORCAGeneralJob(
            molecule=mol,
            settings=self.settings,
            label=label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

    def _prepare_both_qrc_jobs(self):
        """
        Create both QRCr and QRCf jobs for the given TS/OPT file.
        Generates two ORCAGeneralJob objects, one for each QRC job.
        """
        jobs = [self._prepare_qrc_job(direction) for direction in ("f", "r")]
        logger.debug(f"QRC jobs created: {jobs}")
        return jobs

    def _run_both_jobs(self):
        """
        Execute both QRC jobs (forward and reverse) via ``ORCABatchJob``.

        Forward and reverse children run serially, each with the parent
        jobrunner's full resources. Failure policy is run-all-then-raise.
        """
        logger.info("Running QRC jobs using ORCABatchJob")
        run_child_jobs_as_batch(
            batch_cls=ORCABatchJob,
            jobs=self.both_qrc_jobs,
            parent=self,
            label_suffix="_batch",
            fail_fast=False,
        )

    def _run(self, **kwargs):
        """
        Execute the QRC jobs calculation.
        """
        run_nestable_job(self, self._run_both_jobs)

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
