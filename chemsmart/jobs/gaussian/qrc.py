"""
Gaussian quick reaction coordinate (QRC) job implementation.

Provides the GaussianQRCJob for running Gaussian QRC optimizations
from TS structure file.
"""

import logging

from chemsmart.jobs.gaussian.job import GaussianGeneralJob, GaussianJob

logger = logging.getLogger(__name__)


class GaussianQRCJob(GaussianJob):
    """
    Gaussian job class for quick reaction coordinate (QRC) calculations.

    QRC Reference: Tetrahedron Letters, 2003, 44, 45, 8233-8236.
    """

    TYPE = "g16qrc"

    def __init__(
        self,
        molecule,
        settings=None,
        label=None,
        jobrunner=None,
        amplitude=None,
        mode_number=None,
        skip_completed=True,
        **kwargs,
    ):
        """
        Initialize a Gaussian QRC calculation.

        Args:
            molecules (list[Molecule]): List of Molecule objects representing
                conformers to optimize.
            settings (GaussianJobSettings): Calculation configuration settings.
            label (str, optional): Base label for conformer jobs.
            jobrunner (JobRunner, optional): Job execution handler.
            amplitude: displacement amplitude along chosen vibrational mode.
            mode_number: mode number chosen to displace.
            **kwargs: Additional keyword arguments for parent class.

        Raises:
            ValueError: If molecules fails basic list validation.
        """

        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            **kwargs,
        )
        if amplitude is None:
            amplitude = 0.5  # default amplitude if not specified
        if mode_number is None:
            mode_number = 1  # default mode number if not specified

        self.amplitude = amplitude
        self.mode_number = mode_number

    @property
    def num_conformers(self):
        """
        Get the total number of conformers in the ensemble.

        Returns:
            int: Total number of conformer structures loaded.
        """
        return len(self.all_conformers)

    @property
    def last_run_job_index(self):
        """
        Get the index of the first incomplete job.

        Tracks progress through the conformer ensemble by identifying
        the next job to run. Useful for resuming interrupted
        calculations.

        Returns:
            int: Index of the first incomplete job; equals total number
                of conformers if all jobs are complete.
        """
        return self._check_last_finished_job_index()

    @property
    def all_conformers_opt_jobs(self):
        """
        Get all conformer optimization jobs in the ensemble.

        Returns:
            list: List of GaussianGeneralJob objects for all conformers.
        """
        return self._prepare_all_jobs()

    @property
    def incomplete_conformers_opt_jobs(self):
        """
        Get incomplete conformer optimization jobs.

        Filters the job list to return only those jobs that have not
        yet completed successfully. Useful for selective resubmission.

        Returns:
            list: List of incomplete GaussianGeneralJob objects.
        """
        return [
            job
            for job in self.all_conformers_opt_jobs
            if not job.is_complete()
        ]

    def _check_last_finished_job_index(self):
        """
        Find the index of the first incomplete job in the sequence.

        Iterates through all conformer jobs to identify the last
        point of progress. This is used for progress tracking and
        resuming calculations.

        Returns:
            int: Index of the first incomplete job, or total number of
                conformers if all jobs are complete.
        """
        for i, job in enumerate(self.all_conformers_opt_jobs):
            if not job.is_complete():
                return i

        # If all complete
        return self.num_conformers

    def _prepare_all_jobs(self):
        """
        Create optimization jobs for all conformers in the ensemble.

        Generates a list of GaussianGeneralJob objects, one for each
        conformer structure. Each job gets a unique label based on
        the base label and conformer index.

        Returns:
            list: List of GaussianGeneralJob objects for all conformers.
        """
        jobs = []
        for i in range(self.num_conformers):
            label = f"{self.label}_c{i + 1}"  # 1-indexed for conformers
            jobs += [
                GaussianGeneralJob(
                    molecule=self.all_conformers[i],
                    settings=self.settings,
                    label=label,
                    jobrunner=self.jobrunner,
                    skip_completed=self.skip_completed,
                )
            ]
        return jobs

    def _run_all_jobs(self):
        """
        Execute all conformer optimization jobs up to the specified limit.

        Runs the conformer optimization jobs sequentially up to the
        number specified in num_confs_to_opt. This allows for partial
        ensemble processing when needed.
        """
        for job in self.all_conformers_opt_jobs[: self.num_confs_to_opt]:
            job.run()

    def _run(self):
        """
        Execute the CREST conformer ensemble calculation.

        Main execution method that initiates all conformer optimization
        jobs. This is called internally by the job runner framework.
        """
        self._run_all_jobs()

    def is_complete(self):
        """
        Check if all conformer optimization jobs are complete.

        Returns:
            bool: True if all required conformer optimizations have
                completed successfully, False otherwise.
        """
        return self._run_all_crest_opt_jobs_are_complete()

    def _run_all_crest_opt_jobs_are_complete(self):
        """
        Verify completion status of all required conformer jobs.

        Checks the completion status of all conformer optimization
        jobs up to the specified limit (num_confs_to_opt).

        Returns:
            bool: True if all required jobs are complete, False otherwise.
        """
        return all(
            job.is_complete()
            for job in self.all_conformers_opt_jobs[: self.num_confs_to_opt]
        )
