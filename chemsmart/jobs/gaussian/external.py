"""
Gaussian External calculator job for ML potential energy calculations.

This module provides GaussianExternalJob, which drives Gaussian geometry
optimisations and single-point calculations using an ASE-compatible ML
potential via Gaussian's External= keyword.
"""

import logging
import os

from chemsmart.jobs.gaussian.job import GaussianJob

logger = logging.getLogger(__name__)


class GaussianExternalJob(GaussianJob):
    """
    Gaussian job that uses the External= keyword to call an ASE calculator.

    Wraps an ASEExternalCalculatorScript so that the generated Python
    script is written alongside the Gaussian .com input file before the
    job is submitted.

    Attributes:
        TYPE (str): Job type identifier ('g16external').
        external_script (ASEExternalCalculatorScript): Script writer that
            generates the callable Python script Gaussian will invoke.
    """

    TYPE = "g16external"

    def __init__(
        self,
        molecule,
        settings,
        label,
        external_script=None,
        jobrunner=None,
        **kwargs,
    ):
        """
        Initialise a Gaussian External calculator job.

        Args:
            molecule (Molecule): Molecular structure for the calculation.
            settings (GaussianExternalJobSettings): Job configuration
                including the external_script name.
            label (str): Job identifier used for file naming.
            external_script (ASEExternalCalculatorScript, optional):
                Script writer instance. When provided, calling
                write_external_script() writes the callable Python
                script to the job directory.
            jobrunner (JobRunner, optional): Job execution handler.
            **kwargs: Additional keyword arguments for the parent class.
        """
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
        self.external_script = external_script

    def write_external_script(self, directory=None):
        """
        Write the ASE external calculator script to the job directory.

        If no ASEExternalCalculatorScript was provided at construction
        time, this method is a no-op.

        Args:
            directory (str, optional): Target directory. Defaults to
                the job's working folder.

        Returns:
            str or None: Path to the written script file, or None if
                no script writer was configured.
        """
        if self.external_script is None:
            logger.debug(
                "No ASEExternalCalculatorScript attached; "
                "skipping external script write."
            )
            return None

        target_dir = directory or self.folder
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)
            logger.debug(f"Created directory for external script: {target_dir}")

        script_path = self.external_script.write(directory=target_dir)
        logger.info(f"Wrote Gaussian External script to: {script_path}")
        return script_path
