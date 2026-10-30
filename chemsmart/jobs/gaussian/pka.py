"""
Gaussian pKa calculation job implementation.

This module provides the GaussianpKaJob class for performing pKa
calculations using Gaussian with a proper thermodynamic cycle:
1. Gas phase optimization + frequency for both HA and A-
2. Solution phase single point for both HA and A- at the same level of theory

Using the same level of theory ensures proper error cancellation for
solvation free energy calculations. Analyze completed outputs with
``chemsmart run pka``.
"""

from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.settings import GaussianpKaJobSettings
from chemsmart.jobs.gaussian.singlepoint import GaussianSinglePointJob
from chemsmart.jobs.pka import PkaChainMixin


class GaussianpKaJob(PkaChainMixin, GaussianJob):
    """
    Gaussian job class for pKa calculations using a thermodynamic cycle.

    Performs pKa calculations using the following workflow:
    1. Optimize HA in gas phase (opt + freq)
    2. Optimize A- in gas phase (opt + freq)
    3. Run SP on optimized HA in solution
    4. Run SP on optimized A- in solution

    Analyze completed outputs with ``chemsmart run pka``.

    Attributes:
        TYPE (str): Job type identifier ('g16pka').
        molecule (Molecule): Protonated molecular structure (HA).
        settings (GaussianpKaJobSettings): pKa calculation configuration.
        label (str): Base job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the jobs.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "g16pka"
    _opt_job_class = GaussianOptJob
    _sp_job_class = GaussianSinglePointJob

    def __init__(self, molecule, settings=None, **kwargs):
        if not isinstance(settings, GaussianpKaJobSettings):
            raise ValueError(
                f"Settings must be instance of GaussianpKaJobSettings for "
                f"{self.__class__.__name__}, but is {settings} instead!"
            )
        super().__init__(molecule=molecule, settings=settings, **kwargs)

    @property
    def original_mol(self):
        """Original molecule used to initialize the job (usually HA)."""
        return self.molecule

    @property
    def conjugate_base_mol(self):
        """Conjugate base molecule (A-)."""
        if self.conjugate_base_job:
            return self.conjugate_base_job.molecule
        _, conj_mol = self.settings.conjugate_pair_molecules(self.molecule)
        return conj_mol
