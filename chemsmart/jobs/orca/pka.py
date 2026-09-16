"""
ORCA pKa calculation job implementation.

This module provides the ORCApKaJob class for performing pKa
calculations using ORCA with a proper thermodynamic cycle:
1. Gas phase optimization + frequency for both HA and A-
2. Solution phase single point for both HA and A- at the same level of theory

Using the same level of theory ensures proper error cancellation for
solvation free energy calculations. Analyze completed outputs with
``chemsmart run pka``.
"""

import os

from chemsmart.jobs.orca.job import ORCAJob
from chemsmart.jobs.orca.opt import ORCAOptJob
from chemsmart.jobs.orca.settings import ORCApKaJobSettings
from chemsmart.jobs.orca.singlepoint import ORCASinglePointJob
from chemsmart.jobs.pka import PkaChainMixin


class ORCApKaJob(PkaChainMixin, ORCAJob):
    """
    ORCA job class for pKa calculations using the dual-level proton exchange cycle.

    Performs pKa calculations using the following workflow:
    1. Optimize HA in gas phase (opt + freq)
    2. Optimize A- in gas phase (opt + freq)
    3. Run SP on optimized HA in solution
    4. Run SP on optimized A- in solution
    5. (Optional) Same for reference acid Href and Ref-

    Analyze completed outputs with ``chemsmart run pka``.

    Attributes:
        TYPE (str): Job type identifier ('orcapka').
        molecule (Molecule): Protonated molecular structure (HA).
        settings (ORCApKaJobSettings): pKa calculation configuration.
        label (str): Base job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the jobs.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "orcapka"
    _opt_job_class = ORCAOptJob
    _sp_job_class = ORCASinglePointJob

    def __init__(self, molecule, settings=None, **kwargs):
        if not isinstance(settings, ORCApKaJobSettings):
            raise ValueError(
                f"Settings must be instance of ORCApKaJobSettings, "
                f"but got {type(settings).__name__} instead!"
            )

        if settings.proton_index is None:
            raise ValueError(
                "proton_index must be specified in ORCApKaJobSettings "
                "to identify which proton to remove for the conjugate base."
            )

        super().__init__(molecule=molecule, settings=settings, **kwargs)

    @classmethod
    def settings_class(cls):
        return ORCApKaJobSettings

    @property
    def _ref_basename(self):
        """Basename for the reference acid, from the reference geometry file."""
        if not self.settings.has_reference_file:
            return None
        return os.path.splitext(
            os.path.basename(self.settings.reference_file)
        )[0]

    @property
    def _ref_conjugate_base_label(self):
        """Label for the reference conjugate base (Ref⁻)."""
        ref = self._ref_basename
        if ref is None:
            return None
        return f"{ref}_cb"

    def _href_opt_label(self):
        return self._ref_basename

    def _ref_opt_label(self):
        return self._ref_conjugate_base_label

    def _href_sp_label(self):
        ref = self._ref_basename
        if ref is None:
            return None
        return f"{ref}_sp"

    def _ref_sp_label(self):
        ref_cb = self._ref_conjugate_base_label
        if ref_cb is None:
            return None
        return f"{ref_cb}_sp"
