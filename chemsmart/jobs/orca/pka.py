"""
ORCA pKa calculation job implementation.

This module provides the ORCApKaJob class for performing pKa
calculations using ORCA with a proper thermodynamic cycle:
1. Gas phase optimization + frequency for both HA and A-
2. Solution phase single point for both HA and A- at the same level of theory

Using the same level of theory ensures proper error cancellation for
solvation free energy calculations.
"""

import logging

from chemsmart.jobs.orca.job import ORCAJob
from chemsmart.jobs.orca.opt import ORCAOptJob
from chemsmart.jobs.orca.settings import ORCApKaJobSettings
from chemsmart.jobs.orca.singlepoint import ORCASinglePointJob

logger = logging.getLogger(__name__)


class ORCApKaJob(ORCAJob):
    """
    ORCA job class for pKa calculations using the dual-level proton exchange cycle.

    Performs pKa calculations using the following workflow:
    1. Optimize HA in gas phase (opt + freq)
    2. Optimize A- in gas phase (opt + freq)
    3. Run SP on optimized HA in solution
    4. Run SP on optimized A- in solution
    5. (Optional) Same for reference acid HB and B-

    Attributes:
        TYPE (str): Job type identifier ('orcapka').
        molecule (Molecule): Protonated molecular structure (HA).
        settings (ORCApKaJobSettings): pKa calculation configuration.
        label (str): Base job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the jobs.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "orcapka"

    def __init__(
        self,
        molecule,
        settings=None,
        label=None,
        jobrunner=None,
        skip_completed=True,
        **kwargs,
    ):
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

        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            **kwargs,
        )

        self._opt_jobs = None
        self._sp_jobs = None
        self._ref_opt_jobs = None
        self._ref_sp_jobs = None

    @classmethod
    def settings_class(cls):
        return ORCApKaJobSettings

    # ------------------------------------------------------------------
    # Molecule properties
    # ------------------------------------------------------------------

    @property
    def protonated_molecule(self):
        """Get the protonated molecule (HA)."""
        protonated_mol, _ = self.settings.conjugate_pair_molecules(
            self.molecule
        )
        return protonated_mol

    @property
    def conjugate_base_molecule(self):
        """Get the conjugate base molecule (A-)."""
        _, conjugate_base_mol = self.settings.conjugate_pair_molecules(
            self.molecule
        )
        return conjugate_base_mol

    # ------------------------------------------------------------------
    # Optimization jobs
    # ------------------------------------------------------------------

    @property
    def opt_jobs(self):
        """Get gas phase optimization jobs for HA and A-."""
        if self._opt_jobs is None:
            self._opt_jobs = self._prepare_opt_jobs()
        return self._opt_jobs

    @property
    def protonated_job(self):
        """Get gas phase optimization job for HA."""
        return self.opt_jobs[0]

    @property
    def conjugate_base_job(self):
        """Get gas phase optimization job for A-."""
        return self.opt_jobs[1]

    # ------------------------------------------------------------------
    # Single-point jobs
    # ------------------------------------------------------------------

    @property
    def sp_jobs(self):
        """Get solution phase SP jobs for HA and A-."""
        if self._sp_jobs is None:
            self._sp_jobs = self._prepare_sp_jobs()
        return self._sp_jobs

    @property
    def protonated_sp_job(self):
        """Get solution phase SP job for HA."""
        return self.sp_jobs[0]

    @property
    def conjugate_base_sp_job(self):
        """Get solution phase SP job for A-."""
        return self.sp_jobs[1]

    # ------------------------------------------------------------------
    # Reference acid jobs
    # ------------------------------------------------------------------

    @property
    def has_reference_jobs(self):
        """Check if reference acid jobs are configured."""
        return self.settings.has_reference_file

    @property
    def reference_molecule(self):
        """Get the reference acid molecule (HB)."""
        return self.settings.get_reference_molecule()

    @property
    def reference_conjugate_base_molecule(self):
        """Get the reference conjugate base molecule (B-)."""
        return self.settings.get_reference_conjugate_base_molecule()

    @property
    def ref_opt_jobs(self):
        """Get gas phase optimization jobs for HB and B-."""
        if not self.has_reference_jobs:
            return None
        if self._ref_opt_jobs is None:
            self._ref_opt_jobs = self._prepare_ref_opt_jobs()
        return self._ref_opt_jobs

    @property
    def ref_acid_job(self):
        if not self.has_reference_jobs:
            return None
        return self.ref_opt_jobs[0]

    @property
    def ref_conjugate_base_job(self):
        if not self.has_reference_jobs:
            return None
        return self.ref_opt_jobs[1]

    @property
    def ref_sp_jobs(self):
        """Get solution phase SP jobs for HB and B-."""
        if not self.has_reference_jobs:
            return None
        if self._ref_sp_jobs is None:
            self._ref_sp_jobs = self._prepare_ref_sp_jobs()
        return self._ref_sp_jobs

    @property
    def ref_acid_sp_job(self):
        if not self.has_reference_jobs:
            return None
        return self.ref_sp_jobs[0]

    @property
    def ref_conjugate_base_sp_job(self):
        if not self.has_reference_jobs:
            return None
        return self.ref_sp_jobs[1]

    # ------------------------------------------------------------------
    # Job preparation
    # ------------------------------------------------------------------

    def _prepare_opt_jobs(self):
        """Create gas phase optimization jobs for HA and A-."""
        protonated_mol, conjugate_base_mol = (
            self.settings.conjugate_pair_molecules(self.molecule)
        )
        protonated_settings, conjugate_base_settings = (
            self.settings.conjugate_pair_job_settings(self.molecule)
        )

        protonated_job = ORCAOptJob(
            molecule=protonated_mol,
            settings=protonated_settings,
            label=f"{self.label}_HA",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        conjugate_base_job = ORCAOptJob(
            molecule=conjugate_base_mol,
            settings=conjugate_base_settings,
            label=f"{self.label}_A",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        return protonated_job, conjugate_base_job

    def _prepare_sp_jobs(self):
        """Create solution phase SP jobs for HA and A-."""
        protonated_sp_settings, conjugate_base_sp_settings = (
            self.settings._create_solution_phase_sp_settings(self.molecule)
        )

        # Use optimised geometry if available
        prot_out = self.protonated_job._output()
        if prot_out is not None and prot_out.normal_termination:
            protonated_mol = prot_out.molecule
        else:
            protonated_mol = self.protonated_molecule

        cb_out = self.conjugate_base_job._output()
        if cb_out is not None and cb_out.normal_termination:
            conjugate_base_mol = cb_out.molecule
        else:
            conjugate_base_mol = self.conjugate_base_molecule

        protonated_sp_job = ORCASinglePointJob(
            molecule=protonated_mol,
            settings=protonated_sp_settings,
            label=f"{self.label}_HA_sp",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        conjugate_base_sp_job = ORCASinglePointJob(
            molecule=conjugate_base_mol,
            settings=conjugate_base_sp_settings,
            label=f"{self.label}_A_sp",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        return protonated_sp_job, conjugate_base_sp_job

    def _prepare_ref_opt_jobs(self):
        """Create gas phase optimization jobs for HB and B-."""
        ref_acid_mol, ref_cb_mol = self.settings.reference_pair_molecules()
        ref_acid_settings, ref_cb_settings = (
            self.settings.reference_pair_job_settings()
        )

        ref_acid_job = ORCAOptJob(
            molecule=ref_acid_mol,
            settings=ref_acid_settings,
            label=f"{self.label}_HB",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        ref_cb_job = ORCAOptJob(
            molecule=ref_cb_mol,
            settings=ref_cb_settings,
            label=f"{self.label}_B",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        return ref_acid_job, ref_cb_job

    def _prepare_ref_sp_jobs(self):
        """Create solution phase SP jobs for HB and B-."""
        ref_acid_sp_settings, ref_cb_sp_settings = (
            self.settings.reference_pair_sp_job_settings()
        )

        ref_acid_out = self.ref_acid_job._output()
        if ref_acid_out is not None and ref_acid_out.normal_termination:
            ref_acid_mol = ref_acid_out.molecule
        else:
            ref_acid_mol = self.reference_molecule

        ref_cb_out = self.ref_conjugate_base_job._output()
        if ref_cb_out is not None and ref_cb_out.normal_termination:
            ref_cb_mol = ref_cb_out.molecule
        else:
            ref_cb_mol = self.reference_conjugate_base_molecule

        ref_acid_sp_job = ORCASinglePointJob(
            molecule=ref_acid_mol,
            settings=ref_acid_sp_settings,
            label=f"{self.label}_HB_sp",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        ref_cb_sp_job = ORCASinglePointJob(
            molecule=ref_cb_mol,
            settings=ref_cb_sp_settings,
            label=f"{self.label}_B_sp",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        return ref_acid_sp_job, ref_cb_sp_job

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------

    def _run_opt_jobs(self):
        for job in self.opt_jobs:
            logger.info(f"Running gas phase optimization job: {job}")
            job.run()

    def _run_ref_opt_jobs(self):
        if not self.has_reference_jobs:
            return
        for job in self.ref_opt_jobs:
            logger.info(f"Running reference gas phase optimization job: {job}")
            job.run()

    def _run_sp_jobs(self):
        self._sp_jobs = None  # refresh with optimised geometries
        for job in self.sp_jobs:
            logger.info(f"Running solution phase SP job: {job}")
            job.run()

    def _run_ref_sp_jobs(self):
        if not self.has_reference_jobs:
            return
        self._ref_sp_jobs = None
        for job in self.ref_sp_jobs:
            logger.info(f"Running reference solution phase SP job: {job}")
            job.run()

    def _run(self):
        self._run_opt_jobs()
        if self.has_reference_jobs:
            self._run_ref_opt_jobs()
        self._run_sp_jobs()
        if self.has_reference_jobs:
            self._run_ref_sp_jobs()

    # ------------------------------------------------------------------
    # Completion checks
    # ------------------------------------------------------------------

    def is_complete(self):
        if not all(j.is_complete() for j in self.opt_jobs):
            return False
        if self.sp_jobs is None or not all(
            j.is_complete() for j in self.sp_jobs
        ):
            return False
        if self.has_reference_jobs:
            if not all(j.is_complete() for j in self.ref_opt_jobs):
                return False
            if self.ref_sp_jobs is None or not all(
                j.is_complete() for j in self.ref_sp_jobs
            ):
                return False
        return True

    # ------------------------------------------------------------------
    # Output accessors
    # ------------------------------------------------------------------

    @property
    def protonated_output(self):
        return self.protonated_job._output()

    @property
    def conjugate_base_output(self):
        return self.conjugate_base_job._output()

    @property
    def protonated_sp_output(self):
        return self.protonated_sp_job._output()

    @property
    def conjugate_base_sp_output(self):
        return self.conjugate_base_sp_job._output()

    @property
    def ref_acid_output(self):
        if not self.has_reference_jobs:
            return None
        return self.ref_acid_job._output()

    @property
    def ref_conjugate_base_output(self):
        if not self.has_reference_jobs:
            return None
        return self.ref_conjugate_base_job._output()

    @property
    def ref_acid_sp_output(self):
        if not self.has_reference_jobs:
            return None
        return self.ref_acid_sp_job._output()

    @property
    def ref_conjugate_base_sp_output(self):
        if not self.has_reference_jobs:
            return None
        return self.ref_conjugate_base_sp_job._output()

    # ------------------------------------------------------------------
    # Thermochemistry
    # ------------------------------------------------------------------

    def get_pka_outputs(self):
        """Get ORCApKaOutput objects for completed pKa jobs."""
        from chemsmart.io.orca.output import ORCApKaOutput

        ha_file = self.protonated_job.outputfile
        a_file = self.conjugate_base_job.outputfile
        hb_file = None
        b_file = None
        if self.has_reference_jobs:
            hb_file = self.ref_acid_job.outputfile
            b_file = self.ref_conjugate_base_job.outputfile

        return ORCApKaOutput.from_pka_settings(
            settings=self.settings,
            ha_file=ha_file,
            a_file=a_file,
            hb_file=hb_file,
            b_file=b_file,
        )

    def print_thermochemistry(self):
        """Print formatted thermochemistry summary to stdout."""
        from chemsmart.io.orca.output import ORCApKaOutput

        ha_gas = self.protonated_job.outputfile
        a_gas = self.conjugate_base_job.outputfile
        ha_solv = self.protonated_sp_job.outputfile
        a_solv = self.conjugate_base_sp_job.outputfile
        hb_gas = hb_solv = b_gas = b_solv = None
        if self.has_reference_jobs:
            hb_gas = self.ref_acid_job.outputfile
            b_gas = self.ref_conjugate_base_job.outputfile
            hb_solv = self.ref_acid_sp_job.outputfile
            b_solv = self.ref_conjugate_base_sp_job.outputfile

        ORCApKaOutput.print_pka_summary(
            ha_gas_file=ha_gas,
            a_gas_file=a_gas,
            hb_gas_file=hb_gas,
            b_gas_file=b_gas,
            ha_solv_file=ha_solv,
            a_solv_file=a_solv,
            hb_solv_file=hb_solv,
            b_solv_file=b_solv,
            pka_reference=getattr(self.settings, "reference_pka", None),
            temperature=self.settings.temperature,
            concentration=self.settings.concentration,
            pressure=self.settings.pressure,
            cutoff_entropy_grimme=self.settings.cutoff_entropy_grimme,
            cutoff_enthalpy=self.settings.cutoff_enthalpy,
        )
