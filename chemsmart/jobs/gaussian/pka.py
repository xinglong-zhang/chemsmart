"""
Gaussian pKa calculation job implementation.

This module provides the GaussianpKaJob class for performing pKa
calculations using Gaussian with a proper thermodynamic cycle:
1. Gas phase optimization + frequency for both HA and A-
2. Solution phase single point for both HA and A- at the same level of theory

Using the same level of theory ensures proper error cancellation for
solvation free energy calculations.
"""

import logging
import os

from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.settings import GaussianpKaJobSettings
from chemsmart.jobs.gaussian.singlepoint import GaussianSinglePointJob
from chemsmart.jobs.runner import decide_phase_transition, run_phase_jobs

logger = logging.getLogger(__name__)


class GaussianpKaJob(GaussianJob):
    """
    Gaussian job class for pKa calculations using direct thermodynamic cycle.

    Performs pKa calculations using the following workflow:
    1. Optimize HA in gas phase (opt + freq) - get G(HA)_gas
    2. Optimize A- in gas phase (opt + freq) - get G(A-)_gas
    3. Run SP on optimized HA in solution - get E(HA)_aq
    4. Run SP on optimized A- in solution - get E(A-)_aq
    5. Calculate solvation free energies and pKa

    Attributes:
        TYPE (str): Job type identifier ('g16pka').
        molecule (Molecule): Protonated molecular structure (HA).
        settings (GaussianpKaJobSettings): pKa calculation configuration.
        label (str): Base job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the jobs.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "g16pka"
    _shared_reference_molecule_cache = {}

    def _reference_cache_key(cls, settings):
        if settings is None or not settings.has_reference_file:
            return None
        return (
            settings.scheme,
            os.path.abspath(settings.reference_file),
            settings.reference_proton_index,
            settings.reference_charge,
            settings.reference_multiplicity,
            settings.reference_conjugate_base_charge,
            settings.reference_conjugate_base_multiplicity,
        )

    def _get_cached_reference_pair(cls, settings):
        cache_key = cls._reference_cache_key(settings)
        if cache_key is None:
            return None
        if cache_key not in cls._shared_reference_molecule_cache:
            cls._shared_reference_molecule_cache[cache_key] = (
                settings.reference_pair_molecules()
            )
        return cls._shared_reference_molecule_cache[cache_key]

    def __init__(
        self,
        molecule,
        settings=None,
        label=None,
        jobrunner=None,
        skip_completed=True,
        **kwargs,
    ):
        if not isinstance(settings, GaussianpKaJobSettings):
            raise ValueError(
                f"Settings must be instance of GaussianpKaJobSettings for {self.__class__.__name__}, but is {settings} instead!"
            )

        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            **kwargs,
        )

        self.opt_jobs = []
        self.ref_opt_jobs = []
        self.sp_jobs = None
        self.ref_sp_jobs = None

        # Target acid jobs
        self.protonated_job = None
        self.conjugate_base_job = None
        self.protonated_sp_job = None
        self.conjugate_base_sp_job = None

        # Reference acid jobs
        self.ref_acid_job = None
        self.ref_conjugate_base_job = None
        self.ref_acid_sp_job = None
        self.ref_conjugate_base_sp_job = None

        # Check existing reference jobs
        self.has_reference_jobs = (
            self.settings.has_reference_file if self.settings else False
        )

        self._prepare_pka_jobs()

    def _prepare_pka_jobs(self):
        """Prepare optimization jobs for target and reference acids."""
        if self.settings is None:
            return

        # 1. Target Acid (HA / A-)
        prot_opt_settings, conj_opt_settings = (
            self.settings._create_gas_phase_job_settings(self.molecule)
        )
        prot_mol, conj_mol = self.settings.conjugate_pair_molecules(
            self.molecule
        )

        self.protonated_job = GaussianOptJob(
            molecule=prot_mol,
            settings=prot_opt_settings,
            label=f"{self.label}_HA_opt",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        self.conjugate_base_job = GaussianOptJob(
            molecule=conj_mol,
            settings=conj_opt_settings,
            label=f"{self.label}_A_opt",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        self.opt_jobs = [self.protonated_job, self.conjugate_base_job]

        # 2. Reference Acid (HRef / Ref-)
        if self.has_reference_jobs:
            self.ref_opt_jobs = self._prepare_ref_opt_jobs()
            self.ref_acid_job, self.ref_conjugate_base_job = self.ref_opt_jobs

    def _prepare_ref_opt_jobs(self):
        """Prepare gas phase optimization jobs for HRef and Ref-."""
        reference_pair = self._get_cached_reference_pair(self.settings)
        if reference_pair is None:
            ref_acid_mol, ref_conjugate_base_mol = (
                self.settings.reference_pair_molecules()
            )
        else:
            ref_acid_mol, ref_conjugate_base_mol = reference_pair
        ref_acid_settings, ref_conjugate_base_settings = (
            self.settings.reference_pair_job_settings()
        )

        ref_acid_job = GaussianOptJob(
            molecule=ref_acid_mol,
            settings=ref_acid_settings,
            label=f"{self.label}_HRef_opt",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        ref_conjugate_base_job = GaussianOptJob(
            molecule=ref_conjugate_base_mol,
            settings=ref_conjugate_base_settings,
            label=f"{self.label}_Ref_opt",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        return [ref_acid_job, ref_conjugate_base_job]

    def _optimized_molecule_from_job(self, job, fallback_molecule):
        """Return optimized geometry for a finished job, or a fallback molecule."""
        out = job._output()
        if out is not None and out.normal_termination:
            return out.molecule
        return fallback_molecule

    def _run_opt_jobs(self):
        """Run gas phase optimization jobs."""
        run_phase_jobs(
            parent_runner=self.jobrunner,
            jobs=self.opt_jobs,
            stop_on_incomplete=True,
            logger_obj=logger,
            phase_label="gas phase optimization",
        )

    def _run_ref_opt_jobs(self):
        """Run reference gas phase optimization jobs."""
        if not self.has_reference_jobs:
            return
        run_phase_jobs(
            parent_runner=self.jobrunner,
            jobs=self.ref_opt_jobs,
            stop_on_incomplete=True,
            logger_obj=logger,
            phase_label="reference gas phase optimization",
        )

    def _run_sp_jobs(self):
        """Run solution phase single point jobs using optimized geometries."""
        if not self._opt_jobs_are_complete():
            logger.warning(
                "Optimization jobs not complete. Cannot run SP jobs."
            )
            return

        # Create SP jobs if not already created
        if self.sp_jobs is None:
            self._create_sp_jobs()

        run_phase_jobs(
            parent_runner=self.jobrunner,
            jobs=self.sp_jobs,
            stop_on_incomplete=True,
            logger_obj=logger,
            phase_label="solution phase SP",
        )

    def _run_ref_sp_jobs(self):
        """Run reference solution phase single point jobs."""
        if not self.has_reference_jobs:
            return
        if not self._ref_opt_jobs_are_complete():
            logger.warning(
                "Reference optimization jobs not complete. Cannot run reference SP jobs."
            )
            return

        if self.ref_sp_jobs is None:
            self._create_ref_sp_jobs()

        run_phase_jobs(
            parent_runner=self.jobrunner,
            jobs=self.ref_sp_jobs,
            stop_on_incomplete=True,
            logger_obj=logger,
            phase_label="reference solution phase SP",
        )

    def _create_sp_jobs(self):
        """Create solution phase SP jobs from optimized geometries."""
        _, conj_fallback_mol = self.settings.conjugate_pair_molecules(
            self.molecule
        )
        prot_opt_mol = self._optimized_molecule_from_job(
            self.protonated_job, self.molecule
        )
        conj_opt_mol = self._optimized_molecule_from_job(
            self.conjugate_base_job, conj_fallback_mol
        )

        prot_sp_settings, conj_sp_settings = (
            self.settings._create_solution_phase_sp_settings(self.molecule)
        )

        self.protonated_sp_job = GaussianSinglePointJob(
            molecule=prot_opt_mol,
            settings=prot_sp_settings,
            label=f"{self.label}_HA_sp",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        self.conjugate_base_sp_job = GaussianSinglePointJob(
            molecule=conj_opt_mol,
            settings=conj_sp_settings,
            label=f"{self.label}_A_sp",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        self.sp_jobs = [self.protonated_sp_job, self.conjugate_base_sp_job]

    def _create_ref_sp_jobs(self):
        """Create reference solution phase SP jobs from optimized geometries."""
        reference_pair = self._get_cached_reference_pair(self.settings)
        if reference_pair is None:
            ref_acid_fallback_mol, ref_conjugate_base_fallback_mol = (
                self.settings.reference_pair_molecules()
            )
        else:
            ref_acid_fallback_mol, ref_conjugate_base_fallback_mol = (
                reference_pair
            )
        ref_acid_opt_mol = self._optimized_molecule_from_job(
            self.ref_acid_job, ref_acid_fallback_mol
        )
        ref_conjugate_base_opt_mol = self._optimized_molecule_from_job(
            self.ref_conjugate_base_job, ref_conjugate_base_fallback_mol
        )
        ref_acid_sp_settings, ref_conjugate_base_sp_settings = (
            self.settings.reference_pair_sp_job_settings()
        )

        self.ref_acid_sp_job = GaussianSinglePointJob(
            molecule=ref_acid_opt_mol,
            settings=ref_acid_sp_settings,
            label=f"{self.label}_HRef_sp",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        self.ref_conjugate_base_sp_job = GaussianSinglePointJob(
            molecule=ref_conjugate_base_opt_mol,
            settings=ref_conjugate_base_sp_settings,
            label=f"{self.label}_Ref_sp",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        self.ref_sp_jobs = [
            self.ref_acid_sp_job,
            self.ref_conjugate_base_sp_job,
        ]

    def _run(self, **kwargs):
        """
        Execute the pKa calculation.

        Runs gas phase optimization jobs followed by solution phase SP jobs.
        """
        # Default sequential behaviour
        # Run gas phase optimization jobs for target acid (HA, A-)
        self._run_opt_jobs()

        opt_transition = decide_phase_transition(
            phase_name="Opt",
            require_complete=True,
            is_complete=self._opt_jobs_are_complete(),
            stop_message="Opt jobs incomplete, halting serial execution.",
        )
        if not opt_transition.proceed:
            logger.info(opt_transition.message)
            return

        # Run gas phase optimization jobs for reference acid (HRef, Ref-) if provided
        if self.has_reference_jobs:
            self._run_ref_opt_jobs()
            ref_opt_transition = decide_phase_transition(
                phase_name="Ref Opt",
                require_complete=True,
                is_complete=self._ref_opt_jobs_are_complete(),
                stop_message="Ref Opt jobs incomplete, halting serial execution.",
            )
            if not ref_opt_transition.proceed:
                logger.info(ref_opt_transition.message)
                return

        # Run solution phase SP jobs for target acid
        self._run_sp_jobs()

        sp_transition = decide_phase_transition(
            phase_name="SP",
            require_complete=True,
            is_complete=self._sp_jobs_are_complete(),
            stop_message="SP jobs incomplete, halting serial execution.",
        )
        if not sp_transition.proceed:
            logger.info(sp_transition.message)
            return

        # Run solution phase SP jobs for reference acid if provided
        if self.has_reference_jobs:
            self._run_ref_sp_jobs()

    def is_complete(self):
        """
        Check if all pKa jobs are complete.

        Returns:
            bool: True if all optimization jobs and SP jobs
                have completed successfully (including reference jobs if provided).
        """
        # Check target acid optimization jobs
        if not self._opt_jobs_are_complete():
            return False

        # Check target acid SP jobs
        if not self._sp_jobs_are_complete():
            return False

        # Check reference acid jobs if provided
        if self.has_reference_jobs:
            if not self._ref_opt_jobs_are_complete():
                return False
            if not self._ref_sp_jobs_are_complete():
                return False

        return True

    def _opt_jobs_are_complete(self):
        """
        Verify completion status of both gas phase optimization jobs.

        Returns:
            bool: True if all optimization jobs are complete.
        """
        if not self.opt_jobs:
            return False
        return all(job.is_complete() for job in self.opt_jobs)

    def _ref_opt_jobs_are_complete(self):
        """
        Verify completion status of both reference gas phase optimization jobs.

        Returns:
            bool: True if all reference optimization jobs are complete,
                or True if no reference jobs are configured.
        """
        if not self.has_reference_jobs:
            return True
        if not self.ref_opt_jobs:
            return False
        return all(job.is_complete() for job in self.ref_opt_jobs)

    def _sp_jobs_are_complete(self):
        """
        Verify completion status of both solution phase SP jobs.

        Returns:
            bool: True if all SP jobs are complete.
        """
        if not self.sp_jobs:
            return False
        return all(job.is_complete() for job in self.sp_jobs)

    def _ref_sp_jobs_are_complete(self):
        """
        Verify completion status of both reference solution phase SP jobs.

        Returns:
            bool: True if all reference SP jobs are complete,
                or True if no reference jobs are configured.
        """
        if not self.has_reference_jobs:
            return True
        if not self.ref_sp_jobs:
            return False
        return all(job.is_complete() for job in self.ref_sp_jobs)

    @property
    def protonated_output(self):
        """
        Get the output of the protonated gas phase optimization job.

        Returns:
            Gaussian16Output or None: Parsed output for the protonated job.
        """
        return self.protonated_job._output()

    @property
    def conjugate_base_output(self):
        """
        Get the output of the conjugate base gas phase optimization job.

        Returns:
            Gaussian16Output or None: Parsed output for the conjugate base job.
        """
        return self.conjugate_base_job._output()

    @property
    def protonated_sp_output(self):
        """
        Get the output of the protonated solution phase SP job.

        Returns:
            Gaussian16Output or None: Parsed output for the protonated SP job.
        """
        return self.protonated_sp_job._output()

    @property
    def conjugate_base_sp_output(self):
        """
        Get the output of the conjugate base solution phase SP job.

        Returns:
            Gaussian16Output or None: Parsed output for the conjugate base SP job.
        """
        return self.conjugate_base_sp_job._output()

    # =========================================================================
    # Reference acid output properties
    # =========================================================================

    @property
    def ref_acid_output(self):
        """
        Get the output of the reference acid gas phase optimization job.

        Returns:
            Gaussian16Output or None: Parsed output for the reference acid job.
        """
        if not self.has_reference_jobs:
            return None
        return self.ref_acid_job._output()

    @property
    def ref_conjugate_base_output(self):
        """
        Get the output of the reference conjugate base gas phase optimization job.

        Returns:
            Gaussian16Output or None: Parsed output for the reference conjugate base job.
        """
        if not self.has_reference_jobs:
            return None
        return self.ref_conjugate_base_job._output()

    @property
    def ref_acid_sp_output(self):
        """
        Get the output of the reference acid solution phase SP job.

        Returns:
            Gaussian16Output or None: Parsed output for the reference acid SP job.
        """
        if not self.has_reference_jobs:
            return None
        return self.ref_acid_sp_job._output()

    @property
    def ref_conjugate_base_sp_output(self):
        """
        Get the output of the reference conjugate base solution phase SP job.

        Returns:
            Gaussian16Output or None: Parsed output for the reference conjugate base SP job.
        """
        if not self.has_reference_jobs:
            return None
        return self.ref_conjugate_base_sp_job._output()

    # =========================================================================
    # Thermochemistry extraction
    # =========================================================================

    def _pka_output_files(self):
        """Return (ha_file, a_file, href_file, ref_file) output-path tuple.

        Centralises the repeated logic of resolving output file paths from
        the four pKa species jobs.  Reference paths (href_file, ref_file) are
        only populated when reference jobs exist *and* their opt phase is
        complete; otherwise they are ``None``.
        """
        ha_file = (
            self.protonated_job.outputfile if self.protonated_job else None
        )
        a_file = (
            self.conjugate_base_job.outputfile
            if self.conjugate_base_job
            else None
        )
        href_file = None
        ref_file = None
        if self.has_reference_jobs and self._ref_opt_jobs_are_complete():
            href_file = (
                self.ref_acid_job.outputfile if self.ref_acid_job else None
            )
            ref_file = (
                self.ref_conjugate_base_job.outputfile
                if self.ref_conjugate_base_job
                else None
            )
        return ha_file, a_file, href_file, ref_file

    def compute_thermochemistry(self):
        """
        Compute and return thermochemistry results for all species.

        Convenience method that computes thermochemistry for all pKa species
        and returns the results dictionary.

        Returns:
            dict: Dictionary with thermochemistry data for each species.
                See Gaussian16pKaOutput.compute_pka_thermochemistry() for details.
        """
        from chemsmart.io.gaussian.output import Gaussian16pKaOutput

        if not self._opt_jobs_are_complete():
            raise ValueError(
                "Cannot compute thermochemistry: optimization jobs are not complete. "
                "Run the pKa jobs first using job.run()."
            )

        ha_file, a_file, href_file, ref_file = self._pka_output_files()

        return Gaussian16pKaOutput.compute_pka_thermochemistry(
            ha_file=ha_file,
            a_file=a_file,
            href_file=href_file,
            ref_file=ref_file,
            temperature=self.settings.temperature,
            concentration=self.settings.concentration,
            pressure=self.settings.pressure,
            cutoff_entropy_grimme=self.settings.cutoff_entropy_grimme,
            cutoff_enthalpy=self.settings.cutoff_enthalpy,
            energy_units=self.settings.energy_units,
        )

    def print_thermochemistry(self):
        """Print formatted thermochemistry summary to stdout."""
        from chemsmart.io.gaussian.output import Gaussian16pKaOutput

        if not self._opt_jobs_are_complete():
            raise ValueError(
                "Cannot print thermochemistry: optimization jobs are not complete. "
                "Run the pKa jobs first using job.run()."
            )

        ha_gas, a_gas, href_gas, ref_gas = self._pka_output_files()
        ha_solv = (
            self.protonated_sp_job.outputfile
            if self.protonated_sp_job
            else None
        )
        a_solv = (
            self.conjugate_base_sp_job.outputfile
            if self.conjugate_base_sp_job
            else None
        )
        href_solv = ref_solv = None
        if self.has_reference_jobs:
            href_solv = (
                self.ref_acid_sp_job.outputfile
                if self.ref_acid_sp_job
                else None
            )
            ref_solv = (
                self.ref_conjugate_base_sp_job.outputfile
                if self.ref_conjugate_base_sp_job
                else None
            )

        Gaussian16pKaOutput.print_pka_summary(
            ha_gas_file=ha_gas,
            a_gas_file=a_gas,
            href_gas_file=href_gas,
            ref_gas_file=ref_gas,
            ha_solv_file=ha_solv,
            a_solv_file=a_solv,
            href_solv_file=href_solv,
            ref_solv_file=ref_solv,
            pka_reference=getattr(self.settings, "reference_pka", None),
            temperature=self.settings.temperature,
            concentration=self.settings.concentration,
            pressure=self.settings.pressure,
            cutoff_entropy_grimme=self.settings.cutoff_entropy_grimme,
            cutoff_enthalpy=self.settings.cutoff_enthalpy,
            scheme=self.settings.scheme,
            delta_G_proton=getattr(self.settings, "delta_G_proton", None),
        )

    @property
    def original_mol(self):
        """
        Original molecule used to initialize the job (usually HA).
        """
        return self.molecule

    @property
    def conjugate_base_mol(self):
        """
        Conjugate base molecule (A-).
        """
        if self.conjugate_base_job:
            return self.conjugate_base_job.molecule
        # Fallback if jobs not prepared yet (unlikely given __init__)
        _, conj_mol = self.settings.conjugate_pair_molecules(self.molecule)
        return conj_mol


class GaussianpKaAnalyzeJob(GaussianpKaJob):
    """
    Gaussian job class for analyzing pKa calculation results.
    """

    TYPE = "g16pka_analyze"

    def __init__(self, input_file, **kwargs):
        """
        Initialize the analyze job.

        Args:
            input_file (Molecule): The molecule object.
            **kwargs: Additional arguments.
        """
        super().__init__(molecule=input_file, **kwargs)

    def _run(self, **kwargs):
        """Run the analysis (print thermochemistry)."""
        try:
            self.print_thermochemistry()
        except Exception as e:
            logger.error(f"Analysis failed for {self.label}: {e}")


class GaussianpKaThermoJob(GaussianpKaAnalyzeJob):
    """
    Gaussian job class for computing pKa thermochemistry (alias for analyze).
    """

    TYPE = "g16pka_thermo"
