"""
ORCA pKa calculation job implementation.

This module provides the ORCApKaJob class for performing pKa
calculations using ORCA with a proper thermodynamic cycle:
1. Gas phase optimization + frequency for both HA and A-
2. Solution phase single point for both HA and A- at the same level of theory

Using the same level of theory ensures proper error cancellation for
solvation free energy calculations.

Execution parallelism
---------------------
Within one ``ORCApKaJob`` (a single target molecule), sub-jobs in each
workflow phase (HA opt, A- opt, HA SP, A- SP, and reference legs) always run
sequentially via ``run_phase_jobs``. The CLI ``--run-in-parallel`` flag does
not submit HA and A- optimizations concurrently inside one pKa job.

Parallelism applies only across separate pKa target jobs (for example
multi-molecule ``ORCABatchJob`` fan-out), not within the thermodynamic-cycle
sub-jobs of a single molecule.
"""

import logging
import os

from chemsmart.jobs.orca.job import ORCAJob
from chemsmart.jobs.orca.opt import ORCAOptJob
from chemsmart.jobs.orca.settings import ORCApKaJobSettings
from chemsmart.jobs.orca.singlepoint import ORCASinglePointJob
from chemsmart.jobs.runner import decide_phase_transition, run_phase_jobs

logger = logging.getLogger(__name__)


class ORCApKaJob(ORCAJob):
    """
    ORCA job class for pKa calculations using the dual-level proton exchange cycle.

    Performs pKa calculations using the following workflow:
    1. Optimize HA in gas phase (opt + freq)
    2. Optimize A- in gas phase (opt + freq)
    3. Run SP on optimized HA in solution
    4. Run SP on optimized A- in solution
    5. (Optional) Same for reference acid Href and Ref-

    Intra-molecule execution is strictly sequential: HA and A- sub-jobs within
    each phase (gas opt, solvation SP, etc.) never run in parallel. Use
    ``BatchJob`` fan-out for parallelism across multiple target molecules.

    Attributes:
        TYPE (str): Job type identifier ('orcapka').
        molecule (Molecule): Protonated molecular structure (HA).
        settings (ORCApKaJobSettings): pKa calculation configuration.
        label (str): Base job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the jobs.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "orcapka"
    _shared_reference_molecule_cache = {}

    @classmethod
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

    @classmethod
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

    # ------------------------------------------------------------------
    # Basename helpers for label derivation
    # ------------------------------------------------------------------

    @property
    def _ref_basename(self):
        """Basename for the reference acid (Href), derived from the
        reference geometry filename so it stays unique when multiple HA
        share one Href."""
        import os

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
        """Get the reference acid molecule (Href)."""
        reference_pair = self._get_cached_reference_pair(self.settings)
        if reference_pair is None:
            return self.settings.get_reference_molecule()
        return reference_pair[0]

    @property
    def reference_conjugate_base_molecule(self):
        """Get the reference conjugate base molecule (Ref-)."""
        reference_pair = self._get_cached_reference_pair(self.settings)
        if reference_pair is None:
            return self.settings.get_reference_conjugate_base_molecule()
        return reference_pair[1]

    @property
    def ref_opt_jobs(self):
        """Get gas phase optimization jobs for Href and Ref-."""
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
        """Get solution phase SP jobs for Href and Ref-."""
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

    def _subjob_output_paths(self, job, legacy_label=None):
        """Candidate ORCA output files for a pKa sub-job."""
        paths = []
        runner = job.jobrunner
        if runner is not None:
            runner_out = getattr(runner, "job_outputfile", None)
            if runner_out:
                paths.append(runner_out)
        paths.append(job.outputfile)
        if legacy_label is not None:
            paths.append(os.path.join(self.folder, f"{legacy_label}.out"))
        seen = set()
        ordered = []
        for path in paths:
            if path and path not in seen:
                seen.add(path)
                ordered.append(path)
        return ordered

    def _subjob_is_complete(self, job, legacy_label=None):
        from chemsmart.io.orca.output import ORCAOutput

        for path in self._subjob_output_paths(job, legacy_label):
            if not path or not os.path.exists(path):
                continue
            try:
                if ORCAOutput(path).normal_termination:
                    return True
            except Exception:
                continue
        return False

    def _subjob_output(self, job, legacy_label=None):
        from chemsmart.io.orca.output import ORCAOutput

        for path in self._subjob_output_paths(job, legacy_label):
            if not os.path.exists(path):
                continue
            try:
                output = ORCAOutput(path)
            except Exception:
                continue
            if output.normal_termination:
                return output
        return None

    def _bind_subjob(self, job, legacy_label=None):
        """Keep sub-jobs in the parent folder and resolve scratch/legacy outputs."""
        job.folder = self.folder

        def is_complete():
            return self._subjob_is_complete(job, legacy_label)

        job.is_complete = is_complete

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
            label=f"{self.label}_HA_opt",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        conjugate_base_job = ORCAOptJob(
            molecule=conjugate_base_mol,
            settings=conjugate_base_settings,
            label=f"{self.label}_A_opt",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        self._bind_subjob(protonated_job, legacy_label=self.label)
        self._bind_subjob(conjugate_base_job, legacy_label=f"{self.label}_cb")
        return protonated_job, conjugate_base_job

    def _prepare_sp_jobs(self):
        """Create solution phase SP jobs for HA and A-."""
        protonated_sp_settings, conjugate_base_sp_settings = (
            self.settings._create_solution_phase_sp_settings(self.molecule)
        )

        # Use optimised geometry if available
        prot_out = self._subjob_output(
            self.protonated_job, legacy_label=self.label
        )
        if prot_out is not None:
            protonated_mol = prot_out.molecule
        else:
            protonated_mol = self.protonated_molecule

        cb_out = self._subjob_output(
            self.conjugate_base_job, legacy_label=f"{self.label}_cb"
        )
        if cb_out is not None:
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
        self._bind_subjob(protonated_sp_job, legacy_label=f"{self.label}_sp")
        self._bind_subjob(
            conjugate_base_sp_job, legacy_label=f"{self.label}_cb_sp"
        )
        return protonated_sp_job, conjugate_base_sp_job

    def _prepare_ref_opt_jobs(self):
        """Create gas phase optimization jobs for Href and Ref-."""
        reference_pair = self._get_cached_reference_pair(self.settings)
        if reference_pair is None:
            ref_acid_mol, ref_cb_mol = self.settings.reference_pair_molecules()
        else:
            ref_acid_mol, ref_cb_mol = reference_pair
        ref_acid_settings, ref_cb_settings = (
            self.settings.reference_pair_job_settings()
        )

        ref_acid_job = ORCAOptJob(
            molecule=ref_acid_mol,
            settings=ref_acid_settings,
            label=self._ref_basename,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        ref_cb_job = ORCAOptJob(
            molecule=ref_cb_mol,
            settings=ref_cb_settings,
            label=self._ref_conjugate_base_label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        if self._ref_basename is not None:
            self._bind_subjob(ref_acid_job, legacy_label=self._ref_basename)
            self._bind_subjob(
                ref_cb_job, legacy_label=self._ref_conjugate_base_label
            )
        else:
            self._bind_subjob(ref_acid_job)
            self._bind_subjob(ref_cb_job)
        return ref_acid_job, ref_cb_job

    def _prepare_ref_sp_jobs(self):
        """Create solution phase SP jobs for Href and Ref-."""
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
            label=f"{self._ref_basename}_sp",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        ref_cb_sp_job = ORCASinglePointJob(
            molecule=ref_cb_mol,
            settings=ref_cb_sp_settings,
            label=f"{self._ref_conjugate_base_label}_sp",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        if self._ref_basename is not None:
            self._bind_subjob(
                ref_acid_sp_job, legacy_label=f"{self._ref_basename}_sp"
            )
            self._bind_subjob(
                ref_cb_sp_job,
                legacy_label=f"{self._ref_conjugate_base_label}_sp",
            )
        else:
            self._bind_subjob(ref_acid_sp_job)
            self._bind_subjob(ref_cb_sp_job)
        return ref_acid_sp_job, ref_cb_sp_job

    # ------------------------------------------------------------------
    # Execution
    #
    # Intra-molecule phases always run sequentially via run_phase_jobs.
    # --run-in-parallel only affects separate pKa BatchJob targets, not
    # concurrent HA/A submission within one thermodynamic cycle.
    # ------------------------------------------------------------------

    def _run_opt_jobs(self):
        """Run gas phase optimization jobs (HA then A- sequentially)."""
        run_phase_jobs(
            parent_runner=self.jobrunner,
            jobs=self.opt_jobs,
            stop_on_incomplete=False,
            logger_obj=logger,
            phase_label="opt",
        )

    def _run_ref_opt_jobs(self):
        if not self.has_reference_jobs:
            return
        run_phase_jobs(
            parent_runner=self.jobrunner,
            jobs=self.ref_opt_jobs,
            stop_on_incomplete=False,
            logger_obj=logger,
            phase_label="ref opt",
        )

    def _run_sp_jobs(self):
        run_phase_jobs(
            parent_runner=self.jobrunner,
            jobs=None,
            jobs_factory=lambda: self.sp_jobs,
            stop_on_incomplete=True,
            before_run=lambda: setattr(self, "_sp_jobs", None),
            logger_obj=logger,
            phase_label="sp",
        )

    def _run_ref_sp_jobs(self):
        if not self.has_reference_jobs:
            return
        run_phase_jobs(
            parent_runner=self.jobrunner,
            jobs=None,
            jobs_factory=lambda: self.ref_sp_jobs,
            stop_on_incomplete=True,
            before_run=lambda: setattr(self, "_ref_sp_jobs", None),
            logger_obj=logger,
            phase_label="ref sp",
        )

    def _make_sp_job(self, opt_job, fallback_molecule, sp_settings, sp_label):
        """Create SP job using optimized geometry if available."""
        out = opt_job._output()
        if out is not None and out.normal_termination is True:
            mol = out.molecule
        else:
            mol = fallback_molecule

        sp_job = ORCASinglePointJob(
            molecule=mol,
            settings=sp_settings,
            label=sp_label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        return sp_job

    def _run(self):
        self._run_opt_jobs()

        opt_transition = decide_phase_transition(
            phase_name="Opt",
            require_complete=True,
            is_complete=all(j.is_complete() for j in self.opt_jobs),
            stop_message="Opt jobs incomplete, halting serial execution.",
        )
        if not opt_transition.proceed:
            logger.info(opt_transition.message)
            return

        if self.has_reference_jobs:
            self._run_ref_opt_jobs()
            ref_opt_transition = decide_phase_transition(
                phase_name="Ref Opt",
                require_complete=True,
                is_complete=all(j.is_complete() for j in self.ref_opt_jobs),
                stop_message="Ref Opt jobs incomplete, halting serial execution.",
            )
            if not ref_opt_transition.proceed:
                logger.info(ref_opt_transition.message)
                return

        self._run_sp_jobs()

        if self.sp_jobs is None:
            # Should have been created
            return

        sp_transition = decide_phase_transition(
            phase_name="SP",
            require_complete=True,
            is_complete=all(j.is_complete() for j in self.sp_jobs),
            stop_message="SP jobs incomplete, halting serial execution.",
        )
        if not sp_transition.proceed:
            logger.info(sp_transition.message)
            return

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

    def _opt_jobs_are_complete(self):
        """Return True when both target gas-phase optimization jobs finished."""
        if not self.opt_jobs:
            return False
        return all(job.is_complete() for job in self.opt_jobs)

    def _ref_opt_jobs_are_complete(self):
        """Return True when reference opt jobs finished or are not configured."""
        if not self.has_reference_jobs:
            return True
        if not self.ref_opt_jobs:
            return False
        return all(job.is_complete() for job in self.ref_opt_jobs)

    def _pka_output_files(self):
        """Return (ha_file, a_file, href_file, ref_file) output-path tuple."""
        ha_file = self.protonated_job.outputfile
        a_file = self.conjugate_base_job.outputfile
        href_file = None
        ref_file = None
        if self.has_reference_jobs and self._ref_opt_jobs_are_complete():
            href_file = self.ref_acid_job.outputfile
            ref_file = self.ref_conjugate_base_job.outputfile
        return ha_file, a_file, href_file, ref_file

    def compute_thermochemistry(self):
        """Compute and return thermochemistry results for all species."""
        from chemsmart.io.orca.output import ORCApKaOutput

        if not self._opt_jobs_are_complete():
            raise ValueError(
                "Cannot compute thermochemistry: optimization jobs are not complete. "
                "Run the pKa jobs first using job.run()."
            )

        ha_file, a_file, href_file, ref_file = self._pka_output_files()

        return ORCApKaOutput.compute_pka_thermochemistry(
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
        from chemsmart.io.orca.output import ORCApKaOutput

        if not self._opt_jobs_are_complete():
            raise ValueError(
                "Cannot print thermochemistry: optimization jobs are not complete. "
                "Run the pKa jobs first using job.run()."
            )

        ha_gas, a_gas, href_gas, ref_gas = self._pka_output_files()
        ha_solv = self.protonated_sp_job.outputfile
        a_solv = self.conjugate_base_sp_job.outputfile
        href_solv = ref_solv = None
        if self.has_reference_jobs:
            href_solv = self.ref_acid_sp_job.outputfile
            ref_solv = self.ref_conjugate_base_sp_job.outputfile

        ORCApKaOutput.print_pka_summary(
            ha_gas_file=ha_gas,
            a_gas_file=a_gas,
            href_gas_file=href_gas,
            ref_gas_file=ref_gas,
            ha_solv_file=ha_solv,
            a_solv_file=a_solv,
            href_solv_file=href_solv,
            ref_solv_file=ref_solv,
            pka_reference=self.settings.reference_pka,
            temperature=self.settings.temperature,
            concentration=self.settings.concentration,
            pressure=self.settings.pressure,
            cutoff_entropy_grimme=self.settings.cutoff_entropy_grimme,
            cutoff_enthalpy=self.settings.cutoff_enthalpy,
            scheme=self.settings.scheme,
            delta_G_proton=getattr(self.settings, "delta_G_proton", None),
        )
