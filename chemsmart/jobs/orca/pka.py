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
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed

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
    5. (Optional) Same for reference acid Href and Ref-

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
        parallel=False,
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

        # parallel support
        import threading

        self.parallel = bool(parallel)
        if self.jobrunner and getattr(self.jobrunner, "run_in_serial", False):
            if self.parallel:
                logger.info(
                    "Parallel execution disabled due to run_in_serial=True in JobRunner"
                )
            self.parallel = False

        self._pka_lock = threading.Lock()

    # ------------------------------------------------------------------
    # Basename helpers for label derivation
    # ------------------------------------------------------------------

    @property
    def _acid_basename(self):
        """Basename for the target acid (HA) – equals ``self.label``."""
        return self.label

    @property
    def _conjugate_base_label(self):
        """Label for the conjugate base (A⁻)."""
        return f"{self._acid_basename}_cb"

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
        return self.settings.get_reference_molecule()

    @property
    def reference_conjugate_base_molecule(self):
        """Get the reference conjugate base molecule (Ref-)."""
        return self.settings.get_reference_conjugate_base_molecule()

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
            label=self._acid_basename,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        conjugate_base_job = ORCAOptJob(
            molecule=conjugate_base_mol,
            settings=conjugate_base_settings,
            label=self._conjugate_base_label,
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
            label=f"{self._acid_basename}_sp",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        conjugate_base_sp_job = ORCASinglePointJob(
            molecule=conjugate_base_mol,
            settings=conjugate_base_sp_settings,
            label=f"{self._conjugate_base_label}_sp",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        return protonated_sp_job, conjugate_base_sp_job

    def _prepare_ref_opt_jobs(self):
        """Create gas phase optimization jobs for Href and Ref-."""
        ref_acid_mol, ref_cb_mol = self.settings.reference_pair_molecules()
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
        return ref_acid_sp_job, ref_cb_sp_job

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------

    def _run_opt_jobs(self):
        run_in_serial = self.jobrunner and getattr(
            self.jobrunner, "run_in_serial", False
        )
        for job in self.opt_jobs:
            logger.info(f"Running gas phase optimization job: {job}")
            job.run()
            if run_in_serial and not job.is_complete():
                logger.info(f"Job {job} incomplete, breaking opt loop.")
                break

    def _run_ref_opt_jobs(self):
        if not self.has_reference_jobs:
            return

        run_in_serial = self.jobrunner and getattr(
            self.jobrunner, "run_in_serial", False
        )
        for job in self.ref_opt_jobs:
            logger.info(f"Running reference gas phase optimization job: {job}")
            job.run()
            if run_in_serial and not job.is_complete():
                logger.info(f"Job {job} incomplete, breaking ref opt loop.")
                break

    def _run_sp_jobs(self):
        self._sp_jobs = None  # refresh with optimised geometries
        run_in_serial = self.jobrunner and getattr(
            self.jobrunner, "run_in_serial", False
        )
        for job in self.sp_jobs:
            logger.info(f"Running solution phase SP job: {job}")
            job.run()
            if run_in_serial and not job.is_complete():
                logger.info(f"Job {job} incomplete, breaking sp loop.")
                break

    def _run_ref_sp_jobs(self):
        if not self.has_reference_jobs:
            return

        self._ref_sp_jobs = None
        run_in_serial = self.jobrunner and getattr(
            self.jobrunner, "run_in_serial", False
        )
        for job in self.ref_sp_jobs:
            logger.info(f"Running reference solution phase SP job: {job}")
            job.run()
            if run_in_serial and not job.is_complete():
                logger.info(f"Job {job} incomplete, breaking ref sp loop.")
                break

    def _make_sp_job(self, opt_job, fallback_molecule, sp_settings, sp_label):
        # create SP job using optimized geometry if available
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

    def _validate_imaginary_frequencies(self, opt_job, role):
        """Check for imaginary frequencies in a completed optimisation job.

        Returns:
            None if validation passes; an error message string otherwise.
        """
        out = opt_job._output()
        if out is None:
            return (
                f"[{role}] Optimisation produced no output – "
                f"cannot validate frequencies."
            )
        if not out.normal_termination:
            # Abnormal termination is not an imaginary-frequency issue;
            # it will be caught elsewhere.
            return None

        freqs = getattr(out, "vibrational_frequencies", None)
        if freqs is None:
            # No frequency data available (e.g. SP-only job) – skip check.
            return None

        imaginary = [f for f in freqs if f < 0.0]
        if imaginary:
            return (
                f"[{role}] Imaginary frequency check FAILED for "
                f"{opt_job.label}: found {len(imaginary)} imaginary "
                f"mode(s) {imaginary}. The optimised geometry is not a "
                f"true minimum – please re-optimise."
            )
        return None

    def _worker_opt_then_sp(
        self,
        opt_job,
        role,
        sp_index,
        fallback_molecule,
        sp_settings,
        sp_label,
    ):
        """Run opt → validate frequencies → SP for a single species.

        Returns a dict describing the outcome so the caller can collect
        successes and failures without any silent swallowing.
        """
        result = {
            "role": role,
            "label": opt_job.label,
            "success": False,
            "error": None,
        }

        # --- Optimisation ---------------------------------------------------
        logger.info(f"[parallel] Running opt job for role={role}: {opt_job}")
        opt_job.run()

        # --- Imaginary frequency validation ---------------------------------
        freq_error = self._validate_imaginary_frequencies(opt_job, role)
        if freq_error:
            result["error"] = freq_error
            logger.error(freq_error)
            # Store the SP slot as None so downstream code knows it failed.
            with self._pka_lock:
                self._store_sp_slot(role, sp_index, None)
            return result

        # --- Single-point ----------------------------------------------------
        sp_job = self._make_sp_job(
            opt_job=opt_job,
            fallback_molecule=fallback_molecule,
            sp_settings=sp_settings,
            sp_label=sp_label,
        )
        logger.info(f"[parallel] Running SP job for role={role}: {sp_job}")
        sp_job.run()

        with self._pka_lock:
            self._store_sp_slot(role, sp_index, sp_job)

        result["success"] = True
        return result

    def _store_sp_slot(self, role, sp_index, sp_job):
        """Store an SP job (or *None*) into the correct container slot.

        Must be called while holding ``self._pka_lock``.
        """
        if role in ("HA", "A"):
            if self._sp_jobs is None:
                self._sp_jobs = [None, None]
            self._sp_jobs[sp_index] = sp_job
        elif role in ("Href", "Ref"):
            if self._ref_sp_jobs is None:
                self._ref_sp_jobs = [None, None]
            self._ref_sp_jobs[sp_index] = sp_job
        else:
            if self._sp_jobs is None:
                self._sp_jobs = [None, None]
            self._sp_jobs[sp_index] = sp_job

    def _run_parallel(self):
        """Run pKa workflow with ``ThreadPoolExecutor``.

        Each worker runs opt → imaginary-frequency check → SP for one
        species.  Failures are **collected** (not raised) so that every
        worker has a chance to finish.  After all workers complete a
        summary is logged and, if any worker failed, a ``RuntimeError``
        is raised containing every failure message.
        """
        opt_jobs = self.opt_jobs
        protonated_sp_settings, conjugate_base_sp_settings = (
            self.settings.conjugate_pair_sp_job_settings(self.molecule)
        )
        with self._pka_lock:
            self._sp_jobs = [None, None]

        # Build a list of (args-tuple) for each worker submission.
        worker_args = [
            (
                opt_jobs[0],
                "HA",
                0,
                self.protonated_molecule,
                protonated_sp_settings,
                f"{self._acid_basename}_sp",
            ),
            (
                opt_jobs[1],
                "A",
                1,
                self.conjugate_base_molecule,
                conjugate_base_sp_settings,
                f"{self._conjugate_base_label}_sp",
            ),
        ]

        if self.has_reference_jobs:
            ref_opt_jobs = self.ref_opt_jobs
            ref_acid_sp_settings, ref_cb_sp_settings = (
                self.settings.reference_pair_sp_job_settings()
            )
            with self._pka_lock:
                self._ref_sp_jobs = [None, None]
            worker_args.append(
                (
                    ref_opt_jobs[0],
                    "Href",
                    0,
                    self.reference_molecule,
                    ref_acid_sp_settings,
                    f"{self._ref_basename}_sp",
                )
            )
            worker_args.append(
                (
                    ref_opt_jobs[1],
                    "Ref",
                    1,
                    self.reference_conjugate_base_molecule,
                    ref_cb_sp_settings,
                    f"{self._ref_conjugate_base_label}_sp",
                )
            )

        # ----- submit all workers -------------------------------------------
        successes = []
        failures = []

        with ThreadPoolExecutor(max_workers=len(worker_args)) as executor:
            future_to_role = {}
            for args in worker_args:
                fut = executor.submit(self._worker_opt_then_sp, *args)
                future_to_role[fut] = args[1]  # role name

            for fut in as_completed(future_to_role):
                role = future_to_role[fut]
                try:
                    result = fut.result()
                except Exception:
                    tb = traceback.format_exc()
                    msg = f"[{role}] Worker crashed with unhandled exception:\n{tb}"
                    logger.error(msg)
                    failures.append(msg)
                else:
                    if result["success"]:
                        successes.append(result["label"])
                    else:
                        failures.append(result["error"])

        # ----- Final reporting ----------------------------------------------
        total = len(worker_args)
        logger.info(
            f"Parallel pKa run complete: "
            f"{len(successes)}/{total} succeeded, "
            f"{len(failures)}/{total} failed."
        )
        for label in successes:
            logger.info(f"  ✓ {label}")
        for err in failures:
            logger.error(f"  ✗ {err}")

        if failures:
            summary = (
                f"{len(failures)} of {total} parallel pKa worker(s) failed:\n"
                + "\n".join(f"  - {e}" for e in failures)
            )
            raise RuntimeError(summary)

    def _run(self):
        if self.parallel:
            logger.info("Running ORCA pKa calculation in parallel mode")
            self._run_parallel()
            return

        run_in_serial = self.jobrunner and getattr(
            self.jobrunner, "run_in_serial", False
        )

        self._run_opt_jobs()

        if run_in_serial and not all(j.is_complete() for j in self.opt_jobs):
            logger.info("Opt jobs incomplete, halting serial execution.")
            return

        if self.has_reference_jobs:
            self._run_ref_opt_jobs()
            if run_in_serial and not all(
                j.is_complete() for j in self.ref_opt_jobs
            ):
                logger.info(
                    "Ref Opt jobs incomplete, halting serial execution."
                )
                return

        self._run_sp_jobs()

        if self.sp_jobs is None:
            # Should have been created
            return

        if run_in_serial and not all(j.is_complete() for j in self.sp_jobs):
            logger.info("SP jobs incomplete, halting serial execution.")
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
            href_gas_file=hb_gas,
            ref_gas_file=b_gas,
            ha_solv_file=ha_solv,
            a_solv_file=a_solv,
            href_solv_file=hb_solv,
            ref_solv_file=b_solv,
            pka_reference=self.settings.reference_pka,
            temperature=self.settings.temperature,
            concentration=self.settings.concentration,
            pressure=self.settings.pressure,
            cutoff_entropy_grimme=self.settings.cutoff_entropy_grimme,
            cutoff_enthalpy=self.settings.cutoff_enthalpy,
        )
