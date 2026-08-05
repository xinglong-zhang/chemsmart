"""PySCF job runners.

PySCF runs as a **subprocess** even when it is importable in ChemSmart's own
environment. Four reasons, none of which is isolation for its own sake:

1. ``OMP_NUM_THREADS`` and friends must be set *before* an interpreter starts
   to bind libcint and BLAS thread pools. In-process, ChemSmart's numpy has
   already initialised OpenMP and ``lib.num_threads()`` cannot fix BLAS.
2. Peak RSS (measured: 620 MB for an H2O Hessian, GBs for real systems) is
   returned to the OS on process exit instead of fragmenting ChemSmart's heap
   for the rest of a batch.
3. A libcint/libxc segfault or an OOM kill takes down only the child.
4. GPU4PySCF needs cupy + CUDA and can never share ChemSmart's environment.
   Selecting a different interpreter per program is the only mechanism that
   permits it, and it must be the same code path as CPU.
"""

import logging
import os
import shlex
import subprocess
from contextlib import suppress
from functools import lru_cache
from glob import glob
from shutil import copy

from chemsmart.jobs.pyscf.writer import (
    PySCFScriptWriter,
    applied_pyscf_spec,
    write_pyscf_h5,
)
from chemsmart.jobs.runner import JobRunner
from chemsmart.settings.executable import PySCFExecutable

logger = logging.getLogger(__name__)


class PySCFJobRunner(JobRunner):
    """Run a PySCF job by executing its generated driver script.

    Attributes:
        PROGRAM (str): Program identifier ('pyscf').
        JOBTYPES (list): Job types handled by this runner.
        FAKE (bool): False; the fake variant is a separate subclass.
        SCRATCH (bool): False by default -- a PySCF job writes three small
            files and gains nothing from a scratch round trip, unlike ORCA's
            large temporary integral files.
    """

    PROGRAM = "pyscf"
    JOBTYPES = [
        "pyscf_sp",
        "pyscf_opt",
        "pyscf_hess",
        "pyscfjob",
    ]
    FAKE = False
    SCRATCH = False

    def __init__(
        self, server, scratch=None, fake=False, scratch_dir=None, **kwargs
    ):
        if scratch is None:
            scratch = self.SCRATCH
        super().__init__(
            server=server,
            scratch=scratch,
            scratch_dir=scratch_dir,
            fake=fake,
            **kwargs,
        )
        logger.debug(f"PySCF jobrunner num cores: {self.num_cores}")
        logger.debug(f"PySCF jobrunner num gpus: {self.num_gpus}")
        logger.debug(f"PySCF jobrunner mem gb: {self.mem_gb}")
        logger.debug(f"PySCF jobrunner scratch: {self.scratch}")

    @property
    @lru_cache(maxsize=12)
    def executable(self):
        """Return the PySCF executable configuration for this server."""
        try:
            return PySCFExecutable.from_servername(servername=self.server.name)
        except FileNotFoundError as e:
            logger.error(
                f"No server file {self.server} is found: {e}\n"
                f"Available servers are: {PySCFExecutable.available_servers}"
            )
            raise

    def _prerun(self, job):
        self._assign_variables(job)

    def _assign_variables(self, job):
        if self.scratch and self.scratch_dir:
            self._set_up_variables_in_scratch(job)
        else:
            self._set_up_variables_in_job_directory(job)

        if self.executable and self.executable.local_run is not None:
            job.local = self.executable.local_run

    def _set_up_variables_in_scratch(self, job):
        scratch_job_dir = os.path.join(self.scratch_dir, job.label)
        if not os.path.exists(scratch_job_dir):
            with suppress(FileExistsError):
                os.makedirs(scratch_job_dir)
        self.running_directory = scratch_job_dir
        self._set_file_paths(job, scratch_job_dir)

    def _set_up_variables_in_job_directory(self, job):
        self.running_directory = job.folder
        self._set_file_paths(job, job.folder)

    def _set_file_paths(self, job, directory):
        self.job_inputfile = os.path.abspath(
            os.path.join(directory, job.label + ".py")
        )
        self.job_outputfile = os.path.abspath(
            os.path.join(directory, job.label + ".out")
        )
        self.job_resultsfile = os.path.abspath(
            os.path.join(directory, job.label + ".h5")
        )
        self.job_errfile = os.path.abspath(
            os.path.join(directory, job.label + ".err")
        )
        logger.debug(f"Running directory: {self.running_directory}")

    def _write_input(self, job):
        """Generate the driver script in the running directory."""
        PySCFScriptWriter(job=job).write(
            target_directory=self.running_directory
        )

    def _get_command(self, job):
        """Return ``<interpreter> <label>.py``.

        For a library backend the executable is the interpreter that owns
        pyscf, resolved from the server YAML ``PYSCF.EXEFOLDER`` when set and
        falling back to the running interpreter.
        """
        exe = self.executable.get_executable()
        logger.info(f"PySCF interpreter: {exe}")
        return f"{exe} {self.job_inputfile}"

    def _update_os_environ(self, job):
        """Return the child environment with thread and device limits set.

        These must be in place *before* the interpreter starts: BLAS and
        OpenMP size their pools at load time, so setting them afterwards has
        no effect. This is the substantive reason PySCF runs out-of-process.
        """
        env = super()._update_os_environ(job)
        threads = str(self.num_cores or 1)
        for var in (
            "OMP_NUM_THREADS",
            "MKL_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "NUMEXPR_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS",
        ):
            env[var] = threads

        engine = getattr(job.settings, "engine", "cpu")
        if engine == "gpu":
            num_gpus = self.num_gpus or 1
            env.setdefault(
                "CUDA_VISIBLE_DEVICES",
                ",".join(str(i) for i in range(num_gpus)),
            )
        else:
            # Make a CPU run reproducible even on a GPU node: without this a
            # stray gpu4pyscf import could quietly pick up a device.
            env["CUDA_VISIBLE_DEVICES"] = ""
        return env

    def _create_process(self, job, command, env):
        with open(self.job_errfile, "w") as err:
            logger.info(
                f"Command executed: {command}\n"
                f"PySCF log: {self.job_outputfile}\n"
                f"Results: {self.job_resultsfile}\n"
                f"Errors: {self.job_errfile}"
            )
            return subprocess.Popen(
                shlex.split(command),
                stdout=err,
                stderr=subprocess.STDOUT,
                env=env,
                cwd=self.running_directory,
            )

    def _postrun(self, job):
        """Copy artifacts back from scratch, if scratch was used."""
        if not self.scratch:
            return
        for path in glob(f"{self.running_directory}/{job.label}*"):
            try:
                copy(path, job.folder)
            except Exception as e:
                logger.error(f"Failed to copy {path} to {job.folder}: {e}")


class FakePySCFJobRunner(PySCFJobRunner):
    """Fake PySCF runner for tests and dry runs.

    Writes the real driver script and a synthetic results file **without
    importing pyscf**, so ``--fake`` exercises the full
    CLI -> settings -> writer -> output -> parse path in any environment,
    including one where PySCF is not installed at all.

    This is deliberately unlike ``FakeORCA``, which fabricates log text for a
    real regex parser to re-parse. Because the structured results file is the
    only programmatic contract here, the fake can populate it directly.
    """

    FAKE = True

    def __init__(self, server, scratch=None, fake=True, **kwargs):
        super().__init__(server=server, scratch=scratch, fake=fake, **kwargs)

    def run(self, job, **kwargs):
        self._prerun(job)
        self._write_input(job)
        self._write_fake_results(job)
        self._postrun(job)
        self._postrun_cleanup(job)
        return 0

    def _write_fake_results(self, job):
        """Write a synthetic ``label.h5`` and a marker log."""
        import numpy as np

        config = PySCFScriptWriter(job=job).build_config()
        num_atoms = len(config["symbols"])

        spec = applied_pyscf_spec(config)
        spec["num_basis_functions"] = 0
        spec["num_shells"] = 0
        spec["num_electrons"] = 0
        spec["nelec"] = [0, 0]
        spec["fake"] = True

        status = {
            "stages": {
                stage: {"converged": True} for stage in config["stages"]
            },
            "normal_termination": True,
            "failure": None,
        }
        provenance = {
            "pyscf_version": "fake",
            "gpu4pyscf_version": None,
            "engine": config["engine"],
            "num_threads": config["num_threads"],
            "settings_digest": config["settings_digest"],
            "project_yaml_digest": config.get("project_yaml_digest"),
            "wall_seconds": 0.0,
            "core_seconds": 0.0,
            "fake": True,
        }

        results = {
            "energies": np.zeros(1, dtype=float),
            "positions": np.asarray(config["positions"], dtype=float),
            "atomic_numbers": np.zeros(num_atoms, dtype=int),
            "mo_energy": np.zeros(1, dtype=float),
            "mo_occ": np.zeros(1, dtype=float),
            "point_group": "C1",
        }
        if "hess" in config["stages"]:
            num_modes = max(3 * num_atoms - 6, 1)
            results["vibrational_frequencies"] = np.zeros(
                num_modes, dtype=float
            )
            results["normal_modes"] = np.zeros(
                (num_modes, num_atoms, 3), dtype=float
            )

        with open(self.job_outputfile, "w") as handle:
            handle.write(
                "PySCF version 0.0.0-fake\n"
                "PySCF path  (fake run; no calculation was performed)\n"
            )

        write_pyscf_h5(
            self.job_resultsfile,
            spec=spec,
            provenance=provenance,
            status=status,
            results=results,
        )
