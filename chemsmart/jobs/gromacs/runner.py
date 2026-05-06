"""
GROMACS job runner for executing molecular dynamics workflows.
"""

import logging
import os
import shlex
import subprocess

from chemsmart.jobs.runner import JobRunner

logger = logging.getLogger(__name__)


class GromacsJobRunner(JobRunner):
    """
    Job runner for executing GROMACS jobs.
    """

    JOBTYPES = [
        "gmx",
        "gromacs",
        "gmxem",
    ]

    PROGRAM = "gromacs"

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

    def _prerun(self, job):
        self._assign_variables(job)

    def _assign_variables(self, job):
        self.running_directory = job.folder
        self.job_inputfile = os.path.abspath(job.inputfile)
        self.job_outputfile = os.path.abspath(job.outputfile)
        self.job_errfile = os.path.abspath(job.errfile)

    def _write_input(self, job):
        """
        First version: skip actual writer implementation for now.
        Later this should call a GROMACS input writer.
        """
        pass

    def _get_executable(self):
        """
        First version: directly use gmx in PATH.
        Later this can be replaced by executable settings logic.
        """
        return "gmx"

    def _get_command(self, job):
        """
        First version: assume EM-style execution.
        Later this should branch by job type.
        """
        exe = self._get_executable()
        command = f"{exe} mdrun -deffnm {job.label}"
        return command

    def _create_process(self, job, command, env):
        with (
            open(self.job_outputfile, "w") as out,
            open(self.job_errfile, "w") as err,
        ):
            logger.info(f"Command executed: {command}")
            return subprocess.Popen(
                shlex.split(command),
                stdout=out,
                stderr=err,
                env=env,
                cwd=self.running_directory,
            )

    def _postrun(self, job):
        pass
