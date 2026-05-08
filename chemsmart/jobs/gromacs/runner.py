"""
GROMACS job runner for executing molecular dynamics workflows.
"""

import logging
import os
import shlex
import subprocess
from pathlib import Path

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

    @property
    def executable(self):
        """
        First version: directly use gmx in PATH.
        Later this can be replaced by executable settings logic.
        """
        return None

    def _prerun(self, job):
        self._assign_variables(job)
        self._validate_gromacs_inputs(job)

        if not job.has_tpr():
            self._assemble_tpr(job)

    def _assign_variables(self, job):
        self.running_directory = job.folder
        self.job_outputfile = os.path.abspath(
            os.path.join(job.folder, f"{job.label}.out")
        )
        self.job_errfile = os.path.abspath(
            os.path.join(job.folder, f"{job.label}.err")
        )

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
         Return the GROMACS mdrun command.

         The TPR file should already be assembled in _prerun.
         The -deffnm value should match the actual TPR file prefix.
        """
        exe = self._get_executable()
        deffnm = Path(job.tpr_file).with_suffix("")
        command = f"{exe} mdrun -deffnm {deffnm}"
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

    def _assemble_tpr(self, job):
        """
        Assemble a GROMACS TPR file using gmx grompp.
        """
        exe = self._get_executable()

        command = [
            exe,
            "grompp",
            "-f",
            str(job.mdp_file),
            "-c",
            str(job.structure_file),
            "-p",
            str(job.top_file),
            "-o",
            str(job.tpr_file),
        ]

        if job.index_file is not None:
            command.extend(["-n", str(job.index_file)])

        logger.info(f"Assembling TPR with command: {' '.join(command)}")

        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=self.running_directory,
        )
        stdout, stderr = process.communicate()

        if process.returncode != 0:
            raise RuntimeError(
                "Failed to assemble GROMACS TPR file with gmx grompp.\n"
                f"Command: {' '.join(command)}\n"
                f"STDOUT:\n{stdout}\n"
                f"STDERR:\n{stderr}"
            )

    def _validate_gromacs_inputs(self, job):
        """
        Validate required GROMACS input files.
        The first implementation supports prepared inputs:
        mdp_file, structure_file, and top_file.
        """
        required_files = {
            "mdp_file": getattr(job,"mdp_file",None),
            "structure_file": getattr(job,"structure_file",None),
            "top_file": getattr(job,"top_file",None),
        }

        missing = [
            name
            for name, path in required_files.items()
            if path is None or not Path(path).exists()
        ]

        if missing:
            raise FileNotFoundError(
                "Missing required GROMACS input files: "
                + ", ".join(missing)
            )
