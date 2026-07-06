"""GROMACS job runner for executing molecular dynamics workflows."""

from __future__ import annotations

import logging
import os
import shlex
import subprocess
from functools import lru_cache
from pathlib import Path

from chemsmart.jobs.gromacs.state import GromacsWorkflowState
from chemsmart.jobs.runner import JobRunner
from chemsmart.settings.executable import GromacsExecutable

logger = logging.getLogger(__name__)


class GromacsJobRunner(JobRunner):
    """
    Job runner for executing GROMACS jobs.

    The stable prepared workflow is:

        grompp -> mdrun

    The first full setup workflow is:

        pdb2gmx -> editconf -> solvate -> grompp(ions)
        -> genion -> grompp(em) -> mdrun
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
        self,
        server,
        scratch=None,
        fake=False,
        scratch_dir=None,
        gmx_executable=None,
        gmx_modules=None,
        gmx_env=None,
        gmx_source_scripts=None,
        **kwargs,
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

        self.gmx_executable = gmx_executable
        self.gmx_modules = list(gmx_modules or [])
        self.gmx_env = dict(gmx_env or {})
        self.gmx_source_scripts = list(gmx_source_scripts or [])

    @property
    @lru_cache(maxsize=12)
    def executable(self):
        """
        Return the GROMACS executable configuration for the current server.
        """
        return GromacsExecutable.from_servername(servername=self.server.name)

    def _prerun(self, job):
        """
        Prepare inputs before the main mdrun command is launched.
        """
        self._assign_variables(job)

        # Generate missing GROMACS input files before validation.
        # For example, if job.mdp_file is None, this can write em.mdp / nvt.mdp
        # and assign the generated file back to job.mdp_file.
        self._write_input(job)

        workflow = job.workflow

        if workflow == "prepared":
            self._validate_gromacs_inputs(job)

            if not job.has_tpr():
                self._assemble_tpr(job)

            return

        if workflow == "full_setup":
            self._validate_full_setup_inputs(job)
            self._run_full_setup(job)
            return

        raise ValueError(f"Unsupported GROMACS workflow: {workflow}")

    def _assign_variables(self, job):
        """
        Assign paths used by the parent runner and subprocess layer.
        """
        Path(job.folder).mkdir(parents=True, exist_ok=True)

        self.running_directory = str(Path(job.folder).resolve())
        self.job_outputfile = os.path.abspath(
            os.path.join(job.folder, f"{job.label}.out")
        )
        self.job_errfile = os.path.abspath(
            os.path.join(job.folder, f"{job.label}.err")
        )

    def _write_input(self, job):
        """
        Write missing GROMACS input files.

        If the user already supplied an MDP file, it is respected. Otherwise,
        ChemSmart generates a default MDP file from the GROMACS job settings.
        """
        from chemsmart.jobs.gromacs.writer import GromacsInputWriter

        writer = GromacsInputWriter(job=job, jobrunner=self)
        writer.write(target_directory=job.folder)

    def _get_executable(self):
        """
        Return the GROMACS executable.

        The direct runner option has the highest priority. If it is not given,
        fall back to GromacsExecutable.
        """
        if self.gmx_executable is not None:
            return self.gmx_executable

        return self.executable.get_executable()

    def _get_grompp_command(
        self,
        job,
        structure_file=None,
        tpr_file=None,
        mdp_file=None,
    ):
        """
        Build the GROMACS grompp command for assembling a TPR file.
        """
        command = [
            self._get_executable(),
            "grompp",
            "-f",
            str(mdp_file or job.mdp_file),
            "-c",
            str(structure_file or job.structure_file),
            "-p",
            str(job.top_file),
            "-o",
            str(tpr_file or job.tpr_file),
        ]

        if job.index_file is not None:
            command.extend(["-n", str(job.index_file)])

        if job.grompp_maxwarn is not None:
            command.extend(["-maxwarn", str(job.grompp_maxwarn)])

        return command

    def _get_mdrun_command(self, job):
        """
        Build the GROMACS mdrun command.
        """
        command = [
            self._get_executable(),
            "mdrun",
            "-deffnm",
            str(Path(job.tpr_file).with_suffix("")),
        ]

        if job.mdrun_threads is not None:
            command.extend(["-nt", str(job.mdrun_threads)])

        if job.mdrun_ntmpi is not None:
            command.extend(["-ntmpi", str(job.mdrun_ntmpi)])

        if job.mdrun_ntomp is not None:
            command.extend(["-ntomp", str(job.mdrun_ntomp)])

        command.extend(job.mdrun_extra_args)

        return command

    def _get_pdb2gmx_command(self, job):
        """
        Build the GROMACS pdb2gmx command.
        """
        input_pdb = job.input_pdb
        if input_pdb is None:
            input_pdb = job.structure_file

        output_gro = job.processed_structure_file
        if output_gro is None:
            output_gro = job.structure_file

        command = [
            self._get_executable(),
            "pdb2gmx",
            "-f",
            str(input_pdb),
            "-o",
            str(output_gro),
            "-p",
            str(job.top_file),
        ]

        if job.force_field is not None:
            command.extend(["-ff", str(job.force_field)])

        if job.water_model is not None:
            command.extend(["-water", str(job.water_model)])

        return command

    def _get_editconf_command(self, job):
        """
        Build the GROMACS editconf command.
        """
        input_structure = job.processed_structure_file
        if input_structure is None:
            input_structure = job.structure_file

        output_structure = job.boxed_structure_file

        command = [
            self._get_executable(),
            "editconf",
            "-f",
            str(input_structure),
            "-o",
            str(output_structure),
        ]

        if job.box_type is not None:
            command.extend(["-bt", str(job.box_type)])

        if job.box_distance is not None:
            command.extend(["-d", str(job.box_distance)])

        return command

    def _get_solvate_command(self, job):
        """
        Build the GROMACS solvate command.
        """
        input_structure = job.boxed_structure_file
        output_structure = job.solvated_structure_file

        command = [
            self._get_executable(),
            "solvate",
            "-cp",
            str(input_structure),
            "-o",
            str(output_structure),
            "-p",
            str(job.top_file),
        ]

        if job.solvent_file is not None:
            command.extend(["-cs", str(job.solvent_file)])

        return command

    def _get_genion_command(self, job):
        """
        Build the GROMACS genion command.
        """
        input_tpr = job.ions_tpr_file
        if input_tpr is None:
            input_tpr = job.tpr_file

        output_structure = job.ionized_structure_file

        command = [
            self._get_executable(),
            "genion",
            "-s",
            str(input_tpr),
            "-o",
            str(output_structure),
            "-p",
            str(job.top_file),
        ]

        if job.index_file is not None:
            command.extend(["-n", str(job.index_file)])

        if job.positive_ion is not None:
            command.extend(["-pname", str(job.positive_ion)])

        if job.negative_ion is not None:
            command.extend(["-nname", str(job.negative_ion)])

        if job.neutral:
            command.append("-neutral")

        return command

    def _bind_workflow_state(self, job):
        """
        Bind full setup intermediate file names to the job.
        """
        state = GromacsWorkflowState.from_job(job)

        job.processed_structure_file = state.processed_structure_file
        job.boxed_structure_file = state.boxed_structure_file
        job.solvated_structure_file = state.solvated_structure_file
        job.ions_tpr_file = state.ions_tpr_file
        job.ionized_structure_file = state.ionized_structure_file
        job.tpr_file = state.em_tpr_file

        if job.top_file is None:
            job.top_file = Path(job.folder) / "topol.top"

        return state

    def _get_commands(self, job):
        """
        Return the command sequence for a GROMACS workflow.

        This is mainly useful for unit tests and dry structural validation.
        Real execution is handled through _prerun and _get_command.
        """
        workflow = job.workflow

        if workflow == "prepared":
            return [
                self._get_grompp_command(job),
                self._get_mdrun_command(job),
            ]

        if workflow == "pdb2gmx":
            return [
                self._get_pdb2gmx_command(job),
                self._get_grompp_command(job),
                self._get_mdrun_command(job),
            ]

        if workflow == "full_setup":
            self._bind_workflow_state(job)

            ions_mdp_file = job.ions_mdp_file or job.mdp_file

            return [
                self._get_pdb2gmx_command(job),
                self._get_editconf_command(job),
                self._get_solvate_command(job),
                self._get_grompp_command(
                    job,
                    structure_file=job.solvated_structure_file,
                    tpr_file=job.ions_tpr_file,
                    mdp_file=ions_mdp_file,
                ),
                self._get_genion_command(job),
                self._get_grompp_command(
                    job,
                    structure_file=job.ionized_structure_file,
                    tpr_file=job.tpr_file,
                    mdp_file=job.mdp_file,
                ),
                self._get_mdrun_command(job),
            ]

        raise ValueError(f"Unsupported GROMACS workflow: {workflow}")

    def _get_command(self, job):
        """
        Return the final GROMACS mdrun command.

        The TPR file should already be assembled in _prerun.
        """
        return " ".join(self._get_mdrun_command(job))

    def _format_command_for_log(self, command):
        """
        Format a command list for logging.
        """
        return " ".join(shlex.quote(str(part)) for part in command)

    def _prepare_environment(self, env=None):
        """
        Prepare environment variables for subprocess execution.
        """
        run_env = os.environ.copy()
        run_env.update(self.gmx_env)

        if env is not None:
            run_env.update(env)

        return run_env

    def _wrap_command_if_needed(self, command):
        """
        Wrap command with bash when modules or source scripts are required.
        """
        modules = self.gmx_modules
        source_scripts = self.gmx_source_scripts

        if not modules and not source_scripts:
            return list(command)

        shell_steps = [
            "source /etc/profile >/dev/null 2>&1 || true",
        ]

        for script in source_scripts:
            shell_steps.append(f"source {shlex.quote(str(script))}")

        for module in modules:
            shell_steps.append(f"module load {shlex.quote(str(module))}")

        shell_steps.append(self._format_command_for_log(command))

        return [
            "bash",
            "-lc",
            " && ".join(shell_steps),
        ]

    def _append_log(self, path, content):
        """
        Append content to a log file.
        """
        if path is None:
            return

        with open(path, "a", encoding="utf-8") as handle:
            handle.write(content)

            if not content.endswith("\n"):
                handle.write("\n")

    def _run_command(
        self,
        command,
        cwd=None,
        env=None,
        check=True,
        input_text=None,
        stage=None,
    ):
        """
        Run a GROMACS command and return the completed process.
        """
        if cwd is None:
            cwd = self.running_directory

        logger.info(
            "Running GROMACS command: %s",
            self._format_command_for_log(command),
        )

        if self.fake:
            return subprocess.CompletedProcess(
                args=command,
                returncode=0,
                stdout="",
                stderr="",
            )

        prepared_command = self._wrap_command_if_needed(command)
        run_env = self._prepare_environment(env)

        result = subprocess.run(
            prepared_command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=cwd,
            env=run_env,
            input=input_text,
            check=False,
        )

        stage_name = stage or (command[1] if len(command) > 1 else "command")

        self._append_log(
            self.job_outputfile,
            f"\n===== {stage_name} STDOUT =====\n{result.stdout}\n",
        )
        self._append_log(
            self.job_errfile,
            f"\n===== {stage_name} STDERR =====\n{result.stderr}\n",
        )

        if check and result.returncode != 0:
            raise RuntimeError(
                "GROMACS command failed.\n"
                f"Stage: {stage_name}\n"
                f"Command: {self._format_command_for_log(command)}\n"
                f"Return code: {result.returncode}\n"
                f"STDOUT:\n{result.stdout}\n"
                f"STDERR:\n{result.stderr}"
            )

        return result

    def _run_commands(self, commands, cwd=None, env=None, check=True):
        """
        Run a sequence of GROMACS commands in order.
        """
        processes = []

        for command in commands:
            process = self._run_command(
                command,
                cwd=cwd,
                env=env,
                check=check,
            )
            processes.append(process)

        return processes

    def _create_process(self, job, command, env):
        """
        Create the main GROMACS process.

        This method keeps compatibility with the parent runner.
        """
        if isinstance(command, str):
            command = shlex.split(command)

        prepared_command = self._wrap_command_if_needed(command)
        run_env = self._prepare_environment(env)

        with (
            open(self.job_outputfile, "a", encoding="utf-8") as out,
            open(self.job_errfile, "a", encoding="utf-8") as err,
        ):
            logger.info("Command executed: %s", prepared_command)

            return subprocess.Popen(
                prepared_command,
                stdout=out,
                stderr=err,
                env=run_env,
                cwd=self.running_directory,
                text=True,
            )

    def _postrun(self, job):
        pass

    def _assemble_tpr(self, job):
        """
        Assemble a GROMACS TPR file using gmx grompp.
        """
        command = self._get_grompp_command(job)

        self._run_command(
            command,
            cwd=self.running_directory,
            check=True,
            stage="grompp",
        )

    def _run_full_setup(self, job):
        """
        Run the first full setup workflow.

        This first implementation intentionally covers the common non-interactive
        happy path only.
        """
        self._bind_workflow_state(job)

        self._run_command(
            self._get_pdb2gmx_command(job),
            stage="pdb2gmx",
        )

        self._run_command(
            self._get_editconf_command(job),
            stage="editconf",
        )

        self._run_command(
            self._get_solvate_command(job),
            stage="solvate",
        )

        ions_mdp_file = job.ions_mdp_file or job.mdp_file

        self._run_command(
            self._get_grompp_command(
                job,
                structure_file=job.solvated_structure_file,
                tpr_file=job.ions_tpr_file,
                mdp_file=ions_mdp_file,
            ),
            stage="grompp_ions",
        )

        self._run_command(
            self._get_genion_command(job),
            stage="genion",
            input_text=f"{job.genion_group}\n",
        )

        self._run_command(
            self._get_grompp_command(
                job,
                structure_file=job.ionized_structure_file,
                tpr_file=job.tpr_file,
                mdp_file=job.mdp_file,
            ),
            stage="grompp_em",
        )

        job.structure_file = job.ionized_structure_file

    def _validate_gromacs_inputs(self, job):
        """
        Validate required GROMACS input files for prepared workflow.
        """
        required_files = {
            "mdp_file": job.mdp_file,
            "structure_file": job.structure_file,
            "top_file": job.top_file,
        }

        missing = [
            name
            for name, path in required_files.items()
            if path is None or not Path(path).exists()
        ]

        if missing:
            raise FileNotFoundError(
                "Missing required GROMACS input files: " + ", ".join(missing)
            )

    def _validate_full_setup_inputs(self, job):
        """
        Validate required inputs for full setup workflow.
        """
        input_file = job.input_pdb or job.structure_file

        required_files = {
            "input_pdb|structure_file": input_file,
            "mdp_file": job.mdp_file,
        }

        missing = [
            name
            for name, path in required_files.items()
            if path is None or not Path(path).exists()
        ]

        if job.force_field is None:
            missing.append("force_field")

        if job.water_model is None:
            missing.append("water_model")

        if missing:
            raise FileNotFoundError(
                "Missing required GROMACS full_setup inputs: "
                + ", ".join(missing)
            )
