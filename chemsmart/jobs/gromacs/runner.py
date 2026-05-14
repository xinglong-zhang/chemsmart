"""
GROMACS job runner for executing molecular dynamics workflows.
"""

import logging
import os
import shlex
import subprocess
from pathlib import Path

from chemsmart.jobs.runner import JobRunner
from chemsmart.settings.executable import GromacsExecutable

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
        self,
        server,
        scratch=None,
        fake=False,
        scratch_dir=None,
        gmx_executable=None,
        **kwargs
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
        # Default to "gmx" so the runner works when GROMACS is available in PATH.
        # Later this can be replaced by YAML/global settings logic.
        self.gmx_executable = gmx_executable

    @property
    def executable(self):
        """
        First version: directly use gmx in PATH.
        Later this can be replaced by executable settings logic.
        """
        return self._get_executable()

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
           Return the GROMACS executable.

           The executable can be provided directly through gmx_executable.
           If not provided, it falls back to the default GromacsExecutable.
           """
        if getattr(self, "gmx_executable", None) is not None:
            return self.gmx_executable

        return GromacsExecutable().get_executable()

    def _get_grompp_command(self, job):
        """
        Build the GROMACS grompp command for assembling a TPR file.
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

        return command

    def _get_mdrun_command(self, job):
        """
        Build the GROMACS mdrun command.
        """
        exe = self._get_executable()
        deffnm = Path(job.tpr_file).with_suffix("")

        return [
            exe,
            "mdrun",
            "-deffnm",
            str(deffnm),
        ]

    def _get_pdb2gmx_command(self, job):
        """
        Build the GROMACS pdb2gmx command.

        This command is used to generate topology files from an input PDB file.
        It is only needed for workflows that start from a raw structure.
        """
        exe = self._get_executable()

        input_pdb = getattr(job, "input_pdb", None)
        if input_pdb is None:
            input_pdb = getattr(job, "structure_file", None)

        output_gro = getattr(job, "output_gro", None)
        if output_gro is None:
            output_gro = getattr(job, "structure_file", None)

        command = [
            exe,
            "pdb2gmx",
            "-f",
            str(input_pdb),
            "-o",
            str(output_gro),
            "-p",
            str(job.top_file),
        ]

        force_field = getattr(job, "force_field", None)
        if force_field is not None:
            command.extend(["-ff", str(force_field)])

        water_model = getattr(job, "water_model", None)
        if water_model is not None:
            command.extend(["-water", str(water_model)])

        return command

    def _get_editconf_command(self, job):
        """
        Build the GROMACS editconf command.

        This command is commonly used to define the simulation box.
        """
        exe = self._get_executable()

        input_structure = getattr(job, "editconf_input_file", None)
        if input_structure is None:
            input_structure = getattr(job, "structure_file", None)

        output_structure = getattr(job, "editconf_output_file", None)
        if output_structure is None:
            output_structure = getattr(job, "boxed_structure_file", None)

        command = [
            exe,
            "editconf",
            "-f",
            str(input_structure),
            "-o",
            str(output_structure),
        ]

        box_type = getattr(job, "box_type", None)
        if box_type is not None:
            command.extend(["-bt", str(box_type)])

        distance = getattr(job, "box_distance", None)
        if distance is not None:
            command.extend(["-d", str(distance)])

        return command

    def _get_solvate_command(self, job):
        """
        Build the GROMACS solvate command.

        This command is used to add solvent molecules to the simulation box.
        """
        exe = self._get_executable()

        input_structure = getattr(job, "solvate_input_file", None)
        if input_structure is None:
            input_structure = getattr(job, "boxed_structure_file", None)

        output_structure = getattr(job, "solvate_output_file", None)
        if output_structure is None:
            output_structure = getattr(job, "solvated_structure_file", None)

        command = [
            exe,
            "solvate",
            "-cp",
            str(input_structure),
            "-o",
            str(output_structure),
            "-p",
            str(job.top_file),
        ]

        solvent_file = getattr(job, "solvent_file", None)
        if solvent_file is not None:
            command.extend(["-cs", str(solvent_file)])

        return command

    def _get_genion_command(self, job):
        """
        Build the GROMACS genion command.

        This command is used to add ions to the solvated system.
        It usually requires a TPR file generated by grompp.
        """
        exe = self._get_executable()

        input_tpr = getattr(job, "ions_tpr_file", None)
        if input_tpr is None:
            input_tpr = getattr(job, "tpr_file", None)

        output_structure = getattr(job, "ions_output_file", None)
        if output_structure is None:
            output_structure = getattr(job, "ionized_structure_file", None)

        command = [
            exe,
            "genion",
            "-s",
            str(input_tpr),
            "-o",
            str(output_structure),
            "-p",
            str(job.top_file),
        ]

        positive_ion = getattr(job, "positive_ion", None)
        if positive_ion is not None:
            command.extend(["-pname", str(positive_ion)])

        negative_ion = getattr(job, "negative_ion", None)
        if negative_ion is not None:
            command.extend(["-nname", str(negative_ion)])

        neutral = getattr(job, "neutral", None)
        if neutral:
            command.append("-neutral")

        return command

    def _get_commands(self, job):
        """
        Return the command sequence for a GROMACS workflow.

        Different GROMACS workflows may require different command combinations.
        """
        workflow = getattr(job, "workflow", "prepared")

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
            raise NotImplementedError(
                "The full_setup workflow is not fully implemented yet. "
                "Current stable workflow is prepared: grompp -> mdrun."
            )
        raise ValueError(f"Unsupported GROMACS workflow: {workflow}")

    def _get_command(self, job):
        """
         Return the GROMACS mdrun command.

         The TPR file should already be assembled in _prerun.
        """

        return " ".join(self._get_mdrun_command(job))

    def _run_command(self, command, cwd=None, env=None, check=True):
        """
        Run a GROMACS command and return the completed process.

        Parameters
        ----------
        command : list[str]
            Command represented as a list of arguments.
        cwd : str or Path, optional
            Working directory for the command.
        env : dict, optional
            Environment variables used to run the command.
        check : bool
            If True, raise RuntimeError when the command fails.
        """
        if cwd is None:
            cwd = self.running_directory

        logger.info(f"Running GROMACS command: {' '.join(command)}")

        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=cwd,
            env=env,
        )
        stdout, stderr = process.communicate()

        if check and process.returncode != 0:
            raise RuntimeError(
                "GROMACS command failed.\n"
                f"Command: {' '.join(command)}\n"
                f"Return code: {process.returncode}\n"
                f"STDOUT:\n{stdout}\n"
                f"STDERR:\n{stderr}"
            )

        return process

    def _run_commands(self, commands, cwd=None, env=None, check=True):
        """
        Run a sequence of GROMACS commands in order.

        This helper prepares the runner for multi-step workflows, while the
        current stable workflow still uses grompp in _prerun and mdrun as the
        main command.
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

        This method keeps compatibility with the parent runner. It accepts both
        string commands and list commands.
        """
        if isinstance(command, str):
            command = shlex.split(command)

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
        command = self._get_grompp_command(job)
        self._run_command(command, cwd=self.running_directory, check=True)

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
