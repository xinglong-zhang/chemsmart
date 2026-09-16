import inspect
import logging
import os
from abc import abstractmethod
from dataclasses import dataclass
from typing import Optional

from chemsmart.settings.executable import (  # noqa: F401
    Executable,
    GaussianExecutable,
    NCIPLOTExecutable,
    ORCAExecutable,
    PySCFExecutable,
    XTBExecutable,
)
from chemsmart.settings.user import CHEMSMARTUserSettings
from chemsmart.utils.mixins import RegistryMixin

user_settings = CHEMSMARTUserSettings()


logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class SubmissionReceiptV1:
    """What a scheduler answered when a script was submitted.

    The job id is the one handle a later process can use to ask the
    scheduler about the job; the command and script say exactly what was
    submitted, and the timestamp is when.
    """

    scheduler: str
    job_id: str
    submit_command: str
    submit_script: str
    submitted_at: str
    stdout: str = ""


#: How many approved calculations of one cohort may run at the same
#: time, when the operator's profile names no other number. Physical
#: concurrency is the host's and scientific width is the Agent's: a wave
#: of seven is one cohort and one wake, and this only decides how many of
#: its members Slurm runs at once (owner ruling, 2026-09-16).
#:
#: It was opt-in, so a cohort submitted with no explicit argument ran
#: unthrottled, and every test passed the throttle by hand -- proving it
#: was writable, never that it was on.
DEFAULT_MAX_CONCURRENT_TASKS = 4


class RunScript:
    """
    Script generator for computational job execution.

    Creates Python scripts that handle job execution with proper environment
    setup and command line argument passing. Manages the execution context
    for computational chemistry jobs.

    Attributes:
        filename (str): Path to the output script file.
        batch (bool): Whether this is a batch job execution.
        cli_args: Command line arguments to pass to the job.
    """

    def __init__(self, filename, cli_args, batch=False, execution_cwd=None):
        """
        Initialize the run script generator.

        Args:
            filename (str): Path where the script will be written.
            cli_args: Command line arguments for job execution.
            batch (bool): Whether this is a batch job. Defaults to False.
        """
        self.filename = filename
        self.batch = batch
        self.cli_args = cli_args
        self.execution_cwd = (
            None
            if execution_cwd is None
            else os.path.abspath(os.fspath(execution_cwd))
        )

    def write(self):
        """
        Write the run script to the specified file.

        Creates a Python script file that can be executed to run the
        computational job with the specified arguments.
        """
        with open(self.filename, "w") as f:
            self._write(f)

    def _write(self, f):
        """
        Write the script contents to the file handle.

        Generates a Python script with proper environment setup and
        job execution commands.

        Args:
            f: File handle to write the script contents to.
        """
        cwd_statement = (
            f"os.chdir({self.execution_cwd!r})"
            if self.execution_cwd
            else "pass"
        )
        contents = f"""\
        #!/usr/bin/env python
        import os
        import subprocess
        import sys
        os.environ['OMP_NUM_THREADS'] = '1'

        from chemsmart.cli.run import run

        def run_job():
            {cwd_statement}
            run({self.cli_args!r})

        if __name__ == '__main__':
            try:
                run_job()
            except subprocess.CalledProcessError as error:
                sys.exit(error.returncode)
        """

        # Needed to remove leading whitespace in the docstring
        contents = inspect.cleandoc(contents)
        logger.debug(f"{self.cli_args!r}")

        f.write(contents)


class Submitter(RegistryMixin):
    """
    Abstract base class for job submission systems.

    Provides the foundation for scheduler-specific job submitters that handle
    the creation and submission of computational chemistry jobs to various
    cluster management systems.

    Attributes:
        NAME (str): Class-level identifier for the submitter type.
        name (str): Instance identifier for
        this submitter (often same as NAME).
        job (Job): Job instance to be submitted.
        server (Server): Server configuration used for submission.
        kwargs (dict): Additional submission
        parameters passed through to subclasses.
    """

    NAME: Optional[str] = None

    #: The environment variable this scheduler sets to the running array
    #: task's index, and the index the scheduler's own ``--array``/``-J``
    #: range starts from. They are one declaration because they are one
    #: question -- "which task am I?" -- and answering it twice is what
    #: broke here: runscripts were named ``0..N-1`` by ``enumerate``, the
    #: SLURM directive declared ``1-N``, and the dispatch shell added one
    #: again, so every task ran another task's job and two tasks ran a
    #: file that did not exist. A scheduler that declares no variable
    #: declares no array support and is refused by name rather than
    #: silently writing a submit script with no array directive.
    ARRAY_TASK_ID_VARIABLE: Optional[str] = None
    ARRAY_INDEX_BASE: int = 0

    def __init__(self, name, job, server, **kwargs):
        """
        Initialize the job submitter.

        Args:
            name (str): Name identifier for this submitter instance.
            job: Job instance to be submitted.
            server: Server configuration for submission.
            **kwargs: Additional submission parameters.
        """
        self.name = name
        self.job = job
        self.server = server
        self.kwargs = kwargs

    def __str__(self):
        """
        String representation of the submitter.

        Returns:
            str: Human-readable submitter description.
        """
        return f"Submitter: {self.name}"

    def __eq__(self, other):
        """
        Check equality based on submitter name.

        Args:
            other (Submitter): Another submitter instance to compare.

        Returns:
            bool: True if submitter names are equal.
        """
        return self.name == other.name

    def __hash__(self):
        """
        Generate hash based on submitter name.

        Returns:
            int: Hash value for submitter name.
        """
        return hash(self.name)

    def __repr__(self):
        """
        Developer representation of the submitter.

        Returns:
            str: Detailed submitter representation for debugging.
        """
        return f"Submitter(name={self.name})"

    def __call__(self):
        """
        Create a submitter instance based on the name.

        Returns:
            Submitter: Configured submitter instance.

        Raises:
            ValueError: If no submitter is defined for the specified name.
        """
        submitter_cls = [
            s for s in Submitter.subclasses() if self.name == s.NAME
        ]
        if len(submitter_cls) == 0:
            raise ValueError(
                f"No submitter of defined name: {self.name}.\n"
                f"Available submitters: {Submitter.subclasses()}"
            )

        assert len(submitter_cls) == 1
        submitter_cls = submitter_cls[0]
        return submitter_cls(**self.kwargs)

    @classmethod
    def from_dict(cls, d):
        """
        Create a submitter instance from a dictionary.

        Args:
            d (dict): Dictionary containing submitter configuration.

        Returns:
            Submitter: Configured submitter instance.
        """
        return cls(**d)

    @property
    def submit_folder(self):
        """
        Get the submission folder for the job.

        Returns:
            str: Path to the job submission folder.
        """
        return self.job.folder

    @property
    def submit_script(self):
        """
        Get the submission script filename.

        Returns:
            str: Filename for the job submission script.
        """
        if self.job.label is not None:
            return f"chemsmart_sub_{self.job.label}.sh"
        return "chemsmart_sub.sh"

    @property
    def array_submit_script(self):
        """
        Get the array job submission script filename.

        Returns:
            str: Filename for the array job submission script.
        """
        if self.job.label is not None:
            return f"chemsmart_sub_array_{self.job.label}.sh"
        return "chemsmart_sub_array.sh"

    @property
    def scheduler_request(self):
        """The allocation this submission asks the scheduler for.

        Passed in by a caller that holds an approved execution envelope
        (``server.get_submitter(job, scheduler_request=...)``); otherwise
        resolved from the server profile alone, which is the profile's
        own numbers and therefore byte-identical to what this class wrote
        before the envelope had any route here at all.

        The submitter deliberately holds no envelope and no JobRunner:
        one resolved object is the single channel, so there is one place
        that answers "what is this allocation" rather than one per caller.
        """

        from chemsmart.settings.scheduler_request import (
            resolve_scheduler_request,
        )

        request = self.kwargs.get("scheduler_request")
        if request is not None:
            return request
        return resolve_scheduler_request(
            resources=None, server=self.server, sealed=False
        )

    def max_concurrent_tasks(self, requested=None):
        """How many of this cohort's members may run at the same time.

        The caller's number, else the operator's profile, else the host
        default. Never ``None``: an unthrottled cohort is a burst, and
        the default is a safety and fairness choice rather than an
        opt-in.

        Args:
            requested (int | None): An explicit cap from the caller.

        Returns:
            int: A positive concurrency bound.
        """

        for candidate in (
            requested,
            getattr(self.server, "max_concurrent_tasks", None),
        ):
            try:
                value = int(candidate)
            except (TypeError, ValueError):
                continue
            if value > 0:
                return value
        return DEFAULT_MAX_CONCURRENT_TASKS

    def array_run_script(self, task_id):
        """The runscript one array task runs, named by that task's own id.

        The filename *is* the task id, so the dispatch shell needs no
        arithmetic and there is no second place for the convention to be
        restated. ``task_id`` is the scheduler's index, not a list
        position.

        Args:
            task_id (int): Scheduler array task index.

        Returns:
            str: Filename of the runscript for that task.
        """
        return f"chemsmart_run_array_{task_id}.py"

    def array_task_ids(self, count):
        """The scheduler task ids for an array of ``count`` jobs.

        Args:
            count (int): Number of jobs in the array.

        Returns:
            list[int]: Task ids, from this scheduler's own index base.
        """
        base = int(self.ARRAY_INDEX_BASE)
        return list(range(base, base + int(count)))

    def _require_array_support(self):
        """Refuse an array on a scheduler that declares none.

        A scheduler failure must not wear a job's failure word: writing a
        submit script with no array directive would run one job N times
        under one id rather than say that this scheduler is unsupported.
        """
        if not self.ARRAY_TASK_ID_VARIABLE:
            raise ValueError(
                f"{type(self).__name__} declares no array task-id variable, "
                "so CHEMSMART cannot express an array job on this scheduler. "
                "Submit the jobs individually, or declare "
                "ARRAY_TASK_ID_VARIABLE and ARRAY_INDEX_BASE for it."
            )

    @property
    def run_script(self):
        """
        Get the run script filename.

        Returns:
            str: Filename for the job execution script.
        """
        if self.job.label is not None:
            return f"chemsmart_run_{self.job.label}.py"
        return "chemsmart_run.py"

    @property
    def executable(self):
        """
        Get the executable configuration for the job's program.

        Returns:
            Executable: Instance of the ``Executable`` subclass whose
            ``PROGRAM`` matches ``job.PROGRAM``.

        Raises:
            ValueError: If no executable is registered for the job's program.

        Resolved from the ``Executable`` registry rather than a hardcoded
        if/elif chain, so a newly registered program works under ``sub`` as
        soon as it works under ``run``. The previous chain silently limited
        submission to Gaussian, ORCA and NCIPLOT, which broke run/sub parity
        for any program added later.
        """
        program = str(self.job.PROGRAM or "").upper()
        for executable_class in Executable.subclasses():
            if str(executable_class.PROGRAM or "").upper() == program:
                return executable_class.from_servername(self.server.name)
        raise ValueError(
            f"Program {self.job.PROGRAM} not supported: no Executable "
            f"subclass registered with PROGRAM={program!r}. Registered "
            f"programs: "
            f"{sorted(str(c.PROGRAM) for c in Executable.subclasses())}"
        )

    def write(self, cli_args):
        """
        Write the submission and run scripts for the job.

        Creates both the scheduler-specific submission script and the Python
        run script that will execute the computational job.

        Args:
            cli_args: Command line arguments for job execution.
        """
        if self.job.is_complete():
            logger.warning("Submitting an already complete job.")
        os.makedirs(self.submit_folder, exist_ok=True)
        self._write_runscript(cli_args)
        self._write_submitscript()

    def write_array_job(self, jobs, num_nodes=None, cli_args=None):
        """
        Write submission scripts for an array job.

        Creates scripts for submitting multiple independent jobs as a
        scheduler array job, enabling parallel execution across nodes.

        Args:
            jobs (list): List of Job instances to run as an array.
            num_nodes (int): Number of nodes to request.
            cli_args: Command line arguments for job execution.
        """
        if not jobs:
            logger.warning("No jobs provided for array job")
            return

        self._require_array_support()

        # Store job list for array processing
        self.jobs = jobs
        self.num_nodes = num_nodes

        # Write run scripts for each job
        self._write_array_runscripts(jobs, cli_args)

        # Write array submit script
        self._write_array_submitscript(num_nodes)

    def _write_array_runscripts(self, jobs, cli_args):
        """
        Write individual run scripts for each job in the array.

        Args:
            jobs (list): List of jobs in the array.
            cli_args: Command line arguments. This may be either:
                - a single argument list shared by all jobs (backward-compatible), or
                - a sequence (e.g., list or tuple) of per-job argument lists,
                  where ``cli_args[i]`` contains the args for ``jobs[i]``.
        """
        task_ids = self.array_task_ids(len(jobs))
        for i, job in enumerate(jobs):
            # Determine CLI args for this specific job/index.
            # If cli_args looks like a per-job sequence (same length as jobs and
            # elements are themselves sequences), use cli_args[i]. Otherwise,
            # fall back to using cli_args for all jobs (backward-compatible).
            job_cli_args = cli_args
            if isinstance(cli_args, (list, tuple)):
                if len(cli_args) == len(jobs) and isinstance(
                    cli_args[i], (list, tuple)
                ):
                    job_cli_args = cli_args[i]

            # The filename is the scheduler's own task id, from this
            # submitter's declared index base, so the dispatch shell runs
            # ``chemsmart_run_array_${TASK_ID}.py`` with no arithmetic.
            task_id = task_ids[i]
            runscript_name = self.array_run_script(task_id)
            runscript = RunScript(
                os.path.join(self.submit_folder, runscript_name),
                job_cli_args,
                # The single-job writer has always passed this. The array
                # writer did not, so RunScript emitted `pass` instead of
                # os.chdir and every element inherited the submit
                # directory -- an array of N molecules in N directories
                # ran N times in the first molecule's, with N sets of
                # outputs colliding on program-default filenames.
                execution_cwd=getattr(
                    job, "submission_execution_cwd", job.folder
                ),
            )
            logger.debug(
                f"Writing array run script for task {task_id}: "
                f"{runscript_name}"
            )
            runscript.write()

    def _write_array_submitscript(self, num_nodes):
        """
        Write the array job submission script.

        Must be implemented by subclasses to provide scheduler-specific
        array job directives.

        Args:
            num_nodes (int): Number of nodes to request.
        """
        submit_script_path = os.path.join(
            self.submit_folder, self.array_submit_script
        )
        with open(submit_script_path, "w") as f:
            logger.debug(
                f"Writing array submission script: {submit_script_path}"
            )
            self._write_bash_header(f)
            self._write_array_scheduler_options(f, num_nodes)
            self._write_program_specifics(f)
            self._write_extra_commands(f)
            self._write_change_to_job_directory(f)
            self._write_array_job_command(f)

    def _write_array_scheduler_options(self, f, num_nodes):
        """
        Write scheduler options for array job submission.

        Must be implemented by subclasses.

        Args:
            f: File handle.
            num_nodes (int): Number of nodes.
        """
        # Default implementation - subclasses should override
        self._write_scheduler_options(f)

    def _write_array_job_command(self, f):
        """
        Write the command to execute array jobs.

        Args:
            f: File handle.
        """
        self._require_array_support()
        variable = self.ARRAY_TASK_ID_VARIABLE
        # The runscript's name is the task id, so this is an identity and
        # not a mapping. The arithmetic that used to live here is what
        # made every task run another task's job.
        f.write("# Array job execution\n")
        f.write(f'if [ -z "${variable}" ]; then\n')
        f.write(
            f'  echo "Error: {variable} is not set; this script runs as a '
            'scheduler array task." >&2\n'
        )
        f.write("  exit 1\n")
        f.write("fi\n")
        f.write(f"TASK_ID=${variable}\n\n")
        f.write(f"python {self.array_run_script('${TASK_ID}')}\n")

    def _write_runscript(self, cli_args):
        """
        Write the Python run script for job execution.

        Creates a Python script that handles the actual job execution
        with proper environment setup and argument passing.

        Args:
            cli_args: Command line arguments for the job.
        """
        runscript = RunScript(
            os.path.join(self.submit_folder, self.run_script),
            cli_args,
            execution_cwd=getattr(
                self.job, "submission_execution_cwd", self.job.folder
            ),
        )
        logger.debug(f"Writing run script to: {runscript.filename}")
        runscript.write()

    def _write_submitscript(self):
        """
        Write the scheduler submission script.

        Creates a shell script with scheduler directives and job execution
        commands appropriate for the target cluster management system.
        """
        submit_script_path = os.path.join(
            self.submit_folder, self.submit_script
        )
        with open(submit_script_path, "w") as f:
            logger.debug(f"Writing submission script to: {submit_script_path}")
            self._write_bash_header(f)
            self._write_scheduler_options(f)
            self._write_program_specifics(f)
            self._write_extra_commands(f)
            self._write_change_to_job_directory(f)
            self._write_job_command(f)

    @staticmethod
    def _write_bash_header(f):
        """
        Write the bash shebang header to the script.

        Args:
            f: File handle for writing the script.
        """
        f.write("#!/bin/bash\n\n")

    @abstractmethod
    def _write_scheduler_options(self, f):
        """
        Write scheduler-specific options to the submission script.

        This method must be implemented by subclasses to provide
        scheduler-specific directives and resource requests.

        Args:
            f: File handle for writing scheduler options.

        Raises:
            NotImplementedError: If not implemented by subclass.
        """
        raise NotImplementedError

    def _write_program_specifics(self, f):
        """
        Write program-specific environment setup to the script.

        Includes conda environment activation, module loading, script
        sourcing, and environment variable configuration specific to
        the computational program being used.

        Args:
            f: File handle for writing program-specific setup.
        """
        self._write_program_specific_conda_env(f)
        self._write_load_program_specific_modules(f)
        self._write_source_program_specific_script(f)
        self._write_program_specific_environment_variables(f)

    def _write_program_specific_conda_env(self, f):
        """
        Write conda environment activation commands.

        Different computational programs may require different conda
        environments for proper execution. This method writes the
        necessary activation commands.

        Args:
            f: File handle for writing conda environment setup.
        """
        if self.executable.conda_env is not None:
            logger.debug(
                f"Writing conda environment: {self.executable.conda_env}"
            )
            f.write("# conda environment\n")
            for line in self.executable.conda_env:
                f.write(line)
            f.write("\n")

    def _write_load_program_specific_modules(self, f):
        """
        Write module loading commands for program dependencies.

        Different computational programs may require loading different
        environment modules for proper execution.

        Args:
            f: File handle for writing module loading commands.
        """
        if self.executable.modules is not None:
            logger.debug(f"Writing modules: {self.executable.modules}")
            f.write("# modules\n")
            for line in self.executable.modules:
                f.write(line)
            f.write("\n")

    def _write_source_program_specific_script(self, f):
        """
        Write script sourcing commands for program setup.

        Different computational programs may require sourcing specific
        setup scripts for proper environment configuration.

        Args:
            f: File handle for writing script sourcing commands.
        """
        if self.executable.scripts is not None:
            logger.debug(f"Writing scripts: {self.executable.scripts}")
            f.write("# program specific scripts\n")
            for line in self.executable.scripts:
                f.write(line)
            f.write("\n")

    def _write_extra_scheduler_directives(self, f):
        """
        Write additional scheduler directives from server settings.

        Args:
            f: File handle for writing scheduler directives.
        """
        directives = self.server.extra_scheduler_directives
        if directives is None:
            return
        if isinstance(directives, str):
            f.write(directives)
            if directives and not directives.endswith("\n"):
                f.write("\n")
            return
        for line in directives:
            f.write(line)
            if line and not line.endswith("\n"):
                f.write("\n")

    def _write_extra_commands(self, f):
        """
        Write additional server-specific commands.

        Extra commands that may be required for the job execution.
        These commands are needed for all jobs across all programs
        and are specific to the server configuration.

        Args:
            f: File handle for writing extra commands.
        """
        if self.server.extra_commands is not None:
            for line in self.server.extra_commands:
                f.write(line)
            f.write("\n")

    def _write_program_specific_environment_variables(self, f):
        """
        Write program-specific environment variables.

        Different computational programs may require different environment
        variables for proper execution. May need to configure different
        scratch folders for different programs (e.g., Gaussian vs ORCA).

        Args:
            f: File handle for writing environment variable exports.
        """
        if self.executable.envars is not None:
            f.write("# Writing program specific environment variables\n")
            for key, value in self.executable.env.items():
                f.write(f"export {key}={value}\n")
            f.write("\n")

    @abstractmethod
    def _write_change_to_job_directory(self, f):
        """
        Write scheduler-specific directory change command.

        Each scheduler system has different environment variables
        for the job submission directory. This method must be
        implemented by subclasses.

        Args:
            f: File handle for writing directory change command.

        Raises:
            NotImplementedError: If not implemented by subclass.
        """
        raise NotImplementedError

    def _write_job_command(self, f):
        """
        Write the final job execution commands.

        Make the run script executable and replace the scheduler shell with
        that process. ``exec`` preserves the exact child exit status; the old
        background-plus-operandless-``wait`` sequence could report scheduler
        success after the ChemSmart child had failed.

        Args:
            f: File handle for writing job execution commands.
        """
        f.write(f"chmod +x ./{self.run_script}\n")
        f.write(f"exec ./{self.run_script}\n")

    @classmethod
    def from_scheduler_type(cls, scheduler_type, **kwargs):
        """
        Create a submitter instance for the specified scheduler type.

        Factory method that finds and instantiates the appropriate
        submitter subclass based on the scheduler type name.

        Args:
            scheduler_type (str): Name of the scheduler system
                (e.g., "PBS", "SLURM", "SLF", "FUGAKU").
            **kwargs: Additional arguments passed to the submitter constructor.

        Returns:
            Submitter: Instance of the appropriate submitter subclass
            (one of PBSSubmitter, SLURMSubmitter, SLFSubmitter, or
            FUGAKUSubmitter) configured with the provided kwargs
            (e.g., job and server).

        Raises:
            ValueError: If no submitter is found
            for the specified scheduler type.
        """
        submitters = cls.subclasses()
        for submitter in submitters:
            if submitter.NAME == scheduler_type:
                return submitter(**kwargs)
        raise ValueError(
            f"Could not find any submitters for scheduler type: {scheduler_type}."
        )


class PBSSubmitter(Submitter):
    """
    PBS (Portable Batch System) job submitter.

    Handles job submission to PBS/Torque cluster management systems.
    Creates PBS-specific submission scripts with appropriate resource
    requests and scheduler directives.

    Attributes:
        NAME (str): Identifier for PBS scheduler type ('PBS').
        name (str): Inherited; instance identifier (often 'PBS').
        job (Job): Job instance to be submitted.
        server (Server): Server configuration used for submission.
        kwargs (dict): Additional submission
        parameters passed to the base class.
    """

    #: No array directive is written for PBS (no _write_array_scheduler_options
    #: override), so no array task-id variable is declared: an array on this
    #: scheduler is refused by name rather than submitted as one job run N
    #: times under one id. PBS_ARRAYID and a -J range are what wiring it
    #: would need.

    NAME = "PBS"

    def __init__(self, name="PBS", job=None, server=None, **kwargs):
        """
        Initialize PBS submitter.

        Args:
            name (str): Name identifier for this submitter. Defaults to "PBS".
            job: Job instance to be submitted.
            server: Server configuration for PBS submission.
            **kwargs: Additional submission parameters.
        """
        super().__init__(name=name, job=job, server=server, **kwargs)

    def _write_scheduler_options(self, f):
        """
        Write PBS-specific scheduler directives.

        Writes PBS directives for output files, resource requests,
        queue selection, walltime, and user notification settings.

        Args:
            f: File handle for writing PBS directives.
        """
        f.write(f"#PBS -o {self.job.label}.pbsout\n")
        f.write(f"#PBS -e {self.job.label}.pbserr\n")
        request = self.scheduler_request
        if request.gpu_count > 0:
            f.write(f"#PBS -l gpus={request.gpu_count}\n")
        f.write(
            f"#PBS -l select=1:ncpus={request.cores}:"
            f"mpiprocs={request.cores}:mem={request.memory_directive}\n"
        )
        # using only one node here
        if self.server.queue_name:
            f.write(f"#PBS -q {self.server.queue_name}\n")
        if request.hours:
            f.write(f"#PBS -l walltime={request.hours}:00:00\n")
        if user_settings is not None:
            if user_settings.data.get("PROJECT"):
                f.write(f"#PBS -P {user_settings.data['PROJECT']}\n")
            if user_settings.data.get("EMAIL"):
                f.write(f"#PBS -M {user_settings.data['EMAIL']}\n")
                f.write("#PBS -m abe\n")
        self._write_extra_scheduler_directives(f)
        f.write("\n")
        f.write("\n")

    def _write_change_to_job_directory(self, f):
        """
        Write PBS-specific directory change command.

        Uses PBS_O_WORKDIR environment variable to change to the
        job submission directory.

        Args:
            f: File handle for writing directory change command.
        """
        f.write("cd $PBS_O_WORKDIR\n\n")


class SLURMSubmitter(Submitter):
    """
    SLURM (Simple Linux Utility for Resource Management) job submitter.

    Handles job submission to SLURM cluster management systems.
    Creates SLURM-specific submission scripts with appropriate resource
    requests and scheduler directives.

    Attributes:
        NAME (str): Identifier for SLURM scheduler type ('SLURM').
        name (str): Inherited; instance identifier (often 'SLURM').
        job (Job): Job instance to be submitted.
        server (Server): Server configuration used for submission.
        kwargs (dict): Additional submission
        parameters passed to the base class.
    """

    #: SLURM numbers array tasks from whatever the --array range says; the
    #: range is written from this base, so zero keeps task ids and list
    #: positions identical.
    ARRAY_TASK_ID_VARIABLE = "SLURM_ARRAY_TASK_ID"
    ARRAY_INDEX_BASE = 0

    NAME = "SLURM"

    def __init__(self, name="SLURM", job=None, server=None, **kwargs):
        """
        Initialize SLURM submitter.

        Args:
            name (str): Name identifier for
            this submitter. Defaults to "SLURM".
            job: Job instance to be submitted.
            server: Server configuration for SLURM submission.
            **kwargs: Additional submission parameters.
        """
        super().__init__(name=name, job=job, server=server, **kwargs)

    def _write_scheduler_options(self, f):
        """
        Write SLURM-specific scheduler directives.

        Writes SLURM directives for job name, output files, resource
        requests, partition selection, time limits, and user notifications.

        Args:
            f: File handle for writing SLURM directives.
        """
        f.write(f"#SBATCH --job-name={self.job.label}\n")
        f.write(f"#SBATCH --output={self.job.label}.slurmout\n")
        f.write(f"#SBATCH --error={self.job.label}.slurmerr\n")
        request = self.scheduler_request
        if request.gpu_count:
            f.write(f"#SBATCH --gres=gpu:{request.gpu_count}\n")
        f.write(
            f"#SBATCH --nodes=1 --ntasks-per-node={request.cores} --mem={request.memory_directive}\n"
        )
        if self.server.queue_name:
            f.write(f"#SBATCH --partition={self.server.queue_name}\n")
        if request.hours:
            f.write(f"#SBATCH --time={request.hours}:00:00\n")
        if user_settings is not None:
            if user_settings.data.get("PROJECT"):
                f.write(f"#SBATCH --account={user_settings.data['PROJECT']}\n")
            if user_settings.data.get("EMAIL"):
                f.write(f"#SBATCH --mail-user={user_settings.data['EMAIL']}\n")
                f.write("#SBATCH --mail-type=END,FAIL\n")
        self._write_extra_scheduler_directives(f)
        f.write("\n")
        f.write("\n")

    def _write_array_scheduler_options(self, f, num_nodes):
        """
        Write SLURM-specific array job scheduler directives.

        Each array task runs exactly one independent Gaussian job on a single
        node. Gaussian is a shared-memory program and cannot span multiple
        nodes, so ``--nodes=1`` is always correct per task. Parallelism across
        jobs is achieved by SLURM scheduling multiple ``--nodes=1`` tasks
        simultaneously on separate nodes — not by one task using multiple nodes.

        ``num_nodes`` controls the maximum number of concurrently running array
        tasks (the ``%N`` throttle on ``--array``). When *None*, all tasks may
        run at the same time.

        Args:
            f: File handle for writing SLURM directives.
            num_nodes (int): Number of nodes for the array job.
        """
        # Get number of jobs in array
        num_jobs = len(self.jobs) if hasattr(self, "jobs") else 1

        f.write(f"#SBATCH --job-name={self.job.label}_array\n")
        f.write(f"#SBATCH --output={self.job.label}_array_%a.slurmout\n")
        f.write(f"#SBATCH --error={self.job.label}_array_%a.slurmerr\n")

        # Array directive over this submitter's own declared task ids, so
        # the range and the runscript filenames cannot disagree. `%N`
        # bounds simultaneously running tasks of this one array job --
        # not nodes, not jobs, not the wave's scientific width.
        task_ids = self.array_task_ids(num_jobs)
        span = f"{task_ids[0]}-{task_ids[-1]}"
        f.write(
            f"#SBATCH --array={span}%{self.max_concurrent_tasks(num_nodes)}\n"
        )

        request = self.scheduler_request
        if request.gpu_count:
            f.write(f"#SBATCH --gres=gpu:{request.gpu_count}\n")
        # Each array task runs one job → always 1 node per task. A
        # shared-memory program cannot use MPI across nodes, and the
        # resources are per task, so a cohort's footprint is this
        # allocation times the %N throttle above.
        f.write(
            f"#SBATCH --nodes=1 --ntasks-per-node={request.cores} --mem={request.memory_directive}\n"
        )
        if self.server.queue_name:
            f.write(f"#SBATCH --partition={self.server.queue_name}\n")
        if request.hours:
            f.write(f"#SBATCH --time={request.hours}:00:00\n")
        if user_settings is not None:
            if user_settings.data.get("PROJECT"):
                f.write(f"#SBATCH --account={user_settings.data['PROJECT']}\n")
            if user_settings.data.get("EMAIL"):
                f.write(f"#SBATCH --mail-user={user_settings.data['EMAIL']}\n")
                f.write("#SBATCH --mail-type=END,FAIL\n")
        # The operator's own directives -- reservation, QoS, anything the
        # site requires -- reach the array exactly as they reach a single
        # job. The single-job writer has always called this and the array
        # writer did not, so moving the Agent onto an array would have
        # dropped the reservation that is the only route to some nodes.
        self._write_extra_scheduler_directives(f)
        f.write("\n")
        f.write("\n")

    def _write_change_to_job_directory(self, f):
        """
        Write SLURM-specific directory change command.

        Uses SLURM_SUBMIT_DIR environment variable to change to the
        job submission directory.

        Args:
            f: File handle for writing directory change command.
        """
        f.write("cd $SLURM_SUBMIT_DIR\n\n")


class SLFSubmitter(Submitter):
    """
    LSF (Load Sharing Facility) job submitter.

    Handles job submission to IBM LSF cluster management systems.
    Creates LSF-specific submission scripts with appropriate resource
    requests and scheduler directives.

    Note: The class name 'SLFSubmitter' appears to be a typo for 'LSFSubmitter'
    but is maintained for compatibility.

    Attributes:
        NAME (str): Identifier for LSF scheduler type ('SLF').
        name (str): Inherited; instance identifier (often 'SLF').
        job (Job): Job instance to be submitted.
        server (Server): Server configuration used for submission.
        kwargs (dict): Additional submission
        parameters passed to the base class.
    """

    #: No array directive is written for LSF, so none is declared; see the
    #: note on PBSSubmitter. LSB_JOBINDEX starts at one, so wiring it would
    #: set ARRAY_INDEX_BASE = 1 -- which is why the base is a per-scheduler
    #: declaration and not a constant.

    NAME = "SLF"

    def __init__(self, name="SLF", job=None, server=None, **kwargs):
        """
        Initialize LSF submitter.

        Args:
            name (str): Name identifier for this submitter. Defaults to "SLF".
            job: Job instance to be submitted.
            server: Server configuration for LSF submission.
            **kwargs: Additional submission parameters.
        """
        super().__init__(name=name, job=job, server=server, **kwargs)

    def _write_scheduler_options(self, f):
        """
        Write LSF-specific scheduler directives.

        Writes LSF directives for job name, output files, project assignment,
        node requests, GPU allocation, and walltime limits.

        Args:
            f: File handle for writing LSF directives.
        """
        f.write(f"#BSUB -J {self.job.label}\n")
        f.write(f"#BSUB -o {self.job.label}.bsubout\n")
        f.write(f"#BSUB -e {self.job.label}.bsuberr\n")
        if user_settings is not None:
            project_number = user_settings.data.get("PROJECT")
        if project_number is not None:
            f.write(f"#BSUB -P {project_number}\n")
        # One node per job, as #SBATCH --nodes=1 and #PJM -L node=1 both
        # already say: a shared-memory chemistry program cannot span
        # nodes. This read self.server.num_nodes, which Server does not
        # define, so every LSF submission raised AttributeError. Unverified
        # against a live LSF cluster; the directive is the documented one
        # and the crash was certain.
        f.write("#BSUB -nnodes 1\n")
        if self.scheduler_request.gpu_count:
            f.write(f"#BSUB -gpu num={self.scheduler_request.gpu_count}\n")
        f.write(f"#BSUB -W {self.scheduler_request.hours}\n")
        f.write("#BSUB -alloc_flags gpumps\n")
        f.write("\n")
        f.write("\n")

    def _write_change_to_job_directory(self, f):
        """
        Write LSF-specific directory change command.

        Uses LS_SUBCWD environment variable to change to the
        job submission directory.

        Args:
            f: File handle for writing directory change command.
        """
        f.write("cd $LS_SUBCWD\n\n")


class FUGAKUSubmitter(Submitter):
    """
    FUGAKU supercomputer job submitter.

    Handles job submission to the FUGAKU supercomputer system using
    the Fujitsu Job Operation and Management (PJM) scheduler.
    Creates PJM-specific submission scripts with appropriate resource
    requests and scheduler directives.

    Attributes:
        NAME (str): Identifier for FUGAKU scheduler type ('FUGAKU').
        name (str): Inherited; instance identifier (often 'FUGAKU').
        job (Job): Job instance to be submitted.
        server (Server): Server configuration used for submission.
        kwargs (dict): Additional submission
        parameters passed to the base class.
    """

    NAME = "FUGAKU"

    def __init__(self, name="FUGAKU", job=None, server=None, **kwargs):
        """
        Initialize FUGAKU submitter.

        Args:
            name (str): Name identifier for this
            submitter. Defaults to "FUGAKU".
            job: Job instance to be submitted.
            server: Server configuration for FUGAKU submission.
            **kwargs: Additional submission parameters.
        """
        super().__init__(name=name, job=job, server=server, **kwargs)

    def _write_scheduler_options(self, f):
        """
        Write FUGAKU PJM-specific scheduler directives.

        Writes PJM directives for resource group, node allocation,
        elapsed time, MPI processes, project assignment, and output files.
        Includes FUGAKU-specific optimizations like LLIO cache settings.

        Args:
            f: File handle for writing PJM directives.
        """
        # The resource group and the project are user settings, read the
        # way every other submitter reads them. This block indexed
        # ``data["RSCGRP"]`` directly (a bare KeyError on an unconfigured
        # host) and then wrote ``self.project``, an attribute this class
        # has never defined -- so a Fugaku submission raised
        # AttributeError before it could reach the scheduler at all.
        if user_settings is not None:
            resource_group = user_settings.data.get("RSCGRP")
            if resource_group:
                f.write(f"#PJM -L rscgrp={resource_group}\n")
        f.write("#PJM -L node=1\n")  # using one node here
        f.write(f"#PJM -L elapse={self.scheduler_request.hours}\n")
        f.write(f"#PJM --mpi proc={self.scheduler_request.cores}\n")
        if user_settings is not None and user_settings.data.get("PROJECT"):
            f.write(f'#PJM -g {user_settings.data["PROJECT"]}\n')
        f.write("#PJM -o pjm.%j.out\n")
        f.write("#PJM -e pjm.%j.err\n")
        f.write("#PJM -x PJM_LLIO_GFSCACHE=/vol0005:/vol0004\n")
        f.write("#PJM -S\n")
        f.write("\n")
        f.write("\n")

    def _write_change_to_job_directory(self, f):
        """
        Write FUGAKU PJM-specific directory change command.

        Uses PJM_O_WORKDIR environment variable to change to the
        job submission directory.

        Args:
            f: File handle for writing directory change command.
        """
        f.write("cd $PJM_O_WORKDIR\n\n")
