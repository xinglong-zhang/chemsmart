import logging
import os
import shlex
import subprocess
import sys
from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Callable, Mapping, Optional, Sequence

from chemsmart.io.yaml import YAMLFile
from chemsmart.settings.submitters import Submitter
from chemsmart.settings.user import CHEMSMARTUserSettings
from chemsmart.utils.mixins import RegistryMixin, cached_property

user_settings = CHEMSMARTUserSettings()

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class SchedulerArrayPolicy:
    """Concurrency policy for submitting a ``BatchJob`` as a scheduler array.

    Controls the array concurrency throttle ``M`` (SLURM ``--array=1-N%M``,
    PBS ``#PBS -J 1-N%M``, LSF ``#BSUB -J name[1-N%M]``):

    - ``no_run_in_parallel`` → ``M=1``
    - else explicit CLI ``-M`` / ``--max-tasks`` when set and positive
    - else ``min(num_jobs, max_concurrent)`` from env/server max submitters
    - else all tasks at once (``M=num_jobs``)
    """

    no_run_in_parallel: bool = False
    array_concurrency: Optional[int] = None
    max_concurrent: Optional[int] = None

    @classmethod
    def from_jobrunner(cls, jobrunner: Any) -> "SchedulerArrayPolicy":
        """Build a policy from a ``JobRunner`` (CLI ``sub`` resources/flags)."""
        from chemsmart.jobs.runner import (
            get_configured_array_concurrency_limit,
        )

        return cls(
            no_run_in_parallel=bool(jobrunner.no_run_in_parallel),
            array_concurrency=jobrunner.cli_array_concurrency,
            max_concurrent=get_configured_array_concurrency_limit(jobrunner),
        )

    def array_throttle(self, num_jobs: int) -> int:
        """Return the array concurrency throttle ``M`` for *num_jobs* tasks."""
        if num_jobs <= 0:
            return 1
        if self.no_run_in_parallel:
            return 1
        if self.array_concurrency is not None and self.array_concurrency > 0:
            return min(num_jobs, int(self.array_concurrency))
        if self.max_concurrent is not None and self.max_concurrent > 0:
            return min(num_jobs, int(self.max_concurrent))
        return num_jobs


class Server(RegistryMixin):
    """
    Base class for computational server configurations.

    Represents a computational server or cluster environment with associated
    settings for job submission, resource allocation, and queue management.
    Provides methods for server comparison, serialization, and configuration
    management.

    Attributes:
        name (str): Unique identifier for the server.
        kwargs (dict): Additional server configuration parameters.
        _num_hours (int): Default number of hours for job allocation.
        _queue_name (str): Default queue name for job submission.
    """

    def __init__(self, name, **kwargs):
        """
        Initialize a server configuration.

        Args:
            name (str): Unique name identifier for the server.
            **kwargs: Additional configuration parameters including:
                NUM_HOURS (int): Default job time allocation in hours.
                QUEUE_NAME (str): Default queue for job submission.
                SCHEDULER (str): Job scheduler type (e.g., 'slurm', 'pbs').
        """
        self.name = name
        self.kwargs = kwargs
        self._num_hours = self.kwargs.get("NUM_HOURS", None)
        self._queue_name = self.kwargs.get("QUEUE_NAME", None)
        self._num_nodes = self.kwargs.get("NUM_NODES", 1)

    def __str__(self):
        """
        String representation of the server.

        Returns:
            str: Human-readable server description.
        """
        return f"Server: {self.name}"

    def __eq__(self, other):
        """
        Check equality based on server name.

        Args:
            other (Server): Another server instance to compare.

        Returns:
            bool: True if server names are equal.
        """
        return self.name == other.name

    def __hash__(self):
        """
        Generate hash based on server name.

        Returns:
            int: Hash value for server name.
        """
        return hash(self.name)

    def __repr__(self):
        """
        Developer representation of the server.

        Returns:
            str: Detailed server representation for debugging.
        """
        return f"Server(name={self.name})"

    @classmethod
    def from_dict(cls, d):
        """
        Create server instance from dictionary.

        Args:
            d (dict): Dictionary containing server configuration.

        Returns:
            Server: Configured server instance.
        """
        return cls(**d)

    @classmethod
    def from_yaml(cls, name):
        """
        Create server instance from YAML configuration file.

        Args:
            name (str): Path to YAML configuration file.

        Returns:
            Server: Server instance loaded from YAML.

        Raises:
            ValueError: If no YAML file is provided.
        """
        if not name:
            raise ValueError("No yaml file provided.")
        server_yaml = YAMLFile(filename=name)
        return cls(name, **server_yaml.yaml_contents_dict["SERVER"])

    @cached_property
    def scheduler(self):
        """
        Get the job scheduler for this server.

        Returns:
            str or None: Scheduler type (e.g.,
            'slurm', 'pbs') or None if not set.
        """
        return self.kwargs.get("SCHEDULER", None)

    @property
    def queue_name(self):
        """
        Get or set the queue name for job submission.

        Returns:
            str or None: Queue name for job submission.
        """
        return self._queue_name

    @queue_name.setter
    def queue_name(self, value):
        """
        Set the queue name for job submission.

        Args:
            value (str): Queue name to set.
        """
        self._queue_name = value

    @property
    def num_hours(self):
        """
        Get or set the number of hours for job allocation.

        Returns:
            int or None: Number of hours for job time limit.
        """
        return self._num_hours

    @num_hours.setter
    def num_hours(self, value):
        """
        Set the number of hours for job allocation.

        Args:
            value (int): Number of hours for job time limit.
        """
        self._num_hours = value

    @property
    def num_nodes(self):
        """
        Get or set the number of nodes for job allocation.

        Returns:
            int: Number of nodes for job submission.
        """
        return self._num_nodes

    @num_nodes.setter
    def num_nodes(self, value):
        """
        Set the number of nodes for job allocation.

        Args:
            value (int): Number of nodes.
        """
        if value is None:
            self._num_nodes = 1
            return
        self._num_nodes = int(value)

    @cached_property
    def mem_gb(self):
        """
        Get memory allocation in gigabytes.

        Returns:
            int: Memory allocation in GB (default: 64).
        """
        return self.kwargs.get("MEM_GB", 64)

    @cached_property
    def num_cores(self):
        """
        Get number of CPU cores allocation.

        Returns:
            int: Number of CPU cores (default: 16).
        """
        return self.kwargs.get("NUM_CORES", 16)

    @cached_property
    def num_gpus(self):
        """
        Get number of GPU allocation.

        Returns:
            int: Number of GPUs (default: 0).
        """
        return self.kwargs.get("NUM_GPUS", 0)

    @cached_property
    def num_threads(self):
        """
        Get number of threads for parallel execution.

        Returns:
            int: Number of threads (default: 16).
        """
        return self.kwargs.get("NUM_THREADS", 16)

    @cached_property
    def submit_command(self):
        """
        Get the job submission command for this server.

        Returns:
            str: Command used to submit jobs to the scheduler.
        """
        command = self.kwargs.get("SUBMIT_COMMAND")
        if command is None:
            command = self._get_submit_command()
        logger.debug(f"Submit command to submit the job: {command}")
        return command

    @cached_property
    def scratch_dir(self):
        """
        Get the scratch directory path for temporary files.

        Returns:
            str or None: Path to scratch directory or None if not configured.
        """
        return self.kwargs.get("SCRATCH_DIR", None)

    @cached_property
    def scratch(self):
        """
        Check if scratch directory is available.

        Returns:
            bool: True if scratch directory is configured.
        """
        return self.scratch_dir is not None

    @cached_property
    def use_hosts(self):
        """
        Get host specification configuration.

        Returns:
            str or None: Host specification settings or None if not configured.
        """
        return self.kwargs.get("USE_HOSTS", None)

    @cached_property
    def extra_commands(self):
        """
        Get additional commands to execute during job setup.

        Returns:
            list or None: List of extra commands or None if not configured.
        """
        return self.kwargs.get("EXTRA_COMMANDS", None)

    @cached_property
    def extra_scheduler_directives(self):
        """
        Get additional scheduler directives for submission scripts.

        Returns:
            str or list or None: Additional scheduler directives or None.
        """
        return self.kwargs.get("EXTRA_SCHEDULER_DIRECTIVES", None)

    def _get_submit_command(self):
        """
        Obtain job submission command based on scheduler type.

        Maps scheduler types to their corresponding submission commands
        for various cluster management systems.

        Returns:
            str or None: Submission command for
            the scheduler or None if unknown.
        """
        scheduler_submit_commands = {
            "SLURM": "sbatch",
            "PBS": "qsub",
            "LSF": "bsub < ",
            "SGE": "qsub",
            "HTCondor": "condor_q",
        }
        return scheduler_submit_commands.get(self.scheduler, None)

    def register(self):
        """
        Register this server in the global server registry.

        Adds the server to the registry if not already present,
        enabling server lookup and management.

        Returns:
            Server: This server instance.
        """
        # if server already in registry, pass
        if self in Server._REGISTRY:
            return self
        Server._REGISTRY.append(self)
        return self

    @classmethod
    def current(cls):
        """
        Get the current server based on detected scheduler type.

        Returns:
            Server: Server instance for the current environment.
        """
        return cls.from_scheduler_type()

    @classmethod
    @lru_cache(maxsize=12)
    def from_scheduler_type(cls):
        """
        Create a Server instance based on the detected scheduler type.

        Automatically detects the scheduler type in the current environment
        and creates an appropriate server instance. Falls back to local
        server if no scheduler is detected.

        Returns:
            Server: Server instance for the detected scheduler
            (or local fallback), typically a YamlServerSettings.

        Raises:
            ValueError: If no server class is
            defined for the detected scheduler type.
        """
        scheduler_type = cls.detect_server_scheduler()

        if scheduler_type == "Unknown Scheduler":
            logger.info("No scheduler detected. Using local server.")
            return cls.from_servername(servername="local")

        # Match scheduler type with available Server subclasses
        for server_cls in cls.subclasses():
            if getattr(server_cls, "SCHEDULER_TYPE", None) == scheduler_type:
                return server_cls.from_servername(scheduler_type)

        raise ValueError(
            f"No server class defined for scheduler type: {scheduler_type}. "
            f"Available servers: {cls.subclasses()}"
        )

    @staticmethod
    @lru_cache(maxsize=12)
    def detect_server_scheduler():
        """
        Detect the server's job scheduler system.

        Checks for environment variables and available commands to identify
        the type of job scheduler running on the current system.

        Returns:
            str: The detected scheduler type
            (e.g., SLURM, PBS, LSF, SGE, HTCondor)
                 or "Unknown Scheduler" if none detected.
        """
        schedulers = [
            {
                "name": "SLURM",
                "env_vars": ["SLURM_JOB_ID", "SLURM_CLUSTER_NAME"],
                "commands": [["squeue"]],
            },
            {
                "name": "PBS",
                "env_vars": ["PBS_JOBID", "PBS_QUEUE"],
                "commands": [["qstat"]],
            },
            {
                "name": "LSF",
                "env_vars": ["LSB_JOBID", "LSB_MCPU_HOSTS"],
                "commands": [["bjobs"]],
            },
            {
                "name": "SGE",
                "env_vars": [],
                "commands": [["qstat"], ["qstat", "-help"]],
                "check_output": lambda output: "Grid Engine" in output
                or "Sun" in output,
            },
            {
                "name": "HTCondor",
                "env_vars": [],
                "commands": [["condor_q"]],
            },
        ]

        for scheduler in schedulers:
            # Check environment variables
            if any(env in os.environ for env in scheduler.get("env_vars", [])):
                logger.info(f"Detected scheduler: {scheduler['name']}")
                return scheduler["name"]

            # Check commands
            for command in scheduler.get("commands", []):
                try:
                    result = subprocess.run(
                        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                    )
                    if "check_output" in scheduler:
                        output = result.stdout.decode()
                        if scheduler["check_output"](output):
                            logger.info(
                                f"Detected scheduler: {scheduler['name']}"
                            )
                            return scheduler["name"]
                    else:
                        logger.info(f"Detected scheduler: {scheduler['name']}")
                        return scheduler["name"]
                except FileNotFoundError:
                    pass  # Command not found, move to the next scheduler

        # Default case: unknown scheduler
        logger.info("No scheduler detected.")
        return "Unknown Scheduler"

    @classmethod
    def from_servername(cls, servername):
        """
        Obtain server instance from server name.

        Loads server configuration from YAML files based on the server name.
        Falls back to current server if no name is specified.

        Args:
            servername (str): Name of the server configuration to load.

        Returns:
            Server: Configured server instance.
        """
        if servername is None:
            # by default return current server
            logger.warning("No server specified. Using current server.")
            return cls.current()
        return cls._from_server_name(servername)

    @classmethod
    def _from_server_name(cls, server_name):
        """
        Get server settings from YAML file based on server name.

        Searches for server configuration files in the user settings directory.
        Provides detailed error messages if the server is not found.

        Args:
            server_name (str): Name of the server
            (with or without .yaml extension).

        Returns:
            Server: Configured server instance.

        Raises:
            ValueError: If no server configuration is found for the given name.
        """
        if server_name.endswith(".yaml"):
            server_name = server_name
        else:
            server_name = f"{server_name}.yaml"
        server_name_yaml_path = os.path.join(
            user_settings.user_server_dir, server_name
        )
        user_settings_manager = ServerSettingsManager(
            filename=server_name_yaml_path
        )
        server = cls._from_servers_manager(user_settings_manager)

        if server is not None:
            return server

        # could not find server settings
        templates_path = os.path.join(os.path.dirname(__file__), "templates")
        raise ValueError(
            f"No server implemented for {server_name}.\n\n"
            f"Place new server .yaml file in {user_settings.user_server_dir}.\n\n"
            f"Templates for server settings .yaml files are available at {templates_path}\n\n "
            f"Currently available servers: {user_settings.all_available_servers}"
        )

    @classmethod
    def _from_servers_manager(cls, manager):
        """
        Load server configuration using a settings manager.

        Internal method for loading server settings through a manager instance.

        Args:
            manager: Server settings manager instance.

        Returns:
            Server or None: Server if loaded;
            None only when the file is missing.

        Raises:
            ValueError: if the YAML is malformed or invalid.
        """
        try:
            return manager.create()
        except FileNotFoundError:
            return None

    def get_submitter(self, job, **kwargs):
        """
        Get a job submitter for this server.

        Creates an appropriate submitter instance based on the server's
        scheduler type and job requirements.

        Args:
            job: Job instance to be submitted.
            **kwargs: Additional submitter configuration parameters.

        Returns:
            Submitter: Configured job submitter for this server.
        """
        submitter = Submitter.from_scheduler_type(
            scheduler_type=self.scheduler, job=job, server=self, **kwargs
        )
        logger.info(f"Obtained Submitter: {submitter}")
        return submitter

    def submit(self, job, test=False, cli_args=None, **kwargs):
        """
        Submit a computational job to the server.

        Handles the complete job submission process including validation,
        script writing, and actual submission to the scheduler.

        Args:
            job (Job): Job instance to be submitted.
            test (bool): If True, only creates
            scripts without actual submission.
                Defaults to False.
            cli_args: Command line arguments for the job.
            **kwargs: Additional submission parameters.
        """
        # First check that the job to be
        # submitted is not already queued/running
        self._check_running_jobs(job)
        # Then write the submission script
        self._write_submission_script(job=job, cli_args=cli_args, **kwargs)
        # Submit the job
        if not test:
            self._submit_job(job)

    @staticmethod
    def _scheduler_labels_for_duplicate_check(job):
        """
        Return scheduler-visible labels to guard against duplicate submission.

        For ``BatchJob`` containers, only the container label is queued by
        the scheduler; child labels execute on the compute node and are not
        checked here.
        """
        from chemsmart.jobs.batch import BatchJob
        from chemsmart.jobs.gaussian import GaussianJob
        from chemsmart.jobs.orca.job import ORCAJob

        if job.label is None:
            return []

        if isinstance(job, (GaussianJob, ORCAJob, BatchJob)):
            return [job.label]

        return []

    @staticmethod
    def _check_running_jobs(job):
        """
        Check if the job is already running or queued.

        Prevents duplicate job submissions by checking the scheduler queue
        for jobs with the same label.

        Args:
            job: Job instance to check.
        """
        labels = Server._scheduler_labels_for_duplicate_check(job)
        if not labels:
            return

        from chemsmart.utils.cluster import ClusterHelper

        cluster_helper = ClusterHelper()
        running_job_ids, running_job_names = (
            cluster_helper.get_gaussian_running_jobs()
        )

        for label in labels:
            if label in running_job_names:
                logger.info(
                    f"Warning: submitting job with duplicate name: {label}"
                )
                sys.exit(f"Duplicate job NOT submitted: {label}")

    def _write_submission_script(self, job, cli_args, **kwargs):
        """
        Write the submission script for the job.

        Creates the necessary submission scripts using the appropriate
        submitter for the server's scheduler type.

        Args:
            job: Job instance to create submission script for.
            cli_args: Command line arguments for the job.
            **kwargs: Additional script writing parameters.
        """
        submitter = self.get_submitter(job, **kwargs)
        submitter.write(cli_args)

    def _submit_job(self, job):
        """
        Submit the job to the scheduler.

        Executes the submission command to queue the job in the scheduler.
        Handles both simple commands and complex shell commands with operators.

        Args:
            job: Job instance to submit.

        Returns:
            int: Exit code from the submission command.

        Raises:
            ValueError: If no submission command is defined for this server.
        """
        submitter = self.get_submitter(job)
        command = self.submit_command
        if command is None:
            raise ValueError(
                f"Cannot submit job on {self} "
                f"since no submit command is defined."
            )
        command += f" {submitter.submit_script}"
        logger.info(f"Submitting job with command: {command}")
        if "<" in command or ">" in command or "|" in command:
            # Use shell=True if the command has shell operators
            p = subprocess.Popen(command, shell=True)
        else:
            p = subprocess.Popen(shlex.split(command), cwd=job.folder)
        return p.wait()

    def submit_array_job(
        self,
        jobs,
        array_concurrency=None,
        test=False,
        cli_args=None,
        batch_label=None,
        **kwargs,
    ):
        """Submit child jobs as a scheduler array job.

        Writes per-task run scripts and an array submit script, then
        optionally queues the array. Script naming uses *batch_label*
        when provided (typically ``BatchJob.label``).

        Args:
            jobs: Child ``Job`` instances in array order.
            array_concurrency: Optional array concurrency throttle ``M``
                (SLURM ``--array=1-N%M``, PBS ``-J 1-N%M``, LSF ``-J name[1-N%M]``).
            test: If True, write scripts only (do not submit).
            cli_args: Shared or per-job CLI arguments for run scripts.
            batch_label: Label for ``chemsmart_sub_array_<label>.sh``.
            **kwargs: Extra submitter construction parameters.
        """
        if not jobs:
            logger.warning("No jobs to submit as array")
            return

        first_job = jobs[0]
        if batch_label is not None:
            # Duplicate-check the batch container label, not every child.
            from chemsmart.jobs.gaussian.batch import GaussianBatchJob
            from chemsmart.jobs.orca.batch import OrcaBatchJob

            program = first_job.PROGRAM.lower() if first_job.PROGRAM else ""
            if program == "orca":
                check_job = OrcaBatchJob(jobs=jobs, label=batch_label)
            else:
                check_job = GaussianBatchJob(jobs=jobs, label=batch_label)
            self._check_running_jobs(check_job)
        else:
            for job in jobs:
                self._check_running_jobs(job)

        submitter = self.get_submitter(first_job, **kwargs)
        submitter.write_array_job(
            jobs=jobs,
            array_concurrency=array_concurrency,
            cli_args=cli_args,
            batch_label=batch_label,
        )

        if not test:
            self._submit_array_job(first_job, submitter)

    def submit_batch(
        self,
        batch_job,
        *,
        policy: Optional[SchedulerArrayPolicy] = None,
        test: bool = False,
        cli_args: Optional[Sequence[str]] = None,
        rewrite_cli: Optional[
            Callable[
                [Sequence[str], Optional[Mapping[str, Any]]],
                list[str],
            ]
        ] = None,
        **kwargs,
    ):
        """Submit a top-level ``BatchJob`` as a scheduler array.

        Resolves per-task CLI args (optional *rewrite_cli* when children carry
        ``batch_entry``), applies *policy* array concurrency throttle, then
        delegates to ``submit_array_job``.

        Args:
            batch_job: Top-level ``BatchJob`` whose children become array tasks.
            policy: Array concurrency policy; defaults to an empty policy
                (throttle equals number of children).
            test: If True, write scripts only (do not queue).
            cli_args: Shared reconstructed ``chemsmart run`` CLI tokens.
            rewrite_cli: Optional callback to rewrite shared CLI for each
                child ``batch_entry`` (homogeneous multi-molecule or pKa).
            **kwargs: Extra submitter construction parameters.
        """
        from chemsmart.jobs.batch import (
            BatchJob,
            get_job_batch_entry,
            resolve_array_cli_args,
        )

        if not isinstance(batch_job, BatchJob):
            raise TypeError(
                f"submit_batch expects a BatchJob, got {type(batch_job)!r}"
            )
        if not batch_job.jobs:
            raise ValueError(
                f"BatchJob {batch_job} has no child jobs to submit."
            )

        if policy is None:
            policy = SchedulerArrayPolicy()

        shared_cli_args = list(cli_args) if cli_args is not None else []
        has_batch_entries = any(
            get_job_batch_entry(job) is not None for job in batch_job.jobs
        )
        if has_batch_entries and rewrite_cli is None:
            raise ValueError(
                "BatchJob children have batch_entry but no "
                "rewrite_cli callback was provided for per-task CLI args."
            )
        active_rewrite = rewrite_cli if has_batch_entries else None
        array_cli_args = resolve_array_cli_args(
            batch_job.jobs,
            shared_cli_args,
            rewrite_cli=active_rewrite,
        )

        throttle = policy.array_throttle(len(batch_job.jobs))
        logger.info(
            "Submitting BatchJob %r as array with %s task(s), "
            "concurrency throttle %%s=%s",
            batch_job.label,
            len(batch_job.jobs),
            throttle,
        )
        return self.submit_array_job(
            jobs=batch_job.jobs,
            array_concurrency=throttle,
            test=test,
            cli_args=array_cli_args,
            batch_label=batch_job.label,
            **kwargs,
        )

    def _submit_array_job(self, job, submitter):
        """Queue the array submit script for *job*'s folder."""
        command = self.submit_command
        if command is None:
            raise ValueError(
                f"Cannot submit job on {self} "
                f"since no submit command is defined."
            )
        command += f" {submitter.array_submit_script}"
        logger.info(f"Submitting array job with command: {command}")
        if "<" in command or ">" in command or "|" in command:
            p = subprocess.Popen(command, shell=True)
        else:
            p = subprocess.Popen(shlex.split(command), cwd=job.folder)
        return p.wait()


class YamlServerSettings(Server):
    """
    YAML-based server settings configuration.

    Extends the base Server class to provide YAML file-based configuration
    loading for server settings. Allows server configurations to be defined
    in YAML files and loaded dynamically.

    Attributes:
        NAME (str): Identifier for YAML-based server settings.
        name (str): YAML file path used as this settings' identifier.
        kwargs (dict): Parsed configuration
        values from the YAML under "SERVER".
        _num_hours (int or None): Default job time allocation in hours.
        _queue_name (str or None): Default submission queue name.
    """

    NAME = "yaml"

    def __init__(self, name, **kwargs):
        """
        Initialize YAML-based server settings.

        Args:
            name (str): Server identifier; typically the YAML filename when
                loaded via `from_yaml`, or a logical server name when
                instantiated directly.
            **kwargs: Server configuration parameters (usually values parsed
                from the YAML under the "SERVER" key).
        """
        super().__init__(name, **kwargs)

    @classmethod
    def from_yaml(cls, filename):
        """
        Create server settings from YAML configuration file.

        Args:
            filename (str): Path to YAML configuration file.

        Returns:
            YamlServerSettings: Server settings loaded from YAML.
        """
        yaml_file = YAMLFile(filename)
        return cls(name=filename, **yaml_file.yaml_contents_dict["SERVER"])

    def __repr__(self):
        """
        Developer representation of YAML server settings.

        Returns:
            str: Detailed representation for debugging.
        """
        return f"YamlServerSettings(name={self.name})"

    def __str__(self):
        """
        String representation of YAML server settings.

        Returns:
            str: Human-readable server description.
        """
        return f"YamlServerSettings: {self.name}"

    def __eq__(self, other):
        """
        Check equality based on server name.

        Args:
            other: Another server instance to compare.

        Returns:
            bool: True if server names are equal.
        """
        return self.name == other.name

    def __hash__(self):
        """
        Generate hash based on server name.

        Returns:
            int: Hash value for server name.
        """
        return hash(self.name)

    def __call__(self):
        """
        Make the server settings callable.

        Returns:
            YamlServerSettings: This server settings instance.
        """
        return self

    def register(self):
        """
        Register this server in the global registry.

        Returns:
            YamlServerSettings: This server settings instance.
        """
        return self


class ServerSettingsManager:
    """
    Manager for server settings from YAML configuration files.

    Provides management interface for loading server configurations from YAML
    files in a specified folder structure. Handles file validation and server
    settings creation for computational cluster environments.

    Attributes:
        filename (str): Absolute path to the YAML configuration file.
    """

    def __init__(self, filename):
        """
        Initialize the server settings manager.

        Args:
            filename (str): Path to YAML configuration file containing
                server settings.

        Raises:
            ValueError: If filename is None or not specified.
        """
        if filename is None:
            raise ValueError("filename is not specified")
        self.filename = os.path.abspath(filename)

    def create(self):
        """
        Create server settings from the specified YAML file.

        Loads and parses the YAML configuration file to create a complete
        server settings instance with all configurations.

        Returns:
            YamlServerSettings: Configured server settings loaded from YAML.

        Raises:
            FileNotFoundError: If the specified YAML file does not exist.
            ValueError: If the YAML file is malformed or invalid.
        """
        return YamlServerSettings.from_yaml(self.filename)


class SLURMServer(YamlServerSettings):
    """
    SLURM-specific server configuration.

    Specialized server class for SLURM (Simple Linux Utility for Resource
    Management) scheduler environments. Provides SLURM-specific defaults
    and configurations for computational clusters.

    Attributes:
        NAME (str): Server type identifier ('SLURM').
        SCHEDULER_TYPE (str): Scheduler system type ('SLURM').
        name (str): Inherited; effective identifier (typically 'SLURM.yaml').
        kwargs (dict): Inherited; configuration values parsed from YAML.
        _num_hours (int or None): Inherited default job time allocation.
        _queue_name (str or None): Inherited default submission queue.
    """

    NAME = "SLURM"
    SCHEDULER_TYPE = "SLURM"

    def __init__(self, **kwargs):
        """
        Initialize SLURM server configuration.

        Args:
            **kwargs: Additional SLURM-specific configuration parameters.
        """
        super().__init__(filename=f"{self.NAME}.yaml", **kwargs)


class PBSServer(YamlServerSettings):
    """
    PBS-specific server configuration.

    Specialized server class for PBS (Portable Batch System) scheduler
    environments. Provides PBS-specific defaults and configurations.

    Attributes:
        NAME (str): Server type identifier ('PBS').
        SCHEDULER_TYPE (str): Scheduler system type ('PBS').
        name (str): Inherited; effective identifier (typically 'PBS.yaml').
        kwargs (dict): Inherited; configuration values parsed from YAML.
        _num_hours (int or None): Inherited default job time allocation.
        _queue_name (str or None): Inherited default submission queue.
    """

    NAME = "PBS"
    SCHEDULER_TYPE = "PBS"

    def __init__(self, **kwargs):
        """
        Initialize PBS server configuration.

        Args:
            **kwargs: Additional PBS-specific configuration parameters.
        """
        super().__init__(self.NAME, **kwargs)


class LSFServer(YamlServerSettings):
    """
    LSF-specific server configuration.

    Specialized server class for LSF (Load Sharing Facility) scheduler
    environments. Provides LSF-specific defaults and configurations.

    Attributes:
        NAME (str): Server type identifier ('LSF').
        SCHEDULER_TYPE (str): Scheduler system type ('LSF').
        name (str): Inherited; effective identifier (typically 'LSF.yaml').
        kwargs (dict): Inherited; configuration values parsed from YAML.
        _num_hours (int or None): Inherited default job time allocation.
        _queue_name (str or None): Inherited default submission queue.
    """

    NAME = "LSF"
    SCHEDULER_TYPE = "LSF"

    def __init__(self, **kwargs):
        """
        Initialize LSF server configuration.

        Args:
            **kwargs: Additional LSF-specific configuration parameters.
        """
        super().__init__(self.NAME, **kwargs)


class SGE_Server(YamlServerSettings):
    """
    SGE-specific server configuration.

    Specialized server class for SGE (Sun Grid Engine) scheduler
    environments. Provides SGE-specific defaults and configurations.

    Attributes:
        NAME (str): Server type identifier ('SGE').
        SCHEDULER_TYPE (str): Scheduler system type ('SGE').
        name (str): Inherited; effective identifier (typically 'SGE.yaml').
        kwargs (dict): Inherited; configuration values parsed from YAML.
        _num_hours (int or None): Inherited default job time allocation.
        _queue_name (str or None): Inherited default submission queue.
    """

    NAME = "SGE"
    SCHEDULER_TYPE = "SGE"

    def __init__(self, **kwargs):
        """
        Initialize SGE server configuration.

        Args:
            **kwargs: Additional SGE-specific configuration parameters.
        """
        super().__init__(self.NAME, **kwargs)
