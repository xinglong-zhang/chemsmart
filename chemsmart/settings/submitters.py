import inspect
import logging
import os

from chemsmart.settings.user import ChemsmartUserSettings
user_settings = ChemsmartUserSettings()


logger = logging.getLogger(__name__)


class RunScript:
    def __init__(self, filename, cli_args, batch=False):
        self.filename = filename
        self.batch = batch
        self.cli_args = cli_args

    def write(self):
        with open(self.filename, "w") as f:
            self._write(f)

    def _write(self, f):
        contents = f"""\
        #!/usr/bin/env python
        import os
        os.environ['OMP_NUM_THREADS'] = '1'
        
        from chemsmart.cli.run import run

        def run_job():
            run({self.cli_args!r})

        if __name__ == '__main__':
            run_job()
        """

        # Needed to remove leading whitespace in the docstring
        contents = inspect.cleandoc(contents)

        f.write(contents)


class Submitter:
    def __init__(self, name, job, server, **kwargs):
        self.name = name
        self.job = job
        self.server = server
        self.kwargs = kwargs

    def __str__(self):
        return f"Submitter: {self.name}"

    def __eq__(self, other):
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)

    def __repr__(self):
        return f"Submitter(name={self.name})"

    def __call__(self):
        submitter_cls = [
            s for s in Submitter.subclasses() if self.name == s.NAME
        ]
        if len(submitter_cls) == 0:
            raise ValueError(
                f"No submitter of defined name: {self.name}.\nAvailable submitters: {Submitter.subclasses()}"
            )

        assert len(submitter_cls) == 1
        submitter_cls = submitter_cls[0]
        return submitter_cls(**self.kwargs)

    @classmethod
    def from_dict(cls, d):
        return cls(**d)

    @property
    def submitscript(self):
        """Submission script for the job."""
        return f"{self.job.PROGRAM}_submitscript_{self.job.label}.x"

    def write(self):
        """Write the submission script for the job."""
        if self.job.is_complete():
            logger.warning("Submitting an already complete job.")
        self._write_runscript()
        self._write_submitscript()

    def _write_runscript(self):
        """Write the run script for the job."""
        runscript = RunScript(
            f"{self.job.PROGRAM}_runscript_{self.job.label}.x",
            self.job.cli_args,
        )
        runscript.write()

    def _write_submitscript(self):
        """Write the submission script for the job."""
        with open(self.submitscript, "w") as f:
            self._write_bash_header(f)
            self._write_scheduler_options(f)
            self._write_program_specifics(f)
            self._write_change_to_job_directory(f)
            self._write_job_command(f)

    def _write_bash_header(self, f):
        f.write("#!/bin/bash\n\n")

    @abstractmethod
    def _write_scheduler_options(self, f):
        raise NotImplementedError

    def _write_program_specifics(self, f):
        self._write_program_specific_environment_variables(f)
        self._write_load_program_specific_modules(f)
        self._write_source_program_specific_script(f)

    def _write_program_specific_environment_variables(self, f):
        """Different programs may require different environment variables.
        May need to run different programs in different scratch folder."""
        program_specific_enviornment_vars = os.path.expanduser(
            f"~/.chemsmart/{self.job.PROGRAM}/{self.job.PROGRAM}.envars"
        )
        if os.path.exists(program_specific_enviornment_vars):
            with open(program_specific_enviornment_vars, "r") as f2:
                for line in f2:
                    f.write(line)
            f.write("\n")

    def _write_load_program_specific_modules(self, f):
        """Different programs may require loading different modules."""
        program_specific_modules = os.path.expanduser(
            f"~/.chemsmart/{self.job.PROGRAM}/{self.job.PROGRAM}.modules"
        )
        if os.path.exists(program_specific_modules):
            with open(program_specific_modules, "r") as f2:
                for line in f2:
                    f.write(line)
            f.write("\n")

    def _write_source_program_specific_script(self, f):
        program_specific_script = os.path.expanduser(
            f"~/.chemsmart/{self.job.PROGRAM}/{self.job.PROGRAM}.sh"
        )
        if os.path.exists(program_specific_script):
            f.write(f"source {program_specific_script}\n\n")

    @abstractmethod
    def _write_change_to_job_directory(self, f):
        raise NotImplementedError

    @abstractmethod
    def _write_job_command(self, f):
        f.write(f"chmod +x ./{self.job_run_script}\n")
        f.write(f"./{self.job_run_script} &\n")
        f.write("wait\n")


class PBSSubmitter(Submitter):
    def __init__(self, name="PBS", job=None, server=None, **kwargs):
        super().__init__(name=name, job=job, server=server, **kwargs)

    def _write_scheduler_options(self, f):
        f.write(f"#PBS -N {self.job.label}\n")
        f.write(f"#PBS -o {self.job.label}.pbsout\n")
        f.write(f"#PBS -e {self.job.label}.pbserr\n")
        f.write(
            f"#PBS -l select=1:ncpus={self.server.num_cores}:ngpus={self.server.num_gpus}:"
            f"mpiprocs={self.server.num_cores}:mem={self.server.mem_gb}\n"
        )
        # using only one node here
        f.write(f"#PBS -q {self.server.queue_name}\n")
        f.write(f"#PBS -l walltime={self.server.num_hours}\n")
        if user_settings.data.get("PROJECT"):
            f.write(f"#PBS -P {user_settings.data['PROJECT']}\n")

        if user_settings.data.get("EMAIL"):
            f.write(f"#PBS -M {user_settings.data['EMAIL']}\n")
            f.write("#PBS -m abe\n")

    def _write_change_to_job_directory(self, f):
        f.write(f"cd $PBS_O_WORKDIR\n\n")


class SLURMSubmitter(Submitter):
    def __init__(self, name="SLURM", job=None, server=None, **kwargs):
        super().__init__(name=name, job=job, server=server, **kwargs)

    def _write_scheduler_options(self, f):
        f.write(f"#SBATCH --job-name={self.job.label}\n")
        f.write(f"#SBATCH --output={self.job.label}.slurmout\n")
        f.write(f"#SBATCH --error={self.job.label}.slurmerr\n")
        f.write(
            f"#SBATCH --nodes=1 --ntasks-per-node={self.server.num_cores} "
            f"--ntasks-per-gpu={self.server.num_gpus} --mem={self.server.mem_gb}G\n"
        )
        f.write(f"#SBATCH --partition={self.server.queue_name}\n")
        f.write(f"#SBATCH --time={self.server.num_hours}:00:00\n")
        if user_settings.data.get("PROJECT"):
            f.write(f"#SBATCH --account={user_settings.data['PROJECT']}\n")

        if user_settings.data.get("EMAIL"):
            f.write(f"#SBATCH --mail-user={user_settings.data['EMAIL']}\n")
            f.write("#SBATCH --mail-type=abe\n")

    def _write_change_to_job_directory(self, f):
        f.write(f"cd $SLURM_SUBMIT_DIR\n\n")
