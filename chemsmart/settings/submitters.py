import inspect
import logging
from abc import abstractmethod

from chemsmart.settings.executable import GaussianExecutable, ORCAExecutable
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
    def submit_folder(self):
        return self.job.folder

    @property
    def submit_script(self):
        """Submission script for the job."""
        if self.job.label is not None:
            return f"chemsmart_sub_{self.job.label}.x"
        return "chemsmart_sub.x"

    @property
    def run_script(self):
        """Run script for the job."""
        if self.job.label is not None:
            return f"chemsmart_run_{self.job.label}.py"
        return "chemsmart_run.py"

    @property
    def executable(self):
        if self.job.PROGRAM.lower() == "gaussian":
            executable = GaussianExecutable.from_servername(self.server.name)
        elif self.job.PROGRAM.lower() == "orca":
            executable = ORCAExecutable.from_servername(self.server.name)
        else:
            # Need to add programs here to be supported for other types of programs
            raise ValueError(f"Program {self.job.PROGRAM} not supported.")
        return executable

    def write(self, cli_args):
        """Write the submission script for the job."""
        if self.job.is_complete():
            logger.warning("Submitting an already complete job.")
        self._write_runscript(cli_args)
        self._write_submitscript()

    def _write_runscript(self, cli_args):
        """Write the run script for the job."""
        runscript = RunScript(self.run_script, cli_args)
        logger.debug(f"Writing run script to: {runscript.filename}")
        runscript.write()

    def _write_submitscript(self):
        """Write the submission script for the job."""
        with open(self.submit_script, "w") as f:
            logger.debug(f"Writing submission script to: {self.submit_script}")
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
        self._write_program_specific_conda_env(f)
        self._write_load_program_specific_modules(f)
        self._write_source_program_specific_script(f)
        self._write_program_specific_environment_variables(f)

    def _write_program_specific_conda_env(self, f):
        if self.executable.conda_env is not None:
            logger.debug(
                f"Writing conda environment: {self.executable.conda_env}"
            )
            for line in self.executable.conda_env.split("\n"):
                logger.debug(f"Writing line: {line}")
                f.write(line)
            f.write("\n")

    def _write_load_program_specific_modules(self, f):
        """Different programs may require loading different modules."""
        if self.executable.modules is not None:
            logger.debug(f"Writing modules: {self.executable.modules}")
            for line in self.executable.modules.split("\n"):
                logger.debug(f"Writing line: {line}")
                f.write(line)
            f.write("\n")

    def _write_source_program_specific_script(self, f):
        """Different programs may require sourcing different scripts."""
        if self.executable.scripts is not None:
            logger.debug(f"Writing scripts: {self.executable.scripts}")
            for line in self.executable.scripts.split("\n"):
                logger.debug(f"Writing line: {line}")
                f.write(line)
            f.write("\n")

    def _write_program_specific_environment_variables(self, f):
        """Different programs may require different environment variables.
        May need to run different programs in different scratch folder e.g. Gaussian vs ORCA.
        """
        if self.executable.envars is not None:
            logger.debug(
                f"Writing environment variables: {self.executable.envars}"
            )
            for key, value in self.executable.env.items():
                f.write(f"export {key}={value}\n")
            f.write("\n")

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
        f.write("cd $PBS_O_WORKDIR\n\n")


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
        f.write("cd $SLURM_SUBMIT_DIR\n\n")


class SLFSubmitter(Submitter):
    def __init__(self, name="SLF", job=None, server=None, **kwargs):
        super().__init__(name=name, job=job, server=server, **kwargs)

    def _write_scheduler_options(self, f):
        f.write(f"#BSUB -J {self.job.label}\n")
        f.write(f"#BSUB -o {self.job.label}.bsubout\n")
        f.write(f"#BSUB -e {self.job.label}.bsuberr\n")
        f.write(
            f"#BSUB -nnodes {self.server.num_nodes} -P {user_settings.data['PROJECT']}\n"
        )
        f.write(f"#BSUB -W {self.server.num_hours}\n")
        f.write("#BSUB -alloc_flags gpumps\n")

    def _write_change_to_job_directory(self, f):
        f.write("cd $LS_SUBCWD\n\n")


class FUGAKUSubmitter(Submitter):
    def __init__(self, name="FUGAKU", job=None, server=None, **kwargs):
        super().__init__(name=name, job=job, server=server, **kwargs)

    def _write_scheduler_options(self, f):
        f.write(f'#PJM -L rscgrp={user_settings.data["RSCGRP"]}\n')
        f.write("#PJM -L node=1\n")  # using one node here
        f.write(f"#PJM -L elapse={self.server.num_hours}\n")
        f.write(f"#PJM --mpi proc={self.server.num_cores}\n")
        f.write(f"#PJM -g {self.project}\n")
        f.write("#PJM -o pjm.%j.out\n")
        f.write("#PJM -e pjm.%j.err\n")
        f.write("#PJM -x PJM_LLIO_GFSCACHE=/vol0005:/vol0004\n")
        f.write("#PJM -S\n")

    def _write_change_to_job_directory(self, f):
        f.write("cd $PJM_O_WORKDIR\n\n")
