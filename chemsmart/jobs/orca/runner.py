import logging
import os
import shlex
import shutil
import subprocess
from contextlib import suppress
from glob import glob
from shutil import copy, rmtree

from pyatoms.cli.submitters import SubmitscriptWriter
from pyatoms.io.orca.inputs import ORCAInput
from pyatoms.jobs.orca.execution import OrcaExecutables
from pyatoms.jobs.runner import JobRunner
from pyatoms.utils.periodictable import chemical_symbol_to_atomic_number

shutil._USE_CP_SENDFILE = False
# to avoid "BlockingIOError: [Errno 11] Resource temporarily unavailable:" Error when copying

logger = logging.getLogger(__name__)


class ORCAJobRunner(JobRunner):
    # creates job runner process
    # combines information about server and program

    JOBTYPES = [
        "orcainp",
        "orcaopt",
        "orcamodred",
        "orcascan",
        "orcats",
        "orcasp",
        "orcairc",
    ]

    def __init__(self, server, num_nodes, scratch=True, **kwargs):
        super().__init__(
            server=server, num_nodes=num_nodes, scratch=scratch, **kwargs
        )

    @property
    def executables(self):
        return OrcaExecutables.from_servername()

    def _prerun(self, job):
        self._assign_variables(job)
        self._copy_over_files(job)

    def _assign_variables(self, job):
        """Creates input file in scratch (done in settings.apply_on()) if running in scratch directory.

        Set up file paths if running in scratch/not running in scratch.
        """
        # keep job output file in job folder regardless of running in scratch or not
        self.job_outputfile = job.outputfile

        if self.scratch:
            logger.info("Setting up run in scratch folder.")
            # set up files in scratch folder
            scratch_job_dir = os.path.join(self.scratch_dir, job.label)
            self.running_directory = scratch_job_dir

            job_inputfile = job.label + ".inp"
            scratch_job_inputfile = os.path.join(
                scratch_job_dir, job_inputfile
            )
            self.job_inputfile = scratch_job_inputfile

            job_gbwfile = job.label + ".gbw"
            scratch_job_gbwfile = os.path.join(scratch_job_dir, job_gbwfile)
            self.job_gbwfile = scratch_job_gbwfile

            job_errfile = job.label + ".err"
            scratch_job_errfile = os.path.join(scratch_job_dir, job_errfile)
            self.job_errfile = scratch_job_errfile
        else:
            logger.info("Setting up run in job folder (no scratch).")
            # keep files as in runnind directory
            self.running_directory = job.folder
            self.job_inputfile = job.inputfile
            self.job_gbwfile = job.gbwfile
            self.job_errfile = job.errfile

        if self.executables.local_run is not None:
            logger.info("Local run is True.")
            job.local = self.executables.local_run

    def _copy_over_files(self, job):
        if self.scratch:
            logger.info("Copying over files to scratch folder.")
            # copy over files
            scratch_job_dir = os.path.join(self.scratch_dir, job.label)
            with suppress(FileExistsError):
                os.mkdir(scratch_job_dir)
            for file in glob(f"{job.folder}/*xyz"):
                # copy over xyz files if running in scratch (the xyz files can be part of the input specification
                # e.g., for neb-TS calculations)
                logger.info(
                    f"Copying file {file} from {job.folder} \nto {scratch_job_dir}\n"
                )
                copy(file, scratch_job_dir)
            # copy over hessian files
            for file in glob(f"{job.folder}/*hess"):
                logger.info(
                    f"Copying file {file} from {job.folder} \nto {scratch_job_dir}\n"
                )
                copy(file, scratch_job_dir)
            # copy over gbw files
            for file in glob(f"{job.folder}/*gbw"):
                logger.info(
                    f"Copying file {file} from {job.folder} \nto {scratch_job_dir}\n"
                )
                copy(file, scratch_job_dir)

    def _run(self, job, process):
        process.communicate()
        return process.poll()

    def _environment_vars(self, job, *args, **kwargs):
        return self.executables.ENVIRONMENT_VARIABLES

    def _create_process(self, job, command, env):
        with (
            open(self.job_outputfile, "w") as out,
            open(self.job_errfile, "w") as err,
        ):
            logger.info(
                f"Command executed: {command}\n"
                f"Writing output file to: {self.job_outputfile} and err file to: {self.job_errfile}"
            )
            return subprocess.Popen(
                shlex.split(command),
                stdout=out,
                stderr=err,
                env=env,
                cwd=self.running_directory,
            )

    def _get_executable(self):
        """Get executable for orca."""
        exe = self.executables.get_executable()
        logger.info(f"Orca executable: {exe}")
        return exe

    def get_command(self, job, host):
        self._assign_variables(job)
        executable = self._get_executable()
        resources = self._resources(job=job)

        command = self.server.command(
            folder=self.running_directory,
            executable=executable,
            resources=resources,
            local=job.local,
            host=host,
        )

        return f"{command} {self.job_inputfile}"

    def _postrun(self, job):
        if self.scratch:
            # if job was run in scratch, copy files to job folder except files ending with .tmp or .tmp.*
            for file in glob(f"{self.running_directory}/{job.label}*"):
                if file.endswith(".tmp") or ".tmp." in file:
                    continue
                logger.info(
                    f"Copying file {file} from {self.running_directory} \nto {job.folder}\n"
                )
                copy(file, job.folder)

        if job.is_complete():
            # if job is completed, remove scratch directory and submit_script and log.info and log.err files
            if self.scratch:
                logger.info(
                    f"Removing scratch directory: {self.running_directory}."
                )
                rmtree(self.running_directory)

            writer = SubmitscriptWriter(job)
            submit_script = writer.job_submit_script
            run_script = writer.job_run_script
            err_filepath = os.path.join(job.folder, f"{job.errfile}")
            joblogerr_filepath = os.path.join(job.folder, "log.err")
            jobloginfo_filepath = os.path.join(job.folder, "log.info")
            pbs_errfile = os.path.join(job.folder, "pbs.err")
            pbs_infofile = os.path.join(job.folder, "pbs.info")

            files_to_remove = [
                submit_script,
                run_script,
                err_filepath,
                joblogerr_filepath,
                jobloginfo_filepath,
                pbs_errfile,
                pbs_infofile,
            ]

            for f in files_to_remove:
                with suppress(FileNotFoundError):
                    os.remove(f)


class FakeORCA:
    def __init__(self, file_to_run):
        if not os.path.exists(file_to_run):
            raise FileNotFoundError(f"File {file_to_run} not found.")
        self.file_to_run = os.path.abspath(file_to_run)
        self.input_object = ORCAInput(inpfile=self.file_to_run)

    @property
    def file_folder(self):
        return os.path.dirname(self.file_to_run)

    @property
    def filename(self):
        return os.path.basename(self.file_to_run)

    @property
    def input_filepath(self):
        return os.path.join(self.file_folder, self.file_to_run)

    @property
    def output_filepath(self):
        output_file = self.filename.split(".")[0] + ".out"
        return os.path.join(self.file_folder, output_file)

    @property
    def input_contents(self):
        return self.input_object.contents

    @property
    def molecule(self):
        return self.input_object.atoms

    @property
    def charge(self):
        return self.input_object.charge

    @property
    def multiplicity(self):
        return self.input_object.multiplicity

    @property
    def spin(self):
        if self.multiplicity == 1:
            return "R"
        return "U"

    @property
    def num_atoms(self):
        return len(self.molecule.chemical_symbols)

    @property
    def atomic_symbols(self):
        return list(self.molecule.chemical_symbols)

    @property
    def atomic_numbers(self):
        return [
            chemical_symbol_to_atomic_number(s) for s in self.atomic_symbols
        ]

    @property
    def atomic_coordinates(self):
        return self.molecule.coordinates.coordinates

    @property
    def empirical_formula(self):
        return self.molecule.get_chemical_formula(empirical=True)

    def run(self):
        # get input lines
        with open(self.input_filepath) as f:
            input_lines = f.readlines()

        with open(self.output_filepath, "w") as g:
            g.write("\n")
            g.write(
                """                                 *****************
                                 * O   R   C   A *
                                 *****************\n"""
            )
            g.write(
                "                         Fake Orca version 0.0.0 -  RELEASE  -                \n"
            )
            g.write(
                "================================================================================\n"
            )
            g.write(
                "                                       INPUT FILE                               \n"
            )
            g.write(
                "================================================================================\n"
            )
            g.write(f"NAME = {self.input_filepath}\n")
            for i, line in enumerate(input_lines):
                g.write(f"|  {i+1}> {line}")
            g.write(f"|  {len(input_lines)+1}> ")
            g.write(
                f"|  {len(input_lines) + 2}>                          ****END OF INPUT****\n"
            )
            g.write("                       *****************************\n")
            g.write("                       *        Fake Orca Run      *\n")
            g.write("                       *****************************\n")
            g.write(
                f"Number of atoms                         .... {self.num_atoms}\n"
            )
            g.write(
                f"  Total Charge           Charge          ....    {self.input_object.charge}\n"
            )
            g.write(
                f"  Multiplicity           Mult            ....    {self.input_object.multiplicity}\n"
            )
            g.write(
                "                     ***********************HURRAY********************\n"
            )
            g.write(
                "                     ***        THE OPTIMIZATION HAS CONVERGED     ***\n"
            )
            g.write(
                "                     ***********************HURRAY********************\n"
            )
            g.write(
                "                  *******************************************************\n"
            )
            g.write(
                "                  *** FINAL ENERGY EVALUATION AT THE STATIONARY POINT ***\n"
            )
            g.write(
                "                  *******************************************************\n"
            )
            g.write(" ---------------------------------\n")
            g.write(" CARTESIAN COORDINATES (ANGSTROEM)\n")
            g.write(" ---------------------------------\n")
            for line in input_lines[2:-2]:
                g.write(f"{line}")
            g.write(" ----------------\n")
            g.write(" TOTAL SCF ENERGY\n")
            g.write(" ----------------\n")
            g.write(
                " Total Energy       :          -76.32331101 Eh           -2076.86288 eV\n"
            )
            g.write(
                "                              ****ORCA TERMINATED NORMALLY****\n"
            )
            g.write(
                "TOTAL RUN TIME: 0 days 0 hours 0 minutes 0 seconds xx msec\n"
            )


class FakeORCAJobRunner(ORCAJobRunner):
    # creates job runner process
    # combines information about server and program
    FAKE = True
    JOBTYPES = [
        "orcainp",
        "orcaopt",
        "orcamodred",
        "orcascan",
        "orcats",
        "orcasp",
        "orcairc",
    ]

    def __init__(self, server, num_nodes, scratch=False, **kwargs):
        super().__init__(
            server=server, num_nodes=num_nodes, scratch=scratch, **kwargs
        )

    def run(self, job):
        self._prerun(job=job)
        returncode = FakeORCA(self.job_inputfile).run()
        self._postrun(job=job)
        return returncode
