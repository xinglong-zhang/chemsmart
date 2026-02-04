import logging

from chemsmart.io.folder import BaseFolder
from chemsmart.io.orca.output import ORCAOutput
from chemsmart.utils.periodictable import PeriodicTable

p = PeriodicTable()
logger = logging.getLogger(__name__)


class ORCAInpFolder(BaseFolder):
    """
    Input folder containing all ORCA input files for postprocessing.

    This class provides utilities for managing collections of ORCA input files
    within a directory structure.
    """

    def __init__(self, folder):
        """
        Initialize ORCA input folder handler.

        Args:
            folder (str): Parent folder for all input files
        """
        self.folder = folder

    @property
    def all_inpfiles(self):
        """
        Get all ORCA input files in the folder and subfolders.

        Returns:
            list: Paths to all .inp files found recursively
        """
        return self.get_all_files_in_current_folder_and_subfolders_by_suffix(
            filetype="inp"
        )

    @property
    def all_inpfiles_in_current_folder(self):
        """
        Get all ORCA input files in the current folder only.

        Returns:
            list: Paths to all .inp files in current directory
        """
        return self.get_all_files_in_current_folder_by_suffix(filetype="inp")


class ORCAOutFolder(BaseFolder):
    """
    Output folder containing all ORCA output files for postprocessing.

    This class provides utilities for managing collections of ORCA output files,
    including calculation statistics, runtime analysis, and resource usage tracking.
    """

    def __init__(self, folder):
        """
        Initialize ORCA output folder handler.

        Args:
            folder (str): Parent folder for all output files
        """
        self.folder = folder

    @property
    def all_outfiles(self):
        """
        Get all ORCA output files in the folder and subfolders.

        Returns:
            list: Paths to all .out files found recursively
        """
        return self.get_all_files_in_current_folder_and_subfolders_by_suffix(
            filetype="out"
        )

    @property
    def all_outfiles_in_current_folder(self):
        """
        Get all ORCA output files in the current folder only.

        Returns:
            list: Paths to all .out files in current directory
        """
        return self.get_all_files_in_current_folder_by_suffix(filetype="out")

    @property
    def total_service_units(self):
        """
        Calculate total computational service units used across all output files.

        Returns:
            float: Total core-hours consumed by all calculations in folder
        """
        total_service_units = 0
        for file in self.all_outfiles:
            output_file = ORCAOutput(file)
            core_hours = output_file.total_core_hours
            total_service_units += core_hours
        return total_service_units

    def write_job_runtime(self, job_runtime_file="job_runtime.txt"):
        """
        Write job runtime summary for all output files in the folder.

        Creates a text file listing individual job runtimes and total
        computational resource usage.

        Args:
            job_runtime_file (str): Output file name for runtime summary
        """
        with open(job_runtime_file, "w") as f:
            for file in self.all_outfiles:
                output_file = ORCAOutput(file)
                core_hours = output_file.total_core_hours
                f.write(
                    f"Job: {file:<130} Total time: {core_hours:6.1f} core-hours\n"
                )
            f.write(
                f"TOTAL core-hours in folder {self.folder} is: {self.total_service_units}\n"
            )

    # def assemble_database(self, database_file='database.json'):
    #     """Assemble a database from all log files in the folder."""
    #     database = {}
    #     for file in self.all_outfiles:
    #         output_file = ORCAOutput(file)
    #         database[file] = output_file.__dict__
    #     with open(database_file, 'w') as f:
    #         json.dump(database, f, indent=4)
    #     return database
