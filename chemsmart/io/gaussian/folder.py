import logging

from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.utils.mixins import FolderMixin
from chemsmart.utils.periodictable import PeriodicTable

p = PeriodicTable()
logger = logging.getLogger(__name__)


class GaussianComFolder(FolderMixin):
    def __init__(self, folder):
        self.folder = folder

    @property
    def all_comfiles(self):
        return self.get_all_files_in_current_folder_and_subfolders(
            filetype="com"
        )

    @property
    def all_comfiles_in_current_folder(self):
        return self.get_all_files_in_current_folder(filetype="com")


class GaussianLogFolder(FolderMixin):
    """Log folder containing all Gaussian log files for postprocessing."""

    def __init__(self, folder):
        """:param folder: Parent folder for all log files; type of str"""
        self.folder = folder

    @property
    def all_logfiles(self):
        """Get all log files in the folder, including subfolders."""
        return self.get_all_files_in_current_folder_and_subfolders(
            filetype="log"
        )

    @property
    def all_logfiles_in_current_folder(self):
        """Get all log files in the folder."""
        return self.get_all_files_in_current_folder(filetype="log")

    @property
    def total_service_units(self):
        """Get all service units used in all the log files contained in the folder."""
        total_service_units = 0
        for file in self.all_logfiles:
            output_file = Gaussian16Output(file)
            core_hours = output_file.total_core_hours
            total_service_units += core_hours
        return total_service_units

    def write_job_runtime(self, job_runtime_file="job_runtime.txt"):
        """Write job runtime for all log files in the folder."""
        with open(job_runtime_file, "w") as f:
            for file in self.all_logfiles:
                output_file = Gaussian16Output(file)
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
    #     for file in self.all_logfiles:
    #         output_file = Gaussian16Output(file)
    #         database[file] = output_file.__dict__
    #     with open(database_file, 'w') as f:
    #         json.dump(database, f, indent=4)
    #     return database
