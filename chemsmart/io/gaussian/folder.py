import logging

from chemsmart.io.folder import BaseFolder
from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.utils.periodictable import PeriodicTable

p = PeriodicTable()
logger = logging.getLogger(__name__)


class GaussianInputFolder(BaseFolder):
    """
    Manager for directories containing Gaussian input files.
    """

    def __init__(self, folder):
        """
        Initialize Gaussian input folder manager.
        """
        super().__init__(folder=folder)

    @property
    def all_input_files(self):
        """
        Get all Gaussian input files (.com, .gjf) in folder and subfolders.
        """
        com_files = (
            self.get_all_files_in_current_folder_and_subfolders_by_suffix(
                filetype="com"
            )
        )
        gjf_files = (
            self.get_all_files_in_current_folder_and_subfolders_by_suffix(
                filetype="gjf"
            )
        )
        return com_files + gjf_files

    @property
    def all_input_files_in_current_folder(self):
        """
        Get all Gaussian input files (.com, .gjf) in the current folder only.
        """
        com_files = self.get_all_files_in_current_folder_by_suffix(
            filetype="com"
        )
        gjf_files = self.get_all_files_in_current_folder_by_suffix(
            filetype="gjf"
        )
        return com_files + gjf_files

    @property
    def all_com_files(self):
        """
        Get all .com files in folder and subfolders.
        """
        return self.get_all_files_in_current_folder_and_subfolders_by_suffix(
            filetype="com"
        )

    @property
    def all_com_files_in_current_folder(self):
        """
        Get all .com files in the current folder only.
        """
        return self.get_all_files_in_current_folder_by_suffix(filetype="com")

    @property
    def all_gjf_files(self):
        """
        Get all .gjf files in folder and subfolders.
        """
        return self.get_all_files_in_current_folder_and_subfolders_by_suffix(
            filetype="gjf"
        )

    @property
    def all_gjf_files_in_current_folder(self):
        """
        Get all .gjf files in the current folder only.
        """
        return self.get_all_files_in_current_folder_by_suffix(filetype="gjf")


class GaussianOutputFolder(BaseFolder):
    """
    Output folder containing all Gaussian output files for postprocessing.
    """

    def __init__(self, folder):
        """
        Initialize Gaussian output folder manager.
        """
        super().__init__(folder=folder)

    @property
    def all_output_files(self):
        """
        Get all Gaussian output files in the folder, including subfolders.
        """
        return self.get_all_output_files_in_current_folder_and_subfolders_by_program(
            program="gaussian"
        )

    @property
    def all_output_files_in_current_folder(self):
        """
        Get all Gaussian output files in the current folder.
        """
        return self.get_all_output_files_in_current_folder_by_program(
            program="gaussian"
        )

    @property
    def all_log_files(self):
        """
        Get all .log files in the folder, including subfolders.
        Uses suffix-based detection (fast but may miss .out files).
        For content-based detection, use all_output_files instead.
        """
        return self.get_all_files_in_current_folder_and_subfolders_by_suffix(
            filetype="log"
        )

    @property
    def all_log_files_in_current_folder(self):
        """
        Get all .log files in the current folder.
        Uses suffix-based detection (fast but may miss .out files).
        """
        return self.get_all_files_in_current_folder_by_suffix(filetype="log")

    @property
    def all_molecules(self):
        """
        Get all molecules in the folder.
        """
        molecules = []
        for file in self.all_output_files:
            output_file = Gaussian16Output(file)
            molecules.append(output_file.molecule)
        return molecules

    @property
    def total_service_units(self):
        """
        Get all service units used in all the output files contained in the folder.
        """
        total_service_units = 0
        for file in self.all_output_files:
            output_file = Gaussian16Output(file)
            core_hours = output_file.total_core_hours
            total_service_units += core_hours
        return total_service_units

    def write_job_runtime(self, job_runtime_file="job_runtime.txt"):
        """
        Generate runtime report for all output files in the folder.
        """
        with open(job_runtime_file, "w") as f:
            for file in self.all_output_files:
                output_file = Gaussian16Output(file)
                core_hours = output_file.total_core_hours
                f.write(
                    f"Job: {file:<130} Total time: {core_hours:6.1f} "
                    f"core-hours\n"
                )
            f.write(
                f"TOTAL core-hours in folder {self.folder} is: "
                f"{self.total_service_units}\n"
            )
