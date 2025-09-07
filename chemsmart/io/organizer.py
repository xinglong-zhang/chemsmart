import logging
import os
from shutil import copy

import pandas as pd

from chemsmart.utils.logger import create_logger
from chemsmart.utils.utils import search_file

logger = logging.getLogger(__name__)
os.environ["OMP_NUM_THREADS"] = "1"

create_logger()


class FileOrganizer:
    """
    Class for organizing files in an Excel file into folder.
    It sorts a directory of files into ordered folders, at the same time, it
    renames the filename used during calculation to final file name to be used
    in supporting information.
    It reads in an Excel file with the following structure (shown an example):

        column A         | column B             | column C
        folder_pathway1  | finalname_rct        | runtimename_input_rct_wb97xd
                         | finalname_ts         | runtimename_input_ts_wb97xd
                         | finalname_prc        | runtimename_input_prd_wb97xd
                         |                      |
        folder_pathway2  | finalname_rct2       | runtimename_input_rct2_wb97xd
                         | finalname_ts2        | runtimename_input_ts2_wb97xd
                         | finalname_prc2       | runtimename_input_prd2_wb97xd
                         |                      |
                         |                      |
        folder_pathway3  | finalname_rct3       | runtimename_input_rct3_wb97xd
                         | finalname_ts3        | runtimename_input_ts3_wb97xd
                         | finalname_prc3       | runtimename_input_prd3_wb97xd

    Notes: (1) empty rows are okay and will be skipped;
           (2) no filenames should be the same in column C
           (3) no ".log" (or other) extension is needed for each column in Excel document.

    The final result will be that inside the current directory, a number of subdirectories
    named in column A will be created and the files in column C will be named to the corresponding
    names in column B and moved to their corresponding folders in column A.

    Args:
        directory (str): Directory in which the files are searched.
        filename (str): Excel filename to be used for organizing data.
        sheetname (str): Excel sheetname to be used for organizing data.
        type (str): File type to be searched for.
        cols (str): Columns to be read from Excel file.
        skip (int): Number of rows to be skipped in Excel file.
        row (int): Number of rows to be read in Excel file.
        keep_default_na (bool): Keep default na values in Excel file.
    """

    def __init__(
        self,
        directory=None,
        filename=None,
        sheetname=None,
        type="log",
        cols="B:D",
        skip=2,
        row=100,
        keep_default_na=False,
    ):
        if directory is None:
            directory = "."
        self.directory = os.path.abspath(directory)
        assert os.path.exists(filename), f"{filename} does not exist!"
        self.filename = filename
        self.sheetname = sheetname
        self.type = type
        self.cols = cols
        self.skip = skip
        self.row = row
        self.keep_default_na = keep_default_na
        self.data = self.load_data()

    # noinspection PyTypeChecker
    def load_data(self):
        """
        Load data from Excel file.
        
        Reads data from the Excel file using the specified parameters including
        sheet name, columns, rows to skip, and other formatting options.
        
        Returns:
            pd.DataFrame: DataFrame containing the loaded Excel data.
        """
        df = pd.read_excel(
            io=self.filename,
            sheet_name=self.sheetname,
            usecols=self.cols,
            skiprows=self.skip,
            nrows=self.row,
            keep_default_na=self.keep_default_na,
        )
        return df

    def clean_data(self):
        """
        Clean data from Excel file by removing empty structure names in second column.
        
        Processes the loaded Excel data to remove rows where the second column
        (structure names) is empty or contains only whitespace. Converts all
        valid entries to strings and strips whitespace.
        
        Returns:
            tuple: A tuple containing three lists:
                - new_list1 (list): Cleaned folder names from first column
                - new_list2 (list): Cleaned final filenames from second column  
                - new_list3 (list): Cleaned runtime filenames from third column
        """
        columns = self.data.columns

        list1, list2, list3 = (
            self.data[columns[0]],
            self.data[columns[1]],
            self.data[columns[2]],
        )

        new_list1, new_list2, new_list3 = [], [], []
        for i in range(len(list1)):
            if isinstance(list2[i], str) and list2[i].strip():
                new_list1.append(str(list1[i]).strip())
                new_list2.append(str(list2[i]).strip())
                new_list3.append(str(list3[i]).strip())

        return new_list1, new_list2, new_list3

    def create_directories(self, folders):
        """
        Create directories for organizing files.
        
        Creates subdirectories within the working directory based on the 
        provided folder names. Only creates directories that don't already exist.
        
        Args:
            folders (list): List of folder names to create.
            
        Returns:
            None
        """
        for folder in folders:
            if isinstance(folder, str) and not os.path.exists(folder):
                target_folder = os.path.join(self.directory, folder)
                os.makedirs(target_folder, exist_ok=True)
        return None

    def copy_and_rename(self, target_folder, new_filename, old_filename):
        """
        Copy and rename files to target folder.
        
        Searches for the original file, adds file extensions if needed, and 
        copies the file to the target folder with the new filename.
        
        Args:
            target_folder (str): Destination folder for the file.
            new_filename (str): New name for the file.
            old_filename (str): Original name of the file to copy.
            
        Returns:
            None
        """
        if self.type not in old_filename:
            old_filename += f".{self.type}"
        if self.type not in new_filename:
            new_filename += f".{self.type}"

        absolute_file_path, _ = search_file(old_filename)
        new_file_path = os.path.join(
            self.directory, target_folder, new_filename
        )
        if absolute_file_path is not None:
            logger.info(f"Copying {absolute_file_path} to {new_file_path}.")
            copy(absolute_file_path, new_file_path)
        return None

    def organize_files(self):
        """
        Organize files in Excel file into folders.
        
        Main method that orchestrates the file organization process:
        1. Cleans the Excel data to get valid folder names and filenames
        2. Creates necessary directories
        3. Copies and renames files to their target locations
        
        Returns:
            None
        """
        folders, new_filenames, old_filenames = self.clean_data()
        self.create_directories(folders)
        for folder, new_filename, old_filename in zip(
            folders, new_filenames, old_filenames
        ):
            self.copy_and_rename(folder, new_filename, old_filename)
        return None
