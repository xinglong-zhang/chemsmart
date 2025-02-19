import logging
import os
import shutil
from shutil import copyfile
from subprocess import check_output
import pandas as pd

from chemsmart.utils.logger import create_logger

logger = logging.getLogger(__name__)
os.environ["OMP_NUM_THREADS"] = "1"

create_logger()


class FileOrganizer:
    """Class for organizing files in an Excel file into folder.
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

        Notes: (1) the class accepts empty rows
               (2) no filenames should be the same in column C (this should be the case for Gaussian jobs)
               (3) no ".log" extension is needed for each column in Excel document but this script deals
                   with ".log" files only for now.


    The final result will be that inside the current directory, a number of subdirectories named in column A
    will be created and the files in column C will be named to the corresponding names in column B and moved
    to their corresponding folders in column A.


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
        type=".log",
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

    def load_data(self):
        """Load data from Excel file."""
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
        """Clean data from Excel file by removing empty structure names in second column."""
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
        for folder in folders:
            if isinstance(folder, str) and not os.path.exists(folder):
                target_folder = os.path.join(self.directory, folder)
                os.makedirs(target_folder, exist_ok=True)
        return None

def find_absolute_path(self, filename):
    if '.log' not in filename:
        filename += '.log'

    try:
        absFile = check_output(f"find . -name {filename}", shell=True).decode("utf-8").strip()
        absDir = check_output(f"find . -name {filename} -exec dirname {{}} ';'", shell=True).decode("utf-8").strip()
        return absFile.split('\n')[0], absDir.split('\n')[0]
    except FileNotFoundError:
        print(f"{filename} not found! Check your Excel file.")
        return None, None

