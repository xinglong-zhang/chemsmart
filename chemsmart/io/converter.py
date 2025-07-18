import logging
import os

from chemsmart.io.gaussian.folder import GaussianComFolder, GaussianLogFolder
from chemsmart.io.gaussian.input import Gaussian16Input
from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.io.molecules.structure import SDFFile
from chemsmart.io.orca.folder import ORCAInpFolder, ORCAOutFolder
from chemsmart.io.orca.input import ORCAInput
from chemsmart.io.orca.output import ORCAOutput
from chemsmart.io.xyz.file import XYZFile
from chemsmart.io.xyz.folder import XYZFolder
from chemsmart.utils.logger import create_logger
from chemsmart.utils.mixins import BaseFolder

logger = logging.getLogger(__name__)
os.environ["OMP_NUM_THREADS"] = "1"

create_logger()


class FileConverter:
    """Class for converting files in different formats.
    Args:
        directory (str): Directory in which to convert files.
        type (str): Type of file to be converted, if directory is specified.
        filename (str): Input filename to be converted.
        output_filetype (str): Type of files to convert to, defaults to .xzy.
    """

    def __init__(
        self,
        directory=None,
        type=None,
        filename=None,
        output_filetype="xyz",
        include_intermediate_structures=False,
    ):
        self.directory = directory
        self.type = type
        self.filename = filename
        self.output_filetype = output_filetype
        self.include_intermediate_structures = include_intermediate_structures

    def convert_files(self):
        if self.directory is not None:
            logger.info(f"Converting files in directory: {self.directory}")
            assert (
                self.type is not None
            ), "Type of file to be converted must be specified."
            self._convert_all_files(
                self.directory, self.type, self.output_filetype
            )
        else:
            if self.filename is not None:
                # get filetype/extension from filename
                self.type = self.filename.split(".")[-1]
                logger.info(f"Converting file: {self.filename}")
                self._convert_single_file(self.filename, self.output_filetype)
            else:
                raise ValueError(
                    "Either directory or filename must be specified."
                )

    def _convert_all_files(self, directory, type, output_filetype):
        """Convert all files of specified type in the directory."""
        if type == "log":
            g16_folder = GaussianLogFolder(folder=directory)
            all_files = g16_folder.all_logfiles
        elif type == "com":
            g16_folder = GaussianComFolder(folder=directory)
            all_files = g16_folder.all_com_files
        elif type == "gjf":
            g16_folder = GaussianComFolder(folder=directory)
            all_files = g16_folder.all_gjf_files
        elif type == "out":
            orca_folder = ORCAOutFolder(folder=directory)
            all_files = orca_folder.all_outfiles
        elif type == "inp":
            orca_folder = ORCAInpFolder(folder=directory)
            all_files = orca_folder.all_inpfiles
        elif type == "xyz":
            xyz_folder = XYZFolder(folder=directory)
            all_files = xyz_folder.all_xyzfiles
        elif type == "sdf":
            sdf_folder = BaseFolder(folder=directory)
            all_files = sdf_folder.get_all_files_in_current_folder_and_subfolders_by_suffix(
                filetype="sdf"
            )
        else:
            raise ValueError(f"File type {type} is not supported.")

        logger.debug(f"Files to be converted: {all_files}")

        for file in all_files:
            logger.info(f"Converting file: {file}")
            if type == "log":
                outfile = Gaussian16Output(filename=file)
            elif type == "com" or type == "gjf":
                outfile = Gaussian16Input(filename=file)
            elif type == "out":
                outfile = ORCAOutput(filename=file)
            elif type == "inp":
                outfile = ORCAInput(filename=file)
            elif type == "xyz":
                outfile = XYZFile(filename=file)
            elif type == "sdf":
                outfile = SDFFile(filename=file)
            else:
                raise ValueError(f"File type {type} is not supported.")
            if self.include_intermediate_structures:
                mol = outfile.all_structures
            else:
                mol = (
                    outfile.molecule
                )  # for xyz file with multiple mols, only converts the last one
            filedir, filename = os.path.split(file)
            file_basename = os.path.splitext(filename)[0]
            output_filepath = os.path.join(
                filedir, f"{file_basename}.{output_filetype}"
            )
            if isinstance(mol, list):
                for i, m in enumerate(mol):
                    output_filepath = os.path.join(
                        filedir, f"{file_basename}.{output_filetype}"
                    )
                    m.write(output_filepath, format=output_filetype)
            else:
                mol.write(output_filepath, format=output_filetype)

    def _convert_single_file(self, filename, output_filetype):
        """Convert single file to specified format."""
        logger.info(f"Converting file type: {self.type}")
        if self.type == "log":
            outfile = Gaussian16Output(filename=filename)
        elif self.type == "com" or self.type == "gjf":
            outfile = Gaussian16Input(filename=filename)
        elif self.type == "out":
            outfile = ORCAOutput(filename=filename)
        elif self.type == "inp":
            outfile = ORCAInput(filename=filename)
        elif self.type == "xyz":
            outfile = XYZFile(filename=filename)
        elif self.type == "sdf":
            outfile = SDFFile(filename=filename)
        else:
            raise ValueError(f"File type {self.type} is not supported.")
        if self.include_intermediate_structures:
            mol = outfile.all_structures
        else:
            mol = outfile.molecule
        filedir, filename = os.path.split(filename)
        file_basename = os.path.splitext(filename)[0]
        output_filepath = os.path.join(
            filedir, f"{file_basename}.{output_filetype}"
        )
        if isinstance(mol, list):
            for m in mol:
                output_filepath = os.path.join(
                    filedir, f"{file_basename}.{output_filetype}"
                )
                m.write(output_filepath, format=output_filetype)
        else:
            mol.write(output_filepath, format=output_filetype)
