import logging
import os
import tempfile

from chemsmart.io.file import CDXFile, SDFFile
from chemsmart.io.folder import BaseFolder
from chemsmart.io.gaussian.folder import (
    GaussianInputFolder,
    GaussianOutputFolder,
)
from chemsmart.io.gaussian.input import Gaussian16Input
from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.io.orca.folder import ORCAInputFolder, ORCAOutputFolder
from chemsmart.io.orca.input import ORCAInput
from chemsmart.io.orca.output import ORCAOutput
from chemsmart.io.pdb.pdbfile import PDBFile
from chemsmart.io.xyz.folder import XYZFolder
from chemsmart.io.xyz.xyzfile import XYZFile
from chemsmart.utils.io import get_program_type_from_file
from chemsmart.utils.logger import create_logger

logger = logging.getLogger(__name__)
os.environ["OMP_NUM_THREADS"] = "1"

create_logger()


class FileConverter:
    """Class for converting files in different formats.
    Args:
        directory (str): Directory in which to convert files.
        type (str): Type of file to be converted, if directory is specified.
        program (str | None): Computational chemistry program whose output files
            should be converted. Only required when converting files with
            shared extensions.
        filename (str): Input filename to be converted.
        output_filetype (str): Type of files to convert to, defaults to xyz.
        output_filepath (str | None): Explicit output file path. When provided,
            the file is written to this path and the format is inferred from
            the extension. Takes precedence over ``output_filetype``.
        include_intermediate_structures (bool): Include intermediate structures.
    """

    def __init__(
        self,
        directory=None,
        type=None,
        program=None,
        filename=None,
        output_filetype="xyz",
        output_filepath=None,
        include_intermediate_structures=False,
    ):
        self.directory = directory
        self.type = type
        self.program = program
        self.filename = filename
        self.output_filetype = output_filetype
        self.output_filepath = output_filepath
        self.include_intermediate_structures = include_intermediate_structures

    def convert_files(self):
        """
        Convert files based on the specified parameters.

        Converts either all files in a directory (if directory is specified)
        or a single file (if filename is
        specified) to the target output format.
        """
        if self.directory is not None:
            logger.info(f"Converting files in directory: {self.directory}")
            assert (
                self.type is not None
            ), "Type of file (--filetype) to be converted must be specified."
            if self.type == "out" and self.program is None:
                raise ValueError(
                    "Both --filetype out and --program must be specified when "
                    "converting .out files, because both Gaussian and ORCA use "
                    "this extension. Use --program gaussian or --program orca."
                )
            self._convert_all_files(
                self.directory, self.type, self.output_filetype
            )
        else:
            if self.filename is not None:
                # get filetype/extension from filename
                self.type = self.filename.split(".")[-1]
                logger.info(f"Converting file: {self.filename}")
                if self.output_filepath is not None:
                    self.convert_file(self.filename, self.output_filepath)
                else:
                    self._convert_single_file(
                        self.filename, self.output_filetype
                    )
            else:
                raise ValueError(
                    "Either directory or filename must be specified."
                )

    def convert_file(self, input_filepath, output_filepath):
        """
        Convert a single file from one format to another via Molecule object.

        This is the preferred generic method for converting between any two
        supported file formats. It uses :meth:`Molecule.from_filepath` to
        load the input file and :meth:`Molecule.write` to write the output,
        so any file type supported by those methods is accepted.

        Args:
            input_filepath (str): Path to the input file.
            output_filepath (str): Path to the output file. The format is
                inferred from the file extension.

        Raises:
            ValueError: If the input file cannot be read or the output format
                is not supported by :meth:`Molecule.write`.
        """
        from chemsmart.io.molecules.structure import Molecule

        output_format = os.path.splitext(output_filepath)[1].lstrip(".")
        logger.info(f"Converting {input_filepath} to {output_filepath}")

        if self.include_intermediate_structures:
            mols = Molecule.from_filepath(
                input_filepath, index=":", return_list=True
            )
            if mols is None:
                raise ValueError(
                    f"Could not read molecule(s) from {input_filepath}"
                )
            for m in mols:
                m.write(output_filepath, format=output_format)
        else:
            mol = Molecule.from_filepath(input_filepath)
            if mol is None:
                raise ValueError(
                    f"Could not read molecule from {input_filepath}"
                )
            mol.write(output_filepath, format=output_format)
        logger.info(f"Created: {output_filepath}")

    def _convert_all_files(self, directory, type, output_filetype):
        """
        Convert all files of specified type in the directory.

        Args:
            directory (str): Directory containing files to convert.
            type (str): File type to convert
            (log, com, gjf, out, inp, xyz, sdf, pdb).
            output_filetype (str): Target output format.
        """
        if type == "log":
            g16_folder = GaussianOutputFolder(folder=directory)
            all_files = g16_folder.all_log_files
        elif type == "com":
            g16_folder = GaussianInputFolder(folder=directory)
            all_files = g16_folder.all_com_files
        elif type == "gjf":
            g16_folder = GaussianInputFolder(folder=directory)
            all_files = g16_folder.all_gjf_files
        elif type == "out":
            if self.program == "gaussian":
                g16_folder = GaussianOutputFolder(folder=directory)
                all_files = [
                    f
                    for f in g16_folder.all_output_files
                    if f.endswith(f".{type}")
                ]
            else:
                orca_folder = ORCAOutputFolder(folder=directory)
                all_files = [
                    f
                    for f in orca_folder.all_output_files
                    if f.endswith(f".{type}")
                ]
        elif type == "inp":
            orca_folder = ORCAInputFolder(folder=directory)
            all_files = orca_folder.all_inp_files
        elif type == "xyz":
            xyz_folder = XYZFolder(folder=directory)
            all_files = xyz_folder.all_xyzfiles
        elif type == "sdf":
            sdf_folder = BaseFolder(folder=directory)
            all_files = sdf_folder.get_all_files_in_current_folder_and_subfolders_by_suffix(
                filetype="sdf"
            )
        elif type == "pdb":
            pdb_folder = BaseFolder(folder=directory)
            all_files = pdb_folder.get_all_files_in_current_folder_and_subfolders_by_suffix(
                filetype="pdb"
            )
        elif type in ("cdxml", "cdx"):
            cdx_folder = BaseFolder(folder=directory)
            all_files = cdx_folder.get_all_files_in_current_folder_and_subfolders_by_suffix(
                filetype=type
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
                if self.program == "gaussian":
                    outfile = Gaussian16Output(filename=file)
                else:
                    outfile = ORCAOutput(filename=file)
            elif type == "inp":
                outfile = ORCAInput(filename=file)
            elif type == "xyz":
                outfile = XYZFile(filename=file)
            elif type == "sdf":
                outfile = SDFFile(filename=file)
            elif type == "pdb":
                outfile = PDBFile(filename=file)
            elif type in ("cdxml", "cdx"):
                cdxfile = CDXFile(filename=file)
                mols = cdxfile.molecules
                filedir, fname = os.path.split(file)
                file_basename = os.path.splitext(fname)[0]
                if len(mols) == 1:
                    output_path = os.path.join(
                        filedir, f"{file_basename}.{output_filetype}"
                    )
                    mols[0].write(output_path, format=output_filetype)
                    logger.info(f"Created: {output_path}")
                else:
                    for i, m in enumerate(mols, start=1):
                        output_path = os.path.join(
                            filedir,
                            f"{file_basename}_{i}.{output_filetype}",
                        )
                        m.write(output_path, format=output_filetype)
                        logger.info(f"Created: {output_path}")
                continue
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
        """
        Convert single file to specified format.

        Args:
            filename (str): Path to the file to convert.
            output_filetype (str): Target output format.
        """
        logger.info(f"Converting file type: {self.type}")
        if self.type == "log" or self.type == "out":
            detected = get_program_type_from_file(filename)
            if detected == "gaussian":
                outfile = Gaussian16Output(filename=filename)
            elif detected == "orca":
                outfile = ORCAOutput(filename=filename)
            else:
                raise ValueError(
                    f"Could not detect program type for '{filename}'. "
                )
        elif self.type == "com" or self.type == "gjf":
            outfile = Gaussian16Input(filename=filename)
        elif self.type == "inp":
            outfile = ORCAInput(filename=filename)
        elif self.type == "xyz":
            outfile = XYZFile(filename=filename)
        elif self.type == "sdf":
            outfile = SDFFile(filename=filename)
        elif self.type == "pdb":
            outfile = PDBFile(filename=filename)
        elif self.type in ("cdxml", "cdx"):
            cdxfile = CDXFile(filename=filename)
            mols = cdxfile.molecules
            filedir, fname = os.path.split(filename)
            file_basename = os.path.splitext(fname)[0]
            if len(mols) == 1:
                output_path = os.path.join(
                    filedir, f"{file_basename}.{output_filetype}"
                )
                mols[0].write(output_path, format=output_filetype)
                logger.info(f"Created: {output_path}")
            else:
                for i, m in enumerate(mols, start=1):
                    output_path = os.path.join(
                        filedir, f"{file_basename}_{i}.{output_filetype}"
                    )
                    m.write(output_path, format=output_filetype)
                    logger.info(f"Created: {output_path}")
            return
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

    @staticmethod
    def xyz_to_pdb(
        molecule,
        pdb_filename,
        xyz_filename=None,
        mode="w",
        overwrite=True,
        cleanup=True,
    ):
        """
        Convert a molecule's XYZ representation to PDB using Open Babel.

        This is an alternative to ``Molecule.write_pdb`` / ``Molecule.to_pdb``
        for cases where Open Babel's interpretation of connectivity or atom
        typing is preferred over the RDKit-based path.

        Args:
            molecule (Molecule): Source molecule whose coordinates are used.
            pdb_filename (str): Destination PDB file path.
            xyz_filename (str, optional): Path to an existing XYZ file to
                convert. When omitted or the file is absent, a temporary XYZ
                is written from *molecule* via ``write_xyz``.
            mode (str): File mode passed to ``write_xyz`` when the XYZ file
                must be created. Default ``'w'``.
            overwrite (bool): Whether to overwrite *pdb_filename* if it
                already exists. Default ``True``.
            cleanup (bool): Remove any auto-generated temporary XYZ file after
                the conversion. Default ``True``.

        Raises:
            ImportError: If Open Babel (``openbabel``) is not installed.
            ValueError: If the XYZ file cannot be parsed by Open Babel.
        """
        auto_xyz = False
        if xyz_filename is None:
            tmp = tempfile.NamedTemporaryFile(suffix=".xyz", delete=False)
            tmp.close()
            xyz_filename = tmp.name
            auto_xyz = True
            logger.debug(
                f"Created temporary XYZ {xyz_filename} for PDB conversion."
            )
            molecule.write_xyz(xyz_filename, mode=mode)
        elif not os.path.isfile(xyz_filename):
            logger.debug(
                f"XYZ {xyz_filename} missing; writing coordinates before "
                "conversion."
            )
            molecule.write_xyz(xyz_filename, mode=mode)

        try:
            from openbabel import pybel
        except ImportError as exc:  # pragma: no cover
            if auto_xyz and cleanup:
                os.remove(xyz_filename)
            raise ImportError(
                "xyz_to_pdb requires Open Babel. "
                "Use 'conda install -c conda-forge openbabel' to install it."
            ) from exc

        xyz_mol = next(pybel.readfile("xyz", xyz_filename), None)
        if xyz_mol is None:
            if auto_xyz and cleanup:
                os.remove(xyz_filename)
            raise ValueError(f"Unable to read molecule from {xyz_filename}")

        logger.info(
            f"Converting XYZ {xyz_filename} to PDB {pdb_filename} via "
            f"Open Babel (overwrite={overwrite})"
        )
        xyz_mol.write("pdb", pdb_filename, overwrite=overwrite)

        if auto_xyz and cleanup:
            try:
                os.remove(xyz_filename)
                logger.debug(f"Removed temporary XYZ {xyz_filename}")
            except OSError as exc:
                logger.warning(
                    f"Could not remove temporary XYZ {xyz_filename}: {exc}"
                )
