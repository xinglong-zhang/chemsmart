import logging
import os

from chemsmart.io.folder import BaseFolder
from chemsmart.io.gaussian.folder import (
    GaussianInputFolder,
    GaussianOutputFolder,
)
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.orca.folder import ORCAInputFolder, ORCAOutputFolder
from chemsmart.io.xyz.folder import XYZFolder
from chemsmart.utils.logger import create_logger

logger = logging.getLogger(__name__)
os.environ["OMP_NUM_THREADS"] = "1"

create_logger()


class FileConverter:
    """Class for converting files in different formats.

    Single-file and batch conversions both go through
    :meth:`convert_file`, which reads via ``Molecule.from_filepath``
    and writes via ``Molecule.write``. When
    ``include_intermediate_structures`` is ``True`` and the input
    contains multiple structures, each structure is written to a
    separate numbered file (``basename_1.ext``, ``basename_2.ext``, ...).

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
            When ``True``, multi-structure inputs produce numbered output files.
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
                output_path = self.resolve_output_path(
                    self.filename,
                    output=self.output_filepath,
                    output_filetype=self.output_filetype,
                )
                self.convert_file(
                    self.filename,
                    output_path,
                    include_intermediate_structures=(
                        self.include_intermediate_structures
                    ),
                )
            else:
                raise ValueError(
                    "Either directory or filename must be specified."
                )

    def _convert_all_files(self, directory, type, output_filetype):
        """
        Convert all files of specified type in the directory.

        Collects matching files then delegates each conversion to
        :meth:`convert_file`.

        Args:
            directory (str): Directory containing files to convert.
            type (str): File type to convert
            (log, com, gjf, out, inp, xyz, sdf, pdb, cdxml, cdx).
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
            filedir, fname = os.path.split(file)
            file_basename = os.path.splitext(fname)[0]
            output_path = os.path.join(
                filedir, f"{file_basename}.{output_filetype}"
            )
            self.convert_file(
                file,
                output_path,
                include_intermediate_structures=(
                    self.include_intermediate_structures
                ),
            )

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
                convert. When omitted, delegates to
                ``molecule.write_pdb_pybabel``. When the path is absent, XYZ
                is written to that path from *molecule* via ``write_xyz``,
                then converted (caller-supplied paths are never deleted).
            mode (str): File mode passed to ``write_xyz`` when the XYZ file
                must be created. Default ``'w'``.
            overwrite (bool): Whether to overwrite *pdb_filename* if it
                already exists. Default ``True``.
            cleanup (bool): Remove any auto-generated temporary XYZ file after
                the conversion (only when *xyz_filename* is ``None``).
                Default ``True``.

        Raises:
            ImportError: If Open Babel (``openbabel``) is not installed.
            ValueError: If the XYZ file cannot be parsed by Open Babel.
        """
        if xyz_filename is None:
            molecule.write_pdb_pybabel(
                pdb_filename,
                mode=mode,
                overwrite=overwrite,
                cleanup=cleanup,
            )
            return

        if not os.path.isfile(xyz_filename):
            logger.debug(
                f"XYZ {xyz_filename} missing; writing coordinates before "
                "conversion."
            )
            molecule.write_xyz(xyz_filename, mode=mode)

        try:
            from openbabel import pybel
        except ImportError as exc:  # pragma: no cover
            raise ImportError(
                "xyz_to_pdb requires Open Babel. "
                "Use 'conda install -c conda-forge openbabel' to install it."
            ) from exc

        xyz_mol = next(pybel.readfile("xyz", xyz_filename), None)
        if xyz_mol is None:
            raise ValueError(f"Unable to read molecule from {xyz_filename}")

        logger.info(
            f"Converting XYZ {xyz_filename} to PDB {pdb_filename} via "
            f"Open Babel (overwrite={overwrite})"
        )
        xyz_mol.write("pdb", pdb_filename, overwrite=overwrite)

    @staticmethod
    def resolve_output_path(input_path, output=None, output_filetype="xyz"):
        """Resolve a single-file output path from input and output options.

        When *output* is ``None``, returns ``{input_dir}/{input_stem}.{output_filetype}``.
        When *output* is extension-only (no directory separators and either
        starts with ``.`` or contains no ``.``), returns
        ``{input_dir}/{input_stem}.{ext}``.
        Otherwise *output* is treated as a full output path.

        Args:
            input_path (str): Path to the input file.
            output (str | None): Output path or extension shorthand.
            output_filetype (str): Default output extension when *output* is
                ``None``. Defaults to ``xyz``.

        Returns:
            str: Resolved output file path.
        """
        filedir, fname = os.path.split(input_path)
        basename = os.path.splitext(fname)[0]

        if output is None:
            return os.path.join(filedir, f"{basename}.{output_filetype}")

        if os.path.sep not in output and (
            output.startswith(".") or "." not in output
        ):
            ext = output.lstrip(".")
            return os.path.join(filedir, f"{basename}.{ext}")

        return output

    @staticmethod
    def convert_file(
        input_path, output_path, include_intermediate_structures=False
    ):
        """
        Convert a single input file to a single output file via ``Molecule``.

        The input format is inferred from the file extension by
        ``Molecule.from_filepath`` and the output format is inferred from
        *output_path*.

        When *include_intermediate_structures* is ``False`` (default), only
        a single structure is written to *output_path*. When ``True`` and
        the input contains multiple structures, each structure is written
        to a separate file with a ``_1``, ``_2``, ... suffix.

        Args:
            input_path (str): Path to the input file.
            output_path (str): Path to the output file. Used to infer the
                output format from its extension.
            include_intermediate_structures (bool): If ``True``, convert all
                structures in multi-structure inputs to separate files.
                Default ``False``.

        Raises:
            FileNotFoundError: If *input_path* does not exist.
            ValueError: If the output format is not supported by
                ``Molecule.write``.
        """
        input_path = os.path.abspath(input_path)
        output_path = os.path.abspath(output_path)

        if not os.path.exists(input_path):
            raise FileNotFoundError(f"Input file not found: {input_path}")

        logger.info(f"Converting {input_path} -> {output_path}")

        output_dir, output_name = os.path.split(output_path)
        if output_dir and not os.path.isdir(output_dir):
            os.makedirs(output_dir, exist_ok=True)

        output_basename, output_ext = os.path.splitext(output_name)
        output_ext = output_ext.lstrip(".").lower()
        if not output_ext:
            raise ValueError(
                f"Output path must include a file extension: {output_path}"
            )

        if include_intermediate_structures:
            molecules = Molecule.from_filepath(
                input_path, index=":", return_list=True
            )
            if molecules is None:
                raise ValueError(
                    f"No molecule could be read from {input_path}"
                )
            if not isinstance(molecules, list):
                molecules = [molecules]
            if len(molecules) == 0:
                raise ValueError(f"No molecules found in {input_path}")

            if len(molecules) == 1:
                molecules[0].write(output_path, format=output_ext)
                logger.info(f"Created: {output_path}")
            else:
                for i, molecule in enumerate(molecules, start=1):
                    path = os.path.join(
                        output_dir, f"{output_basename}_{i}.{output_ext}"
                    )
                    if os.path.exists(path):
                        logger.warning(
                            "Overwriting existing numbered output file: %s",
                            path,
                        )
                    molecule.write(path, format=output_ext)
                    logger.info(f"Created: {path}")
        else:
            molecule = Molecule.from_filepath(input_path)
            if molecule is None:
                raise ValueError(
                    f"No molecule could be read from {input_path}"
                )
            molecule.write(output_path, format=output_ext)
            logger.info(f"Created: {output_path}")
