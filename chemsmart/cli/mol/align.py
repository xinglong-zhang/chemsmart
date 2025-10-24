import glob
import logging
import os

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.mol.mol import (
    click_file_options,
    click_pymol_align_options,
    click_pymol_visualization_options,
    mol,
)
from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@mol.command("align", cls=MyCommand)
@click_job_options
@click_file_options
@click_pymol_visualization_options
@click_pymol_align_options
def align(
    filenames,
    label,
    append_label,
    index,
    file,
    style,
    trace,
    vdw,
    quiet,
    command_line_only,
    label_offset,
    filetype,
    skip_completed,
    **kwargs,
):
    """CLI for PyMOL alignment of multiple molecule files.
    Example:
        chemsmart run mol align -f a.log -f b.xyz -f c.gjf
        chemsmart run mol align -t log
    """

    if filetype:

        filetype_pattern = "*." + filetype
        matched_files = glob.glob(filetype_pattern)
        if not matched_files:
            logger.warning(f"No files matched pattern: {filetype_pattern}")
            raise click.BadParameter(
                f"No files found matching pattern: {filetype_pattern}"
            )

        molecules = []
        for file_path in matched_files:
            mols = Molecule.from_filepath(
                filepath=file_path,
                index=index,
                return_list=True,
            )
            for mol in mols:
                mol.name = os.path.splitext(os.path.basename(file_path))[0]
            molecules += mols
        logger.debug(
            f"Loaded {len(molecules)} molecules from {len(matched_files)} files using filetype pattern with index={index}"
        )

        base_file = matched_files[0]
        if label is None:
            label = os.path.splitext(os.path.basename(base_file))[0]

    elif filenames:

        logger.debug(f"Received filename parameter: {filenames}")
        logger.debug(f"Type of filename: {type(filenames)}")

        if not filenames or (
            isinstance(filenames, (list, tuple)) and len(filenames) == 0
        ):
            raise click.BadParameter("No valid filenames provided")

        if isinstance(filenames, str):
            filenames = [filenames]

        molecules = []
        for i, file_path in enumerate(filenames):
            logger.debug(
                f"Processing file {i+1}/{len(filenames)}: {file_path}"
            )

            if not file_path or not isinstance(file_path, str):
                logger.warning(f"Skipping invalid file path: {file_path}")
                continue

            if not os.path.exists(file_path):
                logger.error(f"File not found: {file_path}")
                raise FileNotFoundError(f"File not found: {file_path}")

            try:
                mols = Molecule.from_filepath(
                    filepath=file_path,
                    index=index,
                    return_list=True,
                )
                for mol in mols:
                    mol.name = os.path.splitext(os.path.basename(file_path))[0]
                molecules += mols
                logger.debug(
                    f"Successfully loaded {len(mols)} molecules from {file_path}"
                )
            except Exception as e:
                logger.error(f"Error loading molecules from {file_path}: {e}")
                raise

        logger.debug(
            f"Loaded {len(molecules)} molecules from {len(filenames)} files using align-specific filenames with index={index}"
        )

        if filenames:
            base_file = filenames[0]
            if label is None:
                label = os.path.splitext(os.path.basename(base_file))[0]

    if append_label is not None:
        if label:
            label = f"{label}_{append_label}"
        else:
            raise ValueError("Cannot append label without a base label")

    if not isinstance(molecules, list) or len(molecules) < 2:
        raise click.BadParameter("Need at least two molecules for alignment")

    from chemsmart.jobs.mol.align import PyMOLAlignJob

    return PyMOLAlignJob(
        molecule=molecules,
        label=label,
        pymol_script=file,
        style=style,
        trace=trace,
        vdw=vdw,
        quiet_mode=quiet,
        command_line_only=command_line_only,
        label_offset=label_offset,
        skip_completed=skip_completed,
        **kwargs,
    )
