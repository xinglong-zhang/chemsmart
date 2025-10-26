import glob
import logging
import os
from pathlib import Path

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.mol.mol import (
    click_pymol_align_options,
    click_pymol_visualization_options,
    mol,
)
from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@mol.command("align", cls=MyCommand)
@click_job_options
@click_pymol_visualization_options
@click_pymol_align_options
@click.pass_context
def align(
    ctx,
    filenames,
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
        chemsmart run mol align -F a.log -f b.xyz -F c.gjf
        chemsmart run mol align -T log
    """

    index = ctx.obj["index"]
    label = ctx.obj["label"]

    if filetype:

        filetype_pattern = f"*.{filetype}"
        matched_files = glob.glob(filetype_pattern)
        if not matched_files:
            logger.warning(f"No files matched pattern: {filetype_pattern}")
            raise click.BadParameter(
                f"No files found matching pattern: {filetype_pattern}"
            )

        molecules = []
        for file_path in matched_files:
            # Pass index to per-file reader in user string form so each file
            # yields the structure(s) corresponding to that index.
            mols = Molecule.from_filepath(
                filepath=file_path,
                index=index,
                return_list=True,
            )
            # assign unique names per-structure when file contains multiple structures
            base = os.path.splitext(os.path.basename(file_path))[0]
            if isinstance(mols, list) and len(mols) > 1:
                for j, mol in enumerate(mols, start=1):
                    mol.name = f"{base}_{j}"
            else:
                for mol in mols:
                    mol.name = base
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

            if not file_path:
                logger.warning(f"Skipping invalid file path: {file_path}")
                continue

            p = Path(file_path)
            if not p.is_file():
                logger.error(
                    f"File not found or not a regular file: {file_path}"
                )
                raise FileNotFoundError(f"File not found: {file_path}")

            file_path = str(p)

            try:
                mols = Molecule.from_filepath(
                    filepath=file_path,
                    index=index,
                    return_list=True,
                )
                # assign unique names per-structure when file contains multiple structures
                base = os.path.splitext(os.path.basename(file_path))[0]
                if isinstance(mols, list) and len(mols) > 1:
                    for j, mol in enumerate(mols, start=1):
                        mol.name = f"{base}_{j}"
                else:
                    for mol in mols:
                        mol.name = base
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
