import glob
import logging
import os

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.mol.mol import (
    click_pymol_visualization_options,
    mol,
)
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.io import load_molecules_from_paths

logger = logging.getLogger(__name__)


@mol.command("align", cls=MyCommand)
@click_job_options
@click_pymol_visualization_options
@click.pass_context
def align(
    ctx,
    file,
    style,
    trace,
    vdw,
    quiet,
    command_line_only,
    label_offset,
    skip_completed,
    **kwargs,
):
    """CLI subcommand for aligning multiple molecule files in PyMOL.
    
    Examples:
        # Align multiple files
        chemsmart run mol -f a.log -f b.xyz -f c.gjf align
        chemsmart run mol -t log -d . align
        
        # Align all structures from a single multi-structure file
        chemsmart run mol -f crest_conformers.xyz -i : align
        
        # Align specific structures from a single file
        chemsmart run mol -f dft_opt.log -i 1,-1 align
        chemsmart run mol -f trajectory.xyz -i 1,5,10 align
        chemsmart run mol -f scan.log -i 1-5 align
    """

    index = ctx.obj["index"]
    label = ctx.obj["label"]
    filenames = ctx.obj["filenames"]
    directory = ctx.obj["directory"]
    filetype = ctx.obj["filetype"]
    molecules = []  # Initialize molecules list

    # Variable to store the base file for label generation
    base_file_for_label = None

    if directory:
        directory = os.path.abspath(directory)
        logger.info(
            f"Obtaining files in directory: {directory} for alignment."
        )
        if filetype:
            filetype_pattern = os.path.join(directory, f"*.{filetype}")
            matched_files = glob.glob(filetype_pattern)
            if not matched_files:
                logger.warning(f"No files matched pattern: {filetype_pattern}")
                raise click.BadParameter(
                    f"No files found matching pattern: {filetype_pattern}"
                )

            molecules += load_molecules_from_paths(
                matched_files,
                index=index,
                add_index_suffix_for_single=False,
                check_exists=False,
            )
            logger.debug(
                f"Loaded {len(molecules)} molecules from {len(matched_files)} files using filetype pattern with index={index}"
            )

            base_file_for_label = matched_files[0]

    elif filenames:

        logger.debug(f"Received filename parameter: {filenames}")
        logger.debug(f"Type of filename: {type(filenames)}")

        if not filenames or (
            isinstance(filenames, (list, tuple)) and len(filenames) == 0
        ):
            raise click.BadParameter("No valid filenames provided")

        if isinstance(filenames, str):
            filenames = [filenames]

        # Load molecules from files using the index specification
        # For single file: enables chemsmart run mol -f file.xyz -i : align
        #                  or: chemsmart run mol -f file.log -i 1,-1 align
        # For multiple files: loads specified index from each file
        molecules += load_molecules_from_paths(
            filenames,
            index=index,
            add_index_suffix_for_single=True,
            check_exists=True,
        )

        logger.debug(
            f"Loaded {len(molecules)} molecules from {len(filenames)} files using align-specific filenames with index={index}"
        )

        base_file_for_label = filenames[0]

    if not isinstance(molecules, list) or len(molecules) < 2:
        raise click.BadParameter("Need at least two molecules for alignment")

    # Generate align-specific label if user didn't provide one
    if label is not None:
        # User provided label - use it directly
        align_label = label
    else:
        # Generate label based on first file and molecule count
        if base_file_for_label:
            base_label = os.path.splitext(
                os.path.basename(base_file_for_label)
            )[0]
        else:
            base_label = "molecules"

        # Generate align-specific naming
        n_molecules = len(molecules)
        if n_molecules > 2:
            align_label = f"{base_label}_and_{n_molecules-1}_molecules_align"
        else:
            align_label = f"{base_label}_and_1_molecule_align"

    from chemsmart.jobs.mol.align import PyMOLAlignJob

    return PyMOLAlignJob(
        molecule=molecules,
        label=align_label,
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
