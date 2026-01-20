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
        chemsmart run mol -t log -d . -i : align

        # Align all structures from a single multi-structure file
        chemsmart run mol -f crest_conformers.xyz align

        # Align specific structures from a single file
        chemsmart run mol -f dft_opt.log -i 1,-1 align
        chemsmart run mol -f trajectory.xyz -i 1,5-10 align
        chemsmart run mol -f scan.log -i 1-5 align
        chemsmart run mol -f scan.log -i 1:3:10 align
    """

    index = ctx.obj["index"]
    label = ctx.obj["label"]
    filenames = ctx.obj["filenames"]
    directory = ctx.obj["directory"]
    filetype = ctx.obj["filetype"]

    # Validate input parameters first
    if not filenames and not directory:
        raise click.BadParameter(
            "No input files specified. Use -f/--file for files or -d/--directory for directory."
        )

    if directory and not filetype:
        raise click.BadParameter(
            "Directory specified but no filetype provided. Use -t/--type to specify file type (e.g., xyz, log, gjf)."
        )

    # Set default index values for align task
    if index is None:
        if directory is not None and filetype is not None:
            # Using directory - default to last structure (-1)
            index = "-1"
        elif filenames and len(filenames) > 1:
            # Multiple files - default to last structure (-1)
            index = "-1"
        elif filenames and len(filenames) == 1:
            # Single file - default to all structures (:)
            index = ":"

    molecules = []  # Initialize molecules list

    # Variable to store the base file for label generation
    base_file_for_label = None

    def process_files_with_per_file_indexing(
        files_list, add_index_suffix=True, check_exists=True
    ):
        """Helper function to process multiple files with per-file indexing."""
        for file_path in files_list:
            logger.debug(f"Processing file: {file_path}")

            # Load all structures from this file first
            file_molecules = load_molecules_from_paths(
                [file_path],
                index=":",  # Load all structures first
                add_index_suffix_for_single=add_index_suffix,
                check_exists=check_exists,
            )

            # Apply index selection to this file's structures
            if index and index != ":":
                from chemsmart.utils.utils import parse_index_specification

                try:
                    selected_indices = parse_index_specification(
                        index,
                        total_count=len(file_molecules),
                        allow_duplicates=False,  # Don't allow duplicates for alignment
                        allow_out_of_range=False,  # Don't allow out of range indices
                    )

                    # Handle different return types for this file
                    if isinstance(selected_indices, list):
                        molecules.extend(
                            [file_molecules[i] for i in selected_indices]
                        )
                    elif isinstance(selected_indices, int):
                        molecules.append(file_molecules[selected_indices])
                    elif isinstance(selected_indices, slice):
                        molecules.extend(file_molecules[selected_indices])
                    else:
                        raise ValueError(
                            f"Unexpected index type: {type(selected_indices)}"
                        )

                except ValueError as e:
                    raise click.BadParameter(
                        f"Error processing file {file_path}: {str(e)}"
                    )
            else:
                # Use all structures from this file
                molecules.extend(file_molecules)

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

            # Process all matched files with per-file indexing
            process_files_with_per_file_indexing(
                matched_files, add_index_suffix=False, check_exists=False
            )

            logger.debug(
                f"Loaded {len(molecules)} molecules from {len(matched_files)} files using filetype pattern with index={index}"
            )
            base_file_for_label = matched_files[0]
        else:
            # This should not happen due to validation above, but keep as safeguard
            raise click.BadParameter(
                "Directory specified but no filetype provided. Use -t/--type to specify file type."
            )

    elif filenames:
        logger.debug(f"Received filename parameter: {filenames}")
        logger.debug(f"Type of filename: {type(filenames)}")

        if isinstance(filenames, str):
            filenames = [filenames]

        if len(filenames) == 1:
            # Single file: first load all structures, then apply user's index
            all_molecules = load_molecules_from_paths(
                filenames,
                index=":",  # Load all structures first
                add_index_suffix_for_single=True,
                check_exists=True,
            )

            # Now apply user's index selection if specified
            if index and index != ":":
                from chemsmart.utils.utils import parse_index_specification

                try:
                    # Use the enhanced parse_index_specification function with validation
                    selected_indices = parse_index_specification(
                        index,
                        total_count=len(all_molecules),
                        allow_duplicates=False,  # Don't allow duplicates for alignment
                        allow_out_of_range=False,  # Don't allow out of range indices
                    )

                    # Handle different return types
                    if isinstance(selected_indices, list):
                        # List of 0-based indices - use directly
                        molecules.extend(
                            [all_molecules[i] for i in selected_indices]
                        )
                    elif isinstance(selected_indices, int):
                        # Single 0-based index
                        molecules.append(all_molecules[selected_indices])
                    elif isinstance(selected_indices, slice):
                        # Slice object - Python handles bounds automatically
                        molecules.extend(all_molecules[selected_indices])
                    else:
                        raise ValueError(
                            f"Unexpected index type: {type(selected_indices)}"
                        )

                except ValueError as e:
                    raise click.BadParameter(str(e))
            else:
                # Use all structures
                molecules.extend(all_molecules)
        else:
            # Multiple files: use the same per-file processing as directory mode
            process_files_with_per_file_indexing(
                filenames, add_index_suffix=True, check_exists=True
            )

        logger.debug(
            f"Loaded {len(molecules)} molecules from {len(filenames)} files using align-specific filenames with index={index}"
        )

        base_file_for_label = filenames[0]

    else:
        # This should not happen due to validation above, but keep as safeguard
        raise click.BadParameter(
            "No input files specified. This should have been caught earlier."
        )

    # Check if we have enough molecules for alignment
    if not isinstance(molecules, list):
        molecules = list(molecules) if molecules else []

    if len(molecules) < 2:
        error_msg = f"Need at least 2 molecules for alignment, but only loaded {len(molecules)} molecule(s)."
        if filenames and len(filenames) == 1 and index and index != ":":
            error_msg += f" Index '{index}' may not select enough structures."
        raise click.BadParameter(error_msg)

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
        n_structures = len(molecules)

        # Check if it's a single file case
        if filenames and len(filenames) == 1:
            # Single file with multiple structures: filename_n_structures_align
            align_label = f"{n_structures}_structures_in_{base_label}_align"
        else:
            # Multiple files: original naming scheme
            if n_structures > 2:
                align_label = (
                    f"{base_label}_and_{n_structures-1}_structures_align"
                )
            else:
                align_label = f"{base_label}_and_1_structure_align"

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
