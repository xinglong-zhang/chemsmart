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

                def validate_and_add_index(idx, desc=""):
                    """Helper function to validate index and check for overlaps"""
                    if idx < 0 or idx >= len(all_molecules):
                        raise click.BadParameter(
                            f"Index{desc} is out of range. File has {len(all_molecules)} structures."
                        )
                    if idx in selected_positions:
                        raise click.BadParameter(
                            f"Index overlap detected{desc}."
                        )
                    selected_positions.add(idx)
                    selected_molecules.append(all_molecules[idx])

                selected_molecules = []
                selected_positions = set()  # Track actual 0-based positions

                if "-" in index:
                    # Handle indices containing negative numbers (both mixed and pure negative)
                    for part in [p.strip() for p in index.split(",")]:
                        if part.startswith("-"):
                            neg_idx = int(part)
                            if abs(neg_idx) > len(all_molecules):
                                raise click.BadParameter(
                                    f"Index '{neg_idx}' is out of range. File has {len(all_molecules)} structures."
                                )
                            actual_pos = (
                                len(all_molecules) + neg_idx
                            )  # Convert to 0-based
                            validate_and_add_index(actual_pos, f" '{neg_idx}'")
                        elif "-" in part:
                            # Range like '1-3'
                            start, end = map(int, part.split("-"))
                            for i in range(start, end + 1):
                                if i < 1 or i > len(all_molecules):
                                    raise click.BadParameter(
                                        f"Index '{i}' is out of range. File has {len(all_molecules)} structures (1-based indexing)."
                                    )
                                validate_and_add_index(
                                    i - 1, f" '{i}'"
                                )  # Convert to 0-based
                        else:
                            # Single positive index
                            pos_idx = int(part)
                            if pos_idx < 1 or pos_idx > len(all_molecules):
                                raise click.BadParameter(
                                    f"Index '{pos_idx}' is out of range. File has {len(all_molecules)} structures (1-based indexing)."
                                )
                            validate_and_add_index(
                                pos_idx - 1, f" '{pos_idx}'"
                            )  # Convert to 0-based
                    molecules += selected_molecules
                else:
                    # Standard positive indices only
                    from chemsmart.utils.utils import (
                        get_list_from_string_range,
                    )

                    selected_indices = get_list_from_string_range(index)

                    # Check for duplicates and out-of-range
                    if len(selected_indices) != len(set(selected_indices)):
                        duplicates = [
                            idx
                            for idx in set(selected_indices)
                            if selected_indices.count(idx) > 1
                        ]
                        raise click.BadParameter(
                            f"Index overlap detected. Indices {duplicates} are specified multiple times."
                        )

                    for i in selected_indices:
                        if i < 1 or i > len(all_molecules):
                            raise click.BadParameter(
                                f"Index '{i}' is out of range. File has {len(all_molecules)} structures (1-based indexing)."
                            )

                    molecules += [
                        all_molecules[i - 1] for i in selected_indices
                    ]
            else:
                # Use all structures
                molecules += all_molecules
        else:
            # Multiple files: use original logic
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
