#!/usr/bin/env python
"""
Molecular structure filtering script.

This script filters and groups molecular structures based on various
criteria such as similarity, energy, or structural features, helping
to organize and analyze large sets of molecular conformations.

The script supports automatic file type detection based on file content
using Gaussian or ORCA specifications, eliminating the need for manual
file type specification.
"""

import glob
import logging
import os

import click

from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.grouper import StructureGrouperFactory
from chemsmart.utils.io import get_program_type_from_file
from chemsmart.utils.logger import create_logger
from chemsmart.utils.utils import extract_number

logger = logging.getLogger(__name__)
os.environ["OMP_NUM_THREADS"] = "1"


@click.command()
@click.option(
    "-d",
    "--directory",
    required=True,
    type=str,
    help="Directory containing files to filter.",
)
@click.option(
    "-t",
    "--type",
    type=str,
    default=None,
    help="Type of files to filter (e.g., 'log', 'out'). "
    "If not provided, the script will attempt to auto-detect the file type "
    "based on file content for Gaussian and ORCA output files.",
)
@click.option(
    "-g",
    "--grouping-strategy",
    type=click.Choice(
        [
            "rmsd",
            "fingerprint",
            "isomorphism",
            "formula",
            "connectivity",
        ],
        case_sensitive=False,
    ),
    default="rmsd",
    help="Grouping strategy to use for grouping. \n"
    "Available options are 'rmsd', 'fingerprint', 'isomorphism', "
    "'formula', 'connectivity'",
)
@click.option(
    "-v",
    "--value",
    default=None,
    help="Threshold for grouping strategies. "
    "For RMSD, it is the rmsd_threshold value. "
    "For Tanimoto, it is the similarity_threshold value. "
    "For connectivity, it is the bond_cutoff_buffer value.",
)
@click.option(
    "-n",
    "--num-grouper-processors",
    type=int,
    default=4,
    help="Number of processors for grouping.",
)
def entry_point(
    directory, type, grouping_strategy, num_grouper_processors, value, **kwargs
):
    """
    Filter molecular structures to remove duplicates.

    This script automatically detects file types based on content if not
    explicitly provided. Supports Gaussian (.log, .out) and ORCA (.out)
    output files.
    """
    create_logger()
    directory = os.path.abspath(directory)

    # Auto-detect file type if not provided
    if type is None:
        logger.info(
            "File type not provided. Attempting to auto-detect from directory contents..."
        )
        try:
            # Get all potential output files
            candidate_files = []
            for ext in ["*.log", "*.out"]:
                candidate_files.extend(glob.glob(f"{directory}/{ext}"))

            if not candidate_files:
                raise FileNotFoundError(
                    f"No .log or .out files found in directory: {directory}"
                )

            # Try to detect file type from the first file
            first_file = candidate_files[0]
            program = get_program_type_from_file(first_file)

            if program == "gaussian":
                # Determine file extension based on the first file
                type = os.path.splitext(first_file)[1][
                    1:
                ]  # Remove leading dot
                logger.info(
                    f"Auto-detected Gaussian output files with extension '.{type}'"
                )
            elif program == "orca":
                type = "out"
                logger.info(
                    "Auto-detected ORCA output files with extension '.out'"
                )
            elif program == "unknown":
                raise ValueError(
                    f"Could not determine file type for: {first_file}. "
                    "The file does not appear to be a supported Gaussian or ORCA output file. "
                    "Please specify the file type manually using the --type option."
                )
            else:
                # Program detected but not Gaussian or ORCA
                raise ValueError(
                    f"Detected {program} output file, but only Gaussian and ORCA "
                    "are supported for structure filtering. "
                    "Please specify the file type manually using the --type option."
                )

        except FileNotFoundError as e:
            logger.error(str(e))
            raise
        except ValueError as e:
            logger.error(str(e))
            raise
        except Exception as e:
            logger.error(
                f"Error during file type auto-detection: {e}. "
                "Please specify the file type manually using the --type option."
            )
            raise

    # Obtain all structures from files in the directory
    filenames = glob.glob(f"{directory}/*.{type}")

    if not filenames:
        error_msg = f"No files with extension '.{type}' found in directory: {directory}"
        logger.error(error_msg)
        raise FileNotFoundError(error_msg)

    # Sort filenames based on numeric part
    sorted_filenames = sorted(filenames, key=extract_number)

    # remove last 7 characters, i.e., _c1.log
    base_filename = os.path.basename(sorted_filenames[0])[:-7]

    # Load molecular structures from files
    logger.info(f"Loading {len(sorted_filenames)} molecular structures...")
    try:
        molecules = [Molecule.from_filepath(file) for file in sorted_filenames]
    except Exception as e:
        logger.error(f"Error loading molecular structures: {e}")
        raise

    # Create grouper based on the specified grouping strategy
    try:
        grouper = StructureGrouperFactory.create(
            molecules, strategy=grouping_strategy, **kwargs
        )
        logger.info(
            f"Grouper created using {grouping_strategy} grouping strategy."
        )
        logger.info(f"Threshold for grouping strategy: {value}")
    except Exception as e:
        logger.error(f"Error creating grouper: {e}")
        raise e

    # Perform grouping analysis
    groups, group_indices = grouper.group()
    unique_structures = grouper.unique()
    logger.info(f"Identified {len(unique_structures)} unique structures.")

    # Write filtering results to output file
    output_file = os.path.join(
        directory, f"filter_{base_filename}_{type}files.txt"
    )

    f = open(output_file, "w")

    f.write(f"Initial number of structures to filter: {len(molecules)}\n")
    f.write(f"Final unique number of structures to filter: {len(molecules)}\n")

    # Convert to 1-indexed to be consistent with conformer naming
    group_indices_one = []
    for group_index in group_indices:
        group_index_one = []
        for idx in group_index:
            idx_one = int(idx + 1)
            group_index_one.append(idx_one)
        group_indices_one.append(group_index_one)
    f.write(f"1-indexed groups: {group_indices_one}\n")

    f.close()

    # Write unique structure filenames to separate file
    g = open(f"unique_structures_{base_filename}_{type}.txt", "w")
    logger.info("Writing unique structures to file.")
    for group_index in group_indices:
        unique_file = sorted_filenames[group_index[0]]
        unique_file_basename = os.path.basename(unique_file)
        g.write(f"{unique_file_basename}\n")
    g.close()


if __name__ == "__main__":
    entry_point()
