#!/usr/bin/env python
"""
Molecular structure filtering script.

This script filters and groups molecular structures based on various
criteria such as similarity, energy, or structural features, helping
to organize and analyze large sets of molecular conformations.
"""

import glob
import logging
import os

import click

from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.grouper import StructureGrouperFactory
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
    "-t", "--type", type=str, default=None, 
    help="Type of files to filter (e.g., 'log', 'out')."
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
    """
    create_logger()
    directory = os.path.abspath(directory)
    
    # Set default file type if not provided
    if type is None:
        logger.info(
            "Type of files is not provided! Assuming .log file type "
            "for filtering."
        )
        type = "log"

    # Obtain all structures from files in the directory
    filenames = glob.glob(f"{directory}/*.{type}")

    # Sort filenames based on numeric part
    sorted_filenames = sorted(filenames, key=extract_number)

    # remove last 7 characters, i.e., _c1.log
    base_filename = os.path.basename(sorted_filenames[0])[:-7]

    # Load molecular structures from files
    molecules = [Molecule.from_filepath(file) for file in sorted_filenames]

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
