#!/usr/bin/env python
"""XYZ file writer script for molecular structures.

This script converts molecular structure data from quantum chemistry
calculation output files to standard XYZ coordinate format files.
"""

import logging
import os

import click

from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.logger import create_logger

logger = logging.getLogger(__name__)
os.environ["OMP_NUM_THREADS"] = "1"


@click.command()
@click.option(
    "-f", "--filename", default=None, 
    help="Input filename to be converted to XYZ format."
)
@click.option(
    "-i",
    "--index",
    default="-1",
    help="Index of structure to be written to file (1-indexed).",
)
@click.option(
    "-s/",
    "--single-file/--no-single-files",
    type=bool,
    is_flag=True,
    default=True,
    help="Write all structures to a single .xyz file if more than one "
         "structure is present.\nDefault is to write all structures to "
         "a single file.",
)
def entry_point(filename, index, single_file):
    """Convert molecular structures to XYZ coordinate files.
    
    This function reads molecular structure data from quantum chemistry
    calculation output files and writes the coordinates in standard XYZ
    format. It can handle single or multiple structures per file.
    
    The script can write a single structure or multiple structures based
    on 1-indexing. The default is to write the last structure in the file.
    
    Args:
        filename: Path to input file containing molecular structures
        index: Structure index to extract (default: -1 for last structure)
        single_file: Write all structures to one file vs separate files
    """
    create_logger()

    # Load molecular structures from the input file
    try:
        molecules = Molecule.from_filepath(
            filepath=filename, index=index, return_list=True
        )
    except FileNotFoundError as err:
        logger.error(err)
        return

    # Extract base filename for output naming
    file_basename = os.path.splitext(filename)[0]

    # Write structures to XYZ files
    if len(molecules) == 1:
        molecules[0].write_xyz(file_basename + "_single.xyz")
    else:
        # Handle multiple structures
        for i, molecule in enumerate(molecules):
            if single_file:
                molecule.write_xyz(file_basename + "_all.xyz", mode="a")
            else:
                molecule.write_xyz(file_basename + f"_{i+1}.xyz")


if __name__ == "__main__":
    entry_point()
