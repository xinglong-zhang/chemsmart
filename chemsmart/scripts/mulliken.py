#!/usr/bin/env python
import logging
import os

import click

from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.io.orca.output import ORCAOutput
from chemsmart.utils.logger import create_logger
from chemsmart.utils.utils import (
    get_key_by_value_and_number,
    get_value_by_number,
)

logger = logging.getLogger(__name__)
os.environ["OMP_NUM_THREADS"] = "1"


@click.command()
@click.option(
    "-f",
    "--filename",
    required=True,
    default=None,
    type=str,
    help="Gaussian or ORCA output file.",
)
@click.option(
    "-n",
    "--numbers",
    default=None,
    type=int,
    multiple=True,
    help="Atom numbers from which to obtain Mulliken Charges and Spins. "
         "1-indexed.",
)
def entry_point(filename, numbers):
    """
    Extract and display Mulliken charges and spin densities.
    """
    create_logger()
    
    # Parse output file based on extension
    if filename.endswith(".log"):
        outputfile = Gaussian16Output(filename=filename)
    elif filename.endswith(".out"):
        outputfile = ORCAOutput(filename=filename)
    else:
        raise TypeError(f"File {filename} is of unknown filetype.")
        
    # Extract and display Mulliken charges
    mulliken_charges = outputfile.mulliken_atomic_charges
    logger.info("\nMulliken Charges:")
    for hkey, hvalue in mulliken_charges.items():
        logger.info(f"{hkey:<6}  :  {hvalue:>8.3f}")
    logger.info("\n")

    # Extract and display Mulliken spin densities
    mulliken_spins = outputfile.mulliken_spin_densities
    logger.info("\nMulliken Spin densities:")
    for hkey, hvalue in mulliken_spins.items():
        logger.info(f"{hkey:<6}  :  {hvalue:>8.3f}")
    logger.info("\n")

    # Display specific atom charges if requested
    if numbers is not None:
        for n in numbers:
            charge_value = get_value_by_number(n, mulliken_charges)
            hk = get_key_by_value_and_number(
                charge_value, n, mulliken_charges
            )
            logger.info(f"Mulliken Charge at {hk} is {charge_value:.3f}.")
        logger.info("\n")

    if numbers is not None:
        for n in numbers:
            spin_value = get_value_by_number(n, mulliken_spins)
            hk = get_key_by_value_and_number(spin_value, n, mulliken_spins)
            logger.info(
                f"Mulliken Spin densities at {hk} is {mulliken_spins:.3f}."
            )
        logger.info("\n")


if __name__ == "__main__":
    entry_point()
