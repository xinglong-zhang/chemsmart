#!/usr/bin/env python
"""
Wiberg Bond Index (WBI) analysis script.

This script extracts and analyzes Wiberg Bond Indices from Gaussian
output files, providing bond order analysis and molecular bonding
characterization for chemical systems.
"""

import logging
import os

import click

from chemsmart.io.gaussian.output import Gaussian16WBIOutput
from chemsmart.utils.logger import create_logger
from chemsmart.utils.utils import get_value_by_number

logger = logging.getLogger(__name__)
os.environ["OMP_NUM_THREADS"] = "1"


@click.command()
@click.option(
    "-f",
    "--filename",
    required=True,
    default=None,
    type=str,
    help="Gaussian output file for Wiberg Bond Index analysis.",
)
@click.option(
    "-n",
    "--numbers",
    default=None,
    type=int,
    multiple=True,
    help="Atom numbers from which to obtain data from Wiberg Bond Index "
    "output file. 1-indexed.",
)
def entry_point(filename, numbers):
    """
    Analyze Wiberg Bond Index and Natural Population data.
    """
    create_logger()
    outputfile = Gaussian16WBIOutput(filename=filename)

    # Extract Natural Population Analysis data (alternative method commented)
    # npa = outputfile.natural_population_analysis
    # logger.info("\nNatural Population Analysis:")
    # for npkey, npvalue in npa.items():
    #     natural_charge = npvalue["natural_charge"]
    #     logger.info(f"{npkey}  :  {natural_charge:.3f}")
    # logger.info("\n")

    # Extract and display natural charges
    natural_charges = outputfile.natural_charges
    logger.info("\nNatural Charges:")
    for nkey, nvalue in natural_charges.items():
        logger.info(f"{nkey:<6}  :  {nvalue:>8.3f}")
    logger.info("\n")

    # Display specific atom charges if requested
    if numbers is not None:
        for n in numbers:
            # Get natural population analysis value corresponding to the atom number
            charge_value = get_value_by_number(n, natural_charges)
            logger.info(f"Natural Charge at atom {n} is {charge_value:.>3f}.")
            logger.info("\n")
        logger.info("\n")


if __name__ == "__main__":
    entry_point()
