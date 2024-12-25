#!/usr/bin/env python
import logging
import os
import click
from chemsmart.utils.utils import create_logger
from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.io.orca.outputs import ORCAOutput
from chemsmart.utils.utils import (
    get_value_by_number,
    get_key_by_value_and_number,
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
    help="Atom numbers from which to obtain Hirshfeld Charges and Spins. 1-indexed.",
)
def entry_point(filename, numbers):
    create_logger()
    if filename.endswith(".log"):
        outputfile = Gaussian16Output(filename=filename)
    elif filename.endswith(".out"):
        outputfile = ORCAOutput(filename=filename)
    else:
        raise TypeError(f"File {filename} is of unknown filetype.")
    hirshfeld_charges = outputfile.hirshfeld_charges
    logger.info("\nHirshfeld Charges:")
    for hkey, hvalue in hirshfeld_charges.items():
        logger.info(f"{hkey:<6}  :  {hvalue:>8.3f}")
    logger.info("\n")

    hirshfeld_spins = outputfile.hirshfeld_spin_densities
    logger.info("\nHirshfeld Spins:")
    for hkey, hvalue in hirshfeld_spins.items():
        logger.info(f"{hkey:<6}  :  {hvalue:>8.3f}")
    logger.info("\n")

    if numbers is not None:
        for n in numbers:
            charge_value = get_value_by_number(n, hirshfeld_charges)
            hk = get_key_by_value_and_number(
                charge_value, n, hirshfeld_charges
            )
            logger.info(f"Hirshfeld Charge at {hk} is {charge_value:.3f}.")
        logger.info("\n")

    if numbers is not None:
        for n in numbers:
            spin_value = get_value_by_number(n, hirshfeld_spins)
            hk = get_key_by_value_and_number(spin_value, n, hirshfeld_spins)
            logger.info(f"Hirshfeld Spin at {hk} is {spin_value:.3f}.")
        logger.info("\n")


if __name__ == "__main__":
    entry_point()
