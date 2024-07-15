#!/usr/bin/env python
import logging
import os
import click
from chemsmart.utils.utils import create_logger
from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.io.orca.outputs import ORCAOutput
from chemsmart.utils.utils import get_value_by_number, get_key_by_value

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
    "-u",
    "--unit",
    default="eV",
    type=click.Choice(["eV", "kcal/mol"], case_sensitive=False),
    help="Unit of FMO energy.",
)
@click.option(
    "-n",
    "--numbers",
    default=None,
    type=int,
    multiple=True,
    help="Atom numbers from which to obtain Hirshfeld Charges and Spins. 1-indexed.",
)
def entry_point(filename, unit, numbers):
    create_logger()
    if filename.endswith(".log"):
        outputfile = Gaussian16Output(filename=filename)
    elif filename.endswith(".out"):
        outputfile = ORCAOutput(filename=filename)
    else:
        raise TypeError(f"File {filename} is of unknown filetype.")
    hirshfeld_charges = outputfile.hirshfeld_charges
    hirshfeld_spins = outputfile.hirshfeld_spins
    logger.info(f"\nHirshfeld Charges:")
    for hkey, hvalue in hirshfeld_charges.items():
        logger.info(f"{hkey}  :  {hvalue}")
    logger.info(f"\n")
    logger.info(f"\nHirshfeld Spins:")
    for hkey, hvalue in hirshfeld_spins.items():
        logger.info(f"{hkey}  :  {hvalue}")
    logger.info(f"\n")

    hkeys = hirshfeld_charges.keys()
    if numbers is not None:
        # print(numbers)
        # print(type(numbers))
        for n in numbers:
            charge_value = get_value_by_number(n, hirshfeld_charges)
            hk = get_key_by_value(charge_value, hirshfeld_charges)
            logger.info(f"Hirshfeld Charge at {hk} is {charge_value}.")
        logger.info(f"\n")

        for n in numbers:
            spin_value = get_value_by_number(n, hirshfeld_spins)
            hk = get_key_by_value(spin_value, hirshfeld_spins)
            logger.info(f"Hirshfeld Spin at {hk} is {spin_value}.")
        logger.info(f"\n")


if __name__ == "__main__":
    entry_point()
