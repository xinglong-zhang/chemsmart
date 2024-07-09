#!/usr/bin/env python
import logging
import os
import click
from chemsmart.utils.utils import create_logger
from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.io.orca.outputs import ORCAOutput

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
def entry_point(filename, unit):
    create_logger()
    if filename.endswith(".log"):
        outputfile = Gaussian16Output(filename=filename)
    elif filename.endswith(".out"):
        outputfile = ORCAOutput(filename=filename)
    else:
        raise TypeError(f"File {filename} is of unknown filetype.")
    homo_energy = outputfile.homo_energy
    lumo_energy = outputfile.lumo_energy
    fmo_gap = outputfile.fmo_gap
    energy_unit = "eV"
    if unit.lower() == "kcal/mol":
        homo_energy *= 23.06054195
        lumo_energy *= 23.06054195
        fmo_gap *= 23.06054195
        energy_unit = "kcal/mol"
    logger.info(f"HOMO energy: {homo_energy:.4f} {energy_unit}")
    logger.info(f"LUMO energy: {lumo_energy:.4f} {energy_unit}")
    logger.info(f"HOMO-LUMO gap: {fmo_gap:.4f} {energy_unit}")


if __name__ == "__main__":
    entry_point()
