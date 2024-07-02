#!/usr/bin/env python
import logging
import os
import click
from chemsmart.utils.utils import create_logger
from chemsmart.io.gaussian.output import Gaussian16Output

logger = logging.getLogger(__name__)
os.environ["OMP_NUM_THREADS"] = "1"


@click.command()
@click.option(
    "-f",
    "--filename",
    required=True,
    default=None,
    type=str,
    help="Gaussian output file.",
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
    gaussian_output = Gaussian16Output(filename=filename)
    homo_energy = gaussian_output.homo_energy
    lumo_energy = gaussian_output.lumo_energy
    homo_lumo_gap = homo_energy - lumo_energy
    energy_unit = "eV"
    if unit.lower() == "kcal/mol":
        homo_energy *= 23.06054195
        lumo_energy *= 23.06054195
        homo_lumo_gap *= 23.06054195
        energy_unit = "kcal/mol"
    logger.info(f"HOMO energy: {gaussian_output.homo_energy:.4f} {energy_unit}")
    logger.info(f"LUMO energy: {gaussian_output.lumo_energy:.4f} {energy_unit}")
    logger.info(f"HOMO-LUMO gap: {homo_lumo_gap:.4f} {energy_unit}")


if __name__ == "__main__":
    entry_point()
