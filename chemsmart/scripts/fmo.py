#!/usr/bin/env python
"""
Frontier Molecular Orbital (FMO) analysis script.

This script extracts and analyzes Frontier Molecular Orbital data from
Gaussian and ORCA output files.
"""

import logging
import os

import click

from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.io.orca.output import ORCAOutput
from chemsmart.utils.io import get_program_type_from_file
from chemsmart.utils.logger import create_logger

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
    """
    Calculate and display frontier molecular orbital (FMO) properties.
    """
    create_logger()
    program = get_program_type_from_file(filename)
    if program == "gaussian":
        outputfile = Gaussian16Output(filename=filename)
    elif program == "orca":
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

    # obtain chemical potential, μ = 1/2 * (lumo_energy + homo_energy)
    # Chemical hardness, η = 1/2 * (lumo_energy - homo_energy)
    # Electrophilicity index = ω = μ^2/2η

    chemical_potential = 1 / 2 * (lumo_energy + homo_energy)
    chemical_hardness = 1 / 2 * (lumo_energy - homo_energy)
    electrophilicity_index = chemical_potential**2 / (2 * chemical_hardness)

    # print the results
    logger.info(f"HOMO energy: {homo_energy:.4f} {energy_unit}")
    logger.info(f"LUMO energy: {lumo_energy:.4f} {energy_unit}")
    logger.info(f"HOMO-LUMO gap: {fmo_gap:.4f} {energy_unit}")
    logger.info(
        f"Chemical potential, mu: {chemical_potential:.4} {energy_unit}"
    )
    logger.info(
        f"Chemical hardness, eta: {chemical_hardness:.4} {energy_unit}"
    )
    logger.info(
        f"Electrophilicity Index, omega: {electrophilicity_index:.4} "
        f"{energy_unit}"
    )


if __name__ == "__main__":
    entry_point()
