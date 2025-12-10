#!/usr/bin/env python
"""
Frontier Molecular Orbital (FMO) analysis script.

This script extracts and analyzes Frontier Molecular Orbital data from
Gaussian and ORCA output files.
"""

import logging
import os

import click

from chemsmart.cli.logger import logger_options
from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.io.orca.output import ORCAOutput
from chemsmart.utils.io import get_outfile_format
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
@logger_options
def entry_point(filename, unit, debug, stream):
    """
    Calculate and display frontier molecular orbital (FMO) properties.

    Supports both closed-shell and open-shell systems.
    """
    create_logger(debug=debug, stream=stream)
    program = get_outfile_format(filename)
    if program == "gaussian":
        outputfile = Gaussian16Output(filename=filename)
    elif program == "orca":
        outputfile = ORCAOutput(filename=filename)
    else:
        raise TypeError(f"File {filename} is of unknown filetype.")

    multiplicity = outputfile.multiplicity
    if unit.lower() == "kcal/mol":
        conversion_factor = 23.06054195
        energy_unit = "kcal/mol"
    else:
        conversion_factor = 1.0
        energy_unit = "eV"

    logger.info(f"System multiplicity: {multiplicity}")
    logger.info("")

    if multiplicity == 1:
        # Closed-shell system
        logger.info("=" * 60)
        logger.info("CLOSED-SHELL SYSTEM FMO ANALYSIS")
        logger.info("=" * 60)
        logger.info("")

        homo_energy = outputfile.homo_energy
        lumo_energy = outputfile.lumo_energy
        fmo_gap = outputfile.fmo_gap

        if homo_energy is None or lumo_energy is None:
            logger.error(
                "Could not extract HOMO/LUMO energies from output file."
            )
            return

        # Apply unit conversion
        homo_energy *= conversion_factor
        lumo_energy *= conversion_factor
        fmo_gap *= conversion_factor

        # obtain chemical potential, μ = 1/2 * (lumo_energy + homo_energy)
        # Chemical hardness, η = 1/2 * (lumo_energy - homo_energy)
        # Electrophilicity index = ω = μ^2/2η
        chemical_potential = 1 / 2 * (lumo_energy + homo_energy)
        chemical_hardness = 1 / 2 * (lumo_energy - homo_energy)
        electrophilicity_index = chemical_potential**2 / (
            2 * chemical_hardness
        )

        # print the results
        logger.info(f"HOMO energy: {homo_energy:.4f} {energy_unit}")
        logger.info(f"LUMO energy: {lumo_energy:.4f} {energy_unit}")
        logger.info(f"HOMO-LUMO gap: {fmo_gap:.4f} {energy_unit}")
        logger.info(
            f"Chemical potential, μ: {chemical_potential:.4f} {energy_unit}"
        )
        logger.info(
            f"Chemical hardness, η: {chemical_hardness:.4f} {energy_unit}"
        )
        logger.info(
            f"Electrophilicity index, ω: {electrophilicity_index:.4f} {energy_unit}"
        )
        logger.info("")

    else:
        # Open-shell system
        logger.info("=" * 60)
        logger.info("OPEN-SHELL SYSTEM FMO ANALYSIS")
        logger.info("=" * 60)
        logger.info("")

        alpha_homo = outputfile.alpha_homo_energy
        alpha_lumo = outputfile.alpha_lumo_energy
        beta_homo = outputfile.beta_homo_energy
        beta_lumo = outputfile.beta_lumo_energy
        somo_energies = outputfile.somo_energies
        highest_somo = outputfile.highest_somo_energy
        lowest_somo = outputfile.lowest_somo_energy
        alpha_fmo_gap = outputfile.alpha_fmo_gap
        beta_fmo_gap = outputfile.beta_fmo_gap
        fmo_gap = outputfile.fmo_gap

        # Alpha channel analysis
        if alpha_homo is not None and alpha_lumo is not None:
            logger.info("--- Alpha Spin Channel ---")
            alpha_homo *= conversion_factor
            alpha_lumo *= conversion_factor
            alpha_fmo_gap *= conversion_factor

            logger.info(f"α-HOMO energy: {alpha_homo:.4f} {energy_unit}")
            logger.info(f"α-LUMO energy: {alpha_lumo:.4f} {energy_unit}")
            logger.info(f"α-HOMO-LUMO gap: {alpha_fmo_gap:.4f} {energy_unit}")

            # Reactivity descriptors of alpha channel for open-shell systems
            chemical_potential_alpha = 1 / 2 * (alpha_lumo + alpha_homo)
            chemical_hardness_alpha = 1 / 2 * (alpha_lumo - alpha_homo)
            electrophilicity_index_alpha = chemical_potential_alpha**2 / (
                2 * chemical_hardness_alpha
            )

            logger.info(
                f"Chemical potential, μ_α: {chemical_potential_alpha:.4f} {energy_unit}"
            )
            logger.info(
                f"Chemical hardness, η_α: {chemical_hardness_alpha:.4f} {energy_unit}"
            )
            logger.info(
                f"Electrophilicity index, ω_α: {electrophilicity_index_alpha:.4f} {energy_unit}"
            )
            logger.info("")

        # Beta channel analysis
        if beta_homo is not None and beta_lumo is not None:
            logger.info("--- Beta Spin Channel ---")
            beta_homo *= conversion_factor
            beta_lumo *= conversion_factor
            beta_fmo_gap *= conversion_factor

            logger.info(f"β-HOMO energy: {beta_homo:.4f} {energy_unit}")
            logger.info(f"β-LUMO energy: {beta_lumo:.4f} {energy_unit}")
            logger.info(f"β-HOMO-LUMO gap: {beta_fmo_gap:.4f} {energy_unit}")

            # Reactivity descriptors of beta channel for open-shell systems
            chemical_potential_beta = 0.5 * (beta_lumo + beta_homo)
            chemical_hardness_beta = 0.5 * (beta_lumo - beta_homo)
            electrophilicity_index_beta = chemical_potential_beta**2 / (
                2 * chemical_hardness_beta
            )

            logger.info(
                f"Chemical potential, μ_β: {chemical_potential_beta:.4f} {energy_unit}"
            )
            logger.info(
                f"Chemical hardness, η_β: {chemical_hardness_beta:.4f} {energy_unit}"
            )
            logger.info(
                f"Electrophilicity index, ω_β: {electrophilicity_index_beta:.4f} {energy_unit}"
            )
            logger.info("")

        # SOMO analysis
        if somo_energies is not None:
            logger.info("--- Singly Occupied Molecular Orbitals (SOMOs) ---")
            logger.info(
                f"Number of unpaired electrons: {outputfile.num_unpaired_electrons}"
            )
            somo = [e * conversion_factor for e in somo_energies]
            for i, energy in enumerate(somo, 1):
                logger.info(f"SOMO-{i} energy: {energy:.4f} {energy_unit}")

            lowest_somo *= conversion_factor
            highest_somo *= conversion_factor
            logger.info(f"Lowest SOMO energy: {lowest_somo:.4f} {energy_unit}")
            logger.info(
                f"Highest SOMO energy: {highest_somo:.4f} {energy_unit}"
            )
            logger.info("")

        # Overall FMO gap
        if fmo_gap is not None:
            logger.info("--- Overall FMO Gap ---")
            fmo_gap *= conversion_factor
            logger.info(
                f"FMO gap (highest SOMO to lowest LUMO): {fmo_gap:.4f} {energy_unit}"
            )
            logger.info("")


if __name__ == "__main__":
    entry_point()
