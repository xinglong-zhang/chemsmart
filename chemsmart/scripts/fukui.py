#!/usr/bin/env python

"""Script for calculating Fukui reactivity indices based on charges.
Structures for the neutral, radical cationic and radical anionic
 systems should be based on the optimised neural system."""
import logging
import os

import click

from chemsmart.io.gaussian.output import Gaussian16WBIOutput
from chemsmart.io.orca.output import ORCAOutput
from chemsmart.utils.logger import create_logger

logger = logging.getLogger(__name__)
os.environ["OMP_NUM_THREADS"] = "1"


@click.command()
@click.option(
    "-n",
    "--neutral-filename",
    required=True,
    default=None,
    type=str,
    help="Gaussian or ORCA output file for the neutral system.",
)
@click.option(
    "-c",
    "--radical-cation-filename",
    default=None,
    type=str,
    help="Gaussian or ORCA output file for the radical cationic system.",
)
@click.option(
    "-a",
    "--radical-anion-filename",
    default=None,
    type=str,
    help="Gaussian or ORCA output file for the radical anionic system.",
)
@click.option(
    "-m",
    "--mode",
    default="mulliken",
    type=click.Choice(
        ["mulliken", "nbo", "hirshfeld", "cm5"], case_sensitive=False
    ),
    help="Charges to be used for Fukui Indices calculations.",
)
def entry_point(
    neutral_filename,
    radical_cation_filename=None,
    radical_anion_filename=None,
    mode="mulliken",
):
    create_logger()

    if radical_cation_filename is None and radical_anion_filename is None:
        raise ValueError(
            "At least one of radical cation or radical anion files must be provided."
        )

    if neutral_filename.endswith(".log"):
        file = "gaussian"
    elif neutral_filename.endswith(".out"):
        file = "orca"
    else:
        raise TypeError(f"File {neutral_filename} is of unknown filetype.")

    neutral_output = None
    radical_cation_output = None
    radical_anion_output = None
    charge_for_neutral = None
    charge_for_radical_cation = None
    charge_for_radical_anion = None

    if file == "gaussian":
        neutral_output = Gaussian16WBIOutput(neutral_filename)
        if radical_cation_filename is not None:
            radical_cation_output = Gaussian16WBIOutput(
                radical_cation_filename
            )
        if radical_anion_filename is not None:
            radical_anion_output = Gaussian16WBIOutput(radical_anion_filename)
    elif file == "orca":
        neutral_output = ORCAOutput(neutral_filename)
        if radical_cation_filename is not None:
            radical_cation_output = ORCAOutput(radical_cation_filename)
        if radical_anion_filename is not None:
            radical_anion_output = ORCAOutput(radical_anion_filename)

    if mode == "mulliken":
        logger.info(
            "\nUsing Mulliken Charges for computing Fukui Reactivity Indices."
        )
        charge_for_neutral = neutral_output.mulliken_atomic_charges
        if radical_cation_filename is not None:
            charge_for_radical_cation = (
                radical_cation_output.mulliken_atomic_charges
            )
        if radical_anion_filename is not None:
            charge_for_radical_anion = (
                radical_anion_output.mulliken_atomic_charges
            )
    elif mode == "nbo":
        logger.info(
            "\nUsing NBO Charges for computing Fukui Reactivity Indices."
        )
        charge_for_neutral = neutral_output.natural_charges
        if radical_cation_filename is not None:
            charge_for_radical_cation = radical_cation_output.natural_charges
        if radical_anion_filename is not None:
            charge_for_radical_anion = radical_anion_output.natural_charges
    elif mode == "hirshfeld":
        logger.info(
            "\nUsing Hirshfeld Charges for computing Fukui Reactivity Indices."
        )
        charge_for_neutral = neutral_output.hirshfeld_charges
        if radical_cation_filename is not None:
            charge_for_radical_cation = radical_cation_output.hirshfeld_charges
        if radical_anion_filename is not None:
            charge_for_radical_anion = radical_anion_output.hirshfeld_charges
    elif mode == "cm5":
        logger.info(
            "\nUsing CM5 Charges for computing Fukui Reactivity Indices."
        )
        assert (
            file == "gaussian"
        ), "CM5 charges are only available for Gaussian outputs."
        charge_for_neutral = neutral_output.hirshfeld_cm5_charges
        if radical_cation_filename is not None:
            charge_for_radical_cation = (
                radical_cation_output.hirshfeld_cm5_charges
            )
        if radical_anion_filename is not None:
            charge_for_radical_anion = (
                radical_anion_output.hirshfeld_cm5_charges
            )
    else:
        raise ValueError(
            f"Unknown mode {mode}. Supported modes are: mulliken, nbo, hirshfeld, cm5."
        )

    logger.info("\nNeutral System Charges:")
    for key, value in charge_for_neutral.items():
        logger.info(f"{key:<6}  :  {value:>8.3f}")
    logger.info("\n")

    if radical_cation_filename is not None:
        logger.info("\nRadical Cationic System Charges:")
        for key, value in charge_for_radical_cation.items():
            logger.info(f"{key:<6}  :  {value:>8.3f}")
        logger.info("\n")

    if radical_anion_filename is not None:
        logger.info("\nRadical Anionic System Charges:")
        for key, value in charge_for_radical_anion.items():
            logger.info(f"{key:<6}  :  {value:>8.3f}")
        logger.info("\n")

    logger.info(
        "\nAtom        Fukui Minus (f-)   Fukui Plus(f+)    Fukui Zero(f0)    Fukui Dual Descriptor(f(2))"
    )
    for key, value in charge_for_neutral.items():
        if radical_cation_filename is not None:
            fukui_minus = charge_for_radical_cation[key] - value
        else:
            fukui_minus = 0.0
        if radical_anion_filename is not None:
            fukui_plus = value - charge_for_radical_anion[key]
        else:
            fukui_plus = 0.0
        if (
            radical_cation_filename is not None
            and radical_anion_filename is not None
        ):
            fukui_zero = 0.5 * (
                charge_for_radical_cation[key] - charge_for_radical_anion[key]
            )
            fukui_dual = fukui_plus - fukui_minus
        else:
            fukui_zero = 0.0
            fukui_dual = 0.0
        logger.info(
            f"{key:<4}        {fukui_minus:>12.3f}     {fukui_plus:>12.3f}     {fukui_zero:>12.3f}     {fukui_dual:>12.3f}"
        )


if __name__ == "__main__":
    entry_point()
