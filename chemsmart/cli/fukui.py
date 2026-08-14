"""Backend-independent Fukui output analysis.

Registered as ``chemsmart run fukui`` — post-processing only; does not invoke
Gaussian, ORCA, or any other QC backend.

Job submission lives under ``chemsmart sub … gaussian … fukui``.
"""

import logging

import click

from chemsmart.analysis.fukui import (
    FUKUI_MODES,
    analyze_fukui,
    discover_fukui_companion_outputs,
)
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@click.command(name="fukui", cls=MyCommand)
@click.option(
    "-n",
    "--neutral-filename",
    required=True,
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
    show_default=True,
    type=click.Choice(list(FUKUI_MODES), case_sensitive=False),
    help="Charges to be used for Fukui Indices calculations.",
)
def fukui(
    neutral_filename,
    radical_cation_filename=None,
    radical_anion_filename=None,
    mode="mulliken",
):
    """Compute Fukui reactivity indices from existing output files.

    Companion ``_rc`` / ``_ra`` files are auto-discovered from ``-n`` when
    omitted (same labels as Gaussian Fukui job submission).

    \b
    Examples:
      chemsmart run fukui -n mol_n.log -c mol_rc.log -a mol_ra.log

      chemsmart run fukui -n mol_n.log -m nbo
    """
    if radical_cation_filename is None or radical_anion_filename is None:
        discovered = discover_fukui_companion_outputs(neutral_filename)
        if radical_cation_filename is None and discovered["radical_cation"]:
            radical_cation_filename = discovered["radical_cation"]
            logger.info(
                f"Auto-discovered radical cation: {radical_cation_filename}"
            )
        if radical_anion_filename is None and discovered["radical_anion"]:
            radical_anion_filename = discovered["radical_anion"]
            logger.info(
                f"Auto-discovered radical anion: {radical_anion_filename}"
            )

    if radical_cation_filename is None and radical_anion_filename is None:
        raise click.UsageError(
            "At least one of -c/--radical-cation-filename or "
            "-a/--radical-anion-filename must be provided (or auto-discoverable "
            "as <base>_rc / <base>_ra beside -n)."
        )

    analyze_fukui(
        neutral_filename=neutral_filename,
        radical_cation_filename=radical_cation_filename,
        radical_anion_filename=radical_anion_filename,
        mode=mode,
    )
    return None
