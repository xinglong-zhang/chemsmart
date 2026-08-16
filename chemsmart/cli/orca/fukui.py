"""CLI for ORCA Fukui job submission.

For a given molecule, submits population calculations for the neutral,
radical-cation, and radical-anion charge states at the same geometry.

Output analysis is backend-independent:
  chemsmart run fukui -n <label>_n.out ...
"""

import logging

import click

from chemsmart.analysis.fukui import FUKUI_MODES
from chemsmart.cli.job import click_job_options
from chemsmart.cli.orca.orca import orca
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@orca.command("fukui", cls=MyCommand)
@click_job_options
@click.option(
    "-m",
    "--mode",
    default="mulliken",
    show_default=True,
    type=click.Choice(list(FUKUI_MODES), case_sensitive=False),
    help="Charges to be used for Fukui Indices calculations.",
)
@click.option(
    "-rcc",
    "--radical-cation-charge",
    type=int,
    default=None,
    help=(
        "Override charge for the radical-cation job. "
        "Default is derived from the neutral charge."
    ),
)
@click.option(
    "-rcm",
    "--radical-cation-multiplicity",
    type=int,
    default=None,
    help=(
        "Override multiplicity for the radical-cation job. "
        "Default is derived from the neutral multiplicity."
    ),
)
@click.option(
    "-rac",
    "--radical-anion-charge",
    type=int,
    default=None,
    help=(
        "Override charge for the radical-anion job. "
        "Default is derived from the neutral charge."
    ),
)
@click.option(
    "-ram",
    "--radical-anion-multiplicity",
    type=int,
    default=None,
    help=(
        "Override multiplicity for the radical-anion job. "
        "Default is derived from the neutral multiplicity."
    ),
)
@click.pass_context
def fukui(
    ctx,
    skip_completed,
    mode,
    radical_cation_charge,
    radical_cation_multiplicity,
    radical_anion_charge,
    radical_anion_multiplicity,
    **kwargs,
):
    """Submit ORCA Fukui charge-state calculations.

    Always runs neutral, radical-cation, and radical-anion population jobs
    from the parent ``-f`` structure. Analyze completed outputs with
    ``chemsmart run fukui``.
    """
    jobrunner = ctx.obj["jobrunner"]
    project_settings = ctx.obj["project_settings"]
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    if mode.lower() == "nbo":
        pop_settings = project_settings.wbi_settings()
    else:
        pop_settings = project_settings.sp_settings()
    pop_settings = pop_settings.merge(job_settings, keywords=keywords)
    check_charge_and_multiplicity(pop_settings)

    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]
    label = ctx.obj["label"]

    logger.info(f"Creating ORCA Fukui job: mode={mode}")

    from chemsmart.jobs.orca.fukui import ORCAFukuiJob

    return ORCAFukuiJob(
        molecule=molecule,
        settings=pop_settings,
        label=label,
        jobrunner=jobrunner,
        mode=mode,
        radical_cation_charge=radical_cation_charge,
        radical_cation_multiplicity=radical_cation_multiplicity,
        radical_anion_charge=radical_anion_charge,
        radical_anion_multiplicity=radical_anion_multiplicity,
        skip_completed=skip_completed,
        **kwargs,
    )
