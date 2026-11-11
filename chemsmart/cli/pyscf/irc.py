"""PySCF intrinsic reaction coordinate CLI leaf."""

import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.pyscf.common import build_pyscf_jobs
from chemsmart.cli.pyscf.pyscf import pyscf
from chemsmart.jobs.pyscf.settings import PYSCF_IRC_DIRECTIONS
from chemsmart.utils.cli import MyCommand, update_irc_label

logger = logging.getLogger(__name__)


@pyscf.command("irc", cls=MyCommand)
@click_job_options
@click.option(
    "-d",
    "--direction",
    type=click.Choice(PYSCF_IRC_DIRECTIONS, case_sensitive=False),
    default=None,
    help="Branches to follow from the transition state (default: both). "
    "forward and backward are the two signs of the imaginary mode, not a "
    "reactant or product side.",
)
@click.pass_context
def irc(ctx, skip_completed, direction, **kwargs):
    """Follow the intrinsic reaction coordinate from a transition state.

    Uses the ``gas:`` section of the project YAML. The starting structure
    must have exactly one imaginary mode; its analytic Hessian is checked
    before any step is taken.
    """
    from chemsmart.jobs.pyscf.irc import PySCFIRCJob

    settings = ctx.obj["project_settings"].irc_settings()
    settings = settings.merge(
        ctx.obj["job_settings"], keywords=ctx.obj["keywords"]
    )
    if direction is not None:
        settings.irc_direction = direction.lower()
    # Before the label: a direction from the project YAML that PySCF
    # refuses must be reported with PySCF's directions, not the label
    # helper's, which also accepts Gaussian's 'reverse'.
    settings.validate()
    if settings.irc_direction != "both":
        # One branch per label, so the other branch is not skipped as
        # already complete: water_irc -> water_ircf / water_ircr.
        ctx.obj["label"] = update_irc_label(
            ctx.obj["label"], settings.irc_direction, flat_irc=False
        )
    logger.info(f"Final PySCF irc settings: {settings.__dict__}")

    return build_pyscf_jobs(ctx, PySCFIRCJob, settings, skip_completed, kwargs)
