import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.xtb.opt import _build_xtb_jobs
from chemsmart.cli.xtb.xtb import xtb
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@xtb.command("sp", cls=MyCommand)
@click_job_options
@click.pass_context
def sp(ctx, skip_completed, **kwargs):
    """Run xTB single point calculations."""
    if ctx.obj.get("xtb_missing_filename"):
        raise ValueError("xTB jobs require -f/--filename.")
    sp_settings = ctx.obj["project_settings"].sp_settings()
    sp_settings = sp_settings.merge(
        ctx.obj["job_settings"], keywords=ctx.obj["keywords"]
    )
    check_charge_and_multiplicity(sp_settings)
    logger.debug("Final xTB sp settings: %s", sp_settings.__dict__)

    from chemsmart.jobs.xtb.singlepoint import XTBSinglePointJob

    return _build_xtb_jobs(
        ctx, XTBSinglePointJob, sp_settings, skip_completed, kwargs
    )
