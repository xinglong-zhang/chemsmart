import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.xtb.common import build_xtb_jobs
from chemsmart.cli.xtb.xtb import xtb
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@xtb.command("hess", cls=MyCommand)
@click_job_options
@click.pass_context
def hess(ctx, skip_completed, **kwargs):
    """Run xTB Hessian/frequency calculations."""
    hess_settings = ctx.obj["project_settings"].hess_settings()
    hess_settings = hess_settings.merge(
        ctx.obj["job_settings"], keywords=ctx.obj["keywords"]
    )
    check_charge_and_multiplicity(hess_settings)
    logger.info(f"Final xTB hess settings: {hess_settings.__dict__}")

    from chemsmart.jobs.xtb.hess import XTBHessJob

    return build_xtb_jobs(
        ctx, XTBHessJob, hess_settings, skip_completed, kwargs
    )
