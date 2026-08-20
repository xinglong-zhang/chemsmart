"""xTB single-point CLI leaf."""

import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.xtb.common import build_xtb_jobs
from chemsmart.cli.xtb.xtb import require_xtb_filename, xtb
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@xtb.command("sp", cls=MyCommand)
@click_job_options
@click.pass_context
def sp(ctx, skip_completed, **kwargs):
    """Prepare an xTB single-point calculation."""
    require_xtb_filename(ctx)
    settings = (
        ctx.obj["project_settings"]
        .sp_settings()
        .merge(ctx.obj["job_settings"], keywords=ctx.obj["keywords"])
    )
    settings.validate(expected_jobtype="sp")

    from chemsmart.jobs.xtb.singlepoint import XTBSinglePointJob

    return build_xtb_jobs(
        ctx, XTBSinglePointJob, settings, skip_completed, kwargs
    )
