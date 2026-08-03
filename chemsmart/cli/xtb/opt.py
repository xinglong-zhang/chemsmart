"""xTB geometry-optimization CLI leaf."""

import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.xtb.common import build_xtb_jobs
from chemsmart.cli.xtb.xtb import require_xtb_filename, xtb
from chemsmart.io.xtb import XTB_ALL_OPT_LEVELS
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@xtb.command("opt", cls=MyCommand)
@click_job_options
@click.option(
    "--optimization-level",
    type=click.Choice(XTB_ALL_OPT_LEVELS, case_sensitive=False),
    default=None,
    help="xTB optimization convergence level.",
)
@click.pass_context
def opt(ctx, skip_completed, optimization_level, **kwargs):
    """Prepare an xTB geometry optimization."""
    job_settings = ctx.obj.get("job_settings")
    require_xtb_filename(ctx)
    keywords = list(ctx.obj["keywords"])
    if optimization_level is not None:
        job_settings.optimization_level = optimization_level.lower()
        keywords.append("optimization_level")

    settings = ctx.obj["project_settings"].opt_settings().merge(
        job_settings, keywords=tuple(keywords)
    )
    settings.validate(expected_jobtype="opt")

    from chemsmart.jobs.xtb.opt import XTBOptJob

    return build_xtb_jobs(ctx, XTBOptJob, settings, skip_completed, kwargs)
