import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.xtb.common import build_xtb_jobs
from chemsmart.cli.xtb.xtb import require_xtb_filename, xtb
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@xtb.command("opt", cls=MyCommand)
@click_job_options
@click.option(
    "--optimization-level",
    type=click.Choice(
        [
            "crude",
            "sloppy",
            "loose",
            "lax",
            "normal",
            "tight",
            "vtight",
            "extreme",
        ],
        case_sensitive=False,
    ),
    default=None,
    help="xTB optimization convergence level.",
)
@click.pass_context
def opt(ctx, skip_completed, optimization_level, **kwargs):
    """Run xTB geometry optimization calculations."""
    require_xtb_filename(ctx)
    job_settings = ctx.obj["job_settings"]
    keywords = list(ctx.obj["keywords"])
    if optimization_level is not None:
        job_settings.optimization_level = optimization_level.lower()
        keywords.append("optimization_level")

    opt_settings = ctx.obj["project_settings"].opt_settings()
    opt_settings = opt_settings.merge(job_settings, keywords=tuple(keywords))
    check_charge_and_multiplicity(opt_settings)
    logger.info(f"Final xTB opt settings: {opt_settings.__dict__}")

    from chemsmart.jobs.xtb.opt import XTBOptJob

    return build_xtb_jobs(ctx, XTBOptJob, opt_settings, skip_completed, kwargs)
