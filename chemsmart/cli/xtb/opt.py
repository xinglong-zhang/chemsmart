import ast
import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.utils import build_jobs
from chemsmart.cli.xtb.xtb import (
    click_xtb_constrain_options,
    click_xtb_optimization_level_option,
    xtb,
)
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@xtb.command("opt", cls=MyCommand)
@click_job_options
@click_xtb_optimization_level_option
@click_xtb_constrain_options
@click.pass_context
def opt(
    ctx,
    skip_completed,
    optimization_level,
    constrain,
    force_constant,
    **kwargs,
):
    """Run xTB geometry optimization calculations."""
    job_settings = ctx.obj["job_settings"]
    keywords = list(ctx.obj["keywords"])
    if optimization_level is not None:
        job_settings.optimization_level = optimization_level.lower()
        keywords.append("optimization_level")
    if constrain is not None:
        job_settings.constraints = ast.literal_eval(constrain)
        keywords.append("constraints")
        if force_constant is not None:
            job_settings.force_constant = force_constant
            keywords.append("force_constant")

    opt_settings = ctx.obj["project_settings"].opt_settings()
    opt_settings = opt_settings.merge(job_settings, keywords=tuple(keywords))
    check_charge_and_multiplicity(opt_settings)
    logger.info(f"Final xTB opt settings: {opt_settings.__dict__}")

    from chemsmart.jobs.xtb.opt import XTBOptJob

    return build_jobs(ctx, XTBOptJob, opt_settings, skip_completed, kwargs)
