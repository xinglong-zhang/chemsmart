import ast
import logging

import click

from chemsmart.cli.crest.crest import crest
from chemsmart.cli.job import click_job_options
from chemsmart.cli.utils import build_jobs
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@crest.command("conformers", cls=MyCommand)
@click_job_options
@click.option(
    "-c",
    "--constraints",
    type=str,
    default=None,
    help="List of coordinates to be fixed for constrained conformational search. "
    "1-indexed. Example: [[1,2],[2,3],[3,5]].",
)
@click.option(
    "-f",
    "--force-constant",
    type=float,
    default=None,
    help="Force constant for distance constraints.",
)
@click.pass_context
def conformers(ctx, skip_completed, constraints, force_constant, **kwargs):
    """Run CREST conformational search."""
    job_settings = ctx.obj["job_settings"]
    keywords = list(ctx.obj["keywords"])
    if constraints is not None:
        job_settings.constraints = ast.literal_eval(constraints)
        keywords.append("constraints")
        if force_constant is not None:
            job_settings.force_constant = force_constant
            keywords.append("force_constant")

    conformers_settings = ctx.obj["project_settings"].conformer_settings()
    conformers_settings = conformers_settings.merge(
        job_settings, keywords=tuple(keywords)
    )
    check_charge_and_multiplicity(conformers_settings)
    logger.info(
        f"Conformational sampling job settings from project: "
        f"{conformers_settings.__dict__}"
    )

    from chemsmart.jobs.crest.conformers import CRESTConformerSearchJob

    return build_jobs(
        ctx,
        CRESTConformerSearchJob,
        conformers_settings,
        skip_completed,
        kwargs,
    )
