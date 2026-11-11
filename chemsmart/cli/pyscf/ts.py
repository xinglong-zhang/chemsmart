"""PySCF transition-state search CLI leaf."""

import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.pyscf.common import build_pyscf_jobs
from chemsmart.cli.pyscf.pyscf import pyscf
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@pyscf.command("ts", cls=MyCommand)
@click_job_options
@click.pass_context
def ts(ctx, skip_completed, **kwargs):
    """Run a PySCF transition-state search.

    Uses the ``gas:`` section of the project YAML. With ``freq: True`` the
    Hessian at the result must show exactly one imaginary mode, or the job
    fails rather than report a structure that is not a transition state.
    """
    from chemsmart.jobs.pyscf.ts import PySCFTSJob

    settings = ctx.obj["project_settings"].ts_settings()
    settings = settings.merge(
        ctx.obj["job_settings"], keywords=ctx.obj["keywords"]
    )
    logger.info(f"Final PySCF ts settings: {settings.__dict__}")

    return build_pyscf_jobs(ctx, PySCFTSJob, settings, skip_completed, kwargs)
