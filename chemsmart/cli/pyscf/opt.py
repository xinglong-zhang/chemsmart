"""PySCF geometry optimisation CLI leaf."""

import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.pyscf.common import build_pyscf_jobs
from chemsmart.cli.pyscf.pyscf import pyscf
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@pyscf.command("opt", cls=MyCommand)
@click_job_options
@click.pass_context
def opt(ctx, skip_completed, **kwargs):
    """Run a PySCF geometry optimisation.

    Uses the ``opt:`` section of the project YAML. When that section sets
    ``freq: True`` the Hessian runs in the same process against the
    already-converged mean-field object, instead of paying for a second SCF.
    """
    from chemsmart.jobs.pyscf.opt import PySCFOptJob

    settings = ctx.obj["project_settings"].opt_settings()
    settings = settings.merge(
        ctx.obj["job_settings"], keywords=ctx.obj["keywords"]
    )
    logger.info(f"Final PySCF opt settings: {settings.__dict__}")

    return build_pyscf_jobs(ctx, PySCFOptJob, settings, skip_completed, kwargs)
