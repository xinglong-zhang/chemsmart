"""PySCF single point CLI leaf."""

import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.pyscf.common import build_pyscf_jobs
from chemsmart.cli.pyscf.pyscf import pyscf
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@pyscf.command("sp", cls=MyCommand)
@click_job_options
@click.pass_context
def sp(ctx, skip_completed, **kwargs):
    """Run a PySCF single point energy calculation.

    Uses only the complete ``sp:`` section of the project YAML. Solvent and
    every other scientific setting must be declared for that stage; nothing
    is inherited from ``opt`` or ``hess``.
    """
    from chemsmart.jobs.pyscf.singlepoint import PySCFSinglePointJob

    settings = ctx.obj["project_settings"].sp_settings()
    settings = settings.merge(
        ctx.obj["job_settings"], keywords=ctx.obj["keywords"]
    )
    logger.info(f"Final PySCF sp settings: {settings.__dict__}")

    return build_pyscf_jobs(
        ctx, PySCFSinglePointJob, settings, skip_completed, kwargs
    )
