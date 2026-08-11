"""Preview-only PySCF TDA/TDDFT CLI leaf."""

import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.pyscf.common import build_pyscf_jobs
from chemsmart.cli.pyscf.pyscf import pyscf
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@pyscf.command("td", cls=MyCommand)
@click_job_options
@click.pass_context
def td(ctx, skip_completed, **kwargs):
    """Preview a closed-shell gas-phase singlet TDA/TDDFT calculation.

    This experimental leaf deliberately supports ChemSmart fake/test preview
    only.  A non-fake runner rejects the node before launching PySCF.
    """

    from chemsmart.jobs.pyscf.td import PySCFTDJob

    settings = ctx.obj["project_settings"].td_settings()
    settings = settings.merge(
        ctx.obj["job_settings"], keywords=ctx.obj["keywords"]
    )
    logger.info(f"Final PySCF td settings: {settings.__dict__}")

    return build_pyscf_jobs(ctx, PySCFTDJob, settings, skip_completed, kwargs)
