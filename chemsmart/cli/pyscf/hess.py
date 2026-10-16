"""PySCF Hessian / frequency CLI leaf."""

import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.pyscf.common import build_pyscf_jobs
from chemsmart.cli.pyscf.pyscf import pyscf
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@pyscf.command("hess", cls=MyCommand)
@click_job_options
@click.pass_context
def hess(ctx, skip_completed, **kwargs):
    """Run a PySCF analytic Hessian and harmonic frequency analysis.

    The geometry is used as given and is not re-optimised, so it should
    already be a stationary point at the same method and basis. Frequencies
    are reported; free energies are computed by ChemSmart's own
    thermochemistry engine so that they stay comparable with Gaussian and
    ORCA results.
    """
    from chemsmart.jobs.pyscf.hess import PySCFHessJob

    settings = ctx.obj["project_settings"].hess_settings()
    settings = settings.merge(
        ctx.obj["job_settings"], keywords=ctx.obj["keywords"]
    )
    logger.info(f"Final PySCF hess settings: {settings.__dict__}")

    return build_pyscf_jobs(
        ctx, PySCFHessJob, settings, skip_completed, kwargs
    )
