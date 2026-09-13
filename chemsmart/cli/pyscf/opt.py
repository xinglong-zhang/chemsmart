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

    Uses the ``opt:`` section of the project YAML. With ``excited_state_root``
    (and the td fields) the optimisation follows that root of the TDA/TDDFT
    manifold by index and re-evaluates the spectrum at the reached geometry;
    with ``ab_initio: mp2`` or ``ccsd`` it optimises on the correlated
    surface. A Hessian must be a separate ``hess`` node bound to the
    validated optimized-geometry artifact; none exists for an excited or
    correlated surface in this release.
    """
    from chemsmart.jobs.pyscf.opt import PySCFOptJob

    settings = ctx.obj["project_settings"].opt_settings()
    settings = settings.merge(
        ctx.obj["job_settings"], keywords=ctx.obj["keywords"]
    )
    logger.info(f"Final PySCF opt settings: {settings.__dict__}")

    return build_pyscf_jobs(ctx, PySCFOptJob, settings, skip_completed, kwargs)
