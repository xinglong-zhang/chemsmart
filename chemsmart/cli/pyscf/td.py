"""PySCF TDA/TDDFT vertical-excitation CLI leaf."""

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
    """Run TDA/TDDFT vertical excitations of the supplied geometry.

    Uses the ``td:`` section of the project YAML: a Kohn-Sham reference,
    ``response_method`` (tda or tddft), ``state_manifold`` (singlet or
    triplet on a closed shell; unrestricted on an open shell) and
    ``nstates``. Roots are ascending indices within the manifold at this
    geometry, never state identities. Implicit solvent gives energies
    under PySCF's non-equilibrium response, which the artifact records.
    """

    from chemsmart.jobs.pyscf.td import PySCFTDJob

    settings = ctx.obj["project_settings"].td_settings()
    settings = settings.merge(
        ctx.obj["job_settings"], keywords=ctx.obj["keywords"]
    )
    logger.info(f"Final PySCF td settings: {settings.__dict__}")

    return build_pyscf_jobs(ctx, PySCFTDJob, settings, skip_completed, kwargs)
