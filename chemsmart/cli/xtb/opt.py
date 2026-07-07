import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.xtb.xtb import xtb
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@xtb.command("opt", cls=MyCommand)
@click_job_options
@click.pass_context
def opt(ctx, skip_completed, **kwargs):
    """Run xTB geometry optimization calculations."""
    if ctx.obj.get("xtb_missing_filename"):
        raise ValueError("xTB jobs require -f/--filename.")
    opt_settings = ctx.obj["project_settings"].opt_settings()
    opt_settings = opt_settings.merge(
        ctx.obj["job_settings"], keywords=ctx.obj["keywords"]
    )
    check_charge_and_multiplicity(opt_settings)
    logger.debug("Final xTB opt settings: %s", opt_settings.__dict__)

    from chemsmart.jobs.xtb.opt import XTBOptJob

    return _build_xtb_jobs(
        ctx, XTBOptJob, opt_settings, skip_completed, kwargs
    )


def _build_xtb_jobs(ctx, job_cls, settings, skip_completed, kwargs):
    jobrunner = ctx.obj["jobrunner"]
    molecules = ctx.obj["molecules"]
    molecule_indices = ctx.obj["molecule_indices"]
    label = ctx.obj["label"]

    if len(molecules) > 1 and molecule_indices is not None:
        jobs = []
        for molecule, idx in zip(molecules, molecule_indices):
            # Preserve one output directory per selected structure.
            molecule_label = f"{label}_idx{idx}"
            logger.info("Creating xTB job %s", molecule_label)
            jobs.append(
                job_cls(
                    molecule=molecule,
                    settings=settings,
                    label=molecule_label,
                    jobrunner=jobrunner,
                    skip_completed=skip_completed,
                    **kwargs,
                )
            )
        return jobs

    molecule = molecules[-1]
    logger.info("Creating xTB job %s", label)
    return job_cls(
        molecule=molecule,
        settings=settings,
        label=label,
        jobrunner=jobrunner,
        skip_completed=skip_completed,
        **kwargs,
    )
