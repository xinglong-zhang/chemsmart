import logging

import click

from chemsmart.cli.gaussian.gaussian import click_gaussian_td_options, gaussian
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@gaussian.command("td", cls=MyCommand)
@click_job_options
@click_gaussian_td_options
@click.pass_context
def td(ctx, states, root, nstates, eqsolv, **kwargs):
    """CLI subcommand for running Gaussian TDDFT jobs."""

    # get jobrunner for running Gaussian TDDFT jobs
    jobrunner = ctx.obj["jobrunner"]

    from chemsmart.jobs.gaussian.settings import GaussianTDDFTJobSettings

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    td_settings = project_settings.td_settings()

    # job setting from filename or default, with updates from user in cli
    # specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project settings with job settings from cli keywords from
    # cli.gaussian.py subcommands
    td_settings = td_settings.merge(job_settings, keywords=keywords)
    check_charge_and_multiplicity(td_settings)

    # get molecules
    molecules = ctx.obj["molecules"]

    # get label for the job
    label = ctx.obj["label"]
    logger.debug(f"Label for job: {label}")

    # convert from GaussianJobSettings instance to GaussianTDJobSettings
    # instance
    td_settings = GaussianTDDFTJobSettings(**td_settings.__dict__)

    # populate cli options
    td_settings.states = states
    td_settings.root = root
    td_settings.nstates = nstates
    td_settings.eqsolv = eqsolv

    logger.info(f"TDDFT job settings from project: {td_settings.__dict__}")

    from chemsmart.jobs.gaussian.tddft import GaussianTDDFTJob

    # Get the original molecule indices from context
    molecule_indices = ctx.obj.get("molecule_indices", list(range(1, len(molecules) + 1)))

    # Handle multiple molecules: create one job per molecule
    if len(molecules) > 1:
        logger.info(f"Creating {len(molecules)} TDDFT jobs")
        jobs = []
        for molecule, idx in zip(molecules, molecule_indices):
            molecule_label = f"{label}_idx{idx}"
            logger.info(
                f"Running TDDFT for molecule {idx}: {molecule} with label {molecule_label}"
            )

            job = GaussianTDDFTJob(
                molecule=molecule,
                settings=td_settings,
                label=molecule_label,
                jobrunner=jobrunner,
                **kwargs,
            )
            jobs.append(job)
        return jobs
    else:
        # Single molecule case
        molecule = molecules[-1]
        return GaussianTDDFTJob(
            molecule=molecule,
            settings=td_settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
