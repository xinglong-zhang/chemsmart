import logging

import click

from chemsmart.cli.gaussian.gaussian import (
    click_gaussian_jobtype_options,
    gaussian,
)

# Import and register qmmm subcommand
from chemsmart.cli.gaussian.qmmm import create_qmmm_subcommand
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import (
    MyGroup,
    get_setting_from_jobtype_for_gaussian,
)
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@gaussian.group("modred", cls=MyGroup, invoke_without_command=True)
@click_job_options
@click_gaussian_jobtype_options
@click.pass_context
def modred(
    ctx, jobtype, coordinates, step_size, num_steps, skip_completed, **kwargs
):
    """CLI subcommand for running Gaussian modred jobs.

    Can be used standalone for regular modred or with the 'qmmm'
    subcommand for QM/MM modred calculations.

    Examples:
        chemsmart sub gaussian modred              # Regular modred
        chemsmart sub gaussian modred qmmm         # QM/MM modred
    """

    # get jobrunner for running Gaussian modred jobs
    jobrunner = ctx.obj["jobrunner"]

    if jobtype is None:
        jobtype = "modred"

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    modred_settings = get_setting_from_jobtype_for_gaussian(
        project_settings, jobtype, coordinates, step_size, num_steps
    )

    # job setting from filename or default, with updates from user in cli
    # specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`

    # merge project settings with job settings from cli keywords from
    # cli.gaussian.py subcommands
    # Store parent context for potential qmmm subcommand
    ctx.obj["parent_skip_completed"] = skip_completed
    ctx.obj["parent_freeze_atoms"] = None  # modred doesn't have freeze_atoms
    ctx.obj["parent_kwargs"] = kwargs
    ctx.obj["parent_settings"] = modred_settings
    ctx.obj["modred"] = "modred"
    ctx.obj["parent_jobtype"] = jobtype

    if ctx.invoked_subcommand is None:
        check_charge_and_multiplicity(modred_settings)

        # get molecule
        molecules = ctx.obj["molecules"]

        # get label for the job
        label = ctx.obj["label"]
        logger.debug(f"Label for job: {label}")

        logger.info(
            f"Modred settings from project: {modred_settings.__dict__}"
        )

        # If no subcommand invoked, run regular modred
        from chemsmart.jobs.gaussian.modred import GaussianModredJob

        # Get the original molecule indices from context
        molecule_indices = ctx.obj["molecule_indices"]

        # Handle multiple molecules: create one job per molecule
        if len(molecules) > 1 and molecule_indices is not None:
            logger.info(f"Creating {len(molecules)} modred jobs")
            jobs = []
            for molecule, idx in zip(molecules, molecule_indices):
                molecule_label = f"{label}_idx{idx}"
                logger.info(
                    f"Running modred for molecule {idx}: {molecule} with label {molecule_label}"
                )

                job = GaussianModredJob(
                    molecule=molecule,
                    settings=modred_settings,
                    label=molecule_label,
                    jobrunner=jobrunner,
                    **kwargs,
                )
                jobs.append(job)
            return jobs
        else:
            # Single molecule case
            molecule = molecules[-1]
            return GaussianModredJob(
                molecule=molecule,
                settings=modred_settings,
                label=label,
                jobrunner=jobrunner,
                **kwargs,
            )


create_qmmm_subcommand(modred)
