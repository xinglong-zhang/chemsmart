import logging
import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.gaussian.gaussian import click_gaussian_solvent_options
from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@gaussian.command("sp", cls=MyCommand)
@click_job_options
@click_gaussian_solvent_options
@click.pass_context
def sp(ctx, remove_solvent, solvent_model, solvent_id, **kwargs):

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    sp_settings = project_settings.sp_settings()

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project settings with job settings from cli keywords from cli.gaussian.py subcommands
    sp_settings = sp_settings.merge(job_settings, keywords=keywords)
    check_charge_and_multiplicity(sp_settings)

    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]

    # get label for the job
    label = ctx.obj["label"]
    logger.debug(f"Label for job: {label}")

    # cli-supplied solvent model and solvent id
    sp_settings.modify_solvent(
        remove_solvent=remove_solvent,
        solvent_model=solvent_model,
        solvent_id=solvent_id,
    )

    # either supplied or not from cli, would still want label to have model and id, if both given;
    # will not be activated when both are not given, e.g., in the gaussian calculator calling the sp job
    if (
        sp_settings.solvent_model is not None
        and sp_settings.solvent_id is not None
    ):
        # replace , by _ if it occurs in the solvent name, as , in file will cause gaussian run error
        solvent_label = sp_settings.solvent_id.replace(",", "_")
        solvent_label = solvent_label.replace("-", "_")
        label = f"{label}_{sp_settings.solvent_model}_{solvent_label}"
    elif sp_settings.solvent_model is None and sp_settings.solvent_id is None:
        label = (
            f"{label}_gas_phase"
            if sp_settings.custom_solvent is None
            else f"{label}_custom_solvent"
        )

    logger.info(
        f"Single point job settings from project: {sp_settings.__dict__}"
    )

    from chemsmart.jobs.gaussian.singlepoint import GaussianSinglePointJob

    return GaussianSinglePointJob(
        molecule=molecule, settings=sp_settings, label=label, **kwargs
    )
