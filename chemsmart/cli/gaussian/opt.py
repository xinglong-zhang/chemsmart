import click
import logging

from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyCommand
from chemsmart.cli.gaussian.gaussian import gaussian

logger = logging.getLogger(__name__)


@gaussian.command("opt", cls=MyCommand)
@click_job_options
@click.option(
    "-f",
    "--freeze-atoms",
    type=str,
    help="Indices of atoms to freeze for constrained optimization.",
)
@click.pass_context
def opt(ctx, freeze_atoms, skip_completed, **kwargs):
    # folder = ctx.obj["folder"]

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    opt_settings = project_settings.opt_settings()

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project opt settings with job settings from cli keywords from cli.gaussian.py subcommands
    opt_settings = opt_settings.merge(job_settings, keywords=keywords)

    # get atoms
    molecules = ctx.obj["molecules"]
    molecule = molecules[
        -1
    ]  # get last atom from list of atoms from cli.gaussian.py subcommands
    # index = '-1' would access the right structure from the list of atoms returned from cli.gaussian.py subcommands
    # user specified index was used there to return the right atoms and store it as a list of single element/itself

    # get label for the job
    label = ctx.obj["label"]

    # Set atoms to freeze
    from ase.constraints import FixAtoms

    from chemsmart.utils.utils import get_list_from_string_range

    if freeze_atoms is not None:
        frozen_atoms_list = get_list_from_string_range(freeze_atoms)
        # atoms.set_constraint(FixAtoms(frozen_atoms_list))

    logger.info(f"Opt settings from project: {opt_settings.__dict__}")

    from chemsmart.jobs.gaussian import GaussianGeomOptJob

    return GaussianGeomOptJob(
        molecule=molecule,
        settings=opt_settings,
        label=label,
        skip_completed=skip_completed,
        **kwargs,
    )
