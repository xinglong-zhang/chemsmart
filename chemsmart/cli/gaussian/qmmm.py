import ast
import logging

import click

from chemsmart.cli.gaussian.gaussian import (
    gaussian,
)
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import (
    MyCommand,
)
from chemsmart.utils.utils import get_list_from_string_range

logger = logging.getLogger(__name__)


@gaussian.command("qmmm", cls=MyCommand)
@click_job_options
@click.option(
    "-hx",
    "--high-level-functional",
    type=str,
    help="High-level layer functional.",
)
@click.option(
    "-hb",
    "--high-level-basis",
    type=str,
    help="High-level layer basis.",
)
@click.option(
    "-hf",
    "--high-level-force-field",
    type=str,
    help="High-level layer force field.",
)
@click.option(
    "-mx",
    "--medium-level-functional",
    type=str,
    help="Medium-level layer functional.",
)
@click.option(
    "-mb",
    "--medium-level-basis",
    type=str,
    help="Medium-level layer basis.",
)
@click.option(
    "-mf",
    "--Medium-level-force-field",
    type=str,
    help="Medium-level layer force field.",
)
@click.option(
    "-lx",
    "--low-level-functional",
    type=str,
    help="Low level layer functional.",
)
@click.option(
    "-lb",
    "--low-level-basis",
    type=str,
    help="Low level layer basis.",
)
@click.option(
    "-lf",
    "--low-level-force-field",
    type=str,
    help="Low level layer force field.",
)
@click.option(
    "-cr",
    "--real-charge",
    type=int,
    help="Charge of real system.",
)
@click.option(
    "-mr",
    "--real-multiplicity",
    type=int,
    help="Spin multiplicity of real system.",
)
@click.option(
    "-ci",
    "--int-charge",
    type=int,
    help="Charge of intermediate system.",
)
@click.option(
    "-mi",
    "--int-multiplicity",
    type=int,
    help="Spin multiplicity of intermediate system.",
)
@click.option(
    "-cm",
    "--model-charge",
    type=int,
    help="Charge of model system.",
)
@click.option(
    "-mm",
    "--model-multiplicity",
    type=int,
    help="Spin multiplicity of model system.",
)
@click.option(
    "-ha",
    "--high-level-atoms",
    type=str,
    help="Atom indices for high level.",
)
@click.option(
    "-ma",
    "--medium-level-atoms",
    type=str,
    help="Atom indices for medium level.",
)
@click.option(
    "-la",
    "--low-level-atoms",
    type=str,
    help="Atom indices for low level.",
)
@click.option(
    "-b",
    "--bonded-atoms",
    help="List of tuples of the bonds to be cut, specified by two atomic indexes in each tuple.",
)
@click.option(
    "-s",
    "--scale-factors",
    type=dict,
    help="A dictionary of scale factors for QM/MM calculations, where the key is the bonded atom "
    "pair indices and the value is a list of scale factors for (low, medium, high).",
)
@click.pass_context
def qmmm(
    ctx,
    high_level_functional,
    high_level_basis,
    high_level_force_field,
    medium_level_functional,
    medium_level_basis,
    medium_level_force_field,
    low_level_functional,
    low_level_basis,
    low_level_force_field,
    real_charge,
    real_multiplicity,
    int_charge,
    int_multiplicity,
    model_charge,
    model_multiplicity,
    high_level_atoms,
    medium_level_atoms,
    low_level_atoms,
    bonded_atoms,
    scale_factors,
    **kwargs,
):
    from chemsmart.jobs.gaussian.settings import GaussianQMMMJobSettings

    print(high_level_atoms)

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    qmmm_settings = project_settings.qmmm_settings()
    print("Project settings:", ctx.obj["project_settings"].__dict__)

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project settings with job settings from cli keywords from cli.gaussian.py subcommands
    qmmm_settings = qmmm_settings.merge(job_settings, keywords=keywords)

    # get label for the job
    label = ctx.obj["label"]
    logger.debug(f"Label for job: {label}")

    # convert from GaussianJobSettings instance to GaussianQMMMJobSettings instance
    qmmm_settings = GaussianQMMMJobSettings(**qmmm_settings.__dict__)

    # populate cli options
    qmmm_settings.high_level_functional = high_level_functional
    qmmm_settings.high_level_basis = high_level_basis
    qmmm_settings.high_level_force_field = high_level_force_field
    qmmm_settings.medium_level_functional = medium_level_functional
    qmmm_settings.medium_level_basis = medium_level_basis
    qmmm_settings.medium_level_force_field = medium_level_force_field
    qmmm_settings.low_level_functional = low_level_functional
    qmmm_settings.low_level_basis = low_level_basis
    qmmm_settings.low_level_force_field = low_level_force_field
    qmmm_settings.real_charge = real_charge
    qmmm_settings.real_multiplicity = real_multiplicity
    qmmm_settings.int_charge = int_charge
    qmmm_settings.int_multiplicity = int_multiplicity
    qmmm_settings.model_charge = model_charge
    qmmm_settings.model_multiplicity = model_multiplicity

    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]

    # populate cli options by attaching QMMM parameters to the molecule
    if high_level_atoms is not None:
        high_level_atoms = get_list_from_string_range(high_level_atoms)
        molecule.high_level_atoms = high_level_atoms
    if medium_level_atoms is not None:
        medium_level_atoms = get_list_from_string_range(medium_level_atoms)
        molecule.medium_level_atoms = medium_level_atoms
    if low_level_atoms is not None:
        low_level_atoms = get_list_from_string_range(low_level_atoms)
        molecule.low_level_atoms = low_level_atoms
    if bonded_atoms is not None:
        bonded_atoms = ast.literal_eval(bonded_atoms)
        molecule.bonded_atoms = bonded_atoms

    if scale_factors is not None:
        scale_factors = ast.literal_eval(scale_factors)
        molecule.scale_factors = scale_factors

    from chemsmart.jobs.gaussian.qmmm import GaussianQMMMJob

    return GaussianQMMMJob(
        molecule=molecule,
        settings=qmmm_settings,
        label=label,
        **kwargs,
    )
