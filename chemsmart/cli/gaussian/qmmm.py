import logging

import click

from chemsmart.cli.gaussian.gaussian import (
    gaussian,
)
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import (
    MyCommand,
)

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
    "-lf",
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
    "-hc",
    "--high-level-charge",
    type=str,
    help="High level layer charge.",
)
@click.option(
    "-hm",
    "--high-level-multiplicity",
    type=str,
    help="High level layer spin multiplicity.",
)
@click.option(
    "-mc",
    "--medium-level-charge",
    type=str,
    help="Medium level layer charge.",
)
@click.option(
    "-mm",
    "--medium-level-multiplicity",
    type=str,
    help="Medium level layer spin multiplicity.",
)
@click.option(
    "-lc",
    "--low-level-charge",
    type=str,
    help="Low level layer charge.",
)
@click.option(
    "-lm",
    "--low-level-multiplicity",
    type=str,
    help="Low level layer spin multiplicity.",
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
    "-sf1",
    "--scale-factor1",
    type=str,
    help="Scale factor for bonds between QM and MM region,default=1.0.",
)
@click.option(
    "-sf2",
    "--scale-factor2",
    type=str,
    help="Scale factor for angles involving  MM and MM region,default=1.0.",
)
@click.option(
    "-sf3",
    "--scale-factor3",
    type=str,
    help="Scale factor for torsions, default=1.0",
)
@click.option(
    "-na",
    "--num-atoms",
    type=str,
    help="Number of atoms in the system.",
    # can be optional if this can be read from input file.
)
@click.option(
    "-j",
    "--jobtype",
    type=str,
    help="job type, e.g., opt, sp, etc..",
    # can be optional if this can be read from input file.
)

# @click.pass_context
# def cli(ctx, **kwargs):
#     """Main CLI Command"""
#     ctx.ensure_object(dict)
#     ctx.obj["settings"] = kwargs


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
    high_level_charge,
    high_level_multiplicity,
    medium_level_charge,
    medium_level_multiplicity,
    low_level_charge,
    low_level_multiplicity,
    high_level_atoms,
    medium_level_atoms,
    low_level_atoms,
    bonded_atoms,
    scale_factor1,
    scale_factor2,
    scale_factor3,
    num_atoms,
    jobtype,
    **kwargs,
):
    from chemsmart.jobs.gaussian.settings import GaussianQMMMJobSettings

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    qmmm_settings = project_settings.qmmm_settings()

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

    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]

    # populate cli options by attaching QMMM parameters to the molecule
    molecule.high_level_functional = high_level_functional
    molecule.high_level_basis = high_level_basis
    molecule.high_level_force_field = high_level_force_field
    molecule.medium_level_functional = medium_level_functional

    qmmm_settings.high_level_functional = high_level_functional

    # td_settings.states = states
    # td_settings.root = root
    # td_settings.nstates = nstates
    # td_settings.eqsolv = eqsolv
    #
    # # project_settings = ctx.obj["project_settings"]
    # # qmmm_settings = project_settings.qmmm_settings()
    # # if qmmm_settings is None:
    # #     logger.warning("qmmm_settings is None! Using default settings.")
    # qmmm_settings = GaussianQMMMJobSettings()
    #
    # job_settings = ctx.obj.get("job_settings", {})
    # keywords = ctx.obj.get("keywords", {})
    #
    # qmmm_settings = qmmm_settings.merge(job_settings, keywords=keywords)
    # qmmm_settings = GaussianQMMMJobSettings(**qmmm_settings.__dict__)
    #
    # # Assigning QM/MM parameters
    # qmmm_parameters = {
    #     "functional_high": functional_high,
    #     "basis_high": basis_high,
    #     "force_field_high": force_field_high,
    #     "functional_medium": functional_medium,
    #     "basis_medium": basis_medium,
    #     "force_field_medium": force_field_medium,
    #     "functional_low": functional_low,
    #     "basis_low": basis_low,
    #     "force_field_low": force_field_low,
    #     "high_level_charge": high_level_charge,
    #     "high_level_multiplicity": high_level_multiplicity,
    #     "medium_level_charge": medium_level_charge,
    #     "medium_level_multiplicity": medium_level_multiplicity,
    #     "low_level_charge": low_level_charge,
    #     "low_level_multiplicity": low_level_multiplicity,
    #     "high_level_atoms": high_level_atoms,
    #     "medium_level_atoms": medium_level_atoms,
    #     "low_level_atoms": low_level_atoms,
    #     "bonded_atoms": bonded_atoms,
    #     "scale_factor1": scale_factor1,
    #     "scale_factor2": scale_factor2,
    #     "scale_factor3": scale_factor3,
    #     "num_atoms": num_atoms,
    #     "jobtype": jobtype,
    # }
    #
    # if high_level_charge is not None:
    #     qmmm_settings.charge = high_level_charge
    # if high_level_multiplicity is not None:
    #     qmmm_settings.multiplicity = high_level_multiplicity
    #
    # for key, value in qmmm_parameters.items():
    #     setattr(qmmm_settings, key, value)
    #
    # # get molecule
    # molecules = ctx.obj["molecules"]
    # molecule = molecules[-1]
    #
    # # attach QMMM parameters to the molecule
    #
    # logger.info(f"ONIOM calculation of molecule: {molecule}.")
    #
    # # Get label for the job
    # label = ctx.obj.get("label", "default_label")  # ✅ Default added
    # if jobtype:
    #     label = f"{label}_{jobtype}_QM/MM"
    #
    # logger.debug(f"Label for job: {label}")
    # logger.info(f"Running QM/MM job with settings: {qmmm_settings.__dict__}")

    from chemsmart.jobs.gaussian.qmmm import GaussianQMMMJob

    return GaussianQMMMJob(
        molecule=molecule,
        settings=qmmm_settings,
        label=label,
        **kwargs,
    )
