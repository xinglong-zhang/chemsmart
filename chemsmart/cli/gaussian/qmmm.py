import logging

import click

from chemsmart.cli.gaussian.gaussian import gaussian, click_gaussian_jobtype_options
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyCommand, get_setting_from_jobtype_for_gaussian

logger = logging.getLogger(__name__)


@gaussian.command("qmmm", cls=MyCommand)
@click_job_options
@click.option(
    "-fh",
    "--functional-high",
    type=str,
    help="functional of high-level layer.",
)
@click.option(
    "-bh",
    "--basis-high",
    type=str,
    help="functional of high-level layer.",
)

@click.option(
    "-ffh",
    "--force-field-high",
    type=str,
    help="force field of high-level layer.",
)

@click.option(
    "-fm",
    "--functional-medium",
    type=str,
    help="functional of medium-level layer.",
)
@click.option(
    "-bm",
    "--basis-medium",
    type=str,
    help="basis of medium-level layer.",
)

@click.option(
    "-ffm",
    "--force-field-medium",
    type=str,
    help="force field of medium-level layer.",
)


@click.option(
    "-fl",
    "--functional-low",
    type=str,
    help="functional of low-level layer.",
)
@click.option(
    "-bl",
    "--basis-low",
    type=str,
    help="basis of low-level layer.",
)
@click.option(
    "-ffl",
    "--force-field-low",
    type=str,
    help="force field of low-level layer.",
)
@click.option(
    "-hlc",
    "--high-level-charge",
    type=str,
    help="charge of high-level layer.",
)
@click.option(
    "-hlm",
    "--high-level-multiplicity",
    type=str,
    help="multiplicity of high-level layer.",
)

@click.option(
    "-mlc",
    "--medium-level-charge",
    type=str,
    help="charge of medium-level layer.",
)

@click.option(
    "-mlm",
    "--medium-level-multiplicity",
    type=str,
    help="multiplicity of medium-level layer.",
)

@click.option(
    "-llc",
    "--low-level-charge",
    type=str,
    help="charge of low-level layer.",
)
@click.option(
    "-llm",
    "--low-level-multiplicity",
    type=str,
    help="multiplicity of low-level layer.",
)

@click.option(
    "-ha",
    "--high-level-atoms",
    type=str,
    help="atom index for high level.",
)
@click.option(
    "-ma",
    "--medium-level-atoms",
    type=str,
    help="atom index for medium level.",
)
@click.option(
    "-la",
    "--low-level-atoms",
    type=str,
    help="atom index for low level.",
)
@click.option(
    "-b",
    "--bonded-atom",
    type=str,
    help="the bond to be cut, specified by two atomic indexes.",
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
    #can be optional if this can be read from input file.
)

@click.option(
    "-j",
    "--jobtype",
    type=str,
    help="job type, e.g., opt, sp, etc..",
    #can be optional if this can be read from input file.
)

# @click.pass_context
# def cli(ctx, **kwargs):
#     """Main CLI Command"""
#     ctx.ensure_object(dict)
#     ctx.obj["settings"] = kwargs



@click.pass_context
def qmmm(
        ctx,
        functional_high,
        basis_high,
        force_field_high,
        functional_medium,
        basis_medium,
        force_field_medium,
        functional_low,
        basis_low,
        force_field_low,
        high_level_charge,
        high_level_multiplicity,
        medium_level_charge,
        medium_level_multiplicity,
        low_level_charge,
        low_level_multiplicity,
        high_level_atoms,
        medium_level_atoms,
        low_level_atoms,
        bonded_atom,
        scale_factor1,
        scale_factor2,
        scale_factor3,
        num_atoms,
        jobtype,  # ✅ Re-added jobtype
        coordinates=None,  # ✅ Default values added
        step_size=None,
        num_steps=None,
        **kwargs,
):
    # Ensure ctx.obj is initialized
    if ctx.obj is None:
        ctx.obj = {}

    if "settings" not in ctx.obj:
        ctx.obj["settings"] = {}

    if "molecules" not in ctx.obj or not ctx.obj["molecules"]:
        raise ValueError("No molecule found in context!")


    from chemsmart.jobs.gaussian.settings import GaussianQMMMJobSettings

    project_settings = ctx.obj["project_settings"]
    qmmm_settings = project_settings.qmmm_settings()
    if qmmm_settings is None:
        logger.warning("qmmm_settings is None! Using default settings.")
        qmmm_settings = GaussianQMMMJobSettings()  # ✅ Provide a default object

    job_settings = ctx.obj.get("job_settings", {})
    keywords = ctx.obj.get("keywords",  {})

    qmmm_settings = qmmm_settings.merge(job_settings, keywords=keywords)
    qmmm_settings = GaussianQMMMJobSettings(**qmmm_settings.__dict__)

    # Assigning QM/MM parameters
    qmmm_parameters = {
        "functional_high": functional_high,
        "basis_high": basis_high,
        "force_field_high": force_field_high,
        "functional_medium": functional_medium,
        "basis_medium": basis_medium,
        "force_field_medium": force_field_medium,
        "functional_low": functional_low,
        "basis_low": basis_low,
        "force_field_low": force_field_low,
        "high_level_charge": high_level_charge,
        "high_level_multiplicity": high_level_multiplicity,
        "medium_level_charge": medium_level_charge,
        "medium_level_multiplicity": medium_level_multiplicity,
        "low_level_charge": low_level_charge,
        "low_level_multiplicity": low_level_multiplicity,
        "high_level_atoms": high_level_atoms,
        "medium_level_atoms": medium_level_atoms,
        "low_level_atoms": low_level_atoms,
        "bonded_atom": bonded_atom,
        "scale_factor1": scale_factor1,
        "scale_factor2": scale_factor2,
        "scale_factor3": scale_factor3,
        "num_atoms": num_atoms,
        "jobtype": jobtype,
    }

    if high_level_charge is not None:
        qmmm_settings.charge = high_level_charge
    if high_level_multiplicity is not None:
        qmmm_settings.multiplicity = high_level_multiplicity
    print(f"Initial charge: {qmmm_settings.charge}, hlc: {high_level_charge}")
    print(f"Initial multiplicity: {qmmm_settings.multiplicity}, hlm: {high_level_multiplicity}")

    for key, value in qmmm_parameters.items():
        setattr(qmmm_settings, key, value)


    molecule = ctx.obj["molecules"][-1]
    logger.info(f"ONIOM calculation of molecule: {molecule}.")

    # Get label for the job
    label = ctx.obj.get("label", "default_label")  # ✅ Default added
    if jobtype:
        label = f"{label}_{jobtype}_QM/MM"

    logger.debug(f"Label for job: {label}")
    logger.info(f"Running QM/MM job with settings: {qmmm_settings.__dict__}")

    from chemsmart.jobs.gaussian.qmmm import GaussianQMMMJob
    return GaussianQMMMJob(
        molecule=molecule,
        settings=qmmm_settings,
        label=label,
        **kwargs,
    )