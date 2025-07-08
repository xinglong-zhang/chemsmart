import functools
import logging
import os

import click

from chemsmart.jobs.thermochemistry.job import ThermochemistryJob
from chemsmart.jobs.thermochemistry.settings import ThermochemistryJobSettings
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


def click_thermochemistry_options(f):
    """Common click options for Thermochemistry."""

    @click.option(
        "-d",
        "--directory",
        default=None,
        help="Directory in which to compute thermochemistry for all files.",
    )
    @click.option(
        "-t",
        "--filetype",
        default=None,
        help="Type of file to calculate thermochemistry for, if directory is specified.",
    )
    @click.option(
        "-f",
        "--filenames",
        type=str,
        multiple=True,
        default=None,
        help="Gaussian or ORCA output files for parsing thermochemistry.",
    )
    @click.option(
        "-cs",
        "--cutoff-entropy",
        default=None,
        type=float,
        show_default=True,
        help="Cutoff frequency for entropy in wavenumbers",
    )
    @click.option(
        "-ch",
        "--cutoff-enthalpy",
        default=None,
        type=float,
        show_default=True,
        help="Cutoff frequency for enthalpy in wavenumbers",
    )
    @click.option(
        "-c",
        "--concentration",
        default=1.0,
        type=float,
        show_default=True,
        help="Concentration in mol/L",
    )
    @click.option(
        "-p",
        "--pressure",
        default=1.0,
        type=float,
        show_default=True,
        help="Pressure in atm.",
    )
    @click.option(
        "-T",
        "--temperature",
        required=True,
        default=None,
        type=float,
        help="Temperature in Kelvin.",
    )
    @click.option(
        "-a",
        "--alpha",
        default=4,
        type=int,
        show_default=True,
        help="Interpolator exponent used in the quasi-RRHO approximation.",
    )
    @click.option(
        "-w",
        "--weighted",
        is_flag=True,
        default=False,
        show_default=True,
        help="Use natural abundance weighted masses (True) or use most abundant masses (False).\n"
        "Default to False, i.e., use single isotopic mass.",
    )
    @click.option(
        "-u",
        "--units",
        default="hartree",
        show_default=True,
        type=click.Choice(
            ["hartree", "eV", "kcal/mol", "kJ/mol"], case_sensitive=False
        ),
        help="Units of energetic values.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


# use MyGroup to allow potential subcommands in the future
@click.command(cls=MyCommand)
@click_thermochemistry_options
@click.pass_context
def thermochemistry(
    ctx,
    directory,
    filetype,
    filenames,
    cutoff_entropy,
    cutoff_enthalpy,
    concentration,
    pressure,
    temperature,
    alpha,
    weighted,
    units,
):
    """CLI for running thermochemistry jobs using the chemsmart framework."""
    # validate input
    if directory and filenames:
        raise ValueError(
            "Cannot specify both --directory and --filenames. Choose one."
        )
    if directory and not filetype:
        raise ValueError("Must specify --filetype when using --directory.")

    # Create job settings
    job_settings = ThermochemistryJobSettings(
        temperature=temperature,
        concentration=concentration,
        pressure=pressure,
        use_weighted_mass=weighted,
        alpha=alpha,
        s_freq_cutoff=cutoff_entropy,
        h_freq_cutoff=cutoff_enthalpy,
        units=units,
    )

    # Initialize list to store jobs
    jobs = []

    if directory:
        directory = os.path.abspath(directory)
        logger.info(
            f"Obtaining thermochemistry of files in directory: {directory}"
        )
        if filetype == "log":
            from chemsmart.io.gaussian.folder import GaussianLogFolder

            folder = GaussianLogFolder(directory)
            files = folder.all_logfiles
        elif filetype == "out":
            from chemsmart.io.orca.folder import ORCAOutFolder

            folder = ORCAOutFolder(directory)
            files = folder.all_outfiles
        else:
            raise ValueError(
                f"Unsupported filetype '{filetype}'. Use 'log' or 'out'."
            )
        for file in files:
            job = ThermochemistryJob.from_filename(
                filename=file,
                settings=job_settings,
            )
            jobs.append(job)
            logger.info(f"Created thermochemistry job for file: {file}")
    elif filenames:
        for file in filenames:
            if not file.endswith((".log", ".out")):
                raise ValueError(
                    f"Unsupported file extension for '{file}'. Use .log or .out."
                )
            job = ThermochemistryJob.from_filename(
                filename=file,
                settings=job_settings,
            )
            jobs.append(job)
            logger.info(f"Created thermochemistry job for file: {file}")

    # Store objects in context
    ctx.obj["job_settings"] = job_settings
    ctx.obj["jobs"] = jobs
    ctx.obj["filenames"] = (
        filenames if filenames else files if directory else None
    )
    ctx.obj["directory"] = directory
    ctx.obj["filetype"] = filetype

    # if ctx.invoked_subcommand is None:
    #     # If no subcommand is invoked, run the thermochemistry jobs
    #     logger.info("Running thermochemistry calculations on specified jobs.")
    #     for job in jobs:
    #         try:
    #             job.calculate_thermochemistry()
    #             logger.info(f"Thermochemistry calculation completed for {job.label}.")
    #         except Exception as e:
    #             logger.error(f"Error processing job for {job.label}: {e}")

    # @thermochemistry.result_callback()
    # @click.pass_context
    # def thermochemistry_process_pipeline(ctx, *args, **kwargs):
    #     """Process the thermochemistry jobs."""
    #     jobs = ctx.obj.get("jobs", [])
    for job in jobs:
        try:
            job.calculate_thermochemistry()
            logger.info(
                f"Completed thermochemistry job for {job.label}: results saved to {job.outputfile}"
            )
        except Exception as e:
            logger.error(f"Error processing job for {job.label}: {e}")
    # return args[0]

    # if ctx.invoked_subcommand is None:
    #     # If no subcommand is invoked, run the thermochemistry jobs
    #     logger.info("Running thermochemistry calculations on specified jobs.")
    #     for job in ctx.obj["jobs"]:
    #         job.calculate_thermochemistry()
    #         logger.info(f"Thermochemistry calculation completed for {job.label}.")
    #


if __name__ == "__main__":
    thermochemistry()
