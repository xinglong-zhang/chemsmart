import functools
import logging
import os

import click

from chemsmart.cli.job import click_job_options
from chemsmart.jobs.thermochemistry.job import ThermochemistryJob
from chemsmart.jobs.thermochemistry.settings import ThermochemistryJobSettings
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.io import outfile_format

logger = logging.getLogger(__name__)


def click_thermochemistry_options(f):
    """
    Common click options for Thermochemistry.
    """

    @click.option(
        "-d",
        "--directory",
        default=None,
        help="Directory in which to compute thermochemistry for all " "files.",
    )
    @click.option(
        "-t",
        "--filetype",
        default=None,
        type=click.Choice(["gaussian", "orca"], case_sensitive=False),
        help="Type of quantum chemistry output file to calculate thermochemistry for, "
        "if directory is specified.",
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
        "-csg",
        "--cutoff-entropy-grimme",
        default=None,
        type=float,
        show_default=True,
        help="Cutoff frequency for entropy in wavenumbers, using Grimme's "
        "quasi-RRHO method.",
    )
    @click.option(
        "-cst",
        "--cutoff-entropy-truhlar",
        default=None,
        type=float,
        show_default=True,
        help="Cutoff frequency for entropy in wavenumbers, using Truhlar's "
        "quasi-RRHO method.",
    )
    @click.option(
        "-ch",
        "--cutoff-enthalpy",
        default=None,
        type=float,
        show_default=True,
        help="Cutoff frequency for enthalpy in wavenumbers, using "
        "Head-Gordon's quasi-RRHO method.",
    )
    @click.option(
        "-c",
        "--concentration",
        default=None,
        type=float,
        show_default=True,
        help="Concentration in mol/L.",
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
        help="Use natural abundance weighted masses (True) or use most "
        "abundant masses (False).\nDefault to False, i.e., use single "
        "isotopic mass.",
    )
    @click.option(
        "-u",
        "--energy-units",
        default="hartree",
        show_default=True,
        type=click.Choice(
            ["hartree", "eV", "kcal/mol", "kJ/mol"], case_sensitive=False
        ),
        help="Units of energetic values.",
    )
    @click.option(
        "-o",
        "--outputfile",
        default=None,
        type=str,
        help="Output file to save the thermochemistry results. Defaults to "
        "None, which will save results to file_basename.dat.\nIf "
        "specified, it will save all thermochemistry results to this file.",
    )
    @click.option(
        "-O",
        "--overwrite",
        is_flag=True,
        default=False,
        show_default=True,
        help="Overwrite existing output files if they already exist.",
    )
    @click.option(
        "-i",
        "--check-imaginary-frequencies",
        is_flag=True,
        default=True,
        show_default=True,
        help="Check for imaginary frequencies in the calculations.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


# use MyGroup to allow potential subcommands in the future
@click.group(cls=MyGroup, invoke_without_command=True)
@click_thermochemistry_options
@click_job_options
@click.pass_context
def thermochemistry(
    ctx,
    directory,
    filetype,
    filenames,
    cutoff_entropy_grimme,
    cutoff_entropy_truhlar,
    cutoff_enthalpy,
    concentration,
    pressure,
    temperature,
    alpha,
    weighted,
    energy_units,
    outputfile,
    overwrite,
    check_imaginary_frequencies,
    skip_completed,
    **kwargs,
):
    """
    CLI for running thermochemistry jobs using the chemsmart framework.

    This command allows you to compute thermochemistry for Gaussian or ORCA
    output files.

    Examples:
    `chemsmart run thermochemistry -f udc3_mCF3_monomer_c9.log
    -f udc3_mCF3_monomer_c29.log -T 298.15`
    will save results to `udc3_mCF3_monomer_c9.dat` and
    `udc3_mCF3_monomer_c29.dat`.

    `chemsmart run thermochemistry -d /path/to/directory -t gaussian -T 298.15
    -o thermochemistry_results.dat`
    will compute thermochemistry for all Gaussian output files in the specified
    directory and save to `thermochemistry_results.dat`.
    """
    # validate input
    if directory and filenames:
        raise ValueError(
            "Cannot specify both --directory and --filenames. Choose one."
        )
    if directory and not filetype:
        raise ValueError("Must specify --filetype when using --directory.")
    if cutoff_entropy_grimme and cutoff_entropy_truhlar:
        raise ValueError(
            "Cannot specify both --cutoff-entropy-grimme and "
            "--cutoff-entropy-truhlar. Please choose one."
        )

    # choose entropy cutoff
    if cutoff_entropy_grimme is not None:
        cutoff_entropy = cutoff_entropy_grimme
        entropy_method = "grimme"
    elif cutoff_entropy_truhlar is not None:
        cutoff_entropy = cutoff_entropy_truhlar
        entropy_method = "truhlar"
    else:
        cutoff_entropy = None
        entropy_method = None

    # Create job settings
    job_settings = ThermochemistryJobSettings(
        temperature=temperature,
        concentration=concentration,
        pressure=pressure,
        use_weighted_mass=weighted,
        alpha=alpha,
        s_freq_cutoff=cutoff_entropy,
        entropy_method=entropy_method,
        h_freq_cutoff=cutoff_enthalpy,
        energy_units=energy_units,
        outputfile=outputfile,
        overwrite=overwrite,
        check_imaginary_frequencies=check_imaginary_frequencies,
    )

    # Initialize list to store jobs
    jobs = []
    files = []

    if directory:
        directory = os.path.abspath(directory)
        logger.info(
            f"Obtaining thermochemistry of files in directory: {directory}"
        )
        if filetype == "gaussian":
            from chemsmart.io.gaussian.folder import GaussianLogFolder

            folder = GaussianLogFolder(directory)
            logfiles = folder.all_logfiles
            files = [
                file for file in logfiles if outfile_format(file) == "gaussian"
            ]
        elif filetype == "orca":
            from chemsmart.io.orca.folder import ORCAOutFolder

            folder = ORCAOutFolder(directory)
            outfiles = folder.all_outfiles
            files = [
                file for file in outfiles if outfile_format(file) == "orca"
            ]
        else:
            raise ValueError(
                f"Unsupported filetype '{filetype}'. Use 'gaussian' or 'orca'."
            )
        for file in files:
            job = ThermochemistryJob.from_filename(
                filename=file,
                settings=job_settings,
                skip_completed=skip_completed,
            )
            jobs.append(job)
            logger.info(f"Created thermochemistry job for file: {file}")
            logger.debug(f"Job settings: {job_settings.__dict__}")

    elif filenames:
        for file in filenames:
            if outfile_format(file) not in {"gaussian", "orca"}:
                raise ValueError(
                    f"Unsupported output file type for '{file}'. Use Gaussian or "
                    f"ORCA output files."
                )
            job = ThermochemistryJob.from_filename(
                filename=file,
                settings=job_settings,
                skip_completed=skip_completed,
            )
            jobs.append(job)
            logger.info(f"Created thermochemistry job for file: {file}")

    logger.debug(f"Thermochemistry jobs created: {len(jobs)}")

    # Store objects in context
    ctx.obj["job_settings"] = job_settings
    ctx.obj["jobs"] = jobs
    ctx.obj["filenames"] = (
        filenames if filenames else files if directory else None
    )
    ctx.obj["directory"] = directory
    ctx.obj["filetype"] = filetype
    ctx.obj["outputfile"] = outputfile


@thermochemistry.result_callback()
@click.pass_context
def thermochemistry_process_pipeline(ctx, *args, **kwargs):
    """
    Process the thermochemistry jobs.
    """
    logger.debug(f"Context object: {ctx.obj}")
    logger.debug(f"args: {args}")
    logger.debug(f"kwargs: {kwargs}")

    jobs = ctx.obj.get("jobs", [])
    outputfile = ctx.obj.get("outputfile", None)
    logger.debug(f"Jobs to process:{jobs}")
    if ctx.invoked_subcommand is None:
        # If no subcommand is invoked, run the thermochemistry jobs
        logger.info("Running thermochemistry calculations on specified jobs.")
        for job in jobs:
            try:
                job.compute_thermochemistry()
                logger.info(
                    f"Thermochemistry calculation completed for {job.label}."
                )
            except Exception as e:
                logger.error(f"Error processing job for {job.label}: {e}")

        if outputfile is None:
            # If no output file is specified, save results to individual
            # files
            for job in jobs:
                job.show_results()
        else:
            # If output file is specified, save all results to this file
            jobs[0].show_results()
    else:
        # If a subcommand is invoked, just log the jobs
        logger.info(
            "Subcommand invoked. No thermochemistry calculations performed."
        )
        logger.debug(f"Jobs to process: {jobs}")
