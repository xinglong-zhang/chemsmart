import logging

import click

from chemsmart.cli.job import (
    click_file_label_and_index_options,
    click_filenames_options,
    click_folder_options,
    click_job_options,
)
from chemsmart.io.folder import BaseFolder
from chemsmart.jobs.thermochemistry.job import ThermochemistryJob
from chemsmart.jobs.thermochemistry.settings import ThermochemistryJobSettings
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.io import (
    check_program_availability_in_chemsmart,
    get_program_type_from_file,
)

logger = logging.getLogger(__name__)


def thermochemistry_cutoff_options(
    f,
    entropy_grimme_default=None,
    entropy_truhlar_default=None,
    enthalpy_default=None,
):
    """Reusable quasi-RRHO cutoff options."""
    f = click.option(
        "-csg",
        "--cutoff-entropy-grimme",
        default=entropy_grimme_default,
        type=float,
        show_default=True,
        help="Cutoff frequency for entropy in wavenumbers, using Grimme's "
        "quasi-RRHO method.",
    )(f)
    f = click.option(
        "-cst",
        "--cutoff-entropy-truhlar",
        default=entropy_truhlar_default,
        type=float,
        show_default=True,
        help="Cutoff frequency for entropy in wavenumbers, using Truhlar's "
        "quasi-RRHO method.",
    )(f)
    f = click.option(
        "-ch",
        "--cutoff-enthalpy",
        default=enthalpy_default,
        type=float,
        show_default=True,
        help="Cutoff frequency for enthalpy in wavenumbers, using "
        "Head-Gordon's quasi-RRHO method.",
    )(f)
    return f


def resolve_entropy_cutoff(cutoff_entropy_grimme, cutoff_entropy_truhlar):
    """Resolve entropy cutoff and method from CLI options."""
    if (
        cutoff_entropy_grimme is not None
        and cutoff_entropy_truhlar is not None
    ):
        raise ValueError(
            "Cannot specify both --cutoff-entropy-grimme and "
            "--cutoff-entropy-truhlar. Please choose one."
        )
    if cutoff_entropy_truhlar is not None:
        return cutoff_entropy_truhlar, "truhlar"
    if cutoff_entropy_grimme is not None:
        return cutoff_entropy_grimme, "grimme"
    return None, None


def thermochemistry_temp_pressure_conc_options(
    f,
    temperature_required=True,
    temperature_default=None,
    concentration_default=None,
    pressure_default=1.0,
    include_pressure=True,
    concentration_short="-c",
):
    """Reusable temperature, pressure, and concentration options."""
    f = click.option(
        concentration_short,
        "--concentration",
        default=concentration_default,
        type=float,
        show_default=True,
        help="Concentration in mol/L.",
    )(f)
    if include_pressure:
        f = click.option(
            "-P",
            "--pressure",
            default=pressure_default,
            type=float,
            show_default=True,
            help="Pressure in atm.",
        )(f)
    f = click.option(
        "-T",
        "--temperature",
        required=temperature_required,
        default=temperature_default,
        type=float,
        help="Temperature in Kelvin.",
    )(f)
    return f


def click_thermochemistry_options(f):
    """
    Common click options for Thermochemistry.
    """
    f = thermochemistry_temp_pressure_conc_options(f)
    f = thermochemistry_cutoff_options(f)
    f = click.option(
        "-a",
        "--alpha",
        default=4,
        type=int,
        show_default=True,
        help="Interpolator exponent used in the quasi-RRHO approximation.",
    )(f)
    f = click.option(
        "-w/",
        "--weighted/--no-weighted",
        default=True,
        show_default=True,
        help="Use natural abundance weighted masses (True) or use most abundant "
        "masses (False, via --no-weighted).\nDefault to True, i.e., use natural "
        "abundance weighted masses, which is the real world scenario.",
    )(f)
    f = click.option(
        "-u",
        "--energy-units",
        default="hartree",
        show_default=True,
        type=click.Choice(
            ["hartree", "eV", "kcal/mol", "kJ/mol"], case_sensitive=False
        ),
        help="Units of energetic values.",
    )(f)
    f = click.option(
        "-o",
        "--outputfile",
        default=None,
        type=str,
        help="Output file to save the thermochemistry results. Defaults to "
        "None, which will save results to file_basename.dat.\nIf "
        "specified, it will save all thermochemistry results to this file.",
    )(f)
    f = click.option(
        "-O",
        "--overwrite",
        is_flag=True,
        default=False,
        show_default=True,
        help="Overwrite existing output files if they already exist.",
    )(f)
    f = click.option(
        "-i/",
        "--check-imaginary-frequencies/--no-check-imaginary-frequencies",
        default=True,
        show_default=True,
        help="Check for imaginary frequencies in the calculations.",
    )(f)
    return f


# use MyGroup to allow potential subcommands in the future
@click.group(cls=MyGroup, invoke_without_command=True)
@click_thermochemistry_options
@click_job_options
@click_folder_options
@click_filenames_options
@click_file_label_and_index_options
@click.pass_context
def thermochemistry(
    ctx,
    directory,
    filetype,
    program,
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
    CLI subcommand for running thermochemistry
    jobs using the chemsmart framework.

    This command allows you to compute thermochemistry for Gaussian or ORCA
    output files.

    Examples:
    `chemsmart run thermochemistry -f udc3_mCF3_monomer_c9.log
    -f udc3_mCF3_monomer_c29.log -T 298.15`
    will save results to `udc3_mCF3_monomer_c9.dat` and
    `udc3_mCF3_monomer_c29.dat`.

    `chemsmart run thermochemistry -d /path/to/directory -p gaussian -T 298.15
    -o thermochemistry_results.dat`
    will compute thermochemistry for all Gaussian output files in the specified
    directory and save to `thermochemistry_results.dat`.
    """
    # validate input
    if directory and filenames:
        raise ValueError(
            "Cannot specify both --directory and --filenames. Choose one."
        )
    if directory and not program and not filetype:
        raise ValueError(
            "Must specify --program or --filetype when using --directory."
        )
    if program:
        check_program_availability_in_chemsmart(program)

    cutoff_entropy, entropy_method = resolve_entropy_cutoff(
        cutoff_entropy_grimme, cutoff_entropy_truhlar
    )

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
        if program and not filetype:
            # obtain all output files belonging to a program
            files = BaseFolder(
                folder=directory
            ).get_all_output_files_in_current_folder_by_program(
                program=program.lower()
            )
        elif filetype and not program:
            # obtain all files of a specific type, regardless of program
            files = BaseFolder(
                folder=directory
            ).get_all_files_in_current_folder_by_suffix(filetype=filetype)
        elif program and filetype:
            files = BaseFolder(
                folder=directory
            ).get_all_files_in_current_folder_by_program_and_suffix(
                program=program.lower(), filetype=filetype
            )
        else:
            raise ValueError(
                "Must specify either --program or --filetype when using --directory."
            )

        files = sorted(files)
        for file in files:
            job = ThermochemistryJob.from_filename(
                filename=file,
                settings=job_settings,
                skip_completed=skip_completed,
            )
            if outputfile is not None:
                job_settings.overwrite = False
                job_settings.write_header = False
            jobs.append(job)
            logger.info(f"Created thermochemistry job for file: {file}")
            logger.debug(f"Job settings: {job_settings.__dict__}")

    elif filenames:
        for file in filenames:
            if get_program_type_from_file(file) not in {"gaussian", "orca"}:
                raise ValueError(
                    f"Unsupported output file type for '{file}'. Use Gaussian or "
                    f"ORCA output files."
                )
            job = ThermochemistryJob.from_filename(
                filename=file,
                settings=job_settings,
                skip_completed=skip_completed,
            )
            if outputfile is not None:
                job_settings.overwrite = False
                job_settings.write_header = False
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
    ctx.obj["program"] = program
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
