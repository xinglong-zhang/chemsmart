"""
Here is a docstring made on the top of the hill in Kau Sai Chau
"""

import functools
import logging
import os

import click
import tomlkit

from chemsmart.cli.job import click_filename_options
from chemsmart.jobs.iterate.job import IterateJob
from chemsmart.jobs.iterate.runner import IterateJobRunner
from chemsmart.jobs.iterate.settings import IterateJobSettings
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.iterate import generate_template, validate_config

logger = logging.getLogger(__name__)


def click_iterate_options(f):
    """
    Common click options for Iterate.
    """

    @click.option(
        "-g",
        "--generate-template",
        "generate_template_path",
        is_flag=False,
        flag_value="iterate_template.toml",
        default=None,
        type=str,
        help="Generate a template configuration file and exit. "
        "Optionally specify output path (default: iterate_template.toml).",
    )
    @click.option(
        "--separate-outputs/--no-separate-outputs",
        default=False,
        show_default=True,
        help="Save each structure as a separate XYZ file.",
    )
    @click.option(
        "-np",
        "--nprocs",
        default=1,
        type=click.IntRange(min=1),
        show_default=True,
        help="Number of processes for parallel execution.",
    )
    @click.option(
        "-t",
        "--timeout",
        default=120,
        type=click.IntRange(min=1),
        show_default=True,
        help="Timeout in seconds for each worker process.",
    )
    @click.option(
        "-m",
        "--method",
        default="lagrange_multipliers",
        type=click.Choice(
            [
                "lagrange_multipliers",
            ],
            case_sensitive=False,
        ),
        show_default=True,
        help="Mathematical method to use for substituents' position optimization.",
    )
    @click.option(
        "-s",
        "--sphere-direction-samples-number",
        "sphere_direction_samples_num",
        default=96,
        type=int,
        show_default=True,
        help="Number of points to sample on the unit sphere.",
    )
    @click.option(
        "-a",
        "--axial-rotations-sample-number",
        "axial_rotations_sample_num",
        default=6,
        type=int,
        show_default=True,
        help="Number of axial rotations per sphere point.",
    )
    @click.option(
        "-d",
        "--directory",
        default=None,
        type=click.Path(file_okay=False, dir_okay=True),
        help="Directory to save output files. Use only with --separate-outputs.",
    )
    @click.option(
        "-o",
        "--outputfile",
        default="iterate_out",
        type=str,
        show_default=True,
        help="Output filename (without .xyz extension) for generated structures. Use only with --no-separate-outputs.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


# use MyGroup to allow potential subcommands in the future
@click.group(cls=MyGroup, invoke_without_command=True)
@click_iterate_options
@click_filename_options
@click.pass_context
def iterate(
    ctx,
    filename,
    method,
    nprocs,
    timeout,
    outputfile,
    generate_template_path,
    sphere_direction_samples_num,
    axial_rotations_sample_num,
    directory,
    separate_outputs,
    **kwargs,
):
    """
    CLI subcommand for running iterate jobs using the chemsmart framework.

    This command generates new molecular structures by attaching substituents
    to skeleton molecules at specified positions.

    Examples:

    `chemsmart run iterate -f config.toml`
    will generate structures based on the configuration file.

    `chemsmart run iterate -f config.toml -np 4`
    will use 4 processes for parallel execution.

    `chemsmart run iterate -f config.toml -o ./output`
    will save generated structures to ./output directory.

    `chemsmart run iterate -g`
    will generate a template configuration file (iterate_template.toml).

    `chemsmart run iterate -g my_config.toml`
    will generate a template at the specified path.
    """
    # Handle -g option: generate template and exit
    if generate_template_path is not None:
        template_path = generate_template(
            generate_template_path, overwrite=False
        )
        click.echo(f"Generated template: {template_path}")
        ctx.exit(0)

    # Validate arguments based on separate_outputs flag
    source_output = ctx.get_parameter_source("outputfile")
    source_directory = ctx.get_parameter_source("directory")

    if separate_outputs:
        # If writing separate files, -o is forbidden
        if source_output == click.core.ParameterSource.COMMANDLINE:
            raise click.UsageError(
                "Option '-o' / '--outputfile' is not allowed when '--separate-outputs' "
                "is enabled. Please use '-d' / '--directory' to specify the output location."
            )
        # Set default directory if not provided
        if directory is None:
            directory = os.getcwd()
    else:
        # If writing single file (default), -d is forbidden
        if source_directory == click.core.ParameterSource.COMMANDLINE:
            raise click.UsageError(
                "Option '-d' / '--directory' is not allowed when '--no-separate-outputs' "
                "(default) is active. Please use '-o' / '--outputfile' to specify the output file."
            )

    # Validate filename - expect a single TOML config file
    if not filename:
        raise click.BadParameter(
            "A configuration file is required.",
            param_hint="'-f' / '--filename'",
        )

    if not os.path.exists(filename):
        raise click.BadParameter(
            f"File '{filename}' does not exist.",
            param_hint="'-f' / '--filename'",
        )

    if not filename.endswith(".toml"):
        raise click.BadParameter(
            f"File '{filename}' must be a configuration file (ending with .toml).",
            param_hint="'-f' / '--filename'",
        )

    # Load TOML configuration file
    with open(filename, "r") as f:
        # Load and unwrap to get standard Python dict/list types
        raw_config = tomlkit.load(f).unwrap()

    # Handle empty configuration file
    if raw_config is None:
        raw_config = {}

    # Validate and normalize configuration
    config = validate_config(raw_config, filename)

    logger.info(f"Loaded configuration from '{filename}'")
    logger.info(f"  Skeletons: {len(config['skeletons'])}")
    logger.info(f"  Substituents: {len(config['substituents'])}")

    # Create job settings
    job_settings = IterateJobSettings(
        config_file=filename,
        method=method,
        sphere_direction_samples_num=sphere_direction_samples_num,
        axial_rotations_sample_num=axial_rotations_sample_num,
    )
    job_settings.skeleton_list = config["skeletons"]
    job_settings.substituent_list = config["substituents"]

    # Create job runner
    jobrunner = IterateJobRunner()

    # Create job
    job = IterateJob(
        settings=job_settings,
        jobrunner=jobrunner,
        nprocs=nprocs,
        timeout=timeout,
        outputfile=outputfile,
        separate_outputs=separate_outputs,
        output_directory=directory,
    )

    logger.debug(f"Created IterateJob with {nprocs} process(es)")

    # Store objects in context
    ctx.ensure_object(dict)
    ctx.obj["job"] = job
    ctx.obj["job_settings"] = job_settings
    ctx.obj["filename"] = filename


@iterate.result_callback()
@click.pass_context
def iterate_process_pipeline(ctx, *args, **kwargs):
    """
    Process the iterate job after command execution.
    """
    logger.debug(f"Context object: {ctx.obj}")

    job = ctx.obj.get("job")

    if ctx.invoked_subcommand is None and job is not None:
        # If no subcommand is invoked, run the iterate job
        logger.info("Running iterate job to generate molecular structures.")

        try:
            # Run the job - this already writes results to outputfile
            output_path = job.run()

            # Output path is returned from job.run()
            if output_path and os.path.exists(output_path):
                logger.info("Successfully generated structures.")
                click.echo(f"Output saved to: {output_path}")
            else:
                logger.warning("No structures were generated.")

        except Exception as e:
            logger.error(f"Error running iterate job: {e}")
            raise click.ClickException(str(e))
    else:
        # If a subcommand is invoked, just log
        logger.info("Subcommand invoked. No iterate job executed.")
