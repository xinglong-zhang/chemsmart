"""
CLI subcommand for running iterate jobs from YAML configuration files.
"""

import functools
import logging
import os

import click
import yaml

from chemsmart.cli.job import click_filename_options
from chemsmart.jobs.iterate.job import IterateJob
from chemsmart.jobs.iterate.runner import IterateJobRunner
from chemsmart.jobs.iterate.settings import IterateJobSettings
from chemsmart.utils.iterate import generate_yaml_template, validate_yaml_config

logger = logging.getLogger(__name__)


def click_yaml_iterate_options(f):
    """
    Click options for YAML-based iterate jobs.
    """

    @click.option(
        "-g",
        "--generate-template",
        "generate_template_path",
        is_flag=False,
        flag_value="iterate_template.yaml",
        default=None,
        type=str,
        help="Generate a template configuration file and exit. "
        "Optionally specify output path (default: iterate_template.yaml).",
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
        "-cm",
        "--combination-mode",
        "combination_mode",
        default="independent",
        type=click.Choice(
            ["independent", "global"],
            case_sensitive=False,
        ),
        show_default=True,
        help="Combination strategy for skeleton slots. "
        "Each slot specifies a group number; only substituents belonging to "
        "that group are candidates for the slot. "
        "Each position includes a 'None' (keep original) option. "
        "'independent' (default): each slot is expanded separately. "
        "E.g. slot R1 (group 1), slot R2 (group 2); "
        "sub A, B in group 1, sub C in group 2: "
        "R1 → mol(R1=A), mol(R1=B); R2 → mol(R2=C) → 3 structures. "
        "'global': all slots combined via single Cartesian product. "
        "Same setup → mol(R1=A), mol(R1=B), mol(R2=C), "
        "mol(R1=A,R2=C), mol(R1=B,R2=C) → 5 structures.",
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
        help="Output filename (without .xyz extension) for generated structures. "
        "Use only with --no-separate-outputs.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@click.command(name="yaml")
@click_yaml_iterate_options
@click_filename_options
@click.pass_context
def yaml_cmd(
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
    combination_mode,
    **kwargs,
):
    """
    Run iterate jobs from a YAML configuration file.

    All skeletons participate in global contiguous group numbering.
    Skeletons with link_index occupy one implicit group each;
    skeletons with slots occupy one group per slot.

    Examples:

    \b
    chemsmart run iterate yaml -f config.yaml
    chemsmart run iterate yaml -f config.yaml -np 4
    chemsmart run iterate yaml -f config.yaml -o my_output
    chemsmart run iterate yaml -g
    chemsmart run iterate yaml -g my_config.yaml
    """
    # Handle -g option: generate template and exit
    if generate_template_path is not None:
        template_path = generate_yaml_template(
            generate_template_path, overwrite=False
        )
        click.echo(f"Generated template: {template_path}")
        ctx.exit(0)

    # Validate arguments based on separate_outputs flag
    source_output = ctx.get_parameter_source("outputfile")
    source_directory = ctx.get_parameter_source("directory")

    if separate_outputs:
        if source_output == click.core.ParameterSource.COMMANDLINE:
            raise click.UsageError(
                "Option '-o' / '--outputfile' is not allowed when '--separate-outputs' "
                "is enabled. Please use '-d' / '--directory' to specify the output location."
            )
        if directory is None:
            directory = os.getcwd()
    else:
        if source_directory == click.core.ParameterSource.COMMANDLINE:
            raise click.UsageError(
                "Option '-d' / '--directory' is not allowed when '--no-separate-outputs' "
                "(default) is active. Please use '-o' / '--outputfile' to specify the output file."
            )

    # Validate filename
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

    if not filename.endswith((".yaml", ".yml")):
        raise click.BadParameter(
            f"File '{filename}' must be a YAML file "
            f"(ending with .yaml or .yml).",
            param_hint="'-f' / '--filename'",
        )

    # Load YAML configuration file
    with open(filename, "r") as f:
        raw_config = yaml.safe_load(f)

    if raw_config is None:
        raw_config = {}

    # Validate and normalize configuration
    config = validate_yaml_config(raw_config, filename)

    logger.info(f"Loaded YAML configuration from '{filename}'")
    logger.info(f"  Skeletons: {len(config['skeletons'])}")
    logger.info(f"  Substituents: {len(config['substituents'])}")

    logger.info(f"  Combination mode: {combination_mode}")

    # Create job settings
    job_settings = IterateJobSettings(
        config_file=filename,
        method=method,
        sphere_direction_samples_num=sphere_direction_samples_num,
        axial_rotations_sample_num=axial_rotations_sample_num,
        combination_mode=combination_mode,
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

    # Run the job
    logger.info("Running iterate job to generate molecular structures.")

    try:
        output_path = job.run()

        if output_path and os.path.exists(output_path):
            logger.info("Successfully generated structures.")
            click.echo(f"Output saved to: {output_path}")
        else:
            logger.warning("No structures were generated.")

    except Exception as e:
        logger.error(f"Error running iterate job: {e}")
        raise click.ClickException(str(e))
