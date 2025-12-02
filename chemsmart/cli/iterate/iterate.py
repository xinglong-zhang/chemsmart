import functools
import logging
import os
import yaml

import click

from chemsmart.cli.job import click_filename_options
from chemsmart.jobs.iterate.job import IterateJob
from chemsmart.jobs.iterate.runner import IterateJobRunner
from chemsmart.jobs.iterate.settings import IterateJobSettings
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.iterate import generate_template

logger = logging.getLogger(__name__)

# Define allowed keys for YAML configuration validation
ALLOWED_TOP_LEVEL_KEYS = {'skeletons', 'substituents'}
ALLOWED_SKELETON_KEYS = {'file_path', 'smiles', 'pubchem', 'label', 'link_index', 'skeleton_indices'}
ALLOWED_SUBSTITUENT_KEYS = {'file_path', 'smiles', 'pubchem', 'label', 'link_index'}


def click_iterate_options(f):
    """
    Common click options for Iterate.
    """

    @click.option(
        "-g",
        "--generate-template",
        "generate_template_path",
        is_flag=False,
        flag_value="iterate_template.yaml",
        default=None,
        type=str,
        help="Generate a template YAML configuration file and exit. "
             "Optionally specify output path (default: iterate_template.yaml).",
    )
    @click.option(
        "-P",
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
        "-a",
        "--algorithm",
        default="lagrange_multipliers",
        type=click.Choice(
            [
                "lagrange_multipliers", 
            ],
            case_sensitive=False,
        ),
        show_default=True,
        help="Algorithm to use for substituents' position optimization.",
    )
    @click.option(
        "-o",
        "--outputfile",
        default="iterate_out",
        type=str,
        show_default=True,
        help="Output filename (without .xyz extension) for generated structures.",
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
    algorithm,
    nprocs,
    timeout,
    outputfile,
    generate_template_path,
    **kwargs,
):
    """
    Here is a docstring made on the top of the hill in Kau Sai Chau

    CLI subcommand for running iterate jobs using the chemsmart framework.
    
    This command generates new molecular structures by attaching substituents
    to skeleton molecules at specified positions.
    
    Examples:
    
    `chemsmart run iterate -f config.yaml`
    will generate structures based on the YAML configuration file.
    
    `chemsmart run iterate -f config.yaml -P 4`
    will use 4 processes for parallel execution.
    
    `chemsmart run iterate -f config.yaml -o ./output`
    will save generated structures to ./output directory.
    
    `chemsmart run iterate -g`
    will generate a template YAML configuration file (iterate_template.yaml).
    
    `chemsmart run iterate -g my_config.yaml`
    will generate a template at the specified path.
    """
    # Handle -g option: generate template and exit
    if generate_template_path is not None:
        template_path = generate_template(generate_template_path, overwrite=False)
        click.echo(f"Generated template: {template_path}")
        ctx.exit(0)
    
    # Validate filename - expect a single YAML config file
    if not filename:
        raise click.BadParameter(
            "YAML configuration file is required.", 
            param_hint="'-f' / '--filename'"
        )
    
    if not os.path.exists(filename):
        raise click.BadParameter(
            f"File '{filename}' does not exist.", 
            param_hint="'-f' / '--filename'"
        )
    
    if not filename.endswith(('.yaml', '.yml')):
        raise click.BadParameter(
            f"File '{filename}' must be a YAML file (ending with .yaml or .yml).",
            param_hint="'-f' / '--filename'"
        )

    # Load YAML configuration file
    with open(filename, 'r') as f:
        raw_config = yaml.safe_load(f)
    
    # Handle empty YAML file
    if raw_config is None:
        raw_config = {}
    
    # Validate and normalize configuration
    config = validate_yaml_config(raw_config, filename)
    
    logger.info(f"Loaded configuration from '{filename}'")
    logger.info(f"  Skeletons: {len(config['skeletons'])}")
    logger.info(f"  Substituents: {len(config['substituents'])}")

    # Create job settings
    job_settings = IterateJobSettings(
        config_file=filename,
        algorithm=algorithm,
    )
    job_settings.skeleton_list = config['skeletons']
    job_settings.substituent_list = config['substituents']

    # Create job runner
    jobrunner = IterateJobRunner()

    # Create job
    job = IterateJob(
        settings=job_settings,
        jobrunner=jobrunner,
        nprocs=nprocs,
        timeout=timeout,
        outputfile=outputfile,
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
                logger.info(f"Successfully generated structures.")
                click.echo(f"Output saved to: {output_path}")
            else:
                logger.warning("No structures were generated.")
                
        except Exception as e:
            logger.error(f"Error running iterate job: {e}")
            raise click.ClickException(str(e))
    else:
        # If a subcommand is invoked, just log
        logger.info("Subcommand invoked. No iterate job executed.")


def _parse_index_string(value, entry_type: str, idx: int, field_name: str) -> list[int] | None:
    """
    Parse index string into a list of integers.
    
    Supports formats:
    - Single integer: 5 or "5"
    - Comma-separated: "1, 2, 3" or "1,2,3"
    - Ranges: "1-5" (expands to [1, 2, 3, 4, 5])
    - Mixed: "1, 3-5, 8" (expands to [1, 3, 4, 5, 8])
    
    Parameters
    ----------
    value : int, str, or None
        The value to parse
    entry_type : str
        "skeleton" or "substituent" (for error messages)
    idx : int
        Entry index (for error messages)
    field_name : str
        Field name like "link_index" or "skeleton_indices" (for error messages)
    
    Returns
    -------
    list[int] | None
        List of 1-based indices, or None if value is empty/null
    
    Raises
    ------
    click.BadParameter
        If the format is invalid
    """
    # Handle empty/null values
    if value is None or value == '' or value == 'null':
        return None
    
    # Handle single integer
    if isinstance(value, int):
        return [value]
    
    # Handle string
    if not isinstance(value, str):
        raise click.BadParameter(
            f"{entry_type.capitalize()} entry {idx + 1}: '{field_name}' must be an integer or string, "
            f"got {type(value).__name__}",
            param_hint="'-f' / '--filename'"
        )
    
    result = []
    # Split by comma and process each part
    parts = value.split(',')
    for part in parts:
        part = part.strip()
        if not part:
            continue
        
        # Check if it's a range (contains '-')
        if '-' in part:
            # Handle range like "3-7"
            range_parts = part.split('-')
            if len(range_parts) != 2:
                raise click.BadParameter(
                    f"{entry_type.capitalize()} entry {idx + 1}: Invalid range format '{part}' in '{field_name}'. "
                    f"Expected format like '3-7'.",
                    param_hint="'-f' / '--filename'"
                )
            try:
                start = int(range_parts[0].strip())
                end = int(range_parts[1].strip())
            except ValueError:
                raise click.BadParameter(
                    f"{entry_type.capitalize()} entry {idx + 1}: Invalid range '{part}' in '{field_name}'. "
                    f"Expected integers.",
                    param_hint="'-f' / '--filename'"
                )
            
            if start > end:
                raise click.BadParameter(
                    f"{entry_type.capitalize()} entry {idx + 1}: Invalid range '{part}' in '{field_name}'. "
                    f"Start ({start}) must be <= end ({end}).",
                    param_hint="'-f' / '--filename'"
                )
            
            result.extend(range(start, end + 1))
        else:
            # Handle single integer
            try:
                result.append(int(part))
            except ValueError:
                raise click.BadParameter(
                    f"{entry_type.capitalize()} entry {idx + 1}: Invalid value '{part}' in '{field_name}'. "
                    f"Expected an integer.",
                    param_hint="'-f' / '--filename'"
                )
    
    if not result:
        return None
    
    # Remove duplicates and sort
    return sorted(set(result))


def validate_yaml_config(config: dict, filename: str) -> dict:
    """
    Validate YAML configuration and normalize values.
    
    Parameters
    ----------
    config : dict
        Raw configuration dictionary from YAML
    filename : str
        Path to YAML file (for error messages)
    
    Returns
    -------
    dict
        Validated and normalized configuration
    
    Raises
    ------
    click.BadParameter
        If configuration contains invalid keys or structure
    """
    # Check top-level keys
    unknown_top_keys = set(config.keys()) - ALLOWED_TOP_LEVEL_KEYS
    if unknown_top_keys:
        raise click.BadParameter(
            f"Unknown top-level key(s) in YAML: {unknown_top_keys}. "
            f"Allowed keys: {ALLOWED_TOP_LEVEL_KEYS}",
            param_hint="'-f' / '--filename'"
        )
    
    validated_config = {}
    
    # Validate and normalize skeletons
    if 'skeletons' in config:
        skeletons = config['skeletons']
        if skeletons is None:
            validated_config['skeletons'] = []
        elif not isinstance(skeletons, list):
            raise click.BadParameter(
                f"'skeletons' must be a list, got {type(skeletons).__name__}",
                param_hint="'-f' / '--filename'"
            )
        else:
            validated_config['skeletons'] = [
                _validate_skeleton_entry(entry, idx, filename) 
                for idx, entry in enumerate(skeletons)
            ]
    else:
        validated_config['skeletons'] = []
    
    # Validate and normalize substituents
    if 'substituents' in config:
        substituents = config['substituents']
        if substituents is None:
            validated_config['substituents'] = []
        elif not isinstance(substituents, list):
            raise click.BadParameter(
                f"'substituents' must be a list, got {type(substituents).__name__}",
                param_hint="'-f' / '--filename'"
            )
        else:
            validated_config['substituents'] = [
                _validate_substituent_entry(entry, idx, filename)
                for idx, entry in enumerate(substituents)
            ]
    else:
        validated_config['substituents'] = []
    
    return validated_config


def _validate_skeleton_entry(entry: dict, idx: int, filename: str) -> dict:
    """
    Validate a single skeleton entry.
    
    Parameters
    ----------
    entry : dict
        Skeleton entry from YAML
    idx : int
        Index of the entry (for error messages)
    filename : str
        Path to YAML file (for error messages)
    
    Returns
    -------
    dict
        Validated and normalized skeleton entry
    """
    if entry is None:
        raise click.BadParameter(
            f"Skeleton entry {idx + 1} is empty/null",
            param_hint="'-f' / '--filename'"
        )
    
    if not isinstance(entry, dict):
        raise click.BadParameter(
            f"Skeleton entry {idx + 1} must be a dictionary, got {type(entry).__name__}",
            param_hint="'-f' / '--filename'"
        )
    
    # Check for unknown keys
    unknown_keys = set(entry.keys()) - ALLOWED_SKELETON_KEYS
    if unknown_keys:
        raise click.BadParameter(
            f"Unknown key(s) in skeleton entry {idx + 1}: {unknown_keys}. "
            f"Allowed keys: {ALLOWED_SKELETON_KEYS}",
            param_hint="'-f' / '--filename'"
        )
    
    # Normalize entry: convert empty/null values to None
    normalized = {}
    for key in ALLOWED_SKELETON_KEYS:
        value = entry.get(key)
        
        # Special handling for link_index and skeleton_indices: parse to list[int]
        if key == 'link_index':
            normalized[key] = _parse_index_string(value, 'skeleton', idx, 'link_index')
        elif key == 'skeleton_indices':
            normalized[key] = _parse_index_string(value, 'skeleton', idx, 'skeleton_indices')
        else:
            # Convert empty string, null, or missing to None
            if value is None or value == '' or value == 'null':
                normalized[key] = None
            else:
                normalized[key] = value
    
    return normalized


def _validate_substituent_entry(entry: dict, idx: int, filename: str) -> dict:
    """
    Validate a single substituent entry.
    
    Parameters
    ----------
    entry : dict
        Substituent entry from YAML
    idx : int
        Index of the entry (for error messages)
    filename : str
        Path to YAML file (for error messages)
    
    Returns
    -------
    dict
        Validated and normalized substituent entry
    """
    if entry is None:
        raise click.BadParameter(
            f"Substituent entry {idx + 1} is empty/null",
            param_hint="'-f' / '--filename'"
        )
    
    if not isinstance(entry, dict):
        raise click.BadParameter(
            f"Substituent entry {idx + 1} must be a dictionary, got {type(entry).__name__}",
            param_hint="'-f' / '--filename'"
        )
    
    # Check for unknown keys
    unknown_keys = set(entry.keys()) - ALLOWED_SUBSTITUENT_KEYS
    if unknown_keys:
        raise click.BadParameter(
            f"Unknown key(s) in substituent entry {idx + 1}: {unknown_keys}. "
            f"Allowed keys: {ALLOWED_SUBSTITUENT_KEYS}",
            param_hint="'-f' / '--filename'"
        )
    
    # Normalize entry: convert empty/null values to None
    normalized = {}
    for key in ALLOWED_SUBSTITUENT_KEYS:
        value = entry.get(key)
        
        # Special handling for link_index: parse to list[int]
        if key == 'link_index':
            normalized[key] = _parse_index_string(value, 'substituent', idx, 'link_index')
        else:
            # Convert empty string, null, or missing to None
            if value is None or value == '' or value == 'null':
                normalized[key] = None
            else:
                normalized[key] = value
    
    return normalized

