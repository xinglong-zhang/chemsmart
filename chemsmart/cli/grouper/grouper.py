"""
CLI command group for molecular structure grouping/clustering.

This module provides the main `grouper` command that serves as a parent
for various grouping strategy subcommands (irmsd, tfd, tanimoto, etc.).
"""

import functools
import logging
import os
import re

import click

from chemsmart.cli.job import (
    click_file_label_and_index_options,
    click_filenames_options,
    click_folder_options,
)
from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.io import clean_label

logger = logging.getLogger(__name__)


def click_grouper_common_options(f):
    """Common click options for all grouper subcommands."""

    @click.option(
        "-ih",
        "--ignore-hydrogens",
        is_flag=True,
        default=False,
        help="Whether to ignore hydrogens in grouping.",
    )
    @click.option(
        "-np",
        "--num-procs",
        type=int,
        default=1,
        help="Number of processors to use for the grouper.",
    )
    @click.option(
        "-T",
        "--threshold",
        type=float,
        default=None,
        help="Threshold for grouping. If not specified, uses strategy-specific "
        "defaults: RMSD/HRMSD/SpyRMSD/PyMOLRMSD=0.5, iRMSD=0.125, TFD=0.1, "
        "Tanimoto=0.9, Energy=1.0 (kcal/mol)",
    )
    @click.option(
        "-N",
        "--num-groups",
        type=int,
        default=None,
        help="Target number of groups to return. This uses adaptive threshold "
        "finding to return approximately this many unique structures.",
    )
    @click.option(
        "-o",
        "--output-format",
        type=click.Choice(["xlsx", "csv", "txt"], case_sensitive=False),
        default="xlsx",
        help="Output file format for results. Default is xlsx.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@click.group(name="grouper", cls=MyGroup)
@click_filenames_options
@click_file_label_and_index_options
@click_grouper_common_options
@click_folder_options
@click.pass_context
def grouper(
    ctx,
    filenames,
    label,
    append_label,
    ignore_hydrogens,
    num_procs,
    threshold,
    num_groups,
    output_format,
    directory,
    filetype,
    **kwargs,
):
    """
    Group/cluster molecular structures based on various similarity metrics.

    This command provides multiple grouping strategies as subcommands:
    - irmsd: Invariant RMSD (considers molecular symmetry)
    - tfd: Torsion Fingerprint Deviation
    - tanimoto: Tanimoto similarity with various fingerprints
    - rmsd: Simple Kabsch RMSD
    - hrmsd: Hungarian RMSD
    - spyrmsd: spyrmsd-based RMSD
    - pymolrmsd: PyMOL-based RMSD alignment

    Common options (-f, -T, -N, -np, -ih, -o) go BEFORE the subcommand.
    Strategy-specific options go AFTER the subcommand.

    Input modes:
    - Single file: -f conformers.xyz (all structures from one file)
    - Directory of output files: -d . -t gaussian (or -t orca)
      Loads last structure from each output file with conformer pattern.
      Files should have names like 'xxx_c1_opt.log', 'xxx_c12.out', etc.

    Examples:
        chemsmart run grouper -f conformers.xyz irmsd
        chemsmart run grouper -f conformers.xyz -T 0.5 irmsd --inversion on
        chemsmart run grouper -d . -t gaussian -T 0.5 irmsd
        chemsmart run grouper -d . -t orca -T 0.5 rmsd
        chemsmart run grouper -f conformers.xyz -T 0.1 tfd --use-weights
        chemsmart run grouper -f conformers.xyz -N 10 tanimoto --fingerprint-type morgan
        chemsmart run grouper -f conformers.xyz -o csv rmsd
    """
    ctx.ensure_object(dict)

    # Store common options in context for subcommands
    ctx.obj["ignore_hydrogens"] = ignore_hydrogens
    ctx.obj["num_procs"] = num_procs
    ctx.obj["threshold"] = threshold
    ctx.obj["num_groups"] = num_groups
    ctx.obj["output_format"] = output_format
    ctx.obj["conformer_ids"] = None  # Will be set only in directory mode

    # Validate input
    if filenames and directory:
        raise click.BadParameter(
            "Cannot specify both -f/--filenames and -d/--directory. Choose one."
        )

    if directory and not filetype:
        raise click.BadParameter(
            "Must specify -t/--filetype when using -d/--directory."
        )

    # Mode 1: Directory of output files (-d . -t gaussian/orca)
    if directory is not None and filetype is not None:
        if filetype.lower() not in ("gaussian", "orca"):
            raise click.BadParameter(
                f"For grouper with directory mode, only 'gaussian' or 'orca' "
                f"filetypes are supported. Got: {filetype}"
            )

        molecules, conformer_ids, auto_label = _load_molecules_from_directory(
            directory, filetype.lower()
        )
        grouper_label = _get_label(label, append_label, auto_label)
        ctx.obj["conformer_ids"] = conformer_ids
        logger.info(
            f"Loaded {len(molecules)} molecules from directory with conformer IDs: {conformer_ids}"
        )

    # Mode 2: Single file input (-f file.xyz)
    elif filenames:
        if len(filenames) != 1:
            raise click.BadParameter(
                "Currently only single file input is supported for grouper. "
                "Please provide one file with multiple structures."
            )

        filename = filenames[0]
        molecules = Molecule.from_filepath(
            filepath=filename, index=":", return_list=True
        )
        assert (
            molecules is not None
        ), f"Could not obtain molecules from {filename}!"
        logger.debug(f"Obtained {len(molecules)} molecules from {filename}")

        base_label = os.path.splitext(os.path.basename(filename))[0]
        grouper_label = _get_label(label, append_label, base_label)

    # No input specified
    else:
        logger.warning("[filename] or [directory] has not been specified!")
        ctx.obj["molecules"] = None
        ctx.obj["grouper_label"] = None
        return

    grouper_label = clean_label(grouper_label)

    if len(molecules) < 2:
        logger.warning("Less than 2 molecules loaded, grouping may not work.")

    ctx.obj["molecules"] = molecules
    ctx.obj["grouper_label"] = grouper_label

    logger.info(
        f"Loaded {len(molecules)} molecules for grouping with label: {grouper_label}"
    )


def _extract_conformer_id(filename: str) -> str | None:
    """
    Extract conformer ID from filename pattern like 'xxx_cXX_xxx.log' or 'xxx_cXX.log'.

    Args:
        filename: The filename to parse (e.g., 'structure1_c12_opt.log')

    Returns:
        Conformer ID like 'c12', or None if not found
    """
    basename = os.path.basename(filename)
    match = re.search(r"_c(\d+)[_.]", basename)
    return f"c{match.group(1)}" if match else None


def _load_molecules_from_directory(directory: str, filetype: str) -> tuple:
    """
    Load molecules from output files in a directory, extracting last structure from each.

    Supports both Gaussian and ORCA output files. Validation is done via
    Thermochemistry class:
    - Validates normal termination (via file_object)
    - Validates imaginary frequencies (via cleaned_frequencies)

    Energy extraction:
    - Gibbs free energy = SCF Done + Thermal correction to Gibbs Free Energy
    (Extracted directly from output file, no temperature-dependent recalculation)

    Args:
        directory: Path to directory containing output files
        filetype: Type of output files ('gaussian' or 'orca')

    Returns:
        tuple: (list of Molecule, list of conformer_ids, common_label)

    Raises:
        click.BadParameter: If no valid output files found or no molecules loaded
    """
    from chemsmart.analysis.thermochemistry import Thermochemistry
    from chemsmart.utils.io import find_output_files_in_directory

    directory = os.path.abspath(directory)
    output_files = find_output_files_in_directory(directory, filetype)

    if not output_files:
        raise click.BadParameter(
            f"No {filetype} output files found in directory: {directory}"
        )

    # Extract conformer info and sort
    file_info = []

    for f in output_files:
        conf_id = _extract_conformer_id(f)
        basename = os.path.basename(f)
        name_without_ext = os.path.splitext(basename)[0]

        if conf_id:
            # Has _cXX_ pattern: use cXX as ID, extract number for sorting
            num = int(re.search(r"\d+", conf_id).group())
            file_info.append((f, conf_id, num, name_without_ext))
        else:
            # No pattern: use filename as ID, use None for number
            file_info.append((f, name_without_ext, None, name_without_ext))

    if not file_info:
        raise click.BadParameter(
            f"No valid output files found in directory: {directory}"
        )

    # Sort: files with _cXX_ pattern by number, others by filename at the end
    def sort_key(x):
        if x[2] is not None:
            return (0, x[2], x[3])  # Has pattern: sort by number first
        else:
            return (
                1,
                0,
                x[3],
            )  # No pattern: sort by filename, after patterned files

    file_info.sort(key=sort_key)

    molecules = []
    conformer_ids = []

    for filepath, conf_id, _, _ in file_info:
        try:
            # Thermochemistry validates:
            # - normal_termination (via file_object)
            # - imaginary frequencies (via cleaned_frequencies in __init__)
            thermo = Thermochemistry(
                filename=filepath,
                temperature=298.15,  # Required by Thermochemistry but not used for energy
                check_imaginary_frequencies=True,
            )

            # Get molecule and set Gibbs energy
            mol = thermo.molecule
            if mol is not None:
                mol.name = conf_id

                # Extract Gibbs free energy from file_object
                # Both Gaussian16Output and ORCAOutput have gibbs_free_energy property
                gibbs_energy = thermo.file_object.gibbs_free_energy

                if gibbs_energy is not None:
                    mol._energy = gibbs_energy
                    logger.debug(
                        f"Loaded {conf_id} with Gibbs energy: {gibbs_energy:.6f} Hartree"
                    )
                else:
                    logger.debug(
                        f"Loaded {conf_id} with SCF energy: {mol.energy}"
                    )

                molecules.append(mol)
                conformer_ids.append(conf_id)
            else:
                logger.warning(f"Could not load molecule from {filepath}")

        except ValueError as e:
            # Thermochemistry raises ValueError for validation failures
            logger.warning(f"Skipping {conf_id}: {e}")
        except Exception as e:
            logger.warning(f"Error loading {filepath}: {e}")

    if not molecules:
        raise click.BadParameter(
            "No valid molecules could be loaded from log files"
        )

    # Extract common label from first filename (remove _cXX_ part)
    first_file = os.path.basename(file_info[0][0])
    label_match = re.match(r"(.+?)_c\d+", first_file)
    common_label = (
        label_match.group(1)
        if label_match
        else os.path.splitext(first_file)[0]
    )

    logger.info(
        f"Loaded {len(molecules)} valid molecules from {len(file_info)} log files"
    )

    return molecules, conformer_ids, common_label


def _get_label(label, append_label, base_label):
    """Determine final label from user inputs."""
    if label is not None and append_label is not None:
        raise ValueError("Only give label or append_label, but not both!")

    if append_label is not None:
        return f"{base_label}_{append_label}"
    elif label is not None:
        return label
    return base_label


def create_grouper_job_from_context(
    ctx,
    strategy: str,
    **extra_kwargs,
):
    """
    Create a GrouperJob from CLI context.

    This helper function extracts common parameters from the Click context
    and creates a GrouperJob with the specified strategy.

    Args:
        ctx: Click context object
        strategy: Grouping strategy name (e.g., 'irmsd', 'rmsd', 'tanimoto')
        **extra_kwargs: Strategy-specific arguments (e.g., inversion, fingerprint_type)

    Returns:
        GrouperJob: Configured grouper job instance

    Note:
        Default thresholds are set in each Grouper class, not here.
    """
    from chemsmart.jobs.grouper import GrouperJob

    molecules = ctx.obj["molecules"]
    ignore_hydrogens = ctx.obj["ignore_hydrogens"]
    num_procs = ctx.obj["num_procs"]
    label = ctx.obj["grouper_label"]
    num_groups = ctx.obj["num_groups"]
    conformer_ids = ctx.obj.get("conformer_ids")

    # Use threshold from parent command (None if not specified)
    threshold = ctx.obj["threshold"]
    output_format = ctx.obj.get("output_format", "xlsx")

    return GrouperJob(
        molecules=molecules,
        grouping_strategy=strategy,
        threshold=threshold,
        num_groups=num_groups,
        ignore_hydrogens=ignore_hydrogens,
        num_procs=num_procs,
        label=f"{label}_{strategy}",
        conformer_ids=conformer_ids,
        output_format=output_format,
        **extra_kwargs,
    )
