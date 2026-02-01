"""
CLI command group for molecular structure grouping/clustering.

This module provides the main `grouper` command that serves as a parent
for various grouping strategy subcommands (irmsd, tfd, tanimoto, etc.).
"""

import functools
import glob
import logging
import os
import re

import click

from chemsmart.cli.job import (
    click_file_label_options,
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
        "defaults: RMSD=0.5, Tanimoto=0.9, TFD=0.1, Connectivity=0.0.",
    )
    @click.option(
        "-N",
        "--num-groups",
        type=int,
        default=None,
        help="Target number of groups to return. This uses adaptive threshold "
        "finding to return approximately this many unique structures.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@click.group(name="grouper", cls=MyGroup)
@click_filenames_options
@click_file_label_options
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
    directory,
    filetype,
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

    Common options (-f, -T, -N, -np, -ih) go BEFORE the subcommand.
    Strategy-specific options go AFTER the subcommand.

    Input modes:
    - Single file: -f conformers.xyz (all structures from one file)
    - Directory of logs: -d . -t log (last structure from each log file)
      Files should have names like 'xxx_c1_opt.log', 'xxx_c12.log', etc.

    Examples:
        chemsmart run grouper -f conformers.xyz irmsd
        chemsmart run grouper -f conformers.xyz -T 0.5 irmsd --check-stereo on
        chemsmart run grouper -d . -t log -T 0.5 irmsd
        chemsmart run grouper -f conformers.xyz -T 0.1 tfd --use-weights
        chemsmart run grouper -f conformers.xyz -N 10 tanimoto --fingerprint-type morgan
    """
    ctx.ensure_object(dict)

    # Store common options in context for subcommands
    ctx.obj["ignore_hydrogens"] = ignore_hydrogens
    ctx.obj["num_procs"] = num_procs
    ctx.obj["threshold"] = threshold
    ctx.obj["num_groups"] = num_groups
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

    # Mode 1: Directory of log files (-d . -t log)
    if directory is not None and filetype is not None:
        if filetype.lower() != "log":
            raise click.BadParameter(
                f"For grouper with directory mode, only 'log' filetype is supported. "
                f"Got: {filetype}"
            )

        molecules, conformer_ids, auto_label = (
            _load_molecules_from_log_directory(directory)
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


def _load_molecules_from_log_directory(directory: str) -> tuple:
    """
    Load molecules from log files in a directory, extracting last structure from each.

    Only accepts opt or ts job types with frequency calculations:
    - opt: must have normal termination and no imaginary frequencies
    - ts: must have normal termination and exactly one imaginary frequency

    For valid files, the Gibbs free energy (SCF Done + Thermal correction to
    Gibbs Free Energy) is extracted and stored in mol.energy.

    Args:
        directory: Path to directory containing log files

    Returns:
        tuple: (list of Molecule, list of conformer_ids, common_label)

    Raises:
        click.BadParameter: If no valid log files found or no molecules loaded
    """
    from chemsmart.io.gaussian.output import Gaussian16Output

    directory = os.path.abspath(directory)
    log_files = glob.glob(os.path.join(directory, "*.log"))

    if not log_files:
        raise click.BadParameter(
            f"No .log files found in directory: {directory}"
        )

    # Extract conformer info and sort by conformer number
    file_info = []
    for f in log_files:
        conf_id = _extract_conformer_id(f)
        if conf_id:
            num = int(re.search(r"\d+", conf_id).group())
            file_info.append((f, conf_id, num))
        else:
            logger.warning(
                f"Could not extract conformer ID from {f}, skipping"
            )

    if not file_info:
        raise click.BadParameter(
            "No files with conformer pattern (_cXX_) found in directory. "
            "Expected filenames like 'structure_c1_opt.log' or 'mol_c12.log'"
        )

    # Sort by conformer number
    file_info.sort(key=lambda x: x[2])

    molecules = []
    conformer_ids = []

    for filepath, conf_id, _ in file_info:
        try:
            g16_output = Gaussian16Output(filename=filepath)

            # Validate the file
            is_valid, error_msg = _validate_gaussian_output(
                g16_output, conf_id
            )
            if not is_valid:
                logger.warning(f"Skipping {conf_id}: {error_msg}")
                continue

            mol = Molecule.from_filepath(
                filepath=filepath, index="-1", return_list=False
            )
            if mol is not None:
                mol.name = conf_id

                # Extract Gibbs free energy
                gibbs_energy = _extract_gibbs_energy(g16_output)
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


def _validate_gaussian_output(g16_output, conf_id: str) -> tuple[bool, str]:
    """
    Validate Gaussian output file for grouper usage.

    Requirements:
    - Must be opt or ts job type
    - Must have normal termination
    - opt: no imaginary frequencies
    - ts: exactly one imaginary frequency

    Args:
        g16_output: Gaussian16Output object
        conf_id: Conformer ID for error messages

    Returns:
        tuple: (is_valid, error_message)
    """
    # Check normal termination
    if not g16_output.normal_termination:
        return False, "abnormal termination"

    # Check job type
    jobtype = g16_output.jobtype
    if jobtype not in ("opt", "ts"):
        return (
            False,
            f"unsupported job type '{jobtype}' (only 'opt' or 'ts' allowed)",
        )

    # Check for frequency calculation (needed for Gibbs energy)
    freqs = g16_output.vibrational_frequencies
    if freqs is None or len(freqs) == 0:
        return False, "no frequency calculation found"

    # Count imaginary frequencies
    imaginary_freqs = [f for f in freqs if f < 0]
    num_imaginary = len(imaginary_freqs)

    if jobtype == "opt":
        if num_imaginary > 0:
            return (
                False,
                f"opt job has {num_imaginary} imaginary frequency(s), expected 0",
            )
    elif jobtype == "ts":
        if num_imaginary != 1:
            return (
                False,
                f"ts job has {num_imaginary} imaginary frequency(s), expected 1",
            )

    return True, ""


def _extract_gibbs_energy(g16_output) -> float | None:
    """
    Extract Gibbs free energy from Gaussian output.

    Gibbs free energy = SCF Done + Thermal correction to Gibbs Free Energy

    Args:
        g16_output: Gaussian16Output object

    Returns:
        Gibbs free energy in Hartree, or None if not available
    """
    # Get the last SCF energy
    if not g16_output.scf_energies:
        return None
    scf_energy = g16_output.scf_energies[-1]

    # Search for thermal correction to Gibbs free energy
    thermal_correction = None
    for line in g16_output.contents:
        if "Thermal correction to Gibbs Free Energy=" in line:
            try:
                thermal_correction = float(line.split()[-1])
                break
            except (ValueError, IndexError):
                pass

    if thermal_correction is None:
        return None

    return scf_energy + thermal_correction


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
    default_threshold: float = None,
    **extra_kwargs,
):
    """
    Create a GrouperJob from CLI context.

    This helper function extracts common parameters from the Click context
    and creates a GrouperJob with the specified strategy.

    Args:
        ctx: Click context object
        strategy: Grouping strategy name (e.g., 'irmsd', 'rmsd', 'tanimoto')
        default_threshold: Default threshold if not specified by user
        **extra_kwargs: Strategy-specific arguments (e.g., check_stereo, fingerprint_type)

    Returns:
        GrouperJob: Configured grouper job instance
    """
    from chemsmart.jobs.grouper import GrouperJob

    molecules = ctx.obj["molecules"]
    ignore_hydrogens = ctx.obj["ignore_hydrogens"]
    num_procs = ctx.obj["num_procs"]
    label = ctx.obj["grouper_label"]
    num_groups = ctx.obj["num_groups"]
    conformer_ids = ctx.obj.get("conformer_ids")

    # Use threshold from parent command, with strategy-specific default
    threshold = ctx.obj["threshold"]
    if (
        threshold is None
        and num_groups is None
        and default_threshold is not None
    ):
        threshold = default_threshold

    return GrouperJob(
        molecules=molecules,
        grouping_strategy=strategy,
        threshold=threshold,
        num_groups=num_groups,
        ignore_hydrogens=ignore_hydrogens,
        num_procs=num_procs,
        label=f"{label}_{strategy}",
        conformer_ids=conformer_ids,
        **extra_kwargs,
    )
