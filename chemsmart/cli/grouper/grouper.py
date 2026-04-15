"""
CLI command group for molecular structure grouping/clustering.

This module provides the main `grouper` command that serves as a parent
for various grouping strategy subcommands (irmsd, tfd, tanimoto, etc.).
"""

import functools
import logging
import os
import re
from pathlib import Path

import click

from chemsmart.cli.job import (
    click_file_label_and_index_options,
    click_filenames_options,
    click_folder_options,
)
from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.constants import energy_conversion
from chemsmart.utils.io import clean_label, get_program_type_from_file

logger = logging.getLogger(__name__)


THERMOCHEMISTRY_DEFAULT_KWARGS = {
    "temperature": None,
    "concentration": None,
    "pressure": None,
    "use_weighted_mass": True,
    "alpha": None,
    "s_freq_cutoff": None,
    "entropy_method": None,
    "h_freq_cutoff": None,
    "energy_units": "hartree",
    "check_imaginary_frequencies": True,
}


def build_thermochemistry_kwargs(
    *,
    cutoff_entropy_grimme=None,
    cutoff_entropy_truhlar=None,
    cutoff_enthalpy=None,
    concentration=None,
    pressure=None,
    temperature=None,
    alpha=None,
    weighted=None,
    energy_units=None,
    check_imaginary_frequencies=None,
) -> dict:
    """Build Thermochemistry kwargs from defaults plus CLI overrides."""
    if (
        cutoff_entropy_grimme is not None
        and cutoff_entropy_truhlar is not None
    ):
        raise ValueError(
            "Cannot specify both --cutoff-entropy-grimme and "
            "--cutoff-entropy-truhlar. Please choose one."
        )

    kwargs = dict(THERMOCHEMISTRY_DEFAULT_KWARGS)

    if cutoff_entropy_grimme is not None:
        kwargs["s_freq_cutoff"] = cutoff_entropy_grimme
        kwargs["entropy_method"] = "grimme"
    elif cutoff_entropy_truhlar is not None:
        kwargs["s_freq_cutoff"] = cutoff_entropy_truhlar
        kwargs["entropy_method"] = "truhlar"

    overrides = {
        "h_freq_cutoff": cutoff_enthalpy,
        "concentration": concentration,
        "pressure": pressure,
        "temperature": temperature,
        "alpha": alpha,
        "use_weighted_mass": weighted,
        "energy_units": energy_units,
        "check_imaginary_frequencies": check_imaginary_frequencies,
        "cutoff_entropy_grimme": cutoff_entropy_grimme,
        "cutoff_entropy_truhlar": cutoff_entropy_truhlar,
        "cutoff_enthalpy": cutoff_enthalpy,
    }
    for key, value in overrides.items():
        if value is not None:
            kwargs[key] = value

    return kwargs


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
        "-m",
        "--matrix-format",
        type=click.Choice(["xlsx", "csv", "txt"], case_sensitive=False),
        default="xlsx",
        help="Output file format for results. Default is xlsx.",
    )
    @click.option(
        "-E",
        "--energy-type",
        type=click.Choice(
            ["E", "H", "G", "qhH", "qhG", "sp_qhG"], case_sensitive=False
        ),
        default="E",
        help="Energy type to use for every grouping. Options: E (SCF energy, default), H (enthalpy), "
        "G (Gibbs free energy), qhH (quasi-harmonic enthalpy), qhG (quasi-harmonic Gibbs free energy), "
        "sp_qhG (single-point corrected quasi-harmonic Gibbs free energy)."
        "Note: not for xyz file. "
        "When using qhH, qhG and sp_qhG, chemsmart will call thermochemistry for thermo correction, "
        "Users can specify the thermo-correction parameter.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_grouper_thermochemistry_options(f):
    """Thermochemistry options for directory-based grouper energy extraction."""

    @click.option(
        "-csg",
        "--cutoff-entropy-grimme",
        default=100,
        type=float,
        show_default=True,
        help="For-Thermo-Correction: Cutoff frequency for entropy in wavenumbers, using Grimme's "
        "quasi-RRHO method (Default to 100).",
    )
    @click.option(
        "-cst",
        "--cutoff-entropy-truhlar",
        default=None,
        type=float,
        show_default=True,
        help="For-Thermo-Correction: Cutoff frequency for entropy in wavenumbers, using Truhlar's "
        "quasi-RRHO method.",
    )
    @click.option(
        "-ch",
        "--cutoff-enthalpy",
        default=100,
        type=float,
        show_default=True,
        help="For-Thermo-Correction: Cutoff frequency for enthalpy (cm^-1), Head-Gordon qRRHO method (Default to 100).",
    )
    @click.option(
        "-c",
        "--concentration",
        default=1.0,
        type=float,
        show_default=True,
        help="For-Thermo-Correction: Concentration in mol/L (default to 1.0).",
    )
    @click.option(
        "-P",
        "--pressure",
        default=1.0,
        show_default=True,
        type=float,
        help="For-Thermo-Correction: Pressure in atm (default to 1.0).",
    )
    @click.option(
        "-temp",
        "--temperature",
        default=298.15,
        show_default=True,
        type=float,
        help="For-Thermo-Correction: Temperature in Kelvin (default to 298.15).",
    )
    @click.option(
        "--alpha",
        default=4,
        show_default=True,
        type=int,
        help="For-Thermo-Correction: Interpolator exponent used in the quasi-RRHO approximation (default to 4).",
    )
    @click.option(
        "--weighted/--no-weighted",
        default=True,
        show_default=True,
        help="For-Thermo-Correction: Use natural abundance weighted masses (True) or use most abundant "
        "masses (False, via --no-weighted).\nDefault to True, i.e., use natural "
        "abundance weighted masses, which is the real world scenario.",
    )
    @click.option(
        "--energy-units",
        default="hartree",
        show_default=True,
        type=click.Choice(
            ["hartree", "eV", "kcal/mol", "kJ/mol"], case_sensitive=False
        ),
        help="For-Thermo-Correction: Units of energetic values.",
    )
    @click.option(
        "--check-imaginary-frequencies/--no-check-imaginary-frequencies",
        default=True,
        show_default=True,
        help="For-Thermo-Correction: Whether to check for imaginary frequencies.",
    )
    @functools.wraps(f)
    def wrapper_thermochemistry_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_thermochemistry_options


@click.group(name="grouper", cls=MyGroup)
@click_filenames_options
@click_file_label_and_index_options
@click_grouper_common_options
@click_grouper_thermochemistry_options
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
    matrix_format,
    directory,
    filetype,
    program,
    energy_type,
    cutoff_entropy_grimme,
    cutoff_entropy_truhlar,
    cutoff_enthalpy,
    concentration,
    pressure,
    temperature,
    alpha,
    weighted,
    energy_units,
    check_imaginary_frequencies,
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
    - Directory of output files: -d . -p gaussian (or -p orca)
      Loads last structure from each output file with conformer pattern.
      Files should have names like 'xxx_c1_opt.log', 'xxx_c12.out', etc.

    Examples:
        chemsmart run grouper -f conformers.xyz irmsd
        chemsmart run grouper -f conformers.xyz -T 0.5 irmsd --inversion on
        chemsmart run grouper -d . -p gaussian -T 0.5 irmsd
        chemsmart run grouper -d . -p orca -T 0.5 rmsd
        chemsmart run grouper -d . -t log -T 0.5 irmsd
        chemsmart run grouper -f conformers.xyz -T 0.1 tfd --use-weights
        chemsmart run grouper -f conformers.xyz -N 10 tanimoto --fingerprint-type morgan
        chemsmart run grouper -f conformers.xyz -m csv rmsd
    """
    ctx.ensure_object(dict)

    # Store common options in context for subcommands
    ctx.obj["ignore_hydrogens"] = ignore_hydrogens
    ctx.obj["num_procs"] = num_procs
    ctx.obj["threshold"] = threshold
    ctx.obj["num_groups"] = num_groups
    ctx.obj["matrix_format"] = matrix_format
    ctx.obj["conformer_ids"] = None  # Will be set only in directory mode
    ctx.obj["energy_type"] = energy_type

    try:
        thermo_kwargs = build_thermochemistry_kwargs(
            cutoff_entropy_grimme=cutoff_entropy_grimme,
            cutoff_entropy_truhlar=cutoff_entropy_truhlar,
            cutoff_enthalpy=cutoff_enthalpy,
            concentration=concentration,
            pressure=pressure,
            temperature=temperature,
            alpha=alpha,
            weighted=weighted,
            energy_units=energy_units,
            check_imaginary_frequencies=check_imaginary_frequencies,
        )
    except ValueError as exc:
        raise click.BadParameter(str(exc)) from exc

    ctx.obj["thermo_kwargs"] = thermo_kwargs

    # Validate input
    if filenames and directory:
        raise click.BadParameter(
            "Cannot specify both -f/--filenames and -d/--directory. Choose one."
        )

    if directory and not program and not filetype:
        raise click.BadParameter(
            "Must specify -p/--program or -t/--filetype when using -d/--directory."
        )

    # Mode 1: Directory of output files (-d . -p gaussian)
    if directory is not None:
        molecules, conformer_ids, auto_label = _load_molecules_from_directory(
            directory=directory,
            program=program,
            filetype=filetype,
            energy_type=energy_type,
            thermo_kwargs=thermo_kwargs,
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
        if energy_type != "E":
            raise click.BadParameter(
                "Only energy(default) is supported for xyz files."
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


def _extract_energy_based_on_energy_type(
    thermo, energy_type: str
) -> float | None:
    """Extract the requested energy value."""
    output = thermo.file_object
    key = energy_type.upper()

    if key == "E":
        if output.energies:
            return output.energies[-1]
        return None

    if key == "H":
        try:
            enthalpy = output.enthalpy
        except AttributeError:
            return None
        return enthalpy

    if key == "G":
        try:
            gibbs = output.gibbs_free_energy
        except AttributeError:
            return None
        return gibbs

    if key == "QHH":
        try:
            qh_enthalpy = thermo.qrrho_enthalpy
        except AttributeError:
            return _extract_energy_based_on_energy_type(thermo, "H")
        if qh_enthalpy is None:
            return _extract_energy_based_on_energy_type(thermo, "H")
        return energy_conversion("j/mol", "hartree", qh_enthalpy)

    if key == "QHG":
        try:
            qh_gibbs = thermo.qrrho_gibbs_free_energy
        except AttributeError:
            return _extract_energy_based_on_energy_type(thermo, "G")
        if qh_gibbs is None:
            return _extract_energy_based_on_energy_type(thermo, "G")
        return energy_conversion("j/mol", "hartree", qh_gibbs)

    if key == "SP_QHG":
        # SP-corrected quasi-harmonic Gibbs free energy:
        # qhG - Egas + Esolv
        qh_gibbs = _extract_energy_based_on_energy_type(thermo, "QHG")
        egas = _extract_energy_based_on_energy_type(thermo, "E")

        try:
            source_filepath = thermo.filename
        except AttributeError:
            return qh_gibbs

        sp_filepath = _find_matching_sp_file(source_filepath)
        if sp_filepath is None:
            raise click.ClickException(
                f"Could not find matching SP file for '{source_filepath}'"
            )

        esolv = _extract_last_energy_from_output_file(sp_filepath)

        sp_qh_gibbs = None
        if qh_gibbs is not None and egas is not None and esolv is not None:
            sp_qh_gibbs = qh_gibbs - egas + esolv

        if sp_qh_gibbs is None:
            raise click.ClickException(
                f"Could not obtain final sp_qh_gibbs for '{source_filepath}'"
            )

        return sp_qh_gibbs

    raise click.ClickException(f"Unsupported energy type: {energy_type}")


def _load_molecules_from_directory(
    directory: str,
    program: str,
    filetype: str,
    energy_type: str = "E",
    thermo_kwargs: dict | None = None,
) -> tuple:
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
        program: Program used for calculations ('gaussian' or 'orca')
        filetype: Type of output files ('gaussian' or 'orca')

    Returns:
        tuple: (list of Molecule, list of conformer_ids, common_label)

    Raises:
        click.BadParameter: If no valid output files found or no molecules loaded
    """
    from chemsmart.analysis.thermochemistry import Thermochemistry
    from chemsmart.io.folder import BaseFolder

    directory = os.path.expanduser(directory)
    directory = os.path.abspath(directory)
    folder = BaseFolder(directory)

    if program and not filetype:
        output_files = (
            folder.get_all_output_files_in_current_folder_by_program(
                program=program.lower()
            )
        )
    elif filetype and not program:
        output_files = folder.get_all_files_in_current_folder_by_suffix(
            filetype=filetype
        )
    elif program and filetype:
        output_files = (
            folder.get_all_files_in_current_folder_by_program_and_suffix(
                program=program.lower(), filetype=filetype
            )
        )
    else:
        raise click.BadParameter(
            "Must specify either -p/--program or -t/--filetype when using -d/--directory."
        )

    if not output_files:
        raise click.BadParameter(
            f"No matching output files found in directory: {directory}"
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
    thermo_init_kwargs = dict(THERMOCHEMISTRY_DEFAULT_KWARGS)
    if thermo_kwargs:
        thermo_init_kwargs.update(
            {k: v for k, v in thermo_kwargs.items() if v is not None}
        )

    for filepath, conf_id, _, _ in file_info:
        try:
            # Thermochemistry validates:
            # - normal_termination (via file_object)
            # - imaginary frequencies (via cleaned_frequencies in __init__)
            thermo = Thermochemistry(
                filename=filepath,
                **thermo_init_kwargs,
            )

            mol = thermo.molecule
            if mol is not None:
                mol.name = conf_id

                energy_value = _extract_energy_based_on_energy_type(
                    thermo, energy_type
                )
                if energy_value is None:
                    raise click.ClickException(
                        f"Failed to extract requested energy_type '{energy_type}' for {conf_id} in {filepath}. Stopping."
                    )

                mol._energy = energy_value
                logger.debug(
                    "Loaded %s with %s energy: %.8f Hartree",
                    conf_id,
                    energy_type,
                    energy_value,
                )

                molecules.append(mol)
                conformer_ids.append(conf_id)
            else:
                logger.warning(f"Could not load molecule from {filepath}")

        except click.ClickException:
            raise
        except ValueError as e:
            # Thermochemistry raises ValueError for validation failures
            logger.warning(f"Skipping {conf_id}: {e}")
        except Exception as e:
            logger.warning(f"Error loading {filepath}: {e}")

    if not molecules:
        raise click.BadParameter(
            "No valid molecules could be loaded from output files with the requested energy_type"
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
        "Loaded %d valid molecules from %d output files using energy_type=%s",
        len(molecules),
        len(file_info),
        energy_type,
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
        energy_type: Optional energy type
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
    energy_type = ctx.obj["energy_type"]
    thermo_kwargs = ctx.obj.get("thermo_kwargs")

    thermo_parameters = None
    if (
        energy_type
        and energy_type.upper() in {"QHH", "QHG", "SP_QHG"}
        and thermo_kwargs
    ):
        thermo_items = [
            f"{key}={value}"
            for key, value in thermo_kwargs.items()
            if value is not None
        ]
        if thermo_items:
            thermo_parameters = ", ".join(thermo_items)

    # Use threshold from parent command (None if not specified)
    threshold = ctx.obj["threshold"]
    matrix_format = ctx.obj.get("matrix_format", "xlsx")

    return GrouperJob(
        molecules=molecules,
        grouping_strategy=strategy,
        threshold=threshold,
        num_groups=num_groups,
        ignore_hydrogens=ignore_hydrogens,
        num_procs=num_procs,
        label=f"{label}_{strategy}",
        conformer_ids=conformer_ids,
        matrix_format=matrix_format,
        energy_type=energy_type,
        thermo_parameters=thermo_parameters,
        **extra_kwargs,
    )


def _find_matching_sp_file(source_filepath: str) -> str | None:
    """Find a same-name SP output file under a dedicated sp subfolder."""
    source = Path(source_filepath)
    sp_folder = source.parent / "sp"
    if not sp_folder.is_dir():
        raise ValueError("Can't find SP output folder")

    stem = source.stem
    suffix = source.suffix

    candidates: list[Path] = []
    if suffix:
        candidates.extend(sorted(sp_folder.glob(f"{stem}*{suffix}")))
    candidates.extend(sorted(sp_folder.glob(f"{stem}*.log")))
    candidates.extend(sorted(sp_folder.glob(f"{stem}*.out")))

    seen: set[str] = set()
    for candidate in candidates:
        resolved = str(candidate.resolve())
        if resolved in seen:
            continue
        seen.add(resolved)
        return str(candidate)

    return None


def _extract_last_energy_from_output_file(filepath: str) -> float | None:
    """Extract the last electronic energy from a Gaussian/ORCA output file."""
    program = get_program_type_from_file(filepath)
    if program == "gaussian":
        from chemsmart.io.gaussian.output import Gaussian16Output

        output = Gaussian16Output(filepath)
    elif program == "orca":
        from chemsmart.io.orca.output import ORCAOutput

        output = ORCAOutput(filepath)
    else:
        return None

    if not output.normal_termination:
        return None
    if not output.energies:
        return None
    return output.energies[-1]
