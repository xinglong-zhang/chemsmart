import functools
import logging
import os

import click

from chemsmart.database.export import CSV_OPTIONAL_COLUMNS, DatabaseExporter
from chemsmart.utils.cli import MyCommand

from .database import database

logger = logging.getLogger(__name__)


def click_export_options(f):
    """Common click options for database export."""

    @click.option(
        "-f",
        "--file",
        type=str,
        required=True,
        help="Path to the input database file (.db).",
    )
    @click.option(
        "--ri",
        "--record-index",
        "record_index",
        type=int,
        default=None,
        help="Record index (1-based) to export.",
    )
    @click.option(
        "--rid",
        "--record-id",
        "record_id",
        type=str,
        default=None,
        help="Record ID (or prefix) to export.",
    )
    @click.option(
        "--si",
        "--structure-index",
        "structure_index",
        type=str,
        default=None,
        help="Structure index (1-based) within the selected record. "
        "Defaults to the last structure when used with --ri/--rid.",
    )
    @click.option(
        "--sid",
        "--structure-id",
        "structure_id",
        type=str,
        default=None,
        help="Structure ID (or prefix) to export as a single structure.",
    )
    @click.option(
        "--mid",
        "--molecule-id",
        "molecule_id",
        type=str,
        default=None,
        help="Molecule ID (or prefix); exports every conformer of that molecule.",
    )
    @click.option(
        "-k",
        "--keys",
        type=str,
        default=None,
        help=(
            "Comma-separated extra scalar keys for CSV export. "
            f"Supported: {', '.join(sorted(CSV_OPTIONAL_COLUMNS))}"
        ),
    )
    @click.option(
        "-x",
        "--method-basis",
        "method_basis",
        type=str,
        default=None,
        help=(
            "Filter --sid/--mid XYZ/extXYZ export by 'method/basis' "
            "(e.g. 'MN15/def2tzvp'). For XYZ, structures must have energy "
            "at this level; for extXYZ, structures must have both energy and forces."
        ),
    )
    @click.option(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Output file path. Format inferred from extension (.json, .csv, .xyz, .extxyz).",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@database.command(cls=MyCommand)
@click_export_options
@click.pass_context
def export(
    ctx,
    file,
    record_index,
    record_id,
    structure_index,
    structure_id,
    molecule_id,
    keys,
    method_basis,
    output,
):
    """Export records from a chemsmart database.

    The output format is inferred from the file extension of -o/--output:

    \b
      .json  - Full structured database content
      .csv   - Scalar properties table
      .xyz   - Cartesian coordinates of selected structure(s)
      .extxyz - Extended XYZ with per-frame energy and per-atom forces

    \b
    JSON and CSV always export the entire database; selection options
    (--ri/--rid/--si/--sid/--mid) are accepted only for XYZ/extXYZ.

    \b
    Default CSV columns: record_index, record_id, chemical_formula.
    Use -k to add extra scalar columns.

    \b
    Supported CSV keys:
      program, method, basis, charge, multiplicity, smiles,
      total_energy, homo_energy, lumo_energy, fmo_gap,
      zero_point_energy, enthalpy, entropy, gibbs_free_energy

    \b
    Examples:
        chemsmart run database export -f my.db -o data.json
        chemsmart run database export -f my.db -k total_energy,homo_energy -o training.csv
        chemsmart run database export -f my.db --rid a1b2c3d45e6f -o final.xyz
        chemsmart run database export -f my.db --ri 2 --si 3 -o step3.xyz
        chemsmart run database export -f my.db --ri 2 --si ':' -o traj.xyz
        chemsmart run database export -f my.db --sid 0df6b2ea4bdc -o struct.extxyz
        chemsmart run database export -f my.db --mid BLQJIBCZHWBKSL-U -x 'MN15/def2tzvp' -o conformers.extxyz
    """
    # Validate input database
    file = os.path.abspath(file)
    if not os.path.isfile(file):
        raise click.UsageError(f"Database file not found: {file}")

    from chemsmart.database.utils import (
        check_schema_version,
        is_chemsmart_database,
    )

    if not is_chemsmart_database(file):
        raise click.UsageError(
            f"File {file} is not a valid chemsmart database file."
        )
    try:
        check_schema_version(file)
    except RuntimeError as e:
        raise click.UsageError(str(e))

    output = os.path.abspath(output)
    ext = os.path.splitext(output)[1].lower()

    selectors = {
        "--ri/--record-index": record_index,
        "--rid/--record-id": record_id,
        "--si/--structure-index": structure_index,
        "--sid/--structure-id": structure_id,
        "--mid/--molecule-id": molecule_id,
    }
    used_selectors = [k for k, v in selectors.items() if v is not None]

    if ext in (".json", ".csv"):
        if used_selectors:
            raise click.UsageError(
                f"{ext} export always covers the entire database; "
                f"selection options {', '.join(used_selectors)} are not "
                "allowed. They are only valid for .xyz/.extxyz exports."
            )
        if ext == ".json" and keys is not None:
            raise click.UsageError("-k/--keys is only valid for .csv exports.")
        if method_basis is not None:
            raise click.UsageError(
                "-x/--method-basis is only valid for .xyz/.extxyz exports."
            )
    elif ext in (".xyz", ".extxyz"):
        if keys is not None:
            raise click.UsageError("-k/--keys is only valid for .csv exports.")
        # Mutual exclusivity among record/structure/molecule selectors
        primary = [
            ("--ri/--record-index", record_index),
            ("--rid/--record-id", record_id),
            ("--sid/--structure-id", structure_id),
            ("--mid/--molecule-id", molecule_id),
        ]
        primary_used = [name for name, val in primary if val is not None]
        if len(primary_used) == 0:
            raise click.UsageError(
                f"{ext} export requires exactly one of "
                "--ri/--record-index, --rid/--record-id, "
                "--sid/--structure-id, or --mid/--molecule-id."
            )
        if len(primary_used) > 1:
            raise click.UsageError(
                f"Selectors {', '.join(primary_used)} are mutually "
                "exclusive; please specify exactly one."
            )
        # --si only valid with --ri/--rid
        if structure_index is not None and not (
            record_index is not None or record_id is not None
        ):
            raise click.UsageError(
                "--si/--structure-index can only be used together with "
                "--ri/--record-index or --rid/--record-id."
            )
        if method_basis is not None and not (
            structure_id is not None or molecule_id is not None
        ):
            raise click.UsageError(
                "-x/--method-basis can only be used together with "
                "--sid/--structure-id or --mid/--molecule-id."
            )

    # Parse and validate -x/--method-basis against the database.
    method = basis = None
    if method_basis is not None:
        method, basis = _parse_method_basis(file, method_basis)

    exporter = DatabaseExporter(
        db_file=file,
        output=output,
        record_index=record_index,
        record_id=record_id,
        structure_index=structure_index,
        structure_id=structure_id,
        molecule_id=molecule_id,
        keys=keys,
        method=method,
        basis=basis,
    )

    try:
        exporter.export()
    except ValueError as e:
        raise click.ClickException(str(e))
    logger.info(f"Exported to {os.path.basename(output)}.")

    return None


def _parse_method_basis(db_file, raw):
    """Translate a user-supplied 'method/basis' string into a canonical
    (method, basis) tuple resolved against the database."""
    from chemsmart.database.database import Database

    if "/" not in raw:
        raise click.UsageError(
            f"-x/--method-basis must be given as 'method/basis' "
            f"(got '{raw}'). Example: 'MN15/def2tzvp'."
        )
    method_raw, basis_raw = (s.strip() for s in raw.split("/", 1))
    if not method_raw or not basis_raw:
        raise click.UsageError(
            f"-x/--method-basis must be given as 'method/basis' with "
            f"both parts non-empty (got '{raw}')."
        )

    resolved = Database(db_file).resolve_method_basis(method_raw, basis_raw)
    if resolved is None:
        raise click.UsageError(
            f"-x/--method-basis '{raw}' is not present in the database."
        )
    return resolved
