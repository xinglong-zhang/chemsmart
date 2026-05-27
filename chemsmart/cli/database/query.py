import functools
import logging
import os

import click

from chemsmart.database.query import DatabaseQuery
from chemsmart.utils.cli import MyCommand

from .database import database

logger = logging.getLogger(__name__)


def click_query_options(f):
    """Common click options for database query."""

    @click.option(
        "-f",
        "--file",
        type=str,
        required=True,
        help="Path to the input database file (.db).",
    )
    @click.option(
        "-t",
        "--target",
        type=click.Choice(
            ["records", "molecules", "structures"], case_sensitive=False
        ),
        default="records",
        help="Query target: records (default), molecules, or structures.",
    )
    @click.option(
        "-q",
        "--query",
        type=str,
        default=None,
        help=(
            "Query expression. "
            "If omitted, all entities are shown. "
            "Supported operators: <, <=, >, >=, =, !=, ~. "
            "Logical: AND, OR. "
            "Example: \"fmo_gap < 7 AND program = 'gaussian'\""
        ),
    )
    @click.option(
        "-l",
        "--limit",
        type=int,
        default=None,
        help="Maximum number of results to display.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@database.command(cls=MyCommand)
@click_query_options
@click.pass_context
def query(ctx, file, target, query, limit):
    """Query records, molecules, or structures from a chemsmart database.

    Use -t, --target to switch the query perspective:

    \b
      -t records (default):
        Query calculation records.
        Fields: source_file, jobtype, program, method, basis, total_energy,
        solvent_on, solvent_model, normal_termination, homo_energy, lumo_energy,
        fmo_gap, zero_point_energy, enthalpy, entropy, gibbs_free_energy,

    \b
      -t molecules:
        Query unique chemical species.
        Fields: chemical_formula, smiles, number_of_atoms, mass

    \b
      -t structures:
        Query 3D conformers.
        Fields: chemical_formula, number_of_atoms, charge,
        multiplicity

    Supported operators: <, <=, >, >=, =, !=, ~ (substring match)

    Logical: AND, OR

    \b
    Examples:
        chemsmart run database query -f my.db
        chemsmart run database query -f my.db -q "fmo_gap < 7 AND program = 'gaussian'"
        chemsmart run database query -f my.db -q "source_file ~ 'benzene'"
        chemsmart run database query -f my.db -t molecules
        chemsmart run database query -f my.db -t molecules -q "mass > 100"
        chemsmart run database query -f my.db -t structures -l 10
    """
    # Validate input database
    file = os.path.abspath(file)
    if not os.path.isfile(file):
        raise FileNotFoundError(f"Database file not found: {file}")

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

    if limit is not None:
        if limit <= 0:
            raise click.BadParameter("Limit must be a positive integer.")

    dq = DatabaseQuery(file, query, target=target, limit=limit)

    # Run query for summaries (terminal display)
    try:
        summaries = dq.query_summaries()
    except ValueError as e:
        logger.error(f"Invalid query: {e}")
        return None

    # Log query results
    entity = dq._config["entity_name"]
    total_count = dq.count_total()
    logger.info(
        f"Query returned {len(summaries)} of {total_count} {entity}(s) "
        f"in {os.path.basename(file)}."
    )

    # Print formatted summary
    print("\n" + dq.format_summary(summaries) + "\n")

    return None
