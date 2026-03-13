import functools
import logging
import os

import click

from chemsmart.assembler.query import DatabaseQuery
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.io import resolve_output_path

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
        "-q",
        "--query",
        type=str,
        default=None,
        help=(
            "Query expression. "
            "If omitted, all records are shown. "
            "Supported operators: <, <=, >, >=, =, !=, ~. "
            "Logical: AND, OR. "
            "Example: \"fmo_gap < 7 AND program = 'gaussian'\""
        ),
    )
    @click.option(
        "-o",
        "--output",
        type=str,
        default=None,
        help="Output database file (.db) for saving matching records. ",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@database.command(cls=MyCommand)
@click_query_options
@click.pass_context
def query(ctx, file, query, output):
    """Query records from a chemsmart database.

    Filter records using a query expression and optionally export the
    matching subset to a new database file. If no query is provided,
    all records are displayed.

    Supported logical operators: AND, OR

    Supported comparison operators: <, <=, >, >=, =, !=, ~ (substring match)

    Supported fields: record_id, program, functional, basis, jobtype,
    solvent_on, charge, multiplicity, chemical_formula, smiles,
    number_of_atoms, total_energy, homo_energy, lumo_energy, fmo_gap,
    zero_point_energy, enthalpy, entropy, gibbs_free_energy, source_file

    Note: For molecule-level fields (charge, multiplicity, chemical_formula,
    smiles, number_of_atoms), a record matches if ANY molecule in the record
    satisfies the condition.

    \b
    Examples:
        chemsmart run database query -f database.db
        chemsmart run database query -f my.db -q "chemical_formula = 'CO2'" -o co2.db
        chemsmart run database query -f my.db -q "fmo_gap < 7 AND program = 'gaussian'"
        chemsmart run database query -f my.db -q "source_file ~ 'benzene'"
    """
    # Validate input database
    file = os.path.abspath(file)
    if not os.path.isfile(file):
        raise FileNotFoundError(f"Database file not found: {file}")

    # Prepare output file path if exporting
    if output is not None:
        if not output.lower().endswith(".db"):
            output = output + ".db"

        output, renamed = resolve_output_path(file, output)
        if renamed:
            logger.warning(
                f"Input and output files are the same ({os.path.basename(file)}). "
                f"Writing to {os.path.basename(output)} to avoid overwrite."
            )

    dq = DatabaseQuery(file, query, output)

    # Run query for summaries (terminal display)
    try:
        summaries = dq.query_summaries()
    except ValueError as e:
        logger.error(f"Invalid query: {e}")
        return None

    # Log query results
    total_count = dq.count_records()
    logger.info(
        f"Query matched {len(summaries)} of {total_count} record(s) in {os.path.basename(file)}."
    )

    # Print formatted summary
    print("\n" + dq.format_summary(summaries) + "\n")

    # Export to new database if requested
    if output is not None and summaries:
        records = dq.query()
        dq.export_to_db(records, output)

    return None
