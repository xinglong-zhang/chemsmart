import logging
import os
import re
import sqlite3

from chemsmart.assembler.database import Database
from chemsmart.utils.repattern import (
    query_condition_pattern,
    query_logic_split_pattern,
)

logger = logging.getLogger(__name__)

_CONDITION_RE = re.compile(query_condition_pattern)
_LOGIC_SPLIT_RE = re.compile(query_logic_split_pattern, re.IGNORECASE)

TABLE_FIELDS = {
    # records table
    "r": {
        "record_id",
        "program",
        "functional",
        "basis",
        "jobtype",
        "solvent_on",
        "total_energy",
        "homo_energy",
        "lumo_energy",
        "fmo_gap",
        "zero_point_energy",
        "enthalpy",
        "entropy",
        "gibbs_free_energy",
        "source_file",
    },
    # molecules table
    "m": {
        "charge",
        "multiplicity",
        "chemical_formula",
        "smiles",
        "number_of_atoms",
    },
}

QUERYABLE_FIELDS = {
    field: (alias, field)
    for alias, fields in TABLE_FIELDS.items()
    for field in fields
}

SUPPORTED_OPERATORS = {"<", "<=", ">", ">=", "=", "!=", "~"}


class DatabaseQuery:
    """Query and filter records from a chemsmart SQLite database."""

    _SUMMARY_SQL = """
        SELECT DISTINCT
            r.record_index,
            r.record_id,
            r.program,
            r.functional,
            r.basis,
            r.jobtype,
            r.total_energy,
            r.source_file,
            m.chemical_formula
        FROM records r
        LEFT JOIN molecules m ON r.record_id = m.record_id
    """

    _TABLE_COLUMNS = [
        ("Idx", "record_index", 4, ">"),
        ("Record ID", "record_id", 14, "<"),
        ("File", "source_file", 22, "<"),
        ("Formula", "chemical_formula", 16, "<"),
        ("Job", "jobtype", 6, "<"),
        ("Program", "program", 8, "<"),
        ("Functional", "functional", 32, "<"),
        ("Basis", "basis", 12, "<"),
        ("Total Energy (Eh)", "total_energy", 16, ">"),
    ]

    def __init__(self, db_file, query_string, output_file=None):
        self.db_file = db_file
        self.query_string = query_string
        self.output_file = output_file

    def parse_query(self):
        """Translate a user query string into a parameterised SQL WHERE clause."""
        if not self.query_string:
            raise ValueError("No query string provided.")
        tokens = _LOGIC_SPLIT_RE.split(self.query_string.strip())

        clause_parts = []
        params = []

        for token in tokens:
            token = token.strip()
            if token.upper() in ("AND", "OR"):
                clause_parts.append(token.upper())
                continue

            match = _CONDITION_RE.fullmatch(token)
            if not match:
                raise ValueError(
                    f"Invalid condition: '{token}'. "
                    f"Expected format: FIELD OPERATOR VALUE"
                )

            field, operator, raw_value = match.groups()

            if field not in QUERYABLE_FIELDS:
                raise ValueError(
                    f"Unknown query field: '{field}'. "
                    f"Supported fields: {', '.join(sorted(QUERYABLE_FIELDS))}"
                )
            if operator not in SUPPORTED_OPERATORS:
                raise ValueError(
                    f"Unsupported operator: '{operator}'. "
                    f"Supported: {', '.join(sorted(SUPPORTED_OPERATORS))}"
                )

            alias, column = QUERYABLE_FIELDS[field]
            qualified = f"{alias}.{column}"

            # Convert query literal to Python type
            if (raw_value.startswith("'") and raw_value.endswith("'")) or (
                raw_value.startswith('"') and raw_value.endswith('"')
            ):
                value = raw_value[1:-1]
            else:
                try:
                    value = int(raw_value)
                except ValueError:
                    try:
                        value = float(raw_value)
                    except ValueError:
                        value = raw_value

            if operator == "~":
                clause_parts.append(f"{qualified} LIKE ?")
                params.append(f"%{value}%")
            else:
                clause_parts.append(f"{qualified} {operator} ?")
                params.append(value)

        if not clause_parts:
            raise ValueError("Empty query string.")

        return " ".join(clause_parts), tuple(params)

    def query(self):
        """Return full records matching *query_string*, or all records if no query."""
        conn = sqlite3.connect(self.db_file)
        conn.row_factory = sqlite3.Row
        try:
            if self.query_string:
                where_clause, params = self.parse_query()
                sql = (
                    "SELECT DISTINCT r.* FROM records r "
                    "LEFT JOIN molecules m ON r.record_id = m.record_id "
                    f"WHERE {where_clause} ORDER BY r.record_index"
                )
            else:
                params = ()
                sql = (
                    "SELECT DISTINCT r.* FROM records r "
                    "LEFT JOIN molecules m ON r.record_id = m.record_id "
                    "ORDER BY r.record_index"
                )
            rows = conn.execute(sql, params).fetchall()
            db = Database(self.db_file)
            return [db._row_to_full_record(conn, dict(row)) for row in rows]
        finally:
            conn.close()

    def query_summaries(self):
        """Return lightweight summaries matching *query_string*, or all records if no query."""
        conn = sqlite3.connect(self.db_file)
        conn.row_factory = sqlite3.Row
        try:
            if self.query_string:
                where_clause, params = self.parse_query()
                sql = f"{self._SUMMARY_SQL} WHERE {where_clause} ORDER BY r.record_index"
            else:
                params = ()
                sql = f"{self._SUMMARY_SQL} ORDER BY r.record_index"
            rows = conn.execute(sql, params).fetchall()
            return [dict(row) for row in rows]
        finally:
            conn.close()

    def count_records(self):
        """Count the total number of records in the database.

        Returns:
            Total number of records.
        """
        conn = sqlite3.connect(self.db_file)
        try:
            cursor = conn.execute("SELECT COUNT(*) FROM records")
            return cursor.fetchone()[0]
        finally:
            conn.close()

    def export_to_db(self, records, output_file):
        """Write matching records into a new database."""
        out_db = Database(db_file=output_file)
        out_db.create()

        conn = sqlite3.connect(output_file)
        count = 0
        try:
            for record in records:
                try:
                    out_db.insert_record(record, conn=conn)
                    count += 1
                except Exception as e:
                    logger.error(f"Failed to export record: {e}")
            conn.commit()
        finally:
            conn.close()
        return count

    def format_summary(self, summaries):
        """Format complete summary output with header, table, and footer."""
        db_name = os.path.basename(self.db_file)
        matched_count = len(summaries)

        # Build header
        separator1 = "=" * 3 + " Query Summary " + "=" * 120
        separator2 = "=" * 138
        header_lines = [
            separator1,
            f"DB      : {db_name}",
            f"Query   : {self.query_string or '<all>'}",
            f"Matched : {matched_count} record(s)",
        ]
        if self.output_file:
            output_name = os.path.basename(self.output_file)
            header_lines.append(f"Output  : {output_name}")
        header_lines.append(separator2)

        # Build body
        if summaries:
            body_lines = self._format_table(summaries)
        else:
            body_lines = ["No records found."]

        return "\n".join(header_lines + [""] + body_lines)

    def _format_table(self, summaries):
        """Format summaries as a table (without header/footer)."""
        if not summaries:
            return []

        # Build header row
        header_parts, sep_parts = [], []
        for hdr, _key, width, align in self._TABLE_COLUMNS:
            fmt = f"{hdr:<{width}}" if align == "<" else f"{hdr:>{width}}"
            header_parts.append(fmt)
            sep_parts.append("-" * width)

        lines = [" ".join(header_parts), " ".join(sep_parts)]

        # Build data rows
        for row in summaries:
            parts = []
            for _hdr, key, width, align in self._TABLE_COLUMNS:
                val = row.get(key)
                if val is None:
                    text = ""
                elif key == "record_id":
                    text = str(val)[:width]
                elif key == "source_file":
                    text = os.path.basename(str(val))[:width]
                elif key == "total_energy" and isinstance(val, (int, float)):
                    text = f"{val:.6f}"
                else:
                    text = str(val)
                parts.append(
                    f"{text:>{width}}" if align == ">" else f"{text:<{width}}"
                )
            lines.append(" ".join(parts))

        return lines
