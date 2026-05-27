"""
Database query module for querying records, molecules, or structures
from a chemsmart database.

Supports three query targets via the --target option:
- records (default): Query calculation records.
- molecules: Query unique chemical species.
- structures: Query 3D conformers.
"""

import logging
import os
import re
import sqlite3

from chemsmart.database.utils import format_float, open_connection, separator
from chemsmart.utils.repattern import (
    query_condition_pattern,
    query_logic_split_pattern,
)

logger = logging.getLogger(__name__)

_CONDITION_RE = re.compile(query_condition_pattern)
_LOGIC_SPLIT_RE = re.compile(query_logic_split_pattern, re.IGNORECASE)

SUPPORTED_OPERATORS = {"<", "<=", ">", ">=", "=", "!=", "~"}
VALID_TARGETS = {"records", "molecules", "structures"}

_RECORDS_JOIN = """
    FROM records r
    LEFT JOIN record_structures rs ON r.record_id = rs.record_id
    LEFT JOIN structures s ON rs.structure_id = s.structure_id
    LEFT JOIN molecules m ON s.molecule_id = m.molecule_id
"""

_STRUCTURES_JOIN = """
    FROM structures s
    JOIN molecules m ON s.molecule_id = m.molecule_id
"""

TARGET_CONFIG = {
    "records": {
        "entity_name": "record",
        "queryable_fields": {
            "program": "r.program",
            "method": "r.method",
            "basis": "r.basis",
            "jobtype": "r.jobtype",
            "solvent_on": "r.solvent_on",
            "solvent_model": "r.solvent_model",
            "normal_termination": "r.normal_termination",
            "total_energy": "r.total_energy",
            "homo_energy": "r.homo_energy",
            "lumo_energy": "r.lumo_energy",
            "fmo_gap": "r.fmo_gap",
            "zero_point_energy": "r.zero_point_energy",
            "enthalpy": "r.enthalpy",
            "entropy": "r.entropy",
            "gibbs_free_energy": "r.gibbs_free_energy",
            "source_file": "r.source_file",
        },
        "summary_select": f"""
            SELECT DISTINCT
                r.record_index, r.record_id, r.program, r.method,
                r.basis, r.jobtype, r.total_energy, r.source_file,
                r.normal_termination,
                m.chemical_formula, m.molecule_id
            {_RECORDS_JOIN}
        """,
        "count_total_sql": "SELECT COUNT(*) FROM records",
        "count_matched_template": (
            f"SELECT COUNT(DISTINCT r.record_id) {_RECORDS_JOIN}"
            " WHERE {where}"
        ),
        "order_by": "r.record_index",
        "table_columns": [
            ("Idx", "record_index", 4, ">"),
            ("Record ID", "record_id", 12, "<"),
            ("File", "source_file", 22, "<"),
            ("Program", "program", 8, "<"),
            ("Job", "jobtype", 6, "<"),
            ("Status", "normal_termination", 6, "<"),
            ("Method", "method", 12, "<"),
            ("Basis", "basis", 16, "<"),
            ("Formula", "chemical_formula", 20, "<"),
            ("Energy (Eh)", "total_energy", 11, ">"),
        ],
        "value_formatters": {
            "record_id": lambda v, w: str(v)[:w],
            "source_file": lambda v, w: os.path.basename(str(v))[:w],
            "total_energy": lambda v, w: (
                format_float(v, 2) if isinstance(v, (int, float)) else str(v)
            ),
            "normal_termination": lambda v, w: "normal" if v else "failed",
        },
    },
    "molecules": {
        "entity_name": "molecule",
        "queryable_fields": {
            "chemical_formula": "m.chemical_formula",
            "smiles": "m.smiles",
            "inchi": "m.inchi",
            "number_of_atoms": "m.number_of_atoms",
            "mass": "m.mass",
        },
        "summary_select": """
            SELECT
                m.molecule_id, m.chemical_formula, m.number_of_atoms,
                m.mass, m.smiles,
                (SELECT COUNT(*) FROM structures s2
                 WHERE s2.molecule_id = m.molecule_id) AS num_structures,
                (SELECT COUNT(DISTINCT rs2.record_id)
                 FROM record_structures rs2
                 JOIN structures s3 ON rs2.structure_id = s3.structure_id
                 WHERE s3.molecule_id = m.molecule_id) AS num_records
            FROM molecules m
        """,
        "count_total_sql": "SELECT COUNT(*) FROM molecules",
        "count_matched_template": (
            "SELECT COUNT(*) FROM molecules m WHERE {where}"
        ),
        "order_by": "m.chemical_formula",
        "table_columns": [
            ("Molecule ID", "molecule_id", 27, "<"),
            ("SMILES", "smiles", 18, "<"),
            ("Formula", "chemical_formula", 18, "<"),
            ("Atoms", "number_of_atoms", 6, ">"),
            ("Mass", "mass", 8, ">"),
            ("Structures", "num_structures", 11, ">"),
            ("Records", "num_records", 8, ">"),
        ],
        "value_formatters": {
            "molecule_id": lambda v, w: str(v)[:w],
            "smiles": lambda v, w: (
                str(v or "")
                if len(str(v or "")) <= w
                else str(v or "")[: max(0, w - 3)] + "..."
            ),
            "mass": lambda v, w: (
                format_float(v, 2) if isinstance(v, (int, float)) else str(v)
            ),
        },
    },
    "structures": {
        "entity_name": "structure",
        "queryable_fields": {
            "chemical_formula": "m.chemical_formula",
            "number_of_atoms": "m.number_of_atoms",
            "charge": "s.charge",
            "multiplicity": "s.multiplicity",
        },
        "summary_select": f"""
            SELECT
                s.structure_id, s.molecule_id, m.chemical_formula,
                s.charge, s.multiplicity, 
                (SELECT COUNT(DISTINCT rs2.record_id) FROM record_structures rs2
                 WHERE rs2.structure_id = s.structure_id) AS num_records
            {_STRUCTURES_JOIN}
        """,
        "count_total_sql": "SELECT COUNT(*) FROM structures",
        "count_matched_template": (
            f"SELECT COUNT(*) {_STRUCTURES_JOIN} WHERE {{where}}"
        ),
        "order_by": "s.molecule_id, s.structure_id",
        "table_columns": [
            ("Structure ID", "structure_id", 13, "<"),
            ("Charge", "charge", 7, ">"),
            ("Mult", "multiplicity", 5, ">"),
            ("Formula", "chemical_formula", 20, "<"),
            ("Molecule ID", "molecule_id", 16, "<"),
            ("Records", "num_records", 8, ">"),
        ],
        "value_formatters": {
            "structure_id": lambda v, w: str(v)[:w],
            "molecule_id": lambda v, w: str(v)[:w],
        },
    },
}


class DatabaseQuery:
    """Query and filter records, molecules, or structures from a chemsmart database."""

    def __init__(self, db_file, query_string, target="records", limit=None):
        """Initialize a database query.

        Args:
            db_file: Path to the chemsmart database file.
            query_string: Optional query expression to filter results.
            target: Query target — "records", "molecules", or "structures".
            limit: Maximum number of results to return.
        """
        if target not in VALID_TARGETS:
            raise ValueError(
                f"Invalid target '{target}'. "
                f"Must be one of: {', '.join(sorted(VALID_TARGETS))}"
            )
        self.db_file = db_file
        self.query_string = query_string
        self.target = target
        self.limit = limit
        self._config = TARGET_CONFIG[target]

    def parse_query(self):
        """Translate a user query string into a parameterised SQL WHERE clause."""
        if not self.query_string:
            raise ValueError("No query string provided.")

        queryable = self._config["queryable_fields"]
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

            if field not in queryable:
                raise ValueError(
                    f"Unknown field '{field}' for target '{self.target}'. "
                    f"Supported fields for '{self.target}': "
                    f"{', '.join(sorted(queryable))}"
                )
            if operator not in SUPPORTED_OPERATORS:
                raise ValueError(
                    f"Unsupported operator: '{operator}'. "
                    f"Supported: {', '.join(sorted(SUPPORTED_OPERATORS))}"
                )

            qualified = queryable[field]

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

    def query_summaries(self):
        """Return lightweight summaries for the current target.

        Returns:
            List of dictionaries, one per matched entity.
        """
        conn = open_connection(self.db_file)
        conn.row_factory = sqlite3.Row
        try:
            sql = self._config["summary_select"]
            if self.query_string:
                where_clause, params = self.parse_query()
                sql += f" WHERE {where_clause}"
            else:
                params = ()
            sql += f" ORDER BY {self._config['order_by']}"
            if self.limit is not None and self.limit > 0:
                sql += f" LIMIT {self.limit}"
            rows = conn.execute(sql, params).fetchall()
            return [dict(row) for row in rows]
        finally:
            conn.close()

    def count_total(self):
        """Count total entities for the current target.

        Returns:
            Total number of entities in the database.
        """
        conn = open_connection(self.db_file)
        try:
            return conn.execute(self._config["count_total_sql"]).fetchone()[0]
        finally:
            conn.close()

    def count_matched(self):
        """Count matched entities for the current target (ignoring limit).

        Returns:
            Number of entities matching the query, or total if no query.
        """
        conn = open_connection(self.db_file)
        try:
            if self.query_string:
                where_clause, params = self.parse_query()
                sql = self._config["count_matched_template"].format(
                    where=where_clause
                )
                return conn.execute(sql, params).fetchone()[0]
            else:
                return conn.execute(
                    self._config["count_total_sql"]
                ).fetchone()[0]
        finally:
            conn.close()

    def format_summary(self, summaries):
        """Format complete summary output with header, table, and footer.

        Args:
            summaries: List of summary dicts from query_summaries.

        Returns:
            Formatted multi-line string ready for terminal output.
        """
        entity = self._config["entity_name"]
        db_name = os.path.basename(self.db_file)

        title = f"Query Summary ({self.target})"
        header_lines = [
            separator(title),
            f"DB      : {db_name}",
            f"Total   : {self.count_total()} {entity}(s)",
        ]
        if self.query_string:
            header_lines.append(f"Query   : {self.query_string}")
            header_lines.append(
                f"Matched : {self.count_matched()} {entity}(s)"
            )
        if self.limit:
            header_lines.append(f"Limit   : {self.limit}")
        header_lines.append(separator())

        # Build body
        if summaries:
            body_lines = self._format_table(summaries)
        else:
            body_lines = [f"No {entity}s found."]

        return "\n".join(header_lines + [""] + body_lines)

    def _format_table(self, summaries):
        """Format summaries as a table, using target-specific columns and formatters."""
        if not summaries:
            return []

        columns = self._config["table_columns"]
        formatters = self._config.get("value_formatters", {})

        # Build header row
        header_parts, sep_parts = [], []
        for hdr, _key, width, align in columns:
            fmt = f"{hdr:<{width}}" if align == "<" else f"{hdr:>{width}}"
            header_parts.append(fmt)
            sep_parts.append("-" * width)

        lines = [" ".join(header_parts), " ".join(sep_parts)]

        # Build data rows
        for row in summaries:
            parts = []
            for _hdr, key, width, align in columns:
                val = row.get(key)
                if val is None:
                    text = ""
                elif key in formatters:
                    text = formatters[key](val, width)
                else:
                    text = str(val)
                parts.append(
                    f"{text:>{width}}" if align == ">" else f"{text:<{width}}"
                )
            lines.append(" ".join(parts))

        return lines
