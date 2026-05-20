import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

# Width of the separator lines
LINE_WIDTH = 100

# Keywords that indicate a custom (inline) basis set
CUSTOM_BASIS_KEYWORDS = {"gen", "genecp"}

# Keywords that indicate a custom (inline) solvent
CUSTOM_SOLVENT_KEYWORDS = {"generic,read", "generic"}


def is_chemsmart_database(filepath):
    """Check if a .db file is a chemsmart database.

    A valid chemsmart database contains the four tables:
    'records', 'molecules', 'structures', and 'record_structures'.
    """
    import sqlite3

    required_tables = {
        "records",
        "molecules",
        "structures",
        "record_structures",
    }
    try:
        conn = sqlite3.connect(filepath)
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = {row[0] for row in cursor.fetchall()}
        conn.close()
        return required_tables.issubset(tables)
    except Exception:
        return False


def is_custom_basis(basis):
    """Return True if the basis set keyword indicates a custom inline basis."""
    if basis is None:
        return False
    return basis.strip().lower() in CUSTOM_BASIS_KEYWORDS


def is_custom_solvent(solvent_id):
    """Return True if the solvent ID indicates a custom inline solvent."""
    if solvent_id is None:
        return False
    return solvent_id.strip().lower() in CUSTOM_SOLVENT_KEYWORDS


def get_record_id(structure_id, program, method, basis, jobtype):
    """Generate a stable record ID from molecular identity and calculation fields.

    The hash is computed from the molecule's ``structure_id`` (which already
    encodes canonical geometry, charge and multiplicity) together with the
    computational method descriptors.  This ensures that:

    * The same structure computed with the same method and job type always
      produces the same ``record_id``.
    * Different structures, methods or job types always produce different IDs.

    Args:
        structure_id: SHA-256 hex digest of the canonical geometry, charge
            and multiplicity (from ``Molecule.structure_id``).
        program: Computational chemistry program name (e.g. "gaussian").
        method: DFT functional or method.
        basis: Basis set.
        jobtype: Job type specification.

    Returns:
        SHA-256 hex digest string.
    """
    components = [
        str(structure_id),
        str(program),
        str(method),
        str(basis),
        str(jobtype),
    ]
    payload = "|".join(components)
    return hashlib.sha256(payload.encode()).hexdigest()


def utcnow_iso():
    """Return current UTC time in ISO 8601 format."""
    return datetime.now(timezone.utc).isoformat()


def truncate_iso(iso_str):
    """Shorten an ISO datetime to 'YYYY-MM-DD HH:MM'."""
    if iso_str is None:
        return None
    return iso_str[:16].replace("T", " ")


def file_size(filename):
    """Return file size in bytes."""
    return Path(filename).stat().st_size


def human_size(size_bytes):
    """Convert byte count to human-readable string."""
    if size_bytes is None:
        return "-"
    for unit in ("B", "KB", "MB", "GB"):
        if abs(size_bytes) < 1024:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024
    return f"{size_bytes:.1f} TB"


def sha256_content(output):
    """Compute SHA-256 hash of a file."""
    sha256 = hashlib.sha256(output.content_lines_string.encode()).hexdigest()
    return sha256


def convert_numpy(obj):
    """Recursively convert NumPy types to JSON-serializable types."""
    # NumPy arrays -> lists
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    # NumPy scalars -> Python scalars
    elif isinstance(obj, np.generic):
        return obj.item()
    # Mapping types
    elif isinstance(obj, dict):
        return {k: convert_numpy(v) for k, v in obj.items()}
    # Sequence types
    elif isinstance(obj, (list, tuple)):
        return [convert_numpy(x) for x in obj]
    return obj


def to_json(obj):
    """Convert object to JSON string if not None."""
    if obj is None:
        return None
    return json.dumps(obj)


def from_json(json_str):
    """Parse JSON string, return None if input is None."""
    if json_str is None:
        return None
    return json.loads(json_str)


def separator(title=""):
    """Build a separator line with an optional centred title."""
    if title:
        return f"=== {title} " + "=" * max(1, LINE_WIDTH - len(title) - 5)
    return "=" * LINE_WIDTH


def format_kv(key, value, key_width=30):
    """Format a key-value pair for terminal display."""
    if value is None:
        return f"  {key:<{key_width}}: NULL"
    return f"  {key:<{key_width}}: {value}"


def format_energy(val):
    """Format an energy value (Hartree) for display."""
    if val is None:
        return "NULL"
    return f"{val:.10f}"


def format_float(val, decimals=6):
    """Format a generic float value."""
    if val is None:
        return "NULL"
    return f"{val:.{decimals}f}"


def bool_to_str(val):
    """Convert boolean-ish value to Yes/No."""
    if val is None:
        return "NULL"
    return "Yes" if bool(val) else "No"


def standardize_basis_set(basis):
    """
    Standardize basis set names across different quantum chemistry programs.

    Converts ORCA-style basis names (e.g., 'def2-svp', 'def2-tzvp') to a
    standardized Gaussian-like format (e.g., 'def2svp', 'def2tzvp') by
    removing hyphens between 'def2' and the polarization level.
    """
    if basis.startswith("def2-"):
        return basis.replace("def2-", "def2", 1)
    return basis


def resolve_record(db, record_index=None, record_id=None, return_list=True):
    """Resolve a single record from the database by index or ID.

    This utility function consolidates the common pattern of resolving
    a record by either index or ID, with proper validation and error handling.

    Args:
        db: Database instance (from chemsmart.database.database.Database).
        record_index (int, optional): 1-based record index.
        record_id (str, optional): Record ID or unique prefix.
        return_list (bool): If True, return [record] (for consistency with
            get_all_records). If False, return record directly.
    """
    if record_index is not None:
        record = db.get_record(record_index=record_index)
        if record is None:
            raise ValueError(f"No record found at index {record_index}.")
    elif record_id is not None:
        full_id = db.get_record_by_partial_id(record_id)
        record = db.get_record(record_id=full_id)
        if record is None:
            raise ValueError(f"No record found with ID '{record_id}'.")
    else:
        raise ValueError(
            "Either record_index or record_id must be provided to select "
            "a record from the database."
        )

    return [record] if return_list else record
