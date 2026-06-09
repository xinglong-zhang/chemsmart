import hashlib
import json
import logging
import re
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from chemsmart.utils.repattern import route_split_pattern

logger = logging.getLogger(__name__)

SCHEMA_VERSION = 1

# Width of the separator lines
LINE_WIDTH = 100

# Keywords that indicate a custom (inline) basis set
CUSTOM_BASIS_KEYWORDS = {"gen", "genecp"}

# Keywords that indicate a custom (inline) solvent
CUSTOM_SOLVENT_KEYWORDS = {"generic,read", "generic"}

# Tokens that mark the start of a route line
ROUTE_PREFIX_TOKENS = {"#", "#p", "#n", "#t", "!"}


def open_connection(db_file):
    """Open a SQLite connection with chemsmart's standard pragmas."""
    conn = sqlite3.connect(db_file, timeout=30.0)
    conn.execute("PRAGMA foreign_keys = ON")
    conn.execute("PRAGMA synchronous = NORMAL")
    return conn


def check_schema_version(db_file):
    """Raise RuntimeError if the database schema version does not match SCHEMA_VERSION."""
    conn = open_connection(db_file)
    try:
        version = conn.execute("PRAGMA user_version").fetchone()[0]
    finally:
        conn.close()
    if version != SCHEMA_VERSION:
        raise RuntimeError(
            f"Database schema version mismatch: expected {SCHEMA_VERSION}, got {version}. "
            f"Re-assemble the database with the current version of chemsmart."
        )


def is_chemsmart_database(filepath):
    """Check if a .db file is a chemsmart database.

    A valid chemsmart database contains the four tables:
    'records', 'molecules', 'structures', and 'record_structures'.
    """
    if not is_database_file(filepath):
        return False

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


def is_database_file(filename):
    """Check if a file is a non-empty .db path."""
    return filename is not None and filename.endswith(".db")


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


def get_record_id(
    *,
    program,
    method,
    basis,
    jobtype,
    trajectory_id,
    custom_basis_hash=None,
    solvent_model=None,
    solvent_id=None,
    custom_solvent_hash=None,
    route_hash=None,
):
    """Generate a stable record ID for one full QM calculation.
    The record_id uniquely identifies a single source file's worth of work:
    a specific (program, method, basis, jobtype, solvent setup, route)
    applied to a specific trajectory of structures.
    Args:
        program: Computational chemistry program name (e.g. "gaussian").
        method: DFT functional / method (e.g. "b3lyp").
        basis: Basis set name (e.g. "def2svp"); for inline custom basis,
            pass the literal keyword (e.g. "gen") and provide the actual
            basis content hash via custom_basis_hash.
        jobtype: Job type label (e.g. "opt", "ts", "freq", "sp").
        trajectory_id: SHA-256 hash of the ordered list of structure_ids
            visited in this calculation.
        custom_basis_hash: SHA-256 of canonical-JSON of the custom basis
            definition (heavy_elements / ECP / light_elements / ...) or
            None when no custom basis is used.
        solvent_model: Solvent model name (e.g. "pcm", "smd") or None.
        solvent_id: Solvent identifier (e.g. "water", "toluene") or None.
        custom_solvent_hash: SHA-256 of canonical-JSON of the custom-solvent
            parameters block, or None.
        route_hash: SHA-256 hash of the canonicalized route string. Makes
        record identity sensitive to route-level options not otherwise
        captured in structured fields, such as grid, SCF, dispersion, symmetry,
        and opt/freq settings. None if unavailable.
    Returns:
        SHA-256 hex digest string.
    """
    components = [
        program,
        method,
        basis,
        jobtype,
        trajectory_id,
        custom_basis_hash,
        solvent_model,
        solvent_id,
        custom_solvent_hash,
        route_hash,
    ]
    payload = "|".join("" if c is None else str(c) for c in components)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def canonicalize_route_string(route_string, drop_terms=()):
    """Canonicalize a QM-program route string for hashing into record_id."""
    if not route_string:
        return None
    drop_set = {t.lower().replace("-", "") for t in drop_terms if t}
    kept = []
    for token in re.split(route_split_pattern, route_string):
        if not token:
            continue
        norm = token.lower().replace("-", "")
        if norm in ROUTE_PREFIX_TOKENS:
            continue
        if norm in drop_set:
            continue
        kept.append(norm)
    return sorted(kept) or None


def canonical_json_hash(obj):
    """SHA-256 of canonical JSON of a dict/list object."""
    if not obj:
        return None
    payload = json.dumps(
        obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def compute_trajectory_id(structure_ids):
    """Compute a trajectory ID from an ordered list of structure_id."""
    return hashlib.sha256("|".join(structure_ids).encode("utf-8")).hexdigest()


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
    return f"{val:.8f}"


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
    if not basis:
        return basis
    if basis.startswith("def2-"):
        return basis.replace("def2-", "def2", 1)
    return basis


def collect_energies_for_structure(db_file, structure_id):
    """Return list of (method, basis, energy) for every record that references
    the given structure_id, ordered by record_index.
    """
    conn = open_connection(db_file)
    conn.row_factory = sqlite3.Row
    try:
        cursor = conn.execute(
            """
            SELECT r.method AS method, r.basis AS basis,
                   rs.energy AS energy
            FROM record_structures rs
            JOIN records r ON rs.record_id = r.record_id
            WHERE rs.structure_id = ?
            ORDER BY r.record_index
            """,
            (structure_id,),
        )
        rows = cursor.fetchall()
    finally:
        conn.close()
    return [
        (r["method"], r["basis"], r["energy"])
        for r in rows
        if r["energy"] is not None
    ]


def sort_frames_by_energy(frames, primary=None):
    """Sort frames ascending by energy at a chosen (method, basis) key;
    frames missing that key go to the end (sorted by their own lowest
    available energy).
    """
    if primary is None:
        # 1) Count (method, basis) coverage across all frames
        counts = {}
        for frame in frames:
            for method, basis, _ in frame.get("energies", []):
                key = (method, basis)
                counts[key] = counts.get(key, 0) + 1
        if not counts:
            return frames  # nothing to sort by
        # 2) Pick the most-covered (method, basis) as the primary sorting key
        primary = max(counts.items(), key=lambda kv: kv[1])[0]
        logger.info(
            f"Sorting {len(frames)} frame(s) by energy at "
            f"{primary[0]}/{primary[1]}; frames missing it go to the end."
        )
    else:
        logger.info(
            f"Sorting {len(frames)} frame(s) by specified energy at "
            f"{primary[0]}/{primary[1]}; frames missing it go to the end."
        )

    def sort_key(frame):
        energies = frame.get("energies", [])
        primary_e = next(
            (e for m, b, e in energies if (m, b) == primary), None
        )
        sid = str(frame.get("structure_id") or "")
        if primary_e is not None:
            return (0, primary_e, sid)
        fallback = min((e for _, _, e in energies), default=float("inf"))
        return (1, fallback, sid)

    sorted_frames = sorted(frames, key=sort_key)

    # 3) Within each frame, move the primary (method, basis) entry to the
    #    front so comment lines / display shows it first.
    for frame in sorted_frames:
        energies = frame.get("energies", [])
        head = [(m, b, e) for (m, b, e) in energies if (m, b) == primary]
        tail = [(m, b, e) for (m, b, e) in energies if (m, b) != primary]
        frame["energies"] = head + tail

    return sorted_frames


def sort_structure_dicts_by_energy(db_file, struct_dicts):
    """Sort a list of structure dicts ascending by energy at the most-covered
    (method, basis) pair.
    """
    if not struct_dicts:
        return struct_dicts
    frames = [
        {
            "structure_id": s.get("structure_id"),
            "energies": (
                collect_energies_for_structure(db_file, s.get("structure_id"))
                if s.get("structure_id") is not None
                else []
            ),
            "struct_dict": s,
        }
        for s in struct_dicts
    ]
    sorted_frames = sort_frames_by_energy(frames)
    if sorted_frames:
        first_energies = sorted_frames[0].get("energies", [])
        primary = first_energies[0] if first_energies else None
        primary_mb = (primary[0], primary[1]) if primary else None
    else:
        primary_mb = None
    sorted_dicts = []
    for frame in sorted_frames:
        s = dict(frame["struct_dict"])
        energies = frame.get("energies", [])
        if primary_mb is not None:
            primary_e = next(
                (e for m, b, e in energies if (m, b) == primary_mb), None
            )
        else:
            primary_e = None
        s["primary_energy"] = primary_e
        s["primary_method_basis"] = primary_mb
        sorted_dicts.append(s)
    return sorted_dicts


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
