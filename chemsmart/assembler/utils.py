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
CUSTOM_SOLVENT_KEYWORDS = {"generic,read"}


def is_chemsmart_database(filepath):
    """Check if a .db file is a chemsmart database (has both 'records' and 'molecules' tables)."""
    import sqlite3

    try:
        conn = sqlite3.connect(filepath)
        cursor = conn.cursor()
        cursor.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='records'"
        )
        has_records = cursor.fetchone() is not None
        cursor.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='molecules'"
        )
        has_molecules = cursor.fetchone() is not None
        conn.close()
        return has_records and has_molecules
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


def get_record_id(
    canonical_geometry, charge, multiplicity, program, functional, basis
):
    """Generate a stable record ID from molecular identity fields.

    The hash is computed from canonical geometry (sorted symbols + rounded
    positions), charge, multiplicity, program, functional, and basis set.
    This ensures that the same calculation always produces the same ID
    regardless of file path.

    Args:
        canonical_geometry: String representation of the molecular geometry
            (e.g. sorted chemical symbols + rounded positions).
        charge: Molecular charge.
        multiplicity: Spin multiplicity.
        program: Computational chemistry program name.
        functional: DFT functional or method.
        basis: Basis set.

    Returns:
        SHA-256 hex digest string.
    """
    components = [
        str(canonical_geometry),
        str(charge),
        str(multiplicity),
        str(program),
        str(functional),
        str(basis),
    ]
    payload = "|".join(components)
    return hashlib.sha256(payload.encode()).hexdigest()


def canonical_geometry_string(chemical_symbols, positions, decimals=6):
    """Build a canonical string representation of molecular geometry.

    Atoms are sorted by (symbol, x, y, z) so that the representation is
    invariant to input ordering.

    Args:
        chemical_symbols: List of element symbols.
        positions: Nx3 array-like of Cartesian coordinates.
        decimals: Number of decimal places for coordinate rounding.

    Returns:
        Canonical geometry string.
    """
    rounded = np.round(np.asarray(positions, dtype=float), decimals=decimals)
    atoms = sorted(zip(chemical_symbols, rounded.tolist()))
    parts = [
        f"{sym}:{x:.{decimals}f},{y:.{decimals}f},{z:.{decimals}f}"
        for sym, (x, y, z) in atoms
    ]
    return ";".join(parts)


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


def format_kv(key, value, key_width=28):
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
