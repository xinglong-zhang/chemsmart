import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


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


def file_size(filename):
    return Path(filename).stat().st_size


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
