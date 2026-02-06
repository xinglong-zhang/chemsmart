import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


def get_record_id(filename):
    """Generate a stable record ID from filename."""
    path = Path(filename).resolve()
    return hashlib.sha256(str(path).encode()).hexdigest()


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
