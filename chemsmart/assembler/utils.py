import hashlib
from datetime import datetime, timezone
from pathlib import Path


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
