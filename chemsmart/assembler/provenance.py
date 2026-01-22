from datetime import datetime, timezone

from chemsmart import __version__ as chemsmart_version

# from chemsmart.utils.hashing import sha256_file


def utcnow_iso() -> str:
    """Return current UTC time in ISO 8601 format."""
    return datetime.now(timezone.utc).isoformat()


def build_provenance(filename, output) -> dict:
    return {
        "source_file": filename,
        # "source_file_hash": sha256_file(filename),
        # "source_file_size": filename.stat().st_size,
        "program": output.__class__.__module__.split(".")[2],
        "program_version": output.version,
        "parser": output.__class__.__name__,
        "chemsmart_version": chemsmart_version,
        "assembled_at": utcnow_iso(),
    }
