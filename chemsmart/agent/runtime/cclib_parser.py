"""Safe cclib projection for Gaussian and ORCA calculation outputs."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

try:
    from cclib.io import ccopen
    from cclib.parser.utils import convertor
except ImportError:  # Optional ``result-parsers`` extra.
    ccopen = None
    convertor = None

MAX_CCLIB_FILE_SIZE = 32 * 1024 * 1024


def parse_cclib_output(output_path: str | Path) -> dict[str, Any]:
    """Return a compact cclib result without making cclib a success oracle.

    cclib provides normalized electronic-structure values, but support can lag
    new Gaussian/ORCA output formats. Parse errors are therefore evidence for a
    deterministic fallback, not calculation failures.
    """

    path = Path(output_path)
    if not path.is_file():
        return {
            "status": "not_found",
            "warning": "Output file does not exist.",
        }
    if path.stat().st_size > MAX_CCLIB_FILE_SIZE:
        return {
            "status": "skipped",
            "warning": "Output is larger than the bounded cclib parse limit.",
        }
    if ccopen is None or convertor is None:
        return {
            "status": "unavailable",
            "warning": (
                "cclib is not installed; using the built-in result parser. "
                "Install chemsmart[result-parsers] for normalized advanced "
                "properties."
            ),
        }

    try:
        parser = ccopen(
            str(path),
            quiet=True,
            loglevel=logging.CRITICAL,
        )
    except Exception as exc:  # cclib format detection is plugin-driven.
        return _failure("detection_failed", exc)
    if parser is None:
        return {
            "status": "unsupported",
            "warning": "cclib did not identify a supported output format.",
        }

    try:
        data = parser.parse()
    except Exception as exc:  # cclib may lag newly released output formats.
        return _failure("parse_failed", exc)

    result: dict[str, Any] = {
        "status": "ok",
        "parser": type(parser).__name__.lower(),
    }
    scfenergies = getattr(data, "scfenergies", None)
    if scfenergies is not None and len(scfenergies):
        result["energy"] = float(
            convertor(float(scfenergies[-1]), "eV", "hartree")
        )

    scfvalues = getattr(data, "scfvalues", None)
    if scfvalues is not None and len(scfvalues) and len(scfvalues[-1]):
        result["scf_cycles"] = len(scfvalues[-1])

    vibfreqs = getattr(data, "vibfreqs", None)
    if vibfreqs is not None:
        frequencies = [float(value) for value in vibfreqs]
        result["frequencies"] = frequencies
        result["frequency_count"] = len(frequencies)
        result["imag_freqs"] = [value for value in frequencies if value < 0]
        result["imaginary_frequency_count"] = len(result["imag_freqs"])

    atomnos = getattr(data, "atomnos", None)
    if atomnos is not None:
        result["atom_count"] = len(atomnos)
    charge = getattr(data, "charge", None)
    if charge is not None:
        result["charge"] = int(charge)
    mult = getattr(data, "mult", None)
    if mult is not None:
        result["multiplicity"] = int(mult)
    optdone = getattr(data, "optdone", None)
    if isinstance(optdone, bool):
        result["optimization_converged"] = optdone
    return result


def _failure(status: str, exc: Exception) -> dict[str, Any]:
    detail = str(exc).strip()
    suffix = f": {detail}" if detail else ""
    return {
        "status": status,
        "warning": f"cclib {type(exc).__name__}{suffix}",
    }


__all__ = ["MAX_CCLIB_FILE_SIZE", "parse_cclib_output"]
