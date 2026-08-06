"""Small dependency-free primitives for immutable agent contracts.

The v3.1.4 distribution intentionally does not depend on the v2 agent stack.
This module therefore keeps the reusable contract boundary in the standard
library.  It provides canonical JSON and content identities without making a
model response, Python object identity, or display string authoritative.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from dataclasses import asdict, dataclass, is_dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Mapping

_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_IDENTIFIER_RE = re.compile(r"^[a-z][a-z0-9_.-]*$")

#: Appended wherever a blank analysis unit is refused.  Seen in three separate
#: live cases, always on a count: "must not be empty" named neither the output
#: nor the value a dimensionless quantity should carry, leaving a caller to
#: guess between "", "none", "count" and "1".  Two planes raise this, so the
#: sentence lives here rather than being duplicated and kept in step by hand.
DIMENSIONLESS_UNIT_HINT = (
    "a dimensionless quantity such as a count, a population or an "
    "oscillator strength uses '1', not an empty string"
)


class ContractError(ValueError):
    """Raised when an agent contract is internally inconsistent."""


def require_identifier(value: str, field_name: str) -> str:
    """Return a normalized public identifier or raise ``ContractError``."""

    normalized = str(value or "").strip().lower()
    if _IDENTIFIER_RE.fullmatch(normalized) is None:
        # This validator stands behind every identifier the model supplies, so
        # its message is the one seen for most malformed IDs. Naming the field
        # alone leaves the caller guessing which of several values was wrong
        # and what shape was expected, so quote the value and the rule.
        detail = "it is empty" if not normalized else f"got {value!r}"
        raise ContractError(
            f"{field_name} must be a lower-case public identifier "
            "(start with a letter, then letters, digits, '_', '.' or '-'); "
            f"{detail}"
        )
    return normalized


def require_sha256(value: str, field_name: str) -> str:
    """Validate an exact lower-case SHA-256 string."""

    normalized = str(value or "").strip().lower()
    if _SHA256_RE.fullmatch(normalized) is None:
        raise ContractError(f"{field_name} must be a SHA-256 digest")
    return normalized


def _canonical_value(value: Any) -> Any:
    """Convert supported values to a deterministic JSON-compatible shape."""

    if is_dataclass(value):
        return _canonical_value(asdict(value))
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {
            str(key): _canonical_value(item)
            for key, item in sorted(
                value.items(), key=lambda pair: str(pair[0])
            )
        }
    if isinstance(value, (tuple, list)):
        return [_canonical_value(item) for item in value]
    if isinstance(value, (set, frozenset)):
        return sorted(_canonical_value(item) for item in value)
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ContractError(
                "canonical records cannot contain NaN or infinity"
            )
        return value
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise ContractError(
        f"unsupported canonical record value: {type(value).__name__}"
    )


def canonical_data(value: Any) -> Any:
    """Return a detached canonical representation of ``value``."""

    return _canonical_value(value)


def canonical_json(value: Any) -> str:
    """Serialize ``value`` with the repository-wide canonical JSON shape."""

    return json.dumps(
        canonical_data(value),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    )


def canonical_sha256(value: Any) -> str:
    """Return the SHA-256 identity of a canonical record."""

    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def file_sha256(path: str | Path) -> str:
    """Hash a file without parsing or executing it."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass(frozen=True)
class TrustedArtifactRefV1:
    """Host-resolved artifact binding; never populated from model paths."""

    artifact_id: str
    kind: str
    sha256: str
    size_bytes: int
    path: str
    cli_value: str

    def __post_init__(self) -> None:
        if not str(self.artifact_id).strip():
            raise ContractError("artifact_id must not be empty")
        require_identifier(self.kind, "kind")
        require_sha256(self.sha256, "sha256")
        if self.size_bytes < 0:
            raise ContractError("size_bytes must be non-negative")
        if not Path(self.path).is_absolute():
            raise ContractError("trusted artifact path must be absolute")
        if not str(self.cli_value).strip():
            raise ContractError("cli_value must not be empty")


__all__ = [
    "DIMENSIONLESS_UNIT_HINT",
    "ContractError",
    "TrustedArtifactRefV1",
    "canonical_data",
    "canonical_json",
    "canonical_sha256",
    "file_sha256",
    "require_identifier",
    "require_sha256",
]
