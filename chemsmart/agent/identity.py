"""Approved molecular-identity evidence for exact coordinate artifacts.

Molecular identity and electronic state are intentionally separate.  An
approved identity record may name the molecular system represented by exact
coordinate bytes, but it does not establish charge or multiplicity.  Those
remain explicit task facts bound through :class:`ScientificIdentityBindingV1`.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import re
from typing import Any, Iterable

from chemsmart.agent._contracts import (
    ContractError,
    canonical_data,
    canonical_sha256,
    require_identifier,
    require_sha256,
)


_SYMBOL = re.compile(r"^[A-Z][a-z]?$")
_UNITS = frozenset({"angstrom", "bohr"})


@dataclass(frozen=True)
class ApprovedMolecularIdentityV1:
    """Path-free identity evidence bound to one exact coordinate frame."""

    schema_version: str
    identity_id: str
    approved_names: tuple[str, ...]
    geometry_sha256: str
    coordinate_units: str
    atom_order: tuple[str, ...]
    source_locator: str
    source_record_sha256: str
    approval_scope: str
    identity_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.approved-molecular-identity.v1":
            raise ContractError("unsupported approved molecular identity schema")
        require_identifier(self.identity_id, "identity_id")
        require_sha256(self.geometry_sha256, "geometry_sha256")
        require_sha256(self.source_record_sha256, "source_record_sha256")
        if self.coordinate_units not in _UNITS:
            raise ContractError("coordinate_units must be angstrom or bohr")
        if not self.atom_order or any(
            _SYMBOL.fullmatch(symbol) is None for symbol in self.atom_order
        ):
            raise ContractError("atom_order must contain canonical element symbols")
        if not self.source_locator.strip():
            raise ContractError("molecular identity requires a source locator")
        if self.approval_scope != "user-approved-coordinate-identity":
            raise ContractError("molecular identity approval scope is invalid")
        if (
            not self.approved_names
            or self.approved_names
            != tuple(dict.fromkeys(self.approved_names))
            or any(not name.strip() for name in self.approved_names)
        ):
            raise ContractError("approved_names must be non-empty and canonical")
        body = _without_field(self, "identity_sha256")
        if self.identity_sha256 != canonical_sha256(body):
            raise ContractError("approved molecular identity digest mismatch")

    @property
    def evidence_ref(self) -> str:
        return f"molecular_identity:{self.identity_sha256}"

    def public_record(self) -> dict[str, Any]:
        return {
            "record_kind": "approved_molecular_identity",
            **canonical_data(asdict(self)),
            "evidence_ref": self.evidence_ref,
            "electronic_state_status": "not_established_by_identity_record",
        }


def build_approved_molecular_identity(
    *,
    identity_id: str,
    approved_names: Iterable[str],
    geometry_sha256: str,
    coordinate_units: str,
    atom_order: Iterable[str],
    source_locator: str,
    source_record_sha256: str,
) -> ApprovedMolecularIdentityV1:
    """Build explicit identity evidence without inferring electronic state."""

    names = tuple(dict.fromkeys(str(item).strip() for item in approved_names))
    body = {
        "schema_version": "chemsmart.approved-molecular-identity.v1",
        "identity_id": require_identifier(identity_id, "identity_id"),
        "approved_names": names,
        "geometry_sha256": require_sha256(geometry_sha256, "geometry_sha256"),
        "coordinate_units": str(coordinate_units).strip().lower(),
        "atom_order": tuple(str(item).strip() for item in atom_order),
        "source_locator": str(source_locator).strip(),
        "source_record_sha256": require_sha256(
            source_record_sha256, "source_record_sha256"
        ),
        "approval_scope": "user-approved-coordinate-identity",
    }
    return ApprovedMolecularIdentityV1(
        **body, identity_sha256=canonical_sha256(body)
    )


def validate_identity_for_geometry(
    identity: ApprovedMolecularIdentityV1,
    *,
    geometry_sha256: str,
    atom_order: Iterable[str],
) -> None:
    """Reject identity substitution or atom reordering before model access."""

    if identity.geometry_sha256 != require_sha256(
        geometry_sha256, "geometry_sha256"
    ):
        raise ContractError("approved molecular identity targets another geometry")
    if identity.atom_order != tuple(str(item).strip() for item in atom_order):
        raise ContractError("approved molecular identity atom order differs")


def _without_field(value: Any, field: str) -> dict[str, Any]:
    return {key: item for key, item in asdict(value).items() if key != field}


__all__ = [
    "ApprovedMolecularIdentityV1",
    "build_approved_molecular_identity",
    "validate_identity_for_geometry",
]
