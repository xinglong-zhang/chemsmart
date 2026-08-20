"""Approved molecular-identity evidence for exact coordinate artifacts.

Molecular identity and electronic state are intentionally separate.  An
approved identity record may name the molecular system represented by exact
coordinate bytes, but it does not establish charge or multiplicity.  Those
remain explicit task facts bound through :class:`ScientificIdentityBindingV1`.
"""

from __future__ import annotations

import math
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

import yaml

from chemsmart.agent._contracts import (
    ContractError,
    canonical_data,
    canonical_sha256,
    file_sha256,
    require_identifier,
    require_sha256,
)

_SYMBOL = re.compile(r"^[A-Z][a-z]?$")
_UNITS = frozenset({"angstrom", "bohr"})
_MANIFEST_SCHEMA = "chemsmart.approved-molecular-input-manifest.v1"
_INPUT_SCHEMA = "chemsmart.approved-molecular-input.v1"
_MANIFEST_FIELDS = frozenset({"schema_version", "inputs"})
_INPUT_FIELDS = frozenset(
    {
        "input_id",
        "identity_id",
        "approved_names",
        "geometry_file",
        "geometry_sha256",
        "coordinate_units",
        "geometry_role",
        "charge",
        "multiplicity",
        "source_locator",
        "source_record_sha256",
        "state_source_locator",
    }
)


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
            raise ContractError(
                "unsupported approved molecular identity schema"
            )
        require_identifier(self.identity_id, "identity_id")
        require_sha256(self.geometry_sha256, "geometry_sha256")
        require_sha256(self.source_record_sha256, "source_record_sha256")
        if self.coordinate_units not in _UNITS:
            raise ContractError("coordinate_units must be angstrom or bohr")
        if not self.atom_order or any(
            _SYMBOL.fullmatch(symbol) is None for symbol in self.atom_order
        ):
            raise ContractError(
                "atom_order must contain canonical element symbols"
            )
        if not self.source_locator.strip():
            raise ContractError("molecular identity requires a source locator")
        if self.approval_scope != "user-approved-coordinate-identity":
            raise ContractError("molecular identity approval scope is invalid")
        if (
            not self.approved_names
            or self.approved_names != tuple(dict.fromkeys(self.approved_names))
            or any(not name.strip() for name in self.approved_names)
        ):
            raise ContractError(
                "approved_names must be non-empty and canonical"
            )
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


@dataclass(frozen=True)
class ApprovedMolecularInputV1:
    """User-approved role and electronic state for an exact identity frame.

    Molecular identity remains a separate record.  This assignment adds only
    the state and role that the user or source packet explicitly approved for
    those coordinate bytes; it is not permission to infer another state.
    """

    schema_version: str
    input_id: str
    molecular_identity: ApprovedMolecularIdentityV1
    geometry_role: str
    charge: int
    multiplicity: int
    state_source_locator: str
    assignment_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != _INPUT_SCHEMA:
            raise ContractError("unsupported approved molecular input schema")
        require_identifier(self.input_id, "input_id")
        if not isinstance(
            self.molecular_identity, ApprovedMolecularIdentityV1
        ):
            raise ContractError(
                "approved molecular input requires a typed molecular identity"
            )
        if not self.geometry_role.strip():
            raise ContractError(
                "approved molecular input requires a geometry role"
            )
        if isinstance(self.charge, bool) or not isinstance(self.charge, int):
            raise ContractError(
                "approved molecular input charge must be an integer"
            )
        if (
            isinstance(self.multiplicity, bool)
            or not isinstance(self.multiplicity, int)
            or self.multiplicity < 1
        ):
            raise ContractError(
                "approved molecular input multiplicity must be a positive integer"
            )
        if not self.state_source_locator.strip():
            raise ContractError(
                "approved molecular input requires a state source locator"
            )
        body = _approved_input_body(self)
        if self.assignment_sha256 != canonical_sha256(body):
            raise ContractError("approved molecular input digest mismatch")

    @property
    def evidence_ref(self) -> str:
        return f"molecular_input:{self.assignment_sha256}"

    def public_record(self) -> dict[str, Any]:
        """Return the path-free role/state record visible to the model."""

        return {
            "record_kind": "approved_molecular_input",
            **canonical_data(_approved_input_body(self)),
            "approved_names": canonical_data(
                self.molecular_identity.approved_names
            ),
            "coordinate_units": self.molecular_identity.coordinate_units,
            "atom_order": canonical_data(self.molecular_identity.atom_order),
            "source_locator": self.molecular_identity.source_locator,
            "identity_evidence_ref": self.molecular_identity.evidence_ref,
            "state_evidence_ref": self.evidence_ref,
            "electronic_state_status": "user_approved_for_exact_geometry",
            "assignment_sha256": self.assignment_sha256,
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


def build_approved_molecular_input(
    *,
    input_id: str,
    molecular_identity: ApprovedMolecularIdentityV1,
    geometry_role: str,
    charge: int,
    multiplicity: int,
    state_source_locator: str,
) -> ApprovedMolecularInputV1:
    """Bind a source-declared role and state to an approved identity frame."""

    body = {
        "schema_version": _INPUT_SCHEMA,
        "input_id": require_identifier(input_id, "input_id"),
        "molecular_identity_sha256": molecular_identity.identity_sha256,
        "geometry_sha256": molecular_identity.geometry_sha256,
        "geometry_role": str(geometry_role).strip(),
        "charge": charge,
        "multiplicity": multiplicity,
        "state_source_locator": str(state_source_locator).strip(),
    }
    return ApprovedMolecularInputV1(
        schema_version=_INPUT_SCHEMA,
        input_id=body["input_id"],
        molecular_identity=molecular_identity,
        geometry_role=body["geometry_role"],
        charge=charge,
        multiplicity=multiplicity,
        state_source_locator=body["state_source_locator"],
        assignment_sha256=canonical_sha256(body),
    )


def load_approved_molecular_input_manifest(
    path: str | Path,
    *,
    workspace: str | Path,
) -> tuple[ApprovedMolecularInputV1, ...]:
    """Load exact user-approved molecular inputs from YAML or JSON.

    Geometry paths are workspace-relative and may not traverse symlinks.  The
    loader derives coordinate digests and atom order from the exact files, so
    callers cannot authorize a name or state for different bytes by editing a
    manifest digest.
    """

    source = Path(path).expanduser().resolve()
    workspace_root = Path(workspace).expanduser()
    if not workspace_root.is_absolute():
        raise ContractError("identity manifest workspace must be absolute")
    if not workspace_root.is_dir() or workspace_root.is_symlink():
        raise ContractError(
            "identity manifest workspace must be a regular directory"
        )
    workspace_root = workspace_root.resolve()
    try:
        payload = yaml.safe_load(source.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, yaml.YAMLError) as exc:
        raise ContractError(
            "approved molecular input manifest is unreadable"
        ) from exc
    if not isinstance(payload, Mapping):
        raise ContractError(
            "approved molecular input manifest must be a mapping"
        )
    unknown_manifest = set(payload).difference(_MANIFEST_FIELDS)
    if unknown_manifest:
        raise ContractError(
            "approved molecular input manifest has unsupported fields: "
            + ", ".join(sorted(str(item) for item in unknown_manifest))
        )
    if payload.get("schema_version") != _MANIFEST_SCHEMA:
        raise ContractError(
            "unsupported approved molecular input manifest schema"
        )
    raw_inputs = payload.get("inputs")
    if (
        not isinstance(raw_inputs, list)
        or not raw_inputs
        or any(not isinstance(item, Mapping) for item in raw_inputs)
    ):
        raise ContractError(
            "approved molecular input manifest requires a non-empty inputs list"
        )

    approved: list[ApprovedMolecularInputV1] = []
    for ordinal, raw in enumerate(raw_inputs):
        unknown = set(raw).difference(_INPUT_FIELDS)
        if unknown:
            raise ContractError(
                f"approved molecular input {ordinal} has unsupported fields: "
                + ", ".join(sorted(str(item) for item in unknown))
            )
        required = _INPUT_FIELDS
        missing = tuple(
            sorted(field for field in required if field not in raw)
        )
        if missing:
            raise ContractError(
                f"approved molecular input {ordinal} is missing: "
                + ", ".join(missing)
            )
        geometry_file = str(raw["geometry_file"] or "").strip()
        relative = Path(geometry_file)
        if (
            not geometry_file
            or relative.is_absolute()
            or ".." in relative.parts
        ):
            raise ContractError(
                "approved molecular input geometry_file must be workspace-relative"
            )
        geometry = workspace_root.joinpath(relative)
        if geometry.is_symlink() or not geometry.is_file():
            raise ContractError(
                f"approved molecular input geometry is unavailable: {geometry_file}"
            )
        geometry = geometry.resolve()
        try:
            geometry.relative_to(workspace_root)
        except ValueError as exc:
            raise ContractError(
                "approved molecular input geometry escapes the workspace"
            ) from exc
        atom_order = _inspect_xyz_atom_order(geometry)
        geometry_sha256 = file_sha256(geometry)
        if geometry_sha256 != require_sha256(
            str(raw["geometry_sha256"]), "geometry_sha256"
        ):
            raise ContractError(
                f"approved molecular input geometry digest differs: {geometry_file}"
            )
        names = raw["approved_names"]
        if not isinstance(names, list):
            raise ContractError(
                "approved molecular input approved_names must be a list"
            )
        identity = build_approved_molecular_identity(
            identity_id=str(raw["identity_id"]),
            approved_names=tuple(str(item) for item in names),
            geometry_sha256=geometry_sha256,
            coordinate_units=str(raw["coordinate_units"]),
            atom_order=atom_order,
            source_locator=str(raw["source_locator"]),
            source_record_sha256=str(raw["source_record_sha256"]),
        )
        approved.append(
            build_approved_molecular_input(
                input_id=str(raw["input_id"]),
                molecular_identity=identity,
                geometry_role=str(raw["geometry_role"]),
                charge=raw["charge"],
                multiplicity=raw["multiplicity"],
                state_source_locator=str(raw["state_source_locator"]),
            )
        )

    input_ids = tuple(item.input_id for item in approved)
    assignment_sha256s = tuple(item.assignment_sha256 for item in approved)
    identity_ids = tuple(
        item.molecular_identity.identity_id for item in approved
    )
    if len(input_ids) != len(set(input_ids)):
        raise ContractError("approved molecular input IDs must be unique")
    if len(identity_ids) != len(set(identity_ids)):
        raise ContractError("approved molecular identity IDs must be unique")
    if len(assignment_sha256s) != len(set(assignment_sha256s)):
        raise ContractError(
            "approved molecular input assignments must be unique"
        )
    return tuple(sorted(approved, key=lambda item: item.input_id))


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
        raise ContractError(
            "approved molecular identity targets another geometry"
        )
    if identity.atom_order != tuple(str(item).strip() for item in atom_order):
        raise ContractError("approved molecular identity atom order differs")


def _without_field(value: Any, field: str) -> dict[str, Any]:
    return {key: item for key, item in asdict(value).items() if key != field}


def _approved_input_body(value: ApprovedMolecularInputV1) -> dict[str, Any]:
    return {
        "schema_version": value.schema_version,
        "input_id": value.input_id,
        "molecular_identity_sha256": value.molecular_identity.identity_sha256,
        "geometry_sha256": value.molecular_identity.geometry_sha256,
        "geometry_role": value.geometry_role,
        "charge": value.charge,
        "multiplicity": value.multiplicity,
        "state_source_locator": value.state_source_locator,
    }


def _inspect_xyz_atom_order(path: Path) -> tuple[str, ...]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
        count = int(lines[0].strip())
    except (IndexError, OSError, UnicodeDecodeError, ValueError) as exc:
        raise ContractError(
            "approved molecular input XYZ has a malformed atom count"
        ) from exc
    if count < 1 or len(lines) < count + 2:
        raise ContractError("approved molecular input XYZ is truncated")
    symbols: list[str] = []
    for row in lines[2 : count + 2]:
        fields = row.split()
        if len(fields) < 4 or _SYMBOL.fullmatch(fields[0]) is None:
            raise ContractError(
                "approved molecular input XYZ has a malformed atom row"
            )
        try:
            coordinates = tuple(float(item) for item in fields[1:4])
        except ValueError as exc:
            raise ContractError(
                "approved molecular input XYZ coordinates are not numeric"
            ) from exc
        if not all(math.isfinite(item) for item in coordinates):
            raise ContractError(
                "approved molecular input XYZ coordinates must be finite"
            )
        symbols.append(fields[0])
    return tuple(symbols)


__all__ = [
    "ApprovedMolecularInputV1",
    "ApprovedMolecularIdentityV1",
    "build_approved_molecular_input",
    "build_approved_molecular_identity",
    "load_approved_molecular_input_manifest",
    "validate_identity_for_geometry",
]
