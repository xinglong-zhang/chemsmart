"""Public, arm-neutral context for provider program-management campaigns.

The provider needs semantic access to host-bound evidence, but it must never
receive filesystem paths, credentials, executable claims, or host-owned
readiness decisions.  This module therefore separates a small public packet
from the richer objects passed to :class:`CommandCompiledToolHostV1`.

Planning follows progressive assurance: an exact external artifact carries a
digest, while a future producer output remains an unresolved semantic slot.
The unresolved slot is useful for a command DAG, but cannot be compiled,
previewed, approved, or executed until the producer emits an exact artifact
receipt.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, is_dataclass
from enum import Enum
from types import MappingProxyType
from typing import Any, Mapping

from chemsmart.agent._contracts import (
    ContractError,
    TrustedArtifactRefV1,
    canonical_data,
    canonical_json,
    canonical_sha256,
    require_identifier,
    require_sha256,
)


# Public next-action hints must remain executable in both paired arms. The
# program-management inspection tools are deliberately absent: naming them in
# the arm-neutral packet would tell H0 to call tools it cannot see and would
# confound the H0/H1 comparison.
_COMMON_CAMPAIGN_ACTION_TOOLS = frozenset(
    {
        "inspect_calculation_artifact",
        "plan_command_workflow",
        "preflight_program_node",
        "preview_command",
        "read_project_yaml",
        "render_project_yaml",
        "repair_command",
        "synthesize_command",
        "validate_project_yaml",
    }
)


@dataclass(frozen=True)
class CampaignPublicArtifactV1:
    """A path-free view of one exact host-bound artifact."""

    artifact_id: str
    artifact_class: str
    sha256: str
    size_bytes: int
    purpose: str
    provenance_status: str

    def __post_init__(self) -> None:
        require_identifier(self.artifact_id, "artifact_id")
        require_identifier(self.artifact_class, "artifact_class")
        require_sha256(self.sha256, "artifact sha256")
        if self.size_bytes < 0:
            raise ContractError("public artifact size must be non-negative")
        if not self.purpose.strip():
            raise ContractError("public artifact purpose is required")
        if self.provenance_status not in {
            "approved_exact_source",
            "checked_in_non_generated_source",
            "host_bound_project_candidate",
            "validated_project_fixture",
            "seeded_negative_fixture",
        }:
            raise ContractError("invalid campaign artifact provenance status")
        _reject_private_or_path_shaped_text(self.purpose)


@dataclass(frozen=True)
class CampaignArtifactSlotV1:
    """An exact external input or an unresolved producer-output dependency."""

    slot_id: str
    artifact_class: str
    state: str
    artifact_id: str
    artifact_sha256: str
    producer_node_id: str
    producer_output_id: str

    def __post_init__(self) -> None:
        require_identifier(self.slot_id, "slot_id")
        require_identifier(self.artifact_class, "artifact_class")
        if self.state not in {"external_bound", "future_unresolved"}:
            raise ContractError("invalid campaign artifact-slot state")
        if self.state == "external_bound":
            require_identifier(self.artifact_id, "artifact_id")
            require_sha256(self.artifact_sha256, "artifact_sha256")
            if self.producer_node_id or self.producer_output_id:
                raise ContractError("external artifact cannot name a producer")
        else:
            if self.artifact_id or self.artifact_sha256:
                raise ContractError("future artifact slot cannot invent identity")
            require_identifier(self.producer_node_id, "producer_node_id")
            require_identifier(self.producer_output_id, "producer_output_id")


@dataclass(frozen=True)
class CampaignPublicReceiptRefV1:
    """A concise, non-authoritative projection of a host receipt."""

    role: str
    receipt_sha256: str
    status: str
    rule_ids: tuple[str, ...]
    semantic_fields: tuple[tuple[str, Any], ...] = ()

    def __post_init__(self) -> None:
        require_identifier(self.role, "receipt role")
        require_sha256(self.receipt_sha256, "receipt_sha256")
        if not self.status.strip():
            raise ContractError("public receipt status is required")
        if self.rule_ids != tuple(sorted(set(self.rule_ids))):
            raise ContractError("public receipt rules must be sorted and unique")
        names = tuple(item[0] for item in self.semantic_fields)
        if names != tuple(sorted(set(names))):
            raise ContractError("receipt semantic fields must be sorted and unique")
        _reject_private_or_path_shaped_value(self.semantic_fields)
        canonical_data(self.semantic_fields)


@dataclass(frozen=True)
class CampaignToolInputV1:
    """A typed identifier packet for one permitted next host operation."""

    action_id: str
    tool_name: str
    fields: tuple[tuple[str, Any], ...]
    assurance_level: str

    def __post_init__(self) -> None:
        require_identifier(self.action_id, "action_id")
        require_identifier(self.tool_name, "tool_name")
        names = tuple(item[0] for item in self.fields)
        if names != tuple(sorted(set(names))):
            raise ContractError("tool-input fields must be sorted and unique")
        if self.assurance_level not in {
            "planning",
            "inspection",
            "compile_bound",
            "verification",
        }:
            raise ContractError("invalid tool-input assurance level")
        if self.tool_name not in _COMMON_CAMPAIGN_ACTION_TOOLS:
            raise ContractError(
                "public context can name only tools shared by H0 and H1"
            )
        _reject_private_or_path_shaped_value(self.fields)
        canonical_data(self.fields)


@dataclass(frozen=True)
class CampaignPublicContextV1:
    """Canonical context shown byte-for-byte identically to H0 and H1."""

    schema_version: str
    case_id: str
    task_spec_sha256: str
    registry_sha256: str
    live_cli_schema_sha256: str
    joined_capability_sha256: str
    support_overlay_sha256: str
    scientific_facts: tuple[tuple[str, Any], ...]
    artifacts: tuple[CampaignPublicArtifactV1, ...]
    artifact_slots: tuple[CampaignArtifactSlotV1, ...]
    receipts: tuple[CampaignPublicReceiptRefV1, ...]
    next_actions: tuple[CampaignToolInputV1, ...]
    expected_artifact_classes: tuple[str, ...]
    missing_evidence: tuple[str, ...]
    host_state_sha256: str
    packet_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.campaign-public-context.v1":
            raise ContractError("unsupported campaign public-context schema")
        if not self.case_id.startswith("DS-PM-"):
            raise ContractError("campaign context case ID is not canonical")
        for field_name, digest in (
            ("task_spec_sha256", self.task_spec_sha256),
            ("registry_sha256", self.registry_sha256),
            ("live_cli_schema_sha256", self.live_cli_schema_sha256),
            ("joined_capability_sha256", self.joined_capability_sha256),
            ("support_overlay_sha256", self.support_overlay_sha256),
            ("host_state_sha256", self.host_state_sha256),
        ):
            require_sha256(digest, field_name)
        _require_sorted_unique_objects(
            self.artifacts, lambda item: item.artifact_id, "public artifacts"
        )
        _require_sorted_unique_objects(
            self.artifact_slots, lambda item: item.slot_id, "artifact slots"
        )
        _require_sorted_unique_objects(
            self.receipts,
            lambda item: (item.role, item.receipt_sha256),
            "public receipt references",
        )
        _require_sorted_unique_objects(
            self.next_actions, lambda item: item.action_id, "next actions"
        )
        fact_names = tuple(item[0] for item in self.scientific_facts)
        if fact_names != tuple(sorted(set(fact_names))):
            raise ContractError("scientific facts must be sorted and unique")
        for field_name, values in (
            ("expected artifact classes", self.expected_artifact_classes),
            ("missing evidence", self.missing_evidence),
        ):
            if values != tuple(sorted(set(values))):
                raise ContractError(f"{field_name} must be sorted and unique")
        _reject_private_or_path_shaped_value(
            (
                self.scientific_facts,
                self.artifacts,
                self.artifact_slots,
                self.receipts,
                self.next_actions,
                self.missing_evidence,
            )
        )
        body = _public_context_body(self)
        if self.packet_sha256 != canonical_sha256(body):
            raise ContractError("campaign public-context digest mismatch")

    def provider_message(self) -> dict[str, str]:
        """Return the exact public system message included in plan hashes."""

        return {
            "role": "system",
            "content": (
                "CHEMSMART_CAMPAIGN_PUBLIC_CONTEXT_V1\n"
                + canonical_json(self)
            ),
        }


@dataclass(frozen=True)
class CampaignHostFixtureV1:
    """One arm-neutral set of constructor arguments and its public projection."""

    schema_version: str
    case_id: str
    host_inputs: Mapping[str, Any]
    public_context: CampaignPublicContextV1
    host_state_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.campaign-host-fixture.v1":
            raise ContractError("unsupported campaign host-fixture schema")
        if self.case_id != self.public_context.case_id:
            raise ContractError("fixture and public context target different cases")
        observed = host_state_sha256(self.host_inputs)
        if self.host_state_sha256 != observed:
            raise ContractError("campaign host-state digest mismatch")
        if self.public_context.host_state_sha256 != observed:
            raise ContractError("public context does not bind the host state")


def build_campaign_public_context(
    *,
    case_id: str,
    task_spec_sha256: str,
    registry_sha256: str,
    live_cli_schema_sha256: str,
    joined_capability_sha256: str,
    support_overlay_sha256: str,
    scientific_facts: Mapping[str, Any],
    artifacts: tuple[CampaignPublicArtifactV1, ...],
    artifact_slots: tuple[CampaignArtifactSlotV1, ...],
    receipts: tuple[CampaignPublicReceiptRefV1, ...],
    next_actions: tuple[CampaignToolInputV1, ...],
    expected_artifact_classes: tuple[str, ...],
    missing_evidence: tuple[str, ...],
    host_state_sha256: str,
) -> CampaignPublicContextV1:
    body = {
        "schema_version": "chemsmart.campaign-public-context.v1",
        "case_id": case_id,
        "task_spec_sha256": task_spec_sha256,
        "registry_sha256": registry_sha256,
        "live_cli_schema_sha256": live_cli_schema_sha256,
        "joined_capability_sha256": joined_capability_sha256,
        "support_overlay_sha256": support_overlay_sha256,
        "scientific_facts": tuple(sorted(scientific_facts.items())),
        "artifacts": tuple(sorted(artifacts, key=lambda item: item.artifact_id)),
        "artifact_slots": tuple(
            sorted(artifact_slots, key=lambda item: item.slot_id)
        ),
        "receipts": tuple(
            sorted(receipts, key=lambda item: (item.role, item.receipt_sha256))
        ),
        "next_actions": tuple(
            sorted(next_actions, key=lambda item: item.action_id)
        ),
        "expected_artifact_classes": tuple(
            sorted(set(expected_artifact_classes))
        ),
        "missing_evidence": tuple(sorted(set(missing_evidence))),
        "host_state_sha256": host_state_sha256,
    }
    return CampaignPublicContextV1(
        **body, packet_sha256=canonical_sha256(body)
    )


def host_state_sha256(host_inputs: Mapping[str, Any]) -> str:
    """Digest the host constructor state without exposing its private values."""

    forbidden = {"event_store", "support_overlay", "tool_surface"}
    overlap = forbidden.intersection(host_inputs)
    if overlap:
        raise ContractError(
            "campaign fixture cannot replace runtime authority: "
            + ", ".join(sorted(overlap))
        )
    allowed = {
        "artifacts",
        "scientific_identities",
        "environment_targets",
        "compute_environment_receipts",
        "component_conformance_receipts",
        "settings_objects",
        "run_receipts",
        "scientific_claim_evidence",
        "functional_equivalence_receipts",
        "substitution_approvals",
        "capability_receipts",
        "environment_receipts",
        "program_binding_receipts",
        "engine_binding_receipts",
        "project_validation_receipts",
        "registry",
        "live_schema",
    }
    unknown = sorted(set(host_inputs).difference(allowed))
    if unknown:
        raise ContractError(
            "campaign fixture has unknown host constructor inputs: "
            + ", ".join(unknown)
        )
    projection = []
    for name, value in sorted(host_inputs.items()):
        if name == "registry":
            projected = {"registry_sha256": value.registry_sha256}
        elif name == "live_schema":
            projected = {"schema_sha256": value.schema_sha256}
        elif isinstance(value, Mapping):
            projected = tuple(
                (str(key), _object_sha256(item))
                for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
            )
        elif isinstance(value, tuple):
            projected = tuple(_object_sha256(item) for item in value)
        else:
            projected = _object_sha256(value)
        projection.append((name, projected))
    return canonical_sha256(tuple(projection))


def immutable_host_inputs(values: Mapping[str, Any]) -> Mapping[str, Any]:
    """Return a shallow immutable view used by both paired arms."""

    return MappingProxyType(dict(values))


def public_artifact(
    artifact: TrustedArtifactRefV1,
    *,
    purpose: str,
    provenance_status: str,
) -> CampaignPublicArtifactV1:
    return CampaignPublicArtifactV1(
        artifact_id=artifact.artifact_id,
        artifact_class=artifact.kind,
        sha256=artifact.sha256,
        size_bytes=artifact.size_bytes,
        purpose=purpose,
        provenance_status=provenance_status,
    )


def _public_context_body(value: CampaignPublicContextV1) -> dict[str, Any]:
    return {
        key: item
        for key, item in asdict(value).items()
        if key != "packet_sha256"
    }


def _object_sha256(value: Any) -> str:
    if is_dataclass(value) or isinstance(
        value, (Mapping, tuple, list, str, int, float, bool, type(None), Enum)
    ):
        return canonical_sha256(value)
    attributes = getattr(value, "__dict__", None)
    if isinstance(attributes, Mapping):
        public = {
            (key[1:] if key == "_project_yaml_digest" else key): item
            for key, item in attributes.items()
            if not str(key).startswith("_") or key == "_project_yaml_digest"
        }
        return canonical_sha256(public)
    raise ContractError(
        f"cannot bind host fixture object of type {type(value).__name__}"
    )


def _require_sorted_unique_objects(
    values: tuple[Any, ...], key, field_name: str
) -> None:
    observed = tuple(key(item) for item in values)
    if observed != tuple(sorted(set(observed))):
        raise ContractError(f"{field_name} must be sorted and unique")


def _reject_private_or_path_shaped_value(value: Any) -> None:
    if is_dataclass(value):
        value = asdict(value)
    if isinstance(value, Mapping):
        for key, item in value.items():
            lowered = str(key).lower()
            if any(
                marker in lowered
                for marker in (
                    "api_key",
                    "authorization",
                    "credential",
                    "reasoning_content",
                    "secret",
                )
            ):
                raise ContractError("public context contains a private field")
            _reject_private_or_path_shaped_value(item)
        return
    if isinstance(value, (tuple, list)):
        for item in value:
            _reject_private_or_path_shaped_value(item)
        return
    if isinstance(value, str):
        _reject_private_or_path_shaped_text(value)


def _reject_private_or_path_shaped_text(value: str) -> None:
    text = str(value)
    if text.startswith(("/", "~", "file://")) or "\\Users\\" in text:
        raise ContractError("public context cannot expose a filesystem path")
    if "Bearer " in text:
        raise ContractError("public context cannot expose authorization data")


__all__ = [
    "CampaignArtifactSlotV1",
    "CampaignHostFixtureV1",
    "CampaignPublicArtifactV1",
    "CampaignPublicContextV1",
    "CampaignPublicReceiptRefV1",
    "CampaignToolInputV1",
    "build_campaign_public_context",
    "host_state_sha256",
    "immutable_host_inputs",
    "public_artifact",
]
