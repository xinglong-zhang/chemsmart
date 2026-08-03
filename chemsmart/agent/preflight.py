"""One canonical, cross-bound program-node readiness authority."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from chemsmart.agent._contracts import (
    ContractError,
    canonical_sha256,
    require_identifier,
    require_sha256,
)
from chemsmart.agent.capabilities import (
    CapabilityQueryReceiptV1,
    CapabilityQueryStatus,
    ResolvedEngineBindingV1,
    ResolvedProgramBindingV1,
)
from chemsmart.agent.commands import (
    CanonicalCommandInvocationV1,
    CommandInspectionReceiptV1,
    ScientificIdentityBindingV1,
)
from chemsmart.agent.preview import SafePreviewReceiptV1
from chemsmart.agent.projects import ProjectValidationReceiptV1


@dataclass(frozen=True)
class ProgramValidatorReceiptV1:
    """Typed, target-bound deterministic validator result."""

    schema_version: str
    node_id: str
    program: str
    invocation_sha256: str
    scientific_identity_sha256: str
    source_receipt_sha256: str
    validator_id: str
    status: str
    critical_finding_sha256s: tuple[str, ...]
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.program-validator-receipt.v1":
            raise ContractError("unsupported program validator receipt schema")
        require_identifier(self.node_id, "node_id")
        require_identifier(self.program, "program")
        require_identifier(self.validator_id, "validator_id")
        for field_name, digest in (
            ("invocation_sha256", self.invocation_sha256),
            ("scientific_identity_sha256", self.scientific_identity_sha256),
            ("source_receipt_sha256", self.source_receipt_sha256),
        ):
            require_sha256(digest, field_name)
        if self.status not in {"valid", "invalid"}:
            raise ContractError("invalid program validator status")
        _require_digest_set(
            self.critical_finding_sha256s, "critical finding"
        )
        if self.status == "valid" and self.critical_finding_sha256s:
            raise ContractError("a valid validator cannot carry critical findings")
        body = {
            "schema_version": self.schema_version,
            "node_id": self.node_id,
            "program": self.program,
            "invocation_sha256": self.invocation_sha256,
            "scientific_identity_sha256": self.scientific_identity_sha256,
            "source_receipt_sha256": self.source_receipt_sha256,
            "validator_id": self.validator_id,
            "status": self.status,
            "critical_finding_sha256s": self.critical_finding_sha256s,
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("program validator receipt digest mismatch")


def validator_receipt_from_safe_preview(
    *,
    node_id: str,
    program: str,
    scientific_identity_sha256: str,
    safe_preview: SafePreviewReceiptV1,
) -> ProgramValidatorReceiptV1:
    """Project the host-produced preview validator into a preflight gate."""

    critical = set(safe_preview.critical_finding_sha256s)
    if safe_preview.status != "previewed" or (
        safe_preview.program_validation_status != "valid"
    ):
        critical.update(
            canonical_sha256(
                {
                    "safe_preview_receipt_sha256": safe_preview.receipt_sha256,
                    "rule_id": rule_id,
                }
            )
            for rule_id in safe_preview.rule_ids
        )
    body = {
        "schema_version": "chemsmart.program-validator-receipt.v1",
        "node_id": require_identifier(node_id, "node_id"),
        "program": require_identifier(program, "program"),
        "invocation_sha256": safe_preview.invocation_sha256,
        "scientific_identity_sha256": require_sha256(
            scientific_identity_sha256, "scientific_identity_sha256"
        ),
        "source_receipt_sha256": safe_preview.receipt_sha256,
        "validator_id": "safe_preview_validator.v1",
        "status": "valid" if not critical else "invalid",
        "critical_finding_sha256s": tuple(sorted(critical)),
    }
    return ProgramValidatorReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class ProgramNodePreflightRequestV1:
    schema_version: str
    node_id: str
    capability_receipt_sha256: str
    program_binding_sha256: str
    engine_binding_sha256: str
    environment_receipt_sha256: str
    geometry_artifact_sha256: str
    scientific_identity_sha256: str
    charge: int
    multiplicity: int
    project_receipt_sha256: str
    invocation_sha256: str
    command_inspection_sha256: str
    validator_receipt_sha256s: tuple[str, ...]
    request_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.program-node-preflight-request.v1":
            raise ContractError("unsupported node preflight request schema")
        require_identifier(self.node_id, "node_id")
        for field_name, digest in (
            ("capability_receipt_sha256", self.capability_receipt_sha256),
            ("program_binding_sha256", self.program_binding_sha256),
            ("engine_binding_sha256", self.engine_binding_sha256),
            ("environment_receipt_sha256", self.environment_receipt_sha256),
            ("geometry_artifact_sha256", self.geometry_artifact_sha256),
            ("scientific_identity_sha256", self.scientific_identity_sha256),
            ("invocation_sha256", self.invocation_sha256),
            ("command_inspection_sha256", self.command_inspection_sha256),
        ):
            require_sha256(digest, field_name)
        if self.project_receipt_sha256:
            require_sha256(
                self.project_receipt_sha256, "project_receipt_sha256"
            )
        if self.multiplicity < 1:
            raise ContractError("multiplicity must be positive")
        _require_digest_set(
            self.validator_receipt_sha256s, "validator receipt"
        )
        expected = canonical_sha256(_request_body(self))
        if self.request_sha256 != expected:
            raise ContractError("node preflight request digest mismatch")


@dataclass(frozen=True)
class ProgramNodePreflightReceiptV1:
    schema_version: str
    request_sha256: str
    program: str
    jobtype: str
    safe_preview_receipt_sha256: str
    effective_charge: int | None
    effective_multiplicity: int | None
    effective_settings_sha256: str
    validator_receipt_sha256s: tuple[str, ...]
    critical_finding_sha256s: tuple[str, ...]
    plan_state: str
    execution_ready: bool
    rule_ids: tuple[str, ...]
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.program-node-preflight-receipt.v1":
            raise ContractError("unsupported node preflight receipt schema")
        if self.plan_state not in {
            "blocked",
            "ready_for_safe_preview",
            "previewed",
        }:
            raise ContractError("invalid node preflight state")
        if self.execution_ready and self.plan_state != "previewed":
            raise ContractError("execution readiness requires preview evidence")
        if self.plan_state == "previewed" and not self.safe_preview_receipt_sha256:
            raise ContractError("previewed state requires a preview receipt")
        body = {
            "schema_version": self.schema_version,
            "request_sha256": self.request_sha256,
            "program": self.program,
            "jobtype": self.jobtype,
            "safe_preview_receipt_sha256": self.safe_preview_receipt_sha256,
            "effective_charge": self.effective_charge,
            "effective_multiplicity": self.effective_multiplicity,
            "effective_settings_sha256": self.effective_settings_sha256,
            "validator_receipt_sha256s": self.validator_receipt_sha256s,
            "critical_finding_sha256s": self.critical_finding_sha256s,
            "plan_state": self.plan_state,
            "execution_ready": self.execution_ready,
            "rule_ids": self.rule_ids,
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("node preflight receipt digest mismatch")


def build_program_node_preflight_request(
    *,
    node_id: str,
    capability_receipt_sha256: str,
    program_binding_sha256: str,
    engine_binding_sha256: str,
    environment_receipt_sha256: str,
    geometry_artifact_sha256: str,
    scientific_identity_sha256: str,
    charge: int,
    multiplicity: int,
    project_receipt_sha256: str,
    invocation_sha256: str,
    command_inspection_sha256: str,
    validator_receipts: tuple[ProgramValidatorReceiptV1, ...],
) -> ProgramNodePreflightRequestV1:
    body = {
        "schema_version": "chemsmart.program-node-preflight-request.v1",
        "node_id": require_identifier(node_id, "node_id"),
        "capability_receipt_sha256": capability_receipt_sha256,
        "program_binding_sha256": program_binding_sha256,
        "engine_binding_sha256": engine_binding_sha256,
        "environment_receipt_sha256": environment_receipt_sha256,
        "geometry_artifact_sha256": geometry_artifact_sha256,
        "scientific_identity_sha256": scientific_identity_sha256,
        "charge": int(charge),
        "multiplicity": int(multiplicity),
        "project_receipt_sha256": project_receipt_sha256,
        "invocation_sha256": invocation_sha256,
        "command_inspection_sha256": command_inspection_sha256,
        "validator_receipt_sha256s": tuple(
            sorted({item.receipt_sha256 for item in validator_receipts})
        ),
    }
    return ProgramNodePreflightRequestV1(
        **body, request_sha256=canonical_sha256(body)
    )


def evaluate_program_node_preflight(
    request: ProgramNodePreflightRequestV1,
    *,
    capability: CapabilityQueryReceiptV1,
    program_binding: ResolvedProgramBindingV1,
    engine_binding: ResolvedEngineBindingV1,
    project: ProjectValidationReceiptV1 | None,
    invocation: CanonicalCommandInvocationV1,
    command_inspection: CommandInspectionReceiptV1,
    scientific_identity: ScientificIdentityBindingV1,
    validator_receipts: tuple[ProgramValidatorReceiptV1, ...],
    safe_preview: SafePreviewReceiptV1 | None = None,
) -> ProgramNodePreflightReceiptV1:
    """Cross-check every receipt edge and derive effective electronic state."""

    failures = []
    expected_bindings = (
        (request.capability_receipt_sha256, capability.receipt_sha256),
        (request.program_binding_sha256, program_binding.binding_sha256),
        (request.engine_binding_sha256, engine_binding.binding_sha256),
        (
            request.environment_receipt_sha256,
            engine_binding.environment_receipt_sha256,
        ),
        (request.invocation_sha256, invocation.invocation_sha256),
        (
            request.command_inspection_sha256,
            command_inspection.receipt_sha256,
        ),
    )
    if any(expected != observed for expected, observed in expected_bindings):
        failures.append("preflight.request.binding_mismatch")
    if program_binding.capability_receipt_sha256 != capability.receipt_sha256:
        failures.append("preflight.program.capability_mismatch")
    if engine_binding.capability_receipt_sha256 != capability.receipt_sha256:
        failures.append("preflight.engine.capability_mismatch")
    if engine_binding.program_binding_sha256 != program_binding.binding_sha256:
        failures.append("preflight.engine.program_binding_mismatch")
    if invocation.program_engine_binding_sha256 != engine_binding.binding_sha256:
        failures.append("preflight.invocation.engine_binding_mismatch")
    if invocation.joined_capability_sha256 != capability.joined_capability_sha256:
        failures.append("preflight.invocation.capability_mismatch")
    if (program_binding.selected_program, engine_binding.program) != (
        capability.query.program,
        capability.query.program,
    ):
        failures.append("preflight.selected_program_mismatch")
    if engine_binding.selected_engine != capability.query.engine:
        failures.append("preflight.selected_engine_mismatch")

    observed_project = project.receipt_sha256 if project is not None else ""
    if request.project_receipt_sha256 != observed_project:
        failures.append("preflight.request.project_mismatch")
    if project is not None:
        if project.capability_receipt_sha256 != capability.receipt_sha256:
            failures.append("preflight.project.capability_mismatch")
        if project.status != "valid":
            failures.append("preflight.project.red")
    if request.geometry_artifact_sha256 != invocation.input_sha256:
        failures.append("preflight.request.geometry_mismatch")
    if request.scientific_identity_sha256 != scientific_identity.binding_sha256:
        failures.append("preflight.request.scientific_identity_mismatch")
    if invocation.scientific_identity_sha256 != scientific_identity.binding_sha256:
        failures.append("preflight.invocation.scientific_identity_mismatch")
    if scientific_identity.geometry_artifact_sha256 != request.geometry_artifact_sha256:
        failures.append("preflight.identity.geometry_mismatch")
    if (request.charge, request.multiplicity) != (
        scientific_identity.charge,
        scientific_identity.multiplicity,
    ):
        failures.append("preflight.identity.electronic_state_mismatch")
    if capability.status not in {
        CapabilityQueryStatus.SUPPORTED,
        CapabilityQueryStatus.PREVIEW_ONLY,
    }:
        failures.append("preflight.capability.red")
    if program_binding.state == "blocked" or engine_binding.state == "blocked":
        failures.append("preflight.program_or_engine.red")
    requires_project = bool(
        capability.capability is not None
        and capability.capability.requires_project_configuration
    )
    if requires_project and project is None:
        failures.append("preflight.project.required")
    if command_inspection.status != "valid":
        failures.append("preflight.command.red")
    observed_validator_digests = tuple(
        sorted({item.receipt_sha256 for item in validator_receipts})
    )
    if request.validator_receipt_sha256s != observed_validator_digests:
        failures.append("preflight.validator.receipt_mismatch")
    if not validator_receipts:
        failures.append("preflight.validator.receipt_missing")
    for validator in validator_receipts:
        if (
            validator.node_id != request.node_id
            or validator.program != capability.query.program
            or validator.invocation_sha256 != invocation.invocation_sha256
            or validator.scientific_identity_sha256
            != scientific_identity.binding_sha256
        ):
            failures.append("preflight.validator.target_mismatch")
        if validator.status != "valid":
            failures.append("preflight.validator.red")
    critical_finding_sha256s = tuple(
        sorted(
            {
                digest
                for validator in validator_receipts
                for digest in validator.critical_finding_sha256s
            }
        )
    )
    if critical_finding_sha256s:
        failures.append("preflight.validator.critical_findings")

    effective, state_rules = _effective_electronic_state(
        request=request,
        invocation=invocation,
        project=project,
    )
    failures.extend(state_rules)
    effective_charge = effective.get("charge")
    effective_multiplicity = effective.get("multiplicity")
    effective_sha256 = canonical_sha256(effective)

    if safe_preview is not None:
        if safe_preview.invocation_sha256 != invocation.invocation_sha256:
            failures.append("preflight.preview.stale")
        if safe_preview.input_sha256 != invocation.input_sha256:
            failures.append("preflight.preview.input_mismatch")
        if safe_preview.project_sha256 != invocation.project_sha256:
            failures.append("preflight.preview.project_mismatch")
        if safe_preview.status != "previewed":
            failures.append("preflight.safe_preview.failed")

    if failures:
        state = "blocked"
        execution_ready = False
        rules = tuple(sorted(set(failures)))
    elif safe_preview is None:
        state = "ready_for_safe_preview"
        execution_ready = False
        rules = ("preflight.safe_preview.required",)
    else:
        state = "previewed"
        execution_ready = bool(
            capability.status is CapabilityQueryStatus.SUPPORTED
            and engine_binding.execution_ready
        )
        rules = (
            "preflight.execution_environment.observed"
            if execution_ready
            else "preflight.preview_only.execution_not_authorized",
        )
    body = {
        "schema_version": "chemsmart.program-node-preflight-receipt.v1",
        "request_sha256": request.request_sha256,
        "program": capability.query.program,
        "jobtype": capability.query.jobtype,
        "safe_preview_receipt_sha256": (
            safe_preview.receipt_sha256 if safe_preview is not None else ""
        ),
        "effective_charge": effective_charge,
        "effective_multiplicity": effective_multiplicity,
        "effective_settings_sha256": effective_sha256,
        "validator_receipt_sha256s": request.validator_receipt_sha256s,
        "critical_finding_sha256s": critical_finding_sha256s,
        "plan_state": state,
        "execution_ready": execution_ready,
        "rule_ids": rules,
    }
    return ProgramNodePreflightReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def _effective_electronic_state(
    *,
    request: ProgramNodePreflightRequestV1,
    invocation: CanonicalCommandInvocationV1,
    project: ProjectValidationReceiptV1 | None,
) -> tuple[dict[str, int | None], list[str]]:
    options = {
        item.parameter_name: item.values[0]
        for item in invocation.scoped_options
        if item.values
    }
    project_values = dict(project.settings) if project is not None else {}
    effective: dict[str, int | None] = {}
    failures = []
    for field, expected in (
        ("charge", request.charge),
        ("multiplicity", request.multiplicity),
    ):
        raw: Any = options.get(field, project_values.get(field))
        try:
            observed = int(raw) if raw is not None else None
        except (TypeError, ValueError):
            observed = None
        effective[field] = observed
        if observed is None:
            failures.append(f"preflight.effective_{field}.missing")
        elif observed != expected:
            failures.append(f"preflight.effective_{field}.mismatch")
    return effective, failures


def _request_body(request: ProgramNodePreflightRequestV1) -> dict[str, Any]:
    return {
        "schema_version": request.schema_version,
        "node_id": request.node_id,
        "capability_receipt_sha256": request.capability_receipt_sha256,
        "program_binding_sha256": request.program_binding_sha256,
        "engine_binding_sha256": request.engine_binding_sha256,
        "environment_receipt_sha256": request.environment_receipt_sha256,
        "geometry_artifact_sha256": request.geometry_artifact_sha256,
        "scientific_identity_sha256": request.scientific_identity_sha256,
        "charge": request.charge,
        "multiplicity": request.multiplicity,
        "project_receipt_sha256": request.project_receipt_sha256,
        "invocation_sha256": request.invocation_sha256,
        "command_inspection_sha256": request.command_inspection_sha256,
        "validator_receipt_sha256s": request.validator_receipt_sha256s,
    }


def _require_digest_set(values: tuple[str, ...], field_name: str) -> None:
    if values != tuple(sorted(set(values))):
        raise ContractError(f"{field_name} digests must be sorted and unique")
    for digest in values:
        require_sha256(digest, field_name)


__all__ = [
    "ProgramNodePreflightReceiptV1",
    "ProgramNodePreflightRequestV1",
    "ProgramValidatorReceiptV1",
    "build_program_node_preflight_request",
    "evaluate_program_node_preflight",
    "validator_receipt_from_safe_preview",
]
