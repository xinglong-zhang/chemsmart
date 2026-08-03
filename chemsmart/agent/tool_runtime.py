"""Host-owned dispatcher for the command-compiled model tool surface."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
from typing import Any, Mapping

from chemsmart.agent._contracts import (
    ContractError,
    TrustedArtifactRefV1,
    canonical_data,
    canonical_sha256,
    file_sha256,
)
from chemsmart.agent.capabilities import (
    CapabilityQueryReceiptV1,
    CapabilityQueryV1,
    EnvironmentCapabilityReceiptV1,
    EnvironmentTargetV1,
    ProgramComponentConformanceReceiptV1,
    ProgramSupportOverlayV1,
    ProgramCapabilityRegistryV1,
    ResolvedEngineBindingV1,
    ResolvedProgramBindingV1,
    TrustedComputeEnvironmentReceiptV1,
    build_approved_execution_overlay,
    build_command_compiled_preview_overlay,
    load_program_capabilities,
    query_capability,
    query_environment,
    resolve_engine_binding,
    resolve_program_binding,
)
from chemsmart.agent.cli_schema import LiveClickSchemaV1, build_live_click_schema
from chemsmart.agent.commands import (
    CanonicalCommandInvocationV1,
    CommandCounterexampleV1,
    CommandInspectionReceiptV1,
    CommandProposalV1,
    ScientificIdentityBindingV1,
    build_scientific_identity_binding,
    compile_command,
    inspect_command,
)
from chemsmart.agent.inspection import (
    GeneratedArtifactInspectionReceiptV1,
    inspect_generated_artifact,
)
from chemsmart.agent.execution import (
    ExecutionResourceSpecV1,
    OptimizedGeometryHandoffV1,
    ProgramExecutionReceiptV1,
    ProjectArtifactPromotionV1,
    ScientificDecisionRecordV1,
    WorkflowExecutionApprovalV1,
    bind_project_promotion_validation,
    build_program_execution_invocation,
    build_program_execution_receipt,
    build_scientific_decision_record,
    handoff_optimized_pyscf_geometry,
    promote_project_candidate,
)
from chemsmart.agent.knowledge import (
    FunctionalEquivalenceReceiptV1,
    ProgramSubstitutionReceiptV1,
    ProgramSubstitutionApprovalV1,
    ScientificClaimEvidenceV1,
    assess_typed_program_substitution,
    build_program_substitution_request,
)
from chemsmart.agent.preflight import (
    ProgramNodePreflightReceiptV1,
    ProgramValidatorReceiptV1,
    build_program_node_preflight_request,
    evaluate_program_node_preflight,
    validator_receipt_from_safe_preview,
)
from chemsmart.agent.preview import SafePreviewReceiptV1, execute_safe_preview
from chemsmart.agent.program_verifiers import build_preview_expectation
from chemsmart.agent.projects import (
    ProjectDocumentV1,
    ProjectRenderReceiptV1,
    ProjectValidationReceiptV1,
    project_document,
    read_project_yaml,
    render_project_yaml,
    validate_project_yaml,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.runtime.events import EventKind
from chemsmart.agent.tool_specs import (
    AgentToolSurfaceV1,
    build_approved_execution_tool_surface,
    build_command_compiled_tool_surface,
    build_single_agent_baseline_tool_surface,
)
from chemsmart.agent.workflows import (
    ArtifactInputIntentV1,
    ArtifactOutputIntentV1,
    CommandNodeIntentV1,
    CommandWorkflowDraftV1,
    build_command_workflow_draft,
)


@dataclass(frozen=True)
class _CommandContext:
    proposal: CommandProposalV1
    capability: CapabilityQueryReceiptV1
    program_binding: ResolvedProgramBindingV1
    engine_binding: ResolvedEngineBindingV1
    project_artifact: TrustedArtifactRefV1 | None
    project_validation: ProjectValidationReceiptV1 | None
    input_artifact: TrustedArtifactRefV1
    scientific_identity: ScientificIdentityBindingV1


class CommandCompiledToolHostV1:
    """Resolve every model ID against immutable host-held objects."""

    def __init__(
        self,
        *,
        event_store: RuntimeEventStore,
        artifacts: Mapping[str, TrustedArtifactRefV1] = {},
        scientific_identities: Mapping[str, ScientificIdentityBindingV1] = {},
        environment_targets: tuple[EnvironmentTargetV1, ...] = (),
        compute_environment_receipts: tuple[
            TrustedComputeEnvironmentReceiptV1, ...
        ] = (),
        component_conformance_receipts: tuple[
            ProgramComponentConformanceReceiptV1, ...
        ] = (),
        support_overlay: ProgramSupportOverlayV1 | None = None,
        tool_surface: AgentToolSurfaceV1 | None = None,
        settings_objects: Mapping[str, Any] = {},
        run_receipts: Mapping[str, Mapping[str, Any]] = {},
        scientific_claim_evidence: Mapping[
            str, ScientificClaimEvidenceV1
        ] = {},
        functional_equivalence_receipts: Mapping[
            str, FunctionalEquivalenceReceiptV1
        ] = {},
        substitution_approvals: Mapping[
            str, ProgramSubstitutionApprovalV1
        ] = {},
        capability_receipts: Mapping[str, CapabilityQueryReceiptV1] = {},
        environment_receipts: Mapping[
            str, EnvironmentCapabilityReceiptV1
        ] = {},
        program_binding_receipts: Mapping[
            str, ResolvedProgramBindingV1
        ] = {},
        engine_binding_receipts: Mapping[
            str, ResolvedEngineBindingV1
        ] = {},
        project_validation_receipts: Mapping[
            str, ProjectValidationReceiptV1
        ] = {},
        registry: ProgramCapabilityRegistryV1 | None = None,
        live_schema: LiveClickSchemaV1 | None = None,
        task_spec_sha256s: tuple[str, ...] = (),
        approved_workspace: str | Path | None = None,
        execution_resources: ExecutionResourceSpecV1 | None = None,
        workflow_execution_approval: WorkflowExecutionApprovalV1 | None = None,
        execution_server: str = "",
        execution_environment: Mapping[str, str] = {},
    ) -> None:
        self.event_store = event_store
        self.registry = registry or load_program_capabilities()
        self.live_schema = live_schema or build_live_click_schema()
        preview_overlay = build_command_compiled_preview_overlay(
            self.registry,
            conformance_receipts=component_conformance_receipts,
            live_schema=self.live_schema,
        )
        h1_surface = build_command_compiled_tool_surface(self.registry)
        h0_surface = build_single_agent_baseline_tool_surface(self.registry)
        execution_surface = build_approved_execution_tool_surface(self.registry)
        if tool_surface is not None and tool_surface.tool_schema_sha256 not in {
            h0_surface.tool_schema_sha256,
            h1_surface.tool_schema_sha256,
            execution_surface.tool_schema_sha256,
        }:
            raise ContractError("injected tool surface is not a canonical profile")
        self.surface = tool_surface or h1_surface
        self.artifacts = dict(artifacts)
        self.task_spec_sha256s = frozenset(task_spec_sha256s)
        self.approved_workspace = (
            Path(approved_workspace).resolve()
            if approved_workspace is not None
            else None
        )
        self.execution_resources = execution_resources
        self.workflow_execution_approval = workflow_execution_approval
        self.execution_server = str(execution_server)
        self.execution_environment = {
            str(key): str(value) for key, value in execution_environment.items()
        }
        if self.surface.profile == "command_compiled_approved_execution":
            if self.approved_workspace is None:
                raise ContractError("execution profile requires an approved workspace")
            if self.execution_resources is None:
                raise ContractError("execution profile requires host-owned resources")
            if self.workflow_execution_approval is None:
                raise ContractError("execution profile requires workflow approval")
            if (
                Path(self.workflow_execution_approval.workspace).resolve()
                != self.approved_workspace
            ):
                raise ContractError("workflow approval targets another workspace")
            execution_evidence_sha256 = canonical_sha256(
                {
                    "approval_sha256": (
                        self.workflow_execution_approval.approval_sha256
                    ),
                    "resource_sha256": self.execution_resources.resource_sha256,
                    "compute_environment_receipts": tuple(
                        sorted(
                            item.evidence_sha256
                            for item in compute_environment_receipts
                        )
                    ),
                    "execution_server": self.execution_server,
                    "execution_environment": self.execution_environment,
                }
            )
            evidence_overlay = build_approved_execution_overlay(
                registry=self.registry,
                preview_overlay=preview_overlay,
                approved_nodes=(
                    (item.program, item.jobtype, item.engine)
                    for item in self.workflow_execution_approval.node_bindings
                ),
                execution_evidence_sha256=execution_evidence_sha256,
            )
        else:
            evidence_overlay = preview_overlay
        if support_overlay is not None:
            if support_overlay.base_registry_sha256 != self.registry.registry_sha256:
                raise ContractError("injected support overlay uses another registry")
            if support_overlay.overlay_sha256 != evidence_overlay.overlay_sha256:
                raise ContractError(
                    "injected support overlay lacks matching host evidence"
                )
        self.overlay = support_overlay or evidence_overlay
        self.scientific_identities = dict(scientific_identities)
        self.environment_targets = tuple(environment_targets)
        self.compute_environment_receipts = tuple(
            compute_environment_receipts
        )
        self.settings_objects = dict(settings_objects)
        self.run_receipts = {
            key: dict(value) for key, value in run_receipts.items()
        }
        self.scientific_claim_evidence = dict(scientific_claim_evidence)
        self.functional_equivalence_receipts = dict(
            functional_equivalence_receipts
        )
        self.substitution_approvals = dict(substitution_approvals)
        if any(
            key != value.claim_sha256
            for key, value in self.scientific_claim_evidence.items()
        ):
            raise ContractError("claim evidence registry key mismatch")
        if any(
            key != value.receipt_sha256
            for key, value in self.functional_equivalence_receipts.items()
        ):
            raise ContractError("equivalence receipt registry key mismatch")
        if any(
            key != value.substitution_request_sha256
            for key, value in self.substitution_approvals.items()
        ):
            raise ContractError("substitution approval registry key mismatch")
        self.capabilities: dict[str, CapabilityQueryReceiptV1] = dict(
            capability_receipts
        )
        self.environments: dict[str, EnvironmentCapabilityReceiptV1] = dict(
            environment_receipts
        )
        self.program_bindings: dict[str, ResolvedProgramBindingV1] = dict(
            program_binding_receipts
        )
        self.engine_bindings: dict[str, ResolvedEngineBindingV1] = dict(
            engine_binding_receipts
        )
        self.substitutions: dict[str, ProgramSubstitutionReceiptV1] = {}
        self.project_documents: dict[str, ProjectDocumentV1] = {}
        self.project_renders: dict[str, ProjectRenderReceiptV1] = {}
        self.project_validations: dict[str, ProjectValidationReceiptV1] = dict(
            project_validation_receipts
        )
        self.invocations: dict[str, CanonicalCommandInvocationV1] = {}
        self.command_inspections: dict[str, CommandInspectionReceiptV1] = {}
        self.safe_previews: dict[str, SafePreviewReceiptV1] = {}
        self.validators: dict[str, ProgramValidatorReceiptV1] = {}
        self.preflights: dict[str, ProgramNodePreflightReceiptV1] = {}
        self.result_inspections: dict[
            str, GeneratedArtifactInspectionReceiptV1
        ] = {}
        self.counterexamples: dict[str, CommandCounterexampleV1] = {}
        self.workflow_drafts: dict[str, CommandWorkflowDraftV1] = {}
        self.project_promotions: dict[str, ProjectArtifactPromotionV1] = {}
        self.scientific_decisions: dict[str, ScientificDecisionRecordV1] = {}
        self.execution_receipts: dict[str, ProgramExecutionReceiptV1] = {}
        self.handoffs: dict[str, OptimizedGeometryHandoffV1] = {}
        self._command_contexts: dict[str, _CommandContext] = {}
        self._completion_sets: dict[str, tuple[str, ...]] = {}
        self._latest_environment_by_capability: dict[
            str, EnvironmentCapabilityReceiptV1
        ] = {}
        _require_registry_keys(
            self.capabilities, "receipt_sha256", "capability receipt"
        )
        _require_registry_keys(
            self.environments, "receipt_sha256", "environment receipt"
        )
        _require_registry_keys(
            self.program_bindings, "binding_sha256", "program binding"
        )
        _require_registry_keys(
            self.engine_bindings, "binding_sha256", "engine binding"
        )
        _require_registry_keys(
            self.project_validations,
            "receipt_sha256",
            "project validation",
        )
        for receipt in self.capabilities.values():
            if (
                receipt.registry_sha256 != self.registry.registry_sha256
                or receipt.live_cli_schema_sha256 != self.live_schema.schema_sha256
                or receipt.overlay_sha256 != self.overlay.overlay_sha256
            ):
                raise ContractError("seeded capability receipt is stale")
        for environment in self.environments.values():
            self._latest_environment_by_capability[
                environment.capability_receipt_sha256
            ] = environment

    def record_seeded_evidence(self, turn_id: str) -> None:
        """Persist host-prebound H0 evidence before any model action."""

        for receipt in self.capabilities.values():
            self._emit(
                turn_id,
                EventKind.CAPABILITY_QUERIED,
                receipt.receipt_sha256,
                status=receipt.status.value,
                program=receipt.query.program,
                jobtype=receipt.query.jobtype,
                engine=receipt.query.engine,
            )
        for receipt in self.environments.values():
            self._emit(
                turn_id,
                EventKind.ENVIRONMENT_QUERIED,
                receipt.receipt_sha256,
                status=receipt.status.value,
                program=receipt.program,
                engine=receipt.engine,
            )
        for binding in self.program_bindings.values():
            self._emit_binding(turn_id, EventKind.PROGRAM_BOUND, binding)
        for binding in self.engine_bindings.values():
            self._emit_binding(turn_id, EventKind.ENGINE_BOUND, binding)
        for receipt in self.project_validations.values():
            self._emit(
                turn_id,
                EventKind.PROJECT_VALIDATED,
                receipt.receipt_sha256,
                status=receipt.status,
                program=receipt.program,
                jobtype=receipt.jobtype,
            )

    def register_counterexample(self, value: CommandCounterexampleV1) -> None:
        """Register a deterministic counterexample; never model-authored here."""

        self.counterexamples[value.counterexample_id] = value

    def dispatch(
        self, *, turn_id: str, tool_name: str, arguments: Mapping[str, Any]
    ) -> dict[str, Any]:
        """Validate a call and invoke exactly one approved host operation."""

        values = dict(arguments)
        _validate_tool_arguments(self.surface, tool_name, values)
        handlers = {
            "inspect_program_capability": self._inspect_program_capability,
            "inspect_program_environment": self._inspect_program_environment,
            "assess_program_candidate": self._assess_program_candidate,
            "render_project_yaml": self._render_project_yaml,
            "promote_project_yaml": self._promote_project_yaml,
            "bind_scientific_identity": self._bind_scientific_identity,
            "read_project_yaml": self._read_project_yaml,
            "validate_project_yaml": self._validate_project_yaml,
            "plan_command_workflow": self._plan_command_workflow,
            "synthesize_command": self._synthesize_command,
            "repair_command": self._repair_command,
            "preview_command": self._preview_command,
            "preflight_program_node": self._preflight_program_node,
            "inspect_calculation_artifact": self._inspect_calculation_artifact,
            "record_scientific_decision": self._record_scientific_decision,
            "execute_approved_program_node": self._execute_approved_program_node,
        }
        handler = handlers.get(tool_name)
        if handler is None:
            raise ContractError("tool is absent from command-compiled profile")
        result = handler(turn_id, values)
        return {
            "schema_version": "chemsmart.tool-result.v1",
            "tool": tool_name,
            "status": "ok",
            "result": _model_visible_data(canonical_data(result)),
        }

    def _bind_scientific_identity(self, turn_id: str, values: dict) -> Any:
        task_spec_sha256 = values["task_spec_sha256"]
        if task_spec_sha256 not in self.task_spec_sha256s:
            raise ContractError("scientific identity targets an unknown task spec")
        artifact = self._artifact(values["input_artifact_id"])
        binding = build_scientific_identity_binding(
            task_spec_sha256=task_spec_sha256,
            geometry_artifact=artifact,
            charge=values["charge"],
            multiplicity=values["multiplicity"],
        )
        self.scientific_identities[binding.binding_sha256] = binding
        return binding

    def _inspect_program_capability(self, turn_id: str, values: dict) -> Any:
        receipt = query_capability(
            CapabilityQueryV1(
                program=values["program"],
                jobtype=values["jobtype"],
                engine=values["engine"],
            ),
            registry=self.registry,
            overlay=self.overlay,
            live_schema=self.live_schema,
        )
        self.capabilities[receipt.receipt_sha256] = receipt
        self._emit(
            turn_id,
            EventKind.CAPABILITY_QUERIED,
            receipt.receipt_sha256,
            status=receipt.status.value,
            program=receipt.query.program,
            jobtype=receipt.query.jobtype,
            engine=receipt.query.engine,
        )
        return receipt

    def _inspect_program_environment(self, turn_id: str, values: dict) -> Any:
        capability = self._get(
            self.capabilities,
            values["capability_receipt_sha256"],
            "capability receipt",
        )
        environment = query_environment(
            capability,
            targets=self.environment_targets,
            compute_receipts=self.compute_environment_receipts,
        )
        program_binding = resolve_program_binding(capability)
        engine_binding = resolve_engine_binding(program_binding, environment)
        self.environments[environment.receipt_sha256] = environment
        self.program_bindings[program_binding.binding_sha256] = program_binding
        self.engine_bindings[engine_binding.binding_sha256] = engine_binding
        self._latest_environment_by_capability[
            capability.receipt_sha256
        ] = environment
        self._emit(
            turn_id,
            EventKind.ENVIRONMENT_QUERIED,
            environment.receipt_sha256,
            status=environment.status.value,
            program=environment.program,
            engine=environment.engine,
        )
        self._emit_binding(turn_id, EventKind.PROGRAM_BOUND, program_binding)
        self._emit_binding(turn_id, EventKind.ENGINE_BOUND, engine_binding)
        return {
            "environment": environment,
            "program_binding": program_binding,
            "engine_binding": engine_binding,
        }

    def _assess_program_candidate(self, turn_id: str, values: dict) -> Any:
        capability = self._get(
            self.capabilities,
            values["capability_receipt_sha256"],
            "capability receipt",
        )
        claim_evidence = tuple(
            self._get(
                self.scientific_claim_evidence,
                digest,
                "scientific claim evidence",
            )
            for digest in values["source_claim_sha256s"]
        )
        verified_claims = tuple(
            item for item in claim_evidence if item.status == "verified"
        )
        identity = values["requested_program"] == values["selected_program"]
        equivalence_digest = values.get(
            "functional_equivalence_receipt_sha256", ""
        )
        equivalence = (
            self._get(
                self.functional_equivalence_receipts,
                equivalence_digest,
                "functional equivalence receipt",
            )
            if equivalence_digest
            else None
        )
        if not identity and equivalence is not None:
            if (
                equivalence.requested_program != values["requested_program"]
                or equivalence.selected_program != values["selected_program"]
                or equivalence.method_name != values["method_name"]
            ):
                raise ContractError("functional equivalence targets another request")
            observed_claim_receipts = tuple(
                sorted(item.receipt_sha256 for item in verified_claims)
            )
            if (
                equivalence.status == "verified"
                and equivalence.claim_evidence_receipt_sha256s
                != observed_claim_receipts
            ):
                raise ContractError("functional equivalence uses other claim evidence")
        request = build_program_substitution_request(
            request_id=values["request_id"],
            requested_program=values["requested_program"],
            selected_program=values["selected_program"],
            requested_engine=values["requested_engine"],
            selected_engine=values["selected_engine"],
            job_families=values["job_families"],
            method_family=values["method_family"],
            method_name=values["method_name"],
            basis_mode=values["basis_mode"],
            constraint_kinds=values["constraint_kinds"],
            requires_post_hf=values["requires_post_hf"],
            requires_double_hybrid=values["requires_double_hybrid"],
            functional_semantics_confirmed=bool(
                identity or (equivalence and equivalence.status == "verified")
            ),
            source_claim_sha256s=tuple(
                item.claim_sha256 for item in verified_claims
            ),
        )
        approval = self.substitution_approvals.get(request.request_sha256)
        approval_ref = ""
        if approval is not None:
            if (
                approval.substitution_request_sha256 != request.request_sha256
                or approval.decision != "approved"
            ):
                raise ContractError("substitution approval is stale or red")
            approval_ref = approval.receipt_sha256
        receipt = assess_typed_program_substitution(
            request, capability, approval_ref=approval_ref
        )
        self.substitutions[receipt.receipt_sha256] = receipt
        self._emit(
            turn_id,
            EventKind.SUBSTITUTION_ASSESSED,
            receipt.receipt_sha256,
            decision=receipt.decision,
            requested_program=receipt.requested_program,
            selected_program=receipt.selected_program,
        )
        if receipt.decision in {"exact", "approved"}:
            program_binding = resolve_program_binding(
                capability,
                requested_program=request.requested_program,
                substitution_receipt_sha256=receipt.receipt_sha256,
            )
            self.program_bindings[program_binding.binding_sha256] = program_binding
            self._emit_binding(
                turn_id, EventKind.PROGRAM_BOUND, program_binding
            )
            environment = self._latest_environment_by_capability.get(
                capability.receipt_sha256
            )
            if environment is not None:
                engine_binding = resolve_engine_binding(
                    program_binding, environment
                )
                self.engine_bindings[
                    engine_binding.binding_sha256
                ] = engine_binding
                self._emit_binding(
                    turn_id, EventKind.ENGINE_BOUND, engine_binding
                )
        return receipt

    def _render_project_yaml(self, turn_id: str, values: dict) -> Any:
        document = project_document(
            program=values["program"], sections=values["sections"]
        )
        receipt = render_project_yaml(document, registry=self.registry)
        self.project_documents[document.document_sha256] = document
        self.project_renders[receipt.receipt_sha256] = receipt
        return receipt

    def _promote_project_yaml(self, turn_id: str, values: dict) -> Any:
        if self.approved_workspace is None:
            raise ContractError("project promotion requires a task workspace")
        if values["artifact_id"] in self.artifacts:
            raise ContractError(
                "artifact ID is already registered; use a distinct project ID"
            )
        render = self._get(
            self.project_renders,
            values["render_receipt_sha256"],
            "project render receipt",
        )
        artifact, promotion = promote_project_candidate(
            render,
            approved_workspace=self.approved_workspace,
            artifact_id=values["artifact_id"],
        )
        self.artifacts[artifact.artifact_id] = artifact
        self.project_promotions[artifact.artifact_id] = promotion
        self._emit(
            turn_id,
            EventKind.PROJECT_PROMOTED,
            promotion.receipt_sha256,
            status=promotion.status,
            artifact_id=artifact.artifact_id,
        )
        return {"artifact": artifact, "promotion": promotion}

    def _read_project_yaml(self, turn_id: str, values: dict) -> Any:
        artifact = self._artifact(values["project_artifact_id"])
        document = read_project_yaml(artifact, program=values["program"])
        self.project_documents[document.document_sha256] = document
        return document

    def _validate_project_yaml(self, turn_id: str, values: dict) -> Any:
        artifact = self._artifact(values["project_artifact_id"])
        capability = self._get(
            self.capabilities,
            values["capability_receipt_sha256"],
            "capability receipt",
        )
        receipt = validate_project_yaml(artifact, capability=capability)
        self.project_validations[receipt.receipt_sha256] = receipt
        promotion = self.project_promotions.get(artifact.artifact_id)
        if promotion is not None and promotion.validation_status == "pending":
            self.project_promotions[artifact.artifact_id] = (
                bind_project_promotion_validation(
                    promotion, artifact, receipt
                )
            )
            bound = self.project_promotions[artifact.artifact_id]
            self._emit(
                turn_id,
                EventKind.PROJECT_PROMOTED,
                bound.receipt_sha256,
                status=bound.status,
                artifact_id=artifact.artifact_id,
            )
        self._emit(
            turn_id,
            EventKind.PROJECT_VALIDATED,
            receipt.receipt_sha256,
            status=receipt.status,
            program=receipt.program,
            jobtype=receipt.jobtype,
        )
        return receipt

    def _record_scientific_decision(self, turn_id: str, values: dict) -> Any:
        task_spec_sha256 = values["task_spec_sha256"]
        if task_spec_sha256 not in self.task_spec_sha256s:
            raise ContractError("scientific decision targets an unknown task spec")
        record = build_scientific_decision_record(
            decision_id=values["decision_id"],
            task_spec_sha256=task_spec_sha256,
            stage_order=values["stage_order"],
            assumptions=values["assumptions"],
            method_rationale=values["method_rationale"],
            alternatives=values["alternatives"],
            uncertainties=values["uncertainties"],
            diagnostics=values["diagnostics"],
            evidence_refs=values["evidence_refs"],
        )
        self.scientific_decisions[record.record_sha256] = record
        self._emit(
            turn_id,
            EventKind.SCIENTIFIC_DECISION_RECORDED,
            record.record_sha256,
            status="recorded",
            decision_id=record.decision_id,
        )
        return record

    def _plan_command_workflow(self, turn_id: str, values: dict) -> Any:
        """Record a broad DAG before execution-grade evidence is available."""

        nodes = []
        findings: list[dict[str, str]] = []
        declared_programs = {item.program: item for item in self.registry.programs}
        for raw_node in values["nodes"]:
            inputs = tuple(
                sorted(
                    (
                        ArtifactInputIntentV1(
                            binding_id=item["binding_id"],
                            artifact_class=item["artifact_class"],
                            artifact_id=item.get("artifact_id", ""),
                            producer_node_id=item["producer_node_id"],
                            producer_output_id=item["producer_output_id"],
                        )
                        for item in raw_node["inputs"]
                    ),
                    key=lambda item: item.binding_id,
                )
            )
            outputs = tuple(
                sorted(
                    (
                        ArtifactOutputIntentV1(
                            output_id=item["output_id"],
                            artifact_class=item["artifact_class"],
                        )
                        for item in raw_node["expected_outputs"]
                    ),
                    key=lambda item: item.output_id,
                )
            )
            node = CommandNodeIntentV1(
                node_id=raw_node["node_id"],
                program=raw_node["program"],
                jobtype=raw_node["jobtype"],
                project_role=raw_node["project_role"],
                dependencies=tuple(sorted(set(raw_node["dependencies"]))),
                inputs=inputs,
                expected_outputs=outputs,
                unresolved_fields=tuple(
                    sorted(set(raw_node["unresolved_fields"]))
                ),
            )
            capability = declared_programs.get(node.program)
            if capability is None:
                findings.append(
                    {"node_id": node.node_id, "rule_id": "workflow.program.not_declared"}
                )
            elif node.jobtype not in capability.jobtypes:
                findings.append(
                    {"node_id": node.node_id, "rule_id": "workflow.job.not_declared"}
                )
            for item in node.inputs:
                if (
                    not item.producer_node_id
                    and (
                        not item.artifact_id
                        or item.artifact_id not in self.artifacts
                    )
                ):
                    findings.append(
                        {
                            "node_id": node.node_id,
                            "rule_id": "workflow.input.unresolved",
                        }
                    )
            nodes.append(node)
        draft = build_command_workflow_draft(
            workflow_id=values["workflow_id"],
            task_spec_id=values["task_spec_id"],
            nodes=tuple(nodes),
        )
        self.workflow_drafts[draft.draft_sha256] = draft
        unresolved_nodes = {
            item["node_id"] for item in findings
        } | {
            node.node_id
            for node in draft.nodes
            if node.dependencies or node.unresolved_fields
        }
        actionable = tuple(
            node.node_id
            for node in draft.nodes
            if node.node_id not in unresolved_nodes
        )
        unresolved = tuple(
            node.node_id for node in draft.nodes if node.node_id in unresolved_nodes
        )
        self._emit(
            turn_id,
            EventKind.WORKFLOW_PLANNED,
            draft.draft_sha256,
            status="planned",
            actionable_node_ids=actionable,
            unresolved_node_ids=unresolved,
        )
        return {
            "workflow_draft": draft,
            "actionable_node_ids": actionable,
            "unresolved_node_ids": unresolved,
            "findings": tuple(findings),
        }

    def _synthesize_command(self, turn_id: str, values: dict) -> Any:
        capability = self._get(
            self.capabilities,
            values["capability_receipt_sha256"],
            "capability receipt",
        )
        binding = self._get(
            self.engine_bindings,
            values["engine_binding_sha256"],
            "engine binding",
        )
        input_artifact = self._artifact(values["input_artifact_id"])
        identity = self._get(
            self.scientific_identities,
            values["scientific_identity_sha256"],
            "scientific identity",
        )
        project = (
            self._artifact(values["project_artifact_id"])
            if values["project_artifact_id"]
            else None
        )
        validation_digest = values.get(
            "project_validation_receipt_sha256", ""
        )
        validation = (
            self._get(
                self.project_validations,
                validation_digest,
                "project validation receipt",
            )
            if validation_digest
            else None
        )
        proposal = CommandProposalV1(
            node_id=values["node_id"],
            execution_target=values["execution_target"],
            program=values["program"],
            jobtype=values["jobtype"],
            project_artifact_id=values["project_artifact_id"],
            input_artifact_id=values["input_artifact_id"],
            scientific_identity_sha256=values[
                "scientific_identity_sha256"
            ],
            charge=values["charge"],
            multiplicity=values["multiplicity"],
        )
        invocation = compile_command(
            proposal,
            capability=capability,
            binding=binding,
            project=project,
            project_validation=validation,
            input_artifact=input_artifact,
            scientific_identity=identity,
            live_schema=self.live_schema,
        )
        program_binding = self._get(
            self.program_bindings,
            binding.program_binding_sha256,
            "program binding",
        )
        context = _CommandContext(
            proposal=proposal,
            capability=capability,
            program_binding=program_binding,
            engine_binding=binding,
            project_artifact=project,
            project_validation=validation,
            input_artifact=input_artifact,
            scientific_identity=identity,
        )
        return self._record_compiled_command(turn_id, invocation, context)

    def _repair_command(self, turn_id: str, values: dict) -> Any:
        parent = self._get(
            self.invocations,
            values["invocation_sha256"],
            "canonical invocation",
        )
        context = self._get(
            self._command_contexts,
            parent.invocation_sha256,
            "command context",
        )
        counterexample = self._get(
            self.counterexamples,
            values["counterexample_id"],
            "counterexample",
        )
        if counterexample.invocation_sha256 != parent.invocation_sha256:
            raise ContractError("counterexample targets another invocation")
        if (
            counterexample.task_spec_sha256
            != context.scientific_identity.task_spec_sha256
        ):
            raise ContractError("counterexample targets another task spec")
        invocation = compile_command(
            context.proposal,
            capability=context.capability,
            binding=context.engine_binding,
            project=context.project_artifact,
            project_validation=context.project_validation,
            input_artifact=context.input_artifact,
            scientific_identity=context.scientific_identity,
            live_schema=self.live_schema,
            repair_parent_sha256=parent.invocation_sha256,
            counterexample_sha256=counterexample.counterexample_sha256,
            repair_attempt=parent.repair_attempt + 1,
        )
        return self._record_compiled_command(turn_id, invocation, context)

    def _preview_command(self, turn_id: str, values: dict) -> Any:
        invocation = self._get(
            self.invocations,
            values["invocation_sha256"],
            "canonical invocation",
        )
        context = self._get(
            self._command_contexts,
            invocation.invocation_sha256,
            "command context",
        )
        expectation = build_preview_expectation(
            program=context.proposal.program,
            jobtype=context.proposal.jobtype,
            input_artifact=context.input_artifact,
            project=context.project_validation,
            charge=context.scientific_identity.charge,
            multiplicity=context.scientific_identity.multiplicity,
        )
        receipt = execute_safe_preview(
            invocation,
            input_artifact=context.input_artifact,
            project_artifact=context.project_artifact,
            expectation=expectation,
        )
        validator = validator_receipt_from_safe_preview(
            node_id=context.proposal.node_id,
            program=context.proposal.program,
            scientific_identity_sha256=(
                context.scientific_identity.binding_sha256
            ),
            safe_preview=receipt,
        )
        self.safe_previews[receipt.receipt_sha256] = receipt
        self.validators[validator.receipt_sha256] = validator
        self._emit(
            turn_id,
            EventKind.SAFE_PREVIEWED,
            receipt.receipt_sha256,
            status=receipt.status,
            program_validation_status=receipt.program_validation_status,
            critical_finding_count=len(receipt.critical_finding_sha256s),
            invocation_sha256=receipt.invocation_sha256,
        )
        self._emit(
            turn_id,
            EventKind.VALIDATOR_OBSERVED,
            validator.receipt_sha256,
            status=validator.status,
            critical_finding_count=len(
                validator.critical_finding_sha256s
            ),
            source_receipt_sha256=validator.source_receipt_sha256,
        )
        return {"safe_preview": receipt, "validator": validator}

    def _preflight_program_node(self, turn_id: str, values: dict) -> Any:
        capability = self._get(
            self.capabilities,
            values["capability_receipt_sha256"],
            "capability receipt",
        )
        program_binding = self._get(
            self.program_bindings,
            values["program_binding_sha256"],
            "program binding",
        )
        engine_binding = self._get(
            self.engine_bindings,
            values["engine_binding_sha256"],
            "engine binding",
        )
        invocation = self._get(
            self.invocations,
            values["invocation_sha256"],
            "canonical invocation",
        )
        command_inspection = self._get(
            self.command_inspections,
            values["command_inspection_receipt_sha256"],
            "command inspection receipt",
        )
        identity = self._get(
            self.scientific_identities,
            values["scientific_identity_sha256"],
            "scientific identity",
        )
        project_digest = values.get(
            "project_validation_receipt_sha256", ""
        )
        project = (
            self._get(
                self.project_validations,
                project_digest,
                "project validation receipt",
            )
            if project_digest
            else None
        )
        validators = tuple(
            self._get(self.validators, digest, "program validator receipt")
            for digest in values["validator_receipt_sha256s"]
        )
        preview_digest = values.get("safe_preview_receipt_sha256", "")
        safe_preview = (
            self._get(
                self.safe_previews, preview_digest, "safe preview receipt"
            )
            if preview_digest
            else None
        )
        request = build_program_node_preflight_request(
            node_id=values["node_id"],
            capability_receipt_sha256=capability.receipt_sha256,
            program_binding_sha256=program_binding.binding_sha256,
            engine_binding_sha256=engine_binding.binding_sha256,
            environment_receipt_sha256=(
                engine_binding.environment_receipt_sha256
            ),
            geometry_artifact_sha256=values["geometry_artifact_sha256"],
            scientific_identity_sha256=identity.binding_sha256,
            charge=values["charge"],
            multiplicity=values["multiplicity"],
            project_receipt_sha256=project.receipt_sha256 if project else "",
            invocation_sha256=invocation.invocation_sha256,
            command_inspection_sha256=command_inspection.receipt_sha256,
            validator_receipts=validators,
        )
        receipt = evaluate_program_node_preflight(
            request,
            capability=capability,
            program_binding=program_binding,
            engine_binding=engine_binding,
            project=project,
            invocation=invocation,
            command_inspection=command_inspection,
            scientific_identity=identity,
            validator_receipts=validators,
            safe_preview=safe_preview,
        )
        self.preflights[receipt.receipt_sha256] = receipt
        completion = {
            capability.receipt_sha256,
            program_binding.binding_sha256,
            engine_binding.binding_sha256,
            invocation.invocation_sha256,
            command_inspection.receipt_sha256,
            receipt.receipt_sha256,
            *(item.receipt_sha256 for item in validators),
        }
        if project is not None:
            completion.add(project.receipt_sha256)
        if safe_preview is not None:
            completion.add(safe_preview.receipt_sha256)
        self._completion_sets[receipt.receipt_sha256] = tuple(
            sorted(completion)
        )
        self._emit(
            turn_id,
            EventKind.PROGRAM_PREFLIGHTED,
            receipt.receipt_sha256,
            plan_state=receipt.plan_state,
            critical_finding_count=len(receipt.critical_finding_sha256s),
            safe_preview_receipt_sha256=(
                receipt.safe_preview_receipt_sha256
            ),
            execution_ready=receipt.execution_ready,
        )
        return receipt

    def completion_receipts_for_latest_preflight(self) -> tuple[str, ...]:
        """Return host-derived gates only for a green previewed preflight."""

        if not self.preflights:
            raise ContractError("no node preflight has been observed")
        latest = next(reversed(self.preflights.values()))
        if latest.plan_state != "previewed" or latest.critical_finding_sha256s:
            raise ContractError("latest node preflight is not green")
        completion = list(self._completion_sets[latest.receipt_sha256])
        if self.surface.profile == "command_compiled_approved_execution":
            if not self.execution_receipts:
                raise ContractError("approved workflow has not executed a node")
            approval = self.workflow_execution_approval
            required_nodes = tuple(
                item.node_id for item in approval.node_bindings
            )
            missing = [
                node_id
                for node_id in required_nodes
                if node_id not in self.execution_receipts
                or not self.execution_receipts[node_id].validated
            ]
            if missing:
                raise ContractError(
                    "approved workflow execution remains incomplete: "
                    + ", ".join(missing)
                )
            completion.extend(
                self.execution_receipts[node_id].receipt_sha256
                for node_id in required_nodes
            )
        return tuple(sorted(set(completion)))

    def _execute_approved_program_node(self, turn_id: str, values: dict) -> Any:
        """Resolve and execute one approved node without model-authored argv."""

        if self.surface.profile != "command_compiled_approved_execution":
            raise ContractError("execution tool is absent from this profile")
        node_id = values["node_id"]
        existing = self.execution_receipts.get(node_id)
        if existing is not None:
            if existing.execution_state == "ambiguous":
                raise ContractError(
                    "prior launch is ambiguous; inspect local process state first"
                )
            return {
                "execution": existing,
                "idempotent_replay": True,
                "handoff": self.handoffs.get(node_id),
            }
        approval = self.workflow_execution_approval
        approved_node = approval.node(node_id)
        ordered_node_ids = tuple(item.node_id for item in approval.node_bindings)
        node_index = ordered_node_ids.index(node_id)
        incomplete_predecessors = [
            predecessor
            for predecessor in ordered_node_ids[:node_index]
            if predecessor not in self.execution_receipts
            or not self.execution_receipts[predecessor].validated
        ]
        if incomplete_predecessors:
            raise ContractError(
                "approved predecessors are not validated: "
                + ", ".join(incomplete_predecessors)
            )
        invocation, context = self._latest_invocation_for_node(node_id)
        if context.project_artifact is None or context.project_validation is None:
            raise ContractError("execution requires a validated project artifact")
        if context.project_validation.status != "valid":
            raise ContractError("execution project loader gate is red")
        if context.project_validation.settings_sha256 != approved_node.settings_sha256:
            raise ContractError("effective project settings differ from approval")
        if context.project_artifact.sha256 != approved_node.project_artifact_sha256:
            raise ContractError("project bytes differ from workflow approval")
        if not self._invocation_has_green_preflight(invocation.invocation_sha256):
            raise ContractError("node requires a green safe-preview preflight")

        handoff = self.handoffs.get(node_id)
        real_argv = self._real_execution_argv(invocation)
        execution_invocation = build_program_execution_invocation(
            node_id=node_id,
            approval=approval,
            project_artifact=context.project_artifact,
            input_artifact=context.input_artifact,
            scientific_identity_sha256=(
                context.scientific_identity.binding_sha256
            ),
            environment_receipt_sha256=(
                context.engine_binding.environment_receipt_sha256
            ),
            resources=self.execution_resources,
            argv=real_argv,
            handoff=handoff,
        )
        node_workspace = self.approved_workspace / "nodes" / node_id
        _prepare_execution_node_workspace(node_workspace)
        started = datetime.now(timezone.utc).isoformat()
        command = [sys.executable, "-m", "chemsmart", *real_argv[1:]]
        environment = os.environ.copy()
        environment.update(self.execution_environment)
        source_root = str(Path(__file__).resolve().parents[2])
        current_pythonpath = environment.get("PYTHONPATH", "")
        environment["PYTHONPATH"] = (
            source_root
            if not current_pythonpath
            else source_root + os.pathsep + current_pythonpath
        )
        try:
            completed = subprocess.run(
                command,
                cwd=node_workspace,
                env=environment,
                capture_output=True,
                text=True,
                timeout=self.execution_resources.node_timeout_seconds,
                check=False,
            )
            exit_status = int(completed.returncode)
            (node_workspace / "controller.stdout").write_text(
                completed.stdout, encoding="utf-8"
            )
            (node_workspace / "controller.stderr").write_text(
                completed.stderr, encoding="utf-8"
            )
            execution_state = (
                "engine_complete" if exit_status == 0 else "failed"
            )
        except subprocess.TimeoutExpired as exc:
            (node_workspace / "controller.stdout").write_text(
                str(exc.stdout or ""), encoding="utf-8"
            )
            (node_workspace / "controller.stderr").write_text(
                str(exc.stderr or ""), encoding="utf-8"
            )
            exit_status = None
            execution_state = "ambiguous"
        outputs = self._execution_output_artifacts(node_id, node_workspace)
        validated, validator_sha256s, findings = self._validate_execution_outputs(
            program=context.proposal.program,
            jobtype=context.proposal.jobtype,
            charge=context.scientific_identity.charge,
            multiplicity=context.scientific_identity.multiplicity,
            expected_settings=dict(context.project_validation.settings),
            expected_input_artifact=context.input_artifact,
            expected_project_artifact=context.project_artifact,
            output_artifacts=outputs,
            exit_status=exit_status,
        )
        finished = datetime.now(timezone.utc).isoformat()
        receipt = build_program_execution_receipt(
            execution_invocation,
            execution_state=execution_state,
            exit_status=exit_status,
            engine_complete=execution_state == "engine_complete",
            validated=validated,
            output_artifacts=outputs,
            validator_receipt_sha256s=validator_sha256s,
            findings=findings,
            started_at=started,
            finished_at=finished,
        )
        self.execution_receipts[node_id] = receipt
        self._emit(
            turn_id,
            EventKind.PROGRAM_EXECUTED,
            receipt.receipt_sha256,
            execution_state=receipt.execution_state,
            engine_complete=receipt.engine_complete,
            validated=receipt.validated,
            node_id=node_id,
        )
        produced_handoffs = []
        if receipt.validated and context.proposal.program == "pyscf":
            for edge in approval.producer_edges:
                if edge.producer_node_id != node_id:
                    continue
                result = next(
                    (item for item in outputs if item.kind == "pyscf_hdf5"),
                    None,
                )
                if result is None:
                    raise ContractError("validated PySCF OPT lacks an HDF5 result")
                geometry, observed = handoff_optimized_pyscf_geometry(
                    producer_receipt=receipt,
                    result_artifact=result,
                    producer_edge=edge,
                    approved_workspace=self.approved_workspace,
                    geometry_artifact_id=(
                        f"geometry.{edge.producer_node_id}-to-{edge.consumer_node_id}"
                    ),
                    expected_charge=context.scientific_identity.charge,
                    expected_multiplicity=context.scientific_identity.multiplicity,
                )
                identity = build_scientific_identity_binding(
                    task_spec_sha256=approval.task_spec_sha256,
                    geometry_artifact=geometry,
                    charge=observed.charge,
                    multiplicity=observed.multiplicity,
                )
                self.artifacts[geometry.artifact_id] = geometry
                self.scientific_identities[identity.binding_sha256] = identity
                self.handoffs[edge.consumer_node_id] = observed
                produced_handoffs.append(
                    {
                        "handoff": observed,
                        "artifact": geometry,
                        "scientific_identity": identity,
                    }
                )
                self._emit(
                    turn_id,
                    EventKind.OPTIMIZED_GEOMETRY_HANDED_OFF,
                    observed.receipt_sha256,
                    status=observed.status,
                    producer_node_id=edge.producer_node_id,
                    consumer_node_id=edge.consumer_node_id,
                )
        return {
            "execution": receipt,
            "idempotent_replay": False,
            "produced_handoffs": tuple(produced_handoffs),
        }

    def _latest_invocation_for_node(
        self, node_id: str
    ) -> tuple[CanonicalCommandInvocationV1, _CommandContext]:
        for invocation in reversed(tuple(self.invocations.values())):
            context = self._command_contexts[invocation.invocation_sha256]
            if context.proposal.node_id == node_id:
                return invocation, context
        raise ContractError("node has no compiled command invocation")

    def _invocation_has_green_preflight(self, invocation_sha256: str) -> bool:
        return any(
            invocation_sha256 in self._completion_sets.get(
                receipt.receipt_sha256, ()
            )
            and receipt.plan_state == "previewed"
            and receipt.execution_ready
            and not receipt.critical_finding_sha256s
            for receipt in self.preflights.values()
        )

    def _real_execution_argv(
        self, invocation: CanonicalCommandInvocationV1
    ) -> tuple[str, ...]:
        if invocation.command_path[0] != "run":
            raise ContractError("local execution requires a run command")
        program = invocation.command_path[1]
        try:
            program_index = invocation.argv.index(program, 2)
        except ValueError as exc:
            raise ContractError("compiled argv does not contain its program") from exc
        resources = self.execution_resources
        root = ["chemsmart", "run", "--no-fake"]
        root.append("--no-scratch" if resources.scratch_policy == "none" else "--scratch")
        root.extend(("--num-cores", str(resources.cores)))
        root.extend(("--num-gpus", str(resources.gpu_count)))
        memory = (
            str(int(resources.memory_gb))
            if resources.memory_gb.is_integer()
            else str(resources.memory_gb)
        )
        root.extend(("--mem-gb", memory))
        if self.execution_server:
            root.extend(("--server", self.execution_server))
        return tuple(root + list(invocation.argv[program_index:]))

    def _execution_output_artifacts(
        self, node_id: str, workspace: Path
    ) -> tuple[TrustedArtifactRefV1, ...]:
        artifacts = []
        ordinal = 0
        for path in sorted(workspace.rglob("*")):
            if path.is_symlink():
                raise ContractError("execution emitted a symbolic link")
            if not path.is_file() or path.name.startswith("controller."):
                continue
            ordinal += 1
            kind = (
                "pyscf_hdf5"
                if path.suffix.lower() == ".h5"
                else "geometry_xyz"
                if path.suffix.lower() == ".xyz"
                else "json"
                if path.suffix.lower() == ".json"
                else "program_output"
            )
            artifact = TrustedArtifactRefV1(
                artifact_id=f"result.{node_id}.{ordinal}",
                kind=kind,
                sha256=file_sha256(path),
                size_bytes=path.stat().st_size,
                path=str(path.resolve()),
                cli_value=str(path.resolve()),
            )
            artifacts.append(artifact)
            self.artifacts[artifact.artifact_id] = artifact
        return tuple(artifacts)

    @staticmethod
    def _validate_execution_outputs(
        *,
        program: str,
        jobtype: str,
        charge: int,
        multiplicity: int,
        expected_settings: Mapping[str, Any] | None = None,
        expected_input_artifact: TrustedArtifactRefV1 | None = None,
        expected_project_artifact: TrustedArtifactRefV1 | None = None,
        output_artifacts: tuple[TrustedArtifactRefV1, ...],
        exit_status: int | None,
    ) -> tuple[bool, tuple[str, ...], tuple[str, ...]]:
        findings = []
        observation: dict[str, Any] = {
            "program": program,
            "jobtype": jobtype,
            "exit_status": exit_status,
        }
        if exit_status != 0:
            findings.append("execution.process.nonzero_or_unknown")
        if program == "pyscf" and exit_status == 0:
            results = [item for item in output_artifacts if item.kind == "pyscf_hdf5"]
            if len(results) != 1:
                findings.append("pyscf.result.hdf5_count")
            else:
                try:
                    from chemsmart.io.pyscf.output import read_pyscf_h5
                    from chemsmart.jobs.pyscf.validation import frequency_validation_receipt

                    spec, provenance, status, values = read_pyscf_h5(results[0].path)
                    stages = tuple(str(item) for item in spec.get("stages") or ())
                    observation.update(
                        {
                            "stages": stages,
                            "normal_termination": bool(status.get("normal_termination")),
                            "charge": spec.get("charge"),
                            "multiplicity": spec.get("multiplicity"),
                            "requested_settings_sha256": spec.get("requested_settings_sha256"),
                            "applied_settings_sha256": spec.get("applied_settings_sha256"),
                            "environment_receipt_sha256": provenance.get("environment_receipt_sha256"),
                        }
                    )
                    if not status.get("normal_termination"):
                        findings.append("pyscf.result.normal_termination")
                    if (spec.get("charge"), spec.get("multiplicity")) != (
                        charge,
                        multiplicity,
                    ):
                        findings.append("pyscf.result.electronic_state")
                    if spec.get("requested_settings_sha256") != spec.get(
                        "applied_settings_sha256"
                    ):
                        findings.append("pyscf.result.settings_provenance")
                    energies = values.get("energies")
                    try:
                        finite_energies = tuple(float(value) for value in energies)
                    except (TypeError, ValueError):
                        finite_energies = ()
                    if not finite_energies or not all(
                        math.isfinite(value) for value in finite_energies
                    ):
                        findings.append("pyscf.result.finite_energy")
                    if expected_input_artifact is not None:
                        try:
                            input_lines = Path(
                                expected_input_artifact.path
                            ).read_text(encoding="utf-8").splitlines()
                            atom_count = int(input_lines[0].strip())
                            input_symbols = tuple(
                                line.split()[0]
                                for line in input_lines[2 : atom_count + 2]
                            )
                        except (OSError, UnicodeDecodeError, ValueError, IndexError):
                            input_symbols = ()
                        if not input_symbols or tuple(
                            str(value) for value in spec.get("symbols") or ()
                        ) != input_symbols:
                            findings.append("pyscf.result.atom_identity_order")
                    required_stage = "scf" if jobtype == "sp" else jobtype
                    stage = (status.get("stages") or {}).get(required_stage) or {}
                    if stage.get("converged") is not True:
                        findings.append(f"pyscf.result.{required_stage}_not_converged")
                    if jobtype == "hess":
                        frequency = frequency_validation_receipt(
                            symbols=spec.get("symbols"),
                            positions=values.get("positions"),
                            frequencies=values.get("vibrational_frequencies"),
                            expected_imaginary_modes=0,
                        )
                        observation["frequency_validation"] = frequency
                        if frequency.get("state") != "validated":
                            findings.append("pyscf.result.frequency_validation")
                except Exception as exc:
                    findings.append(
                        "pyscf.result.unreadable:" + type(exc).__name__
                    )
        elif program == "xtb" and exit_status == 0:
            receipts: list[Path] = []
            for artifact in output_artifacts:
                if (
                    artifact.kind != "json"
                    or "xtb-result-receipt"
                    not in Path(artifact.path).name
                ):
                    continue
                try:
                    json.loads(Path(artifact.path).read_text(encoding="utf-8"))
                    receipts.append(Path(artifact.path))
                except (OSError, json.JSONDecodeError):
                    pass
            if len(receipts) != 1:
                findings.append("xtb.result.receipt_count")
            else:
                from chemsmart.jobs.xtb.validation import (
                    audit_xtb_result_receipt,
                )

                xtb_observation, xtb_findings = audit_xtb_result_receipt(
                    receipts[0],
                    expected_jobtype=jobtype,
                    expected_charge=charge,
                    expected_multiplicity=multiplicity,
                    expected_settings=expected_settings,
                    expected_source_sha256=(
                        expected_input_artifact.sha256
                        if expected_input_artifact is not None
                        else ""
                    ),
                    expected_project_sha256=(
                        expected_project_artifact.sha256
                        if expected_project_artifact is not None
                        else ""
                    ),
                )
                observation["xtb"] = xtb_observation
                findings.extend(xtb_findings)
        elif program not in {"pyscf", "xtb"}:
            findings.append("execution.program.validator_unavailable")
        validator = canonical_sha256(
            {"observation": observation, "findings": tuple(findings)}
        )
        return not findings, (validator,), tuple(findings)

    def latest_workflow_draft_receipt(self) -> str:
        """Return the latest useful plan identity without claiming readiness."""

        if not self.workflow_drafts:
            raise ContractError("no command workflow draft has been observed")
        return next(reversed(self.workflow_drafts))

    def _inspect_calculation_artifact(self, turn_id: str, values: dict) -> Any:
        artifact = self._artifact(values["artifact_id"])
        project = self._artifact(values["project_artifact_id"])
        settings = self._get(
            self.settings_objects, values["settings_id"], "settings object"
        )
        run_receipt = self._get(
            self.run_receipts, values["run_receipt_id"], "run receipt"
        )
        receipt = inspect_generated_artifact(
            program=values["program"],
            settings=settings,
            artifact=artifact,
            project_artifact=project,
            expected_receipt=run_receipt,
        )
        self.result_inspections[receipt.receipt_sha256] = receipt
        self._emit(
            turn_id,
            EventKind.RESULT_VERIFIED,
            receipt.receipt_sha256,
            status=receipt.status,
            artifact_sha256=receipt.artifact_sha256,
            expected_receipt_sha256=receipt.expected_receipt_sha256,
        )
        return receipt

    def _record_compiled_command(
        self,
        turn_id: str,
        invocation: CanonicalCommandInvocationV1,
        context: _CommandContext,
    ) -> dict[str, Any]:
        inspection = inspect_command(
            invocation, live_schema=self.live_schema
        )
        self.invocations[invocation.invocation_sha256] = invocation
        self.command_inspections[inspection.receipt_sha256] = inspection
        self._command_contexts[invocation.invocation_sha256] = context
        self._emit(
            turn_id,
            EventKind.COMMAND_COMPILED,
            invocation.invocation_sha256,
            status=invocation.status,
            input_sha256=invocation.input_sha256,
            project_sha256=invocation.project_sha256,
            scientific_identity_sha256=(
                invocation.scientific_identity_sha256
            ),
        )
        self._emit(
            turn_id,
            EventKind.COMMAND_INSPECTED,
            inspection.receipt_sha256,
            status=inspection.status,
            invocation_sha256=inspection.invocation_sha256,
        )
        return {"invocation": invocation, "inspection": inspection}

    def _emit_binding(self, turn_id: str, kind: EventKind, value: Any) -> None:
        self.event_store.append(
            turn_id=turn_id,
            kind=kind.value,
            payload={
                "binding_sha256": value.binding_sha256,
                "state": value.state,
                "program": getattr(value, "program", ""),
            },
            idempotency_key=f"{kind.value}:{value.binding_sha256}",
        )

    def _emit(
        self,
        turn_id: str,
        kind: EventKind,
        receipt_sha256: str,
        **payload: Any,
    ) -> None:
        body = {"receipt_sha256": receipt_sha256, **payload}
        self.event_store.append(
            turn_id=turn_id,
            kind=kind.value,
            payload=body,
            idempotency_key=f"{kind.value}:{receipt_sha256}",
        )

    def _artifact(self, artifact_id: str) -> TrustedArtifactRefV1:
        return self._get(self.artifacts, artifact_id, "trusted artifact")

    @staticmethod
    def _get(values: Mapping[str, Any], key: str, label: str) -> Any:
        try:
            return values[key]
        except KeyError as exc:
            raise ContractError(f"unknown {label} ID") from exc


def _validate_tool_arguments(
    surface: AgentToolSurfaceV1, tool_name: str, arguments: dict[str, Any]
) -> None:
    definition = next(
        (
            item["function"]
            for item in surface.tool_definitions
            if item["function"]["name"] == tool_name
        ),
        None,
    )
    if definition is None:
        raise ContractError("tool is not exposed by this profile")
    schema = definition["parameters"]
    required = set(schema.get("required", ()))
    missing = sorted(required.difference(arguments))
    if missing:
        raise ContractError("tool arguments missing: " + ", ".join(missing))
    properties = schema.get("properties", {})
    unknown = sorted(set(arguments).difference(properties))
    if unknown:
        raise ContractError("tool arguments not allowed: " + ", ".join(unknown))
    for name, value in arguments.items():
        _validate_json_value(name, value, properties[name])


def _validate_json_value(name: str, value: Any, schema: Mapping[str, Any]) -> None:
    expected = schema.get("type")
    type_ok = {
        "string": isinstance(value, str),
        "integer": isinstance(value, int) and not isinstance(value, bool),
        "boolean": isinstance(value, bool),
        "array": isinstance(value, list),
        "object": isinstance(value, dict),
    }.get(expected, True)
    if not type_ok:
        raise ContractError(f"tool argument {name} has the wrong type")
    if "enum" in schema and value not in schema["enum"]:
        raise ContractError(f"tool argument {name} is outside its enum")
    if isinstance(value, str) and schema.get("pattern"):
        if re.fullmatch(str(schema["pattern"]), value) is None:
            raise ContractError(f"tool argument {name} does not match its pattern")
    if isinstance(value, int) and "minimum" in schema:
        if value < int(schema["minimum"]):
            raise ContractError(f"tool argument {name} is below its minimum")
    if isinstance(value, list) and isinstance(schema.get("items"), Mapping):
        for item in value:
            _validate_json_value(name + "[]", item, schema["items"])
    if isinstance(value, dict):
        properties = schema.get("properties", {})
        additional = schema.get("additionalProperties")
        required = set(schema.get("required", ()))
        missing = sorted(required.difference(value))
        if missing:
            raise ContractError(
                f"tool argument {name} missing: " + ", ".join(missing)
            )
        if schema.get("additionalProperties") is False:
            unknown = sorted(set(value).difference(properties))
            if unknown:
                raise ContractError(
                    f"tool argument {name} not allowed: " + ", ".join(unknown)
                )
        for child_name, child_value in value.items():
            child_schema = properties.get(child_name)
            if child_schema is None and isinstance(additional, Mapping):
                child_schema = additional
            if child_schema is not None:
                _validate_json_value(
                    f"{name}.{child_name}", child_value, child_schema
                )


def _require_registry_keys(
    values: Mapping[str, Any], attribute: str, label: str
) -> None:
    if any(key != getattr(value, attribute) for key, value in values.items()):
        raise ContractError(f"{label} registry key mismatch")


def _prepare_execution_node_workspace(path: Path) -> Path:
    """Create one empty node directory without deleting prior evidence."""

    if path.is_symlink():
        raise ContractError("execution node workspace cannot be a symlink")
    if path.exists():
        if not path.is_dir():
            raise ContractError("execution node workspace is not a directory")
        if any(path.iterdir()):
            raise ContractError(
                "execution workspace already contains outputs; import or "
                "inspect them before authorizing another launch"
            )
    else:
        path.mkdir(parents=True)
    return path


def _model_visible_data(value: Any) -> Any:
    """Remove host-local resolution details from model-visible tool results."""

    if isinstance(value, dict):
        return {
            key: _model_visible_data(item)
            for key, item in value.items()
            if key not in {"path", "cli_value", "workspace"}
        }
    if isinstance(value, (list, tuple)):
        return [_model_visible_data(item) for item in value]
    return value


__all__ = ["CommandCompiledToolHostV1"]
