"""Capability-driven program support and environment observations.

This module is a compatibility *view* of :mod:`chemsmart.settings.capabilities`.
It does not carry a second list of chemistry programs.  Agent support overlays
may only narrow the checked-out CLI declaration; they cannot add an executable
program, job type, engine, or project capability.
"""

from __future__ import annotations

import shutil
import hashlib
import json
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping

from chemsmart.agent._contracts import (
    ContractError,
    canonical_sha256,
    require_identifier,
    require_sha256,
)


class SupportLevel(str, Enum):
    AVAILABLE = "available"
    PREVIEW_ONLY = "preview_only"
    REFERENCE_ONLY = "reference_only"
    DISABLED = "disabled"


class CapabilityQueryStatus(str, Enum):
    SUPPORTED = "supported"
    PREVIEW_ONLY = "preview_only"
    REFERENCE_ONLY = "reference_only"
    DISABLED = "disabled"
    UNKNOWN_PROGRAM = "unknown_program"
    UNSUPPORTED_JOBTYPE = "unsupported_jobtype"
    UNSUPPORTED_ENGINE = "unsupported_engine"
    UNSUPPORTED_ENGINE_JOB_COMBINATION = "unsupported_engine_job_combination"
    CLI_SCHEMA_MISMATCH = "cli_schema_mismatch"


class EnvironmentStatus(str, Enum):
    AVAILABLE = "available"
    MISSING = "missing"
    NOT_DECLARED = "not_declared"
    NOT_QUERIED = "not_queried"


@dataclass(frozen=True, order=True)
class EngineJobCapabilityV1:
    """Exact engine/job support projected from the canonical registry."""

    engine: str
    jobtype: str
    preview_supported: bool = True
    execution_supported: bool = True

    def __post_init__(self) -> None:
        require_identifier(self.engine, "engine")
        require_identifier(self.jobtype, "jobtype")
        if self.execution_supported and not self.preview_supported:
            raise ContractError("execution support requires preview support")

    @property
    def pair(self) -> tuple[str, str]:
        return (self.engine, self.jobtype)


@dataclass(frozen=True)
class ProgramCapabilityV1:
    program: str
    requires_project_configuration: bool
    supports_project_configuration: bool
    jobtypes: tuple[str, ...]
    project_owned_parameters: tuple[str, ...]
    engines: tuple[str, ...]
    project_parameter_domains: tuple[tuple[str, tuple[str, ...]], ...] = ()
    #: Which YAML sections this program's loader reads.  Declared in
    #: settings/capabilities.py since the section-shape gate was added, but
    #: never projected here, so a model authoring a project had to guess the
    #: one thing the loader is strictest about.  Across this campaign sessions
    #: guessed a td: section for a phase-keyed program, guessed gas: where the
    #: loader wanted solv:, and learned the rule only from a rejection.
    project_section_names: tuple[str, ...] = ()
    engine_job_capabilities: tuple[EngineJobCapabilityV1, ...] = ()

    def __post_init__(self) -> None:
        require_identifier(self.program, "program")
        if self.requires_project_configuration and not (
            self.supports_project_configuration
        ):
            raise ContractError(
                "required project configuration is unsupported"
            )
        _require_sorted_unique(self.jobtypes, "jobtypes")
        _require_sorted_unique(
            self.project_owned_parameters, "project_owned_parameters"
        )
        _require_sorted_unique(self.engines, "engines")
        domain_names = tuple(
            name for name, _values in self.project_parameter_domains
        )
        _require_sorted_unique(domain_names, "project parameter domains")
        for name, values in self.project_parameter_domains:
            require_identifier(name, "project parameter domain")
            if name not in self.project_owned_parameters:
                raise ContractError(
                    "a parameter domain must target a project-owned parameter"
                )
            if not values or values != tuple(sorted(set(values))):
                raise ContractError(
                    "project parameter domain values must be sorted and unique"
                )
            if any(not value or value != value.lower() for value in values):
                raise ContractError(
                    "project parameter domain values must be lower-case"
                )
        _require_engine_job_capabilities(
            self.engine_job_capabilities,
            engines=self.engines,
            jobtypes=self.jobtypes,
            allow_empty=True,
        )

    @property
    def resolved_engine_job_capabilities(
        self,
    ) -> tuple[EngineJobCapabilityV1, ...]:
        if self.engine_job_capabilities:
            return self.engine_job_capabilities
        return tuple(
            EngineJobCapabilityV1(engine=engine, jobtype=jobtype)
            for engine in self.engines
            for jobtype in self.jobtypes
        )

    @property
    def preview_engine_job_pairs(self) -> tuple[tuple[str, str], ...]:
        return tuple(
            item.pair
            for item in self.resolved_engine_job_capabilities
            if item.preview_supported
        )

    @property
    def execution_engine_job_pairs(self) -> tuple[tuple[str, str], ...]:
        return tuple(
            item.pair
            for item in self.resolved_engine_job_capabilities
            if item.execution_supported
        )


@dataclass(frozen=True)
class ProgramCapabilityRegistryV1:
    schema_version: str
    programs: tuple[ProgramCapabilityV1, ...]
    registry_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.program-capability-registry.v1":
            raise ContractError("unsupported capability registry schema")
        names = tuple(item.program for item in self.programs)
        _require_sorted_unique(names, "programs")
        expected = canonical_sha256(
            {
                "schema_version": self.schema_version,
                "programs": self.programs,
            }
        )
        legacy_expected = ""
        if all(not item.engine_job_capabilities for item in self.programs):
            legacy_expected = canonical_sha256(
                {
                    "schema_version": self.schema_version,
                    "programs": tuple(
                        {
                            "program": item.program,
                            "requires_project_configuration": (
                                item.requires_project_configuration
                            ),
                            "supports_project_configuration": (
                                item.supports_project_configuration
                            ),
                            "jobtypes": item.jobtypes,
                            "project_owned_parameters": (
                                item.project_owned_parameters
                            ),
                            "engines": item.engines,
                            "project_parameter_domains": (
                                item.project_parameter_domains
                            ),
                        }
                        for item in self.programs
                    ),
                }
            )
        if self.registry_sha256 not in {expected, legacy_expected}:
            raise ContractError("capability registry digest mismatch")

    def get(self, program: str) -> ProgramCapabilityV1 | None:
        normalized = str(program or "").strip().lower()
        return next(
            (item for item in self.programs if item.program == normalized),
            None,
        )


@dataclass(frozen=True)
class ProgramSupportRuleV1:
    program: str
    support_level: SupportLevel = SupportLevel.AVAILABLE
    allowed_jobtypes: tuple[str, ...] = ()
    allowed_engines: tuple[str, ...] = ()
    allowed_engine_job_pairs: tuple[tuple[str, str], ...] = ()
    compiler_evidence_sha256: str = ""
    preview_evidence_sha256: str = ""
    preflight_evidence_sha256: str = ""
    verifier_evidence_sha256: str = ""
    execution_evidence_sha256: str = ""
    rule_ids: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        require_identifier(self.program, "program")
        _require_sorted_unique(self.allowed_jobtypes, "allowed_jobtypes")
        _require_sorted_unique(self.allowed_engines, "allowed_engines")
        _require_engine_job_pairs(
            self.allowed_engine_job_pairs,
            engines=self.allowed_engines,
            jobtypes=self.allowed_jobtypes,
            allow_empty=True,
        )
        _require_sorted_unique(self.rule_ids, "rule_ids")
        evidence = (
            self.compiler_evidence_sha256,
            self.preview_evidence_sha256,
            self.preflight_evidence_sha256,
            self.verifier_evidence_sha256,
            self.execution_evidence_sha256,
        )
        for ordinal, digest in enumerate(evidence):
            if digest:
                require_sha256(digest, f"support evidence {ordinal}")
        preview_evidence = evidence[:4]
        if self.support_level in {
            SupportLevel.AVAILABLE,
            SupportLevel.PREVIEW_ONLY,
        } and not all(preview_evidence):
            raise ContractError(
                "compiler, preview, preflight, and verifier evidence are "
                "required for an operable support rule"
            )
        if (
            self.support_level is SupportLevel.AVAILABLE
            and not self.execution_evidence_sha256
        ):
            raise ContractError(
                "available support requires execution evidence"
            )
        if (
            self.support_level is SupportLevel.PREVIEW_ONLY
            and self.execution_evidence_sha256
        ):
            raise ContractError(
                "preview-only support must not carry execution evidence"
            )


@dataclass(frozen=True)
class ProgramSupportOverlayV1:
    schema_version: str
    overlay_id: str
    base_registry_sha256: str
    rules: tuple[ProgramSupportRuleV1, ...]
    overlay_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.program-support-overlay.v1":
            raise ContractError("unsupported support overlay schema")
        require_identifier(self.overlay_id, "overlay_id")
        names = tuple(item.program for item in self.rules)
        _require_sorted_unique(names, "overlay programs")
        expected = canonical_sha256(
            {
                "schema_version": self.schema_version,
                "overlay_id": self.overlay_id,
                "base_registry_sha256": self.base_registry_sha256,
                "rules": self.rules,
            }
        )
        legacy_expected = ""
        if all(not item.allowed_engine_job_pairs for item in self.rules):
            legacy_expected = canonical_sha256(
                {
                    "schema_version": self.schema_version,
                    "overlay_id": self.overlay_id,
                    "base_registry_sha256": self.base_registry_sha256,
                    "rules": tuple(
                        {
                            "program": item.program,
                            "support_level": item.support_level,
                            "allowed_jobtypes": item.allowed_jobtypes,
                            "allowed_engines": item.allowed_engines,
                            "compiler_evidence_sha256": (
                                item.compiler_evidence_sha256
                            ),
                            "preview_evidence_sha256": (
                                item.preview_evidence_sha256
                            ),
                            "preflight_evidence_sha256": (
                                item.preflight_evidence_sha256
                            ),
                            "verifier_evidence_sha256": (
                                item.verifier_evidence_sha256
                            ),
                            "execution_evidence_sha256": (
                                item.execution_evidence_sha256
                            ),
                            "rule_ids": item.rule_ids,
                        }
                        for item in self.rules
                    ),
                }
            )
        if self.overlay_sha256 not in {expected, legacy_expected}:
            raise ContractError("support overlay digest mismatch")

    def get(self, program: str) -> ProgramSupportRuleV1 | None:
        normalized = str(program or "").strip().lower()
        return next(
            (item for item in self.rules if item.program == normalized), None
        )


@dataclass(frozen=True)
class ProgramComponentConformanceReceiptV1:
    """Observed deterministic component conformance for one program."""

    schema_version: str
    program: str
    registry_sha256: str
    live_cli_schema_sha256: str
    fixture_bundle_sha256: str
    covered_jobtypes: tuple[str, ...]
    covered_engines: tuple[str, ...]
    compiler_receipt_sha256: str
    preview_receipt_sha256: str
    preflight_receipt_sha256: str
    verifier_receipt_sha256: str
    execution_receipt_sha256: str
    compiler_status: str
    preview_status: str
    preflight_status: str
    verifier_status: str
    execution_status: str
    receipt_sha256: str
    covered_engine_job_pairs: tuple[tuple[str, str], ...] = ()

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.program-component-conformance.v1":
            raise ContractError("unsupported component conformance schema")
        require_identifier(self.program, "program")
        for field_name, digest in (
            ("registry_sha256", self.registry_sha256),
            ("live_cli_schema_sha256", self.live_cli_schema_sha256),
            ("fixture_bundle_sha256", self.fixture_bundle_sha256),
        ):
            require_sha256(digest, field_name)
        _require_sorted_unique(self.covered_jobtypes, "covered jobtypes")
        _require_sorted_unique(self.covered_engines, "covered engines")
        _require_engine_job_pairs(
            self.covered_engine_job_pairs,
            engines=self.covered_engines,
            jobtypes=self.covered_jobtypes,
            allow_empty=True,
        )
        if self.compiler_status == "passed" and (
            not self.covered_jobtypes or not self.covered_engines
        ):
            raise ContractError(
                "passed conformance requires explicit jobtype and engine coverage"
            )
        statuses = (
            ("compiler", self.compiler_status, self.compiler_receipt_sha256),
            ("preview", self.preview_status, self.preview_receipt_sha256),
            (
                "preflight",
                self.preflight_status,
                self.preflight_receipt_sha256,
            ),
            ("verifier", self.verifier_status, self.verifier_receipt_sha256),
            (
                "execution",
                self.execution_status,
                self.execution_receipt_sha256,
            ),
        )
        for name, status, digest in statuses:
            if status not in {"passed", "failed", "not_observed"}:
                raise ContractError(f"invalid {name} conformance status")
            if status == "passed":
                require_sha256(digest, f"{name}_receipt_sha256")
            elif digest:
                require_sha256(digest, f"{name}_receipt_sha256")
        body = {
            "schema_version": self.schema_version,
            "program": self.program,
            "registry_sha256": self.registry_sha256,
            "live_cli_schema_sha256": self.live_cli_schema_sha256,
            "fixture_bundle_sha256": self.fixture_bundle_sha256,
            "covered_jobtypes": self.covered_jobtypes,
            "covered_engines": self.covered_engines,
            "compiler_receipt_sha256": self.compiler_receipt_sha256,
            "preview_receipt_sha256": self.preview_receipt_sha256,
            "preflight_receipt_sha256": self.preflight_receipt_sha256,
            "verifier_receipt_sha256": self.verifier_receipt_sha256,
            "execution_receipt_sha256": self.execution_receipt_sha256,
            "compiler_status": self.compiler_status,
            "preview_status": self.preview_status,
            "preflight_status": self.preflight_status,
            "verifier_status": self.verifier_status,
            "execution_status": self.execution_status,
        }
        if self.covered_engine_job_pairs:
            body["covered_engine_job_pairs"] = self.covered_engine_job_pairs
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError(
                "component conformance receipt digest mismatch"
            )

    @property
    def effective_engine_job_pairs(self) -> tuple[tuple[str, str], ...]:
        if self.covered_engine_job_pairs:
            return self.covered_engine_job_pairs
        return tuple(
            (engine, jobtype)
            for engine in self.covered_engines
            for jobtype in self.covered_jobtypes
        )


@dataclass(frozen=True)
class CapabilityQueryV1:
    program: str
    jobtype: str = ""
    engine: str = ""

    def __post_init__(self) -> None:
        require_identifier(self.program, "program")
        if self.jobtype:
            require_identifier(self.jobtype, "jobtype")
        if self.engine:
            require_identifier(self.engine, "engine")


@dataclass(frozen=True)
class CapabilityQueryReceiptV1:
    schema_version: str
    query: CapabilityQueryV1
    registry_sha256: str
    live_cli_schema_sha256: str
    joined_capability_sha256: str
    overlay_sha256: str
    status: CapabilityQueryStatus
    capability: ProgramCapabilityV1 | None
    effective_jobtypes: tuple[str, ...]
    effective_engines: tuple[str, ...]
    rule_ids: tuple[str, ...]
    receipt_sha256: str
    effective_engine_job_pairs: tuple[tuple[str, str], ...] = ()

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.capability-query-receipt.v1":
            raise ContractError("unsupported capability query receipt schema")
        _require_engine_job_pairs(
            self.effective_engine_job_pairs,
            engines=self.effective_engines,
            jobtypes=self.effective_jobtypes,
            allow_empty=True,
        )
        body = {
            "schema_version": self.schema_version,
            "query": self.query,
            "registry_sha256": self.registry_sha256,
            "live_cli_schema_sha256": self.live_cli_schema_sha256,
            "joined_capability_sha256": self.joined_capability_sha256,
            "overlay_sha256": self.overlay_sha256,
            "status": self.status,
            "capability": self.capability,
            "effective_jobtypes": self.effective_jobtypes,
            "effective_engines": self.effective_engines,
            "rule_ids": self.rule_ids,
        }
        if self.effective_engine_job_pairs:
            body["effective_engine_job_pairs"] = (
                self.effective_engine_job_pairs
            )
        expected = canonical_sha256(body)
        if self.receipt_sha256 != expected:
            raise ContractError("capability query receipt digest mismatch")


@dataclass(frozen=True)
class EnvironmentTargetV1:
    """Trusted host configuration for a non-executing environment probe."""

    program: str
    engine: str
    target_kind: str
    locator: str
    distribution: str = ""
    required_dependencies: tuple[str, ...] = ()
    required_dependency_versions: tuple[tuple[str, str], ...] = ()
    required_gpu_facts: tuple[str, ...] = ()
    required_gpu_values: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        require_identifier(self.program, "program")
        require_identifier(self.engine, "engine")
        if self.target_kind not in {
            "python_module",
            "executable",
            "compute_receipt",
        }:
            raise ContractError(
                "target_kind must be python_module, executable, or compute_receipt"
            )
        if not str(self.locator).strip():
            raise ContractError("environment target locator must not be empty")
        _require_sorted_unique(
            self.required_dependencies, "required_dependencies"
        )
        version_names = tuple(
            item[0] for item in self.required_dependency_versions
        )
        _require_sorted_unique(version_names, "required_dependency_versions")
        _require_sorted_unique(self.required_gpu_facts, "required_gpu_facts")
        _require_sorted_unique(self.required_gpu_values, "required_gpu_values")


@dataclass(frozen=True)
class ProgramEnvironmentQueryV1:
    capability_receipt_sha256: str
    program: str
    engine: str

    def __post_init__(self) -> None:
        require_identifier(self.program, "program")
        if self.engine:
            require_identifier(self.engine, "engine")


@dataclass(frozen=True)
class TrustedComputeEnvironmentReceiptV1:
    """Evidence produced for the exact interpreter that will run a job."""

    schema_version: str
    program: str
    engine: str
    compute_interpreter_sha256: str
    dependency_versions: tuple[tuple[str, str], ...]
    solver_evidence: tuple[tuple[str, bool], ...]
    gpu_evidence: tuple[tuple[str, Any], ...]
    source_probe_id: str
    evidence_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.compute-environment-receipt.v1":
            raise ContractError(
                "unsupported compute environment receipt schema"
            )
        require_identifier(self.program, "program")
        require_identifier(self.engine, "engine")
        require_sha256(
            self.compute_interpreter_sha256,
            "compute_interpreter_sha256",
        )
        if not self.source_probe_id:
            raise ContractError("compute environment source probe is required")
        names = tuple(item[0] for item in self.dependency_versions)
        if names != tuple(sorted(set(names))):
            raise ContractError(
                "dependency evidence must be sorted and unique"
            )
        solver_names = tuple(item[0] for item in self.solver_evidence)
        if solver_names != tuple(sorted(set(solver_names))):
            raise ContractError("solver evidence must be sorted and unique")
        gpu_names = tuple(item[0] for item in self.gpu_evidence)
        if gpu_names != tuple(sorted(set(gpu_names))):
            raise ContractError("GPU evidence must be sorted and unique")
        expected = canonical_sha256(
            {
                "schema_version": self.schema_version,
                "program": self.program,
                "engine": self.engine,
                "compute_interpreter_sha256": self.compute_interpreter_sha256,
                "dependency_versions": self.dependency_versions,
                "solver_evidence": self.solver_evidence,
                "gpu_evidence": self.gpu_evidence,
                "source_probe_id": self.source_probe_id,
            }
        )
        if self.evidence_sha256 != expected:
            raise ContractError("compute environment evidence digest mismatch")


@dataclass(frozen=True)
class EnvironmentCapabilityReceiptV1:
    schema_version: str
    query: ProgramEnvironmentQueryV1
    capability_receipt_sha256: str
    program: str
    engine: str
    status: EnvironmentStatus
    target_kind: str
    locator: str
    observed_version: str
    observed_location_sha256: str
    compute_interpreter_sha256: str
    compute_evidence_sha256: str
    dependency_versions: tuple[tuple[str, str], ...]
    solver_evidence: tuple[tuple[str, bool], ...]
    gpu_evidence: tuple[tuple[str, Any], ...]
    observation_method: str
    rule_ids: tuple[str, ...]
    receipt_sha256: str

    def __post_init__(self) -> None:
        if (
            self.schema_version
            != "chemsmart.environment-capability-receipt.v1"
        ):
            raise ContractError("unsupported environment receipt schema")
        if self.status is EnvironmentStatus.AVAILABLE and not (
            self.compute_interpreter_sha256 or self.observed_location_sha256
        ):
            raise ContractError(
                "available environment requires target identity evidence"
            )
        expected = canonical_sha256(
            {
                "schema_version": self.schema_version,
                "query": self.query,
                "capability_receipt_sha256": self.capability_receipt_sha256,
                "program": self.program,
                "engine": self.engine,
                "status": self.status,
                "target_kind": self.target_kind,
                "locator": self.locator,
                "observed_version": self.observed_version,
                "observed_location_sha256": self.observed_location_sha256,
                "compute_interpreter_sha256": (
                    self.compute_interpreter_sha256
                ),
                "compute_evidence_sha256": self.compute_evidence_sha256,
                "dependency_versions": self.dependency_versions,
                "solver_evidence": self.solver_evidence,
                "gpu_evidence": self.gpu_evidence,
                "observation_method": self.observation_method,
                "rule_ids": self.rule_ids,
            }
        )
        if self.receipt_sha256 != expected:
            raise ContractError("environment receipt digest mismatch")


@dataclass(frozen=True)
class ProgramCandidateProposalV1:
    schema_version: str
    requested_program: str
    candidate_program: str
    candidate_engine: str
    source_ids: tuple[str, ...]
    proposal_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.program-candidate-proposal.v1":
            raise ContractError("unsupported program candidate schema")
        require_identifier(self.requested_program, "requested_program")
        require_identifier(self.candidate_program, "candidate_program")
        if self.candidate_engine:
            require_identifier(self.candidate_engine, "candidate_engine")
        expected = canonical_sha256(
            {
                "schema_version": self.schema_version,
                "requested_program": self.requested_program,
                "candidate_program": self.candidate_program,
                "candidate_engine": self.candidate_engine,
                "source_ids": self.source_ids,
            }
        )
        if self.proposal_sha256 != expected:
            raise ContractError("program candidate proposal digest mismatch")


@dataclass(frozen=True)
class ResolvedProgramBindingV1:
    schema_version: str
    requested_program: str
    selected_program: str
    capability_receipt_sha256: str
    substitution_receipt_sha256: str
    joined_capability_sha256: str
    state: str
    rule_ids: tuple[str, ...]
    binding_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.resolved-program-binding.v1":
            raise ContractError("unsupported resolved program binding schema")
        require_identifier(self.requested_program, "requested_program")
        require_identifier(self.selected_program, "selected_program")
        if self.state not in {"resolved", "preview_only", "blocked"}:
            raise ContractError("invalid resolved program binding state")
        expected = canonical_sha256(
            {
                "schema_version": self.schema_version,
                "requested_program": self.requested_program,
                "selected_program": self.selected_program,
                "capability_receipt_sha256": self.capability_receipt_sha256,
                "substitution_receipt_sha256": (
                    self.substitution_receipt_sha256
                ),
                "joined_capability_sha256": self.joined_capability_sha256,
                "state": self.state,
                "rule_ids": self.rule_ids,
            }
        )
        if self.binding_sha256 != expected:
            raise ContractError("resolved program binding digest mismatch")

    @property
    def program(self) -> str:
        return self.selected_program


@dataclass(frozen=True)
class ResolvedEngineBindingV1:
    schema_version: str
    program: str
    requested_engine: str
    selected_engine: str
    program_binding_sha256: str
    capability_receipt_sha256: str
    environment_receipt_sha256: str
    state: str
    execution_ready: bool
    rule_ids: tuple[str, ...]
    binding_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.resolved-engine-binding.v1":
            raise ContractError("unsupported resolved engine binding schema")
        require_identifier(self.program, "program")
        require_identifier(self.requested_engine, "requested_engine")
        require_identifier(self.selected_engine, "selected_engine")
        if self.state not in {"resolved", "preview_only", "blocked"}:
            raise ContractError("invalid program-engine binding state")
        expected = canonical_sha256(
            {
                "schema_version": self.schema_version,
                "program": self.program,
                "requested_engine": self.requested_engine,
                "selected_engine": self.selected_engine,
                "program_binding_sha256": self.program_binding_sha256,
                "capability_receipt_sha256": (self.capability_receipt_sha256),
                "environment_receipt_sha256": self.environment_receipt_sha256,
                "state": self.state,
                "execution_ready": self.execution_ready,
                "rule_ids": self.rule_ids,
            }
        )
        if self.binding_sha256 != expected:
            raise ContractError("resolved engine binding digest mismatch")

    @property
    def engine(self) -> str:
        return self.selected_engine


def load_program_capabilities(
    registry_module: Any | None = None,
) -> ProgramCapabilityRegistryV1:
    """Project the checked-out settings registry into the agent contract."""

    if registry_module is None:
        try:
            from chemsmart.settings import capabilities as registry_module
        except ImportError as exc:  # pragma: no cover - integration boundary
            raise ContractError(
                "chemsmart.settings.capabilities is not available"
            ) from exc
    raw_registry = getattr(registry_module, "PROGRAM_CAPABILITIES", None)
    if raw_registry is None:
        raw_registry = getattr(registry_module, "ENGINE_CAPABILITIES", None)
    if not isinstance(raw_registry, Mapping):
        raise ContractError("settings capability registry is not a mapping")

    programs = []
    for raw_name, raw in sorted(
        raw_registry.items(), key=lambda item: str(item[0])
    ):
        program = require_identifier(str(raw_name), "program")
        declared = require_identifier(
            str(getattr(raw, "program", program)), "declared program"
        )
        if declared != program:
            raise ContractError("capability mapping key differs from program")
        jobtypes = _normalized_tuple(getattr(raw, "jobtypes", ()), "jobtypes")
        engines = _normalized_tuple(getattr(raw, "engines", ()), "engines")
        raw_matrix = tuple(getattr(raw, "engine_job_capabilities", ()) or ())
        if raw_matrix:
            engine_job_capabilities = tuple(
                sorted(
                    EngineJobCapabilityV1(
                        engine=require_identifier(
                            str(getattr(item, "engine", "")), "engine"
                        ),
                        jobtype=require_identifier(
                            str(getattr(item, "jobtype", "")), "jobtype"
                        ),
                        preview_supported=bool(
                            getattr(item, "preview_supported", True)
                        ),
                        execution_supported=bool(
                            getattr(item, "execution_supported", True)
                        ),
                    )
                    for item in raw_matrix
                )
            )
        else:
            # Legacy settings registries declared independent lists. Project
            # them once into an explicit matrix so every downstream decision
            # is pair-aware without changing old registry authorship APIs.
            engine_job_capabilities = tuple(
                EngineJobCapabilityV1(engine=engine, jobtype=jobtype)
                for engine in engines
                for jobtype in jobtypes
            )
        programs.append(
            ProgramCapabilityV1(
                program=program,
                requires_project_configuration=bool(
                    getattr(raw, "requires_project_configuration", False)
                ),
                supports_project_configuration=bool(
                    getattr(raw, "supports_project_configuration", False)
                ),
                jobtypes=jobtypes,
                project_owned_parameters=_normalized_tuple(
                    getattr(raw, "project_owned_parameters", ()),
                    "project_owned_parameters",
                ),
                engines=engines,
                project_section_names=_normalized_tuple(
                    getattr(raw, "project_section_names", ()),
                    "project_section_names",
                ),
                project_parameter_domains=tuple(
                    (
                        require_identifier(
                            str(name), "project parameter domain"
                        ),
                        tuple(
                            sorted(
                                {
                                    str(value).strip().lower()
                                    for value in values
                                    if str(value).strip()
                                }
                            )
                        ),
                    )
                    for name, values in getattr(
                        raw, "project_parameter_domains", ()
                    )
                ),
                engine_job_capabilities=engine_job_capabilities,
            )
        )
    body = {
        "schema_version": "chemsmart.program-capability-registry.v1",
        "programs": tuple(programs),
    }
    return ProgramCapabilityRegistryV1(
        **body, registry_sha256=canonical_sha256(body)
    )


def build_support_overlay(
    *,
    overlay_id: str,
    registry: ProgramCapabilityRegistryV1,
    rules: Iterable[ProgramSupportRuleV1],
) -> ProgramSupportOverlayV1:
    ordered = tuple(sorted(rules, key=lambda item: item.program))
    for rule in ordered:
        capability = registry.get(rule.program)
        if capability is None:
            raise ContractError(
                "support overlays cannot add an undeclared executable program"
            )
        if rule.allowed_jobtypes and not set(rule.allowed_jobtypes).issubset(
            capability.jobtypes
        ):
            raise ContractError("support overlay broadens declared jobtypes")
        if rule.allowed_engines and not set(rule.allowed_engines).issubset(
            capability.engines
        ):
            raise ContractError("support overlay broadens declared engines")
        allowed_pairs = _effective_rule_engine_job_pairs(rule)
        if not set(allowed_pairs).issubset(
            capability.preview_engine_job_pairs
        ):
            raise ContractError(
                "support overlay broadens declared engine-job capability"
            )
        if rule.support_level is SupportLevel.AVAILABLE and not set(
            allowed_pairs
        ).issubset(capability.execution_engine_job_pairs):
            raise ContractError(
                "execution overlay includes a preview-only engine-job pair"
            )
    body = {
        "schema_version": "chemsmart.program-support-overlay.v1",
        "overlay_id": require_identifier(overlay_id, "overlay_id"),
        "base_registry_sha256": registry.registry_sha256,
        "rules": ordered,
    }
    return ProgramSupportOverlayV1(
        **body, overlay_sha256=canonical_sha256(body)
    )


def build_program_component_conformance_receipt(
    *,
    program: str,
    registry_sha256: str,
    live_cli_schema_sha256: str,
    fixture_bundle_sha256: str,
    covered_jobtypes: Iterable[str],
    covered_engines: Iterable[str],
    compiler_receipt_sha256: str,
    preview_receipt_sha256: str,
    preflight_receipt_sha256: str,
    verifier_receipt_sha256: str,
    compiler_status: str,
    preview_status: str,
    preflight_status: str,
    verifier_status: str,
    execution_receipt_sha256: str = "",
    execution_status: str = "not_observed",
    covered_engine_job_pairs: Iterable[tuple[str, str]] | None = None,
) -> ProgramComponentConformanceReceiptV1:
    normalized_jobtypes = tuple(sorted(set(covered_jobtypes)))
    normalized_engines = tuple(sorted(set(covered_engines)))
    normalized_pairs = (
        tuple(
            (engine, jobtype)
            for engine in normalized_engines
            for jobtype in normalized_jobtypes
        )
        if covered_engine_job_pairs is None
        else tuple(sorted(set(covered_engine_job_pairs)))
    )
    body = {
        "schema_version": "chemsmart.program-component-conformance.v1",
        "program": require_identifier(program, "program"),
        "registry_sha256": registry_sha256,
        "live_cli_schema_sha256": live_cli_schema_sha256,
        "fixture_bundle_sha256": fixture_bundle_sha256,
        "covered_jobtypes": normalized_jobtypes,
        "covered_engines": normalized_engines,
        "compiler_receipt_sha256": compiler_receipt_sha256,
        "preview_receipt_sha256": preview_receipt_sha256,
        "preflight_receipt_sha256": preflight_receipt_sha256,
        "verifier_receipt_sha256": verifier_receipt_sha256,
        "execution_receipt_sha256": execution_receipt_sha256,
        "compiler_status": compiler_status,
        "preview_status": preview_status,
        "preflight_status": preflight_status,
        "verifier_status": verifier_status,
        "execution_status": execution_status,
    }
    if normalized_pairs:
        body["covered_engine_job_pairs"] = normalized_pairs
    return ProgramComponentConformanceReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def build_command_compiled_preview_overlay(
    registry: ProgramCapabilityRegistryV1 | None = None,
    *,
    conformance_receipts: Iterable[ProgramComponentConformanceReceiptV1] = (),
    live_schema: Any | None = None,
) -> ProgramSupportOverlayV1:
    """Build an observation-bound, fail-closed preview support declaration.

    Source presence is deliberately irrelevant. A program becomes previewable
    only after compiler, preview, preflight, and verifier conformance receipts
    passed against the exact registry, Click schema, and fixture bundle.
    """

    registry = registry or load_program_capabilities()
    if live_schema is None:
        from chemsmart.agent.cli_schema import build_live_click_schema

        live_schema = build_live_click_schema()
    conformance_receipts = tuple(conformance_receipts)
    evidence_by_program = {item.program: item for item in conformance_receipts}
    if len(evidence_by_program) != len(conformance_receipts):
        raise ContractError("program conformance receipts must be unique")
    rules = []
    for capability in registry.programs:
        evidence = evidence_by_program.get(capability.program)
        current = bool(
            evidence is not None
            and evidence.registry_sha256 == registry.registry_sha256
            and evidence.live_cli_schema_sha256 == live_schema.schema_sha256
        )
        if evidence is not None and (
            not set(evidence.covered_jobtypes).issubset(capability.jobtypes)
            or not set(evidence.covered_engines).issubset(capability.engines)
            or not set(evidence.effective_engine_job_pairs).issubset(
                capability.preview_engine_job_pairs
            )
        ):
            raise ContractError(
                "program conformance receipt broadens declared capability"
            )
        operable = bool(
            current
            and evidence is not None
            and evidence.compiler_status == "passed"
            and evidence.preview_status == "passed"
            and evidence.preflight_status == "passed"
            and evidence.verifier_status == "passed"
            and bool(evidence.effective_engine_job_pairs)
        )
        rules.append(
            ProgramSupportRuleV1(
                program=capability.program,
                support_level=(
                    SupportLevel.PREVIEW_ONLY
                    if operable
                    else SupportLevel.REFERENCE_ONLY
                ),
                allowed_jobtypes=(
                    evidence.covered_jobtypes if operable else ()
                ),
                allowed_engines=(evidence.covered_engines if operable else ()),
                allowed_engine_job_pairs=(
                    evidence.effective_engine_job_pairs if operable else ()
                ),
                compiler_evidence_sha256=(
                    evidence.compiler_receipt_sha256 if operable else ""
                ),
                preview_evidence_sha256=(
                    evidence.preview_receipt_sha256 if operable else ""
                ),
                preflight_evidence_sha256=(
                    evidence.preflight_receipt_sha256 if operable else ""
                ),
                verifier_evidence_sha256=(
                    evidence.verifier_receipt_sha256 if operable else ""
                ),
                rule_ids=(
                    ("agent.support.conformance_receipts_passed",)
                    if operable
                    else (
                        (
                            "agent.support.conformance_stale"
                            if evidence is not None and not current
                            else "agent.support.conformance_missing_or_red"
                        ),
                    )
                ),
            )
        )
    return build_support_overlay(
        overlay_id="command_compiled_preview",
        registry=registry,
        rules=rules,
    )


def environment_identity_sha256(receipt: Any) -> str:
    """Identify the *machine*, independently of what was asked of it.

    ``EnvironmentCapabilityReceiptV1.receipt_sha256`` folds in
    ``capability_receipt_sha256``, and a capability receipt changes with the
    active overlay.  A plan session and an execution session therefore observe
    the same interpreter and produce different environment digests: measured on
    one PySCF node, ``compute_evidence_sha256``,
    ``compute_interpreter_sha256``, ``locator`` and every dependency version
    were byte-identical while the receipt digests were ``8aed9113…`` and
    ``cc554e28…``.

    An approval that pinned the receipt digest therefore rejected the very
    machine it had approved, and no reviewer could have supplied the right
    digest, because it does not exist until execution is already authorised.
    This is the environment's own identity: the fields that say which
    interpreter, which versions, and which accelerator, and nothing about who
    was asking.
    """

    return canonical_sha256(
        {
            "schema_version": "chemsmart.environment-identity.v1",
            "program": getattr(receipt, "program", ""),
            "engine": getattr(receipt, "engine", ""),
            "target_kind": getattr(receipt, "target_kind", ""),
            "locator": getattr(receipt, "locator", ""),
            "observation_method": getattr(receipt, "observation_method", ""),
            "observed_version": getattr(receipt, "observed_version", ""),
            "observed_location_sha256": getattr(
                receipt, "observed_location_sha256", ""
            ),
            "compute_interpreter_sha256": getattr(
                receipt, "compute_interpreter_sha256", ""
            ),
            "compute_evidence_sha256": getattr(
                receipt, "compute_evidence_sha256", ""
            ),
            "dependency_versions": tuple(
                tuple(item)
                for item in getattr(receipt, "dependency_versions", ())
            ),
            "solver_evidence": tuple(
                tuple(item) for item in getattr(receipt, "solver_evidence", ())
            ),
            "gpu_evidence": tuple(
                tuple(item) for item in getattr(receipt, "gpu_evidence", ())
            ),
        }
    )


def build_approved_execution_overlay(
    *,
    registry: ProgramCapabilityRegistryV1,
    preview_overlay: ProgramSupportOverlayV1,
    approved_nodes: Iterable[tuple[str, str, str]],
    execution_evidence_sha256: str,
) -> ProgramSupportOverlayV1:
    """Upgrade only explicitly approved, preview-conformant nodes.

    This is the commissioning transition used before a program has historical
    execution receipts.  It does not infer availability from source presence:
    every upgraded program must already have a green preview rule, and the
    normal environment query still decides whether an engine binding is
    execution-ready for a concrete node.
    """

    if preview_overlay.base_registry_sha256 != registry.registry_sha256:
        raise ContractError("preview overlay is stale for this registry")
    require_sha256(execution_evidence_sha256, "execution_evidence_sha256")
    approved: dict[str, set[tuple[str, str]]] = {}
    for program, jobtype, engine in approved_nodes:
        normalized_program = require_identifier(program, "program")
        normalized_jobtype = require_identifier(jobtype, "jobtype")
        normalized_engine = require_identifier(engine, "engine")
        capability = registry.get(normalized_program)
        if capability is None:
            raise ContractError("approved execution uses an unknown program")
        if normalized_jobtype not in capability.jobtypes:
            raise ContractError("approved execution uses an unknown jobtype")
        if normalized_engine not in capability.engines:
            raise ContractError("approved execution uses an unknown engine")
        pair = (normalized_engine, normalized_jobtype)
        if pair not in capability.execution_engine_job_pairs:
            raise ContractError(
                "approved execution uses a preview-only engine-job pair"
            )
        approved.setdefault(normalized_program, set()).add(pair)

    rules = []
    for preview_rule in preview_overlay.rules:
        selection = approved.get(preview_rule.program)
        if selection is None:
            rules.append(preview_rule)
            continue
        if preview_rule.support_level is not SupportLevel.PREVIEW_ONLY:
            raise ContractError(
                "approved execution requires prior preview conformance"
            )
        pairs = tuple(sorted(selection))
        jobtypes = tuple(sorted({jobtype for _engine, jobtype in pairs}))
        engines = tuple(sorted({engine for engine, _jobtype in pairs}))
        if not set(jobtypes).issubset(preview_rule.allowed_jobtypes):
            raise ContractError("approved jobtype lacks preview conformance")
        if not set(engines).issubset(preview_rule.allowed_engines):
            raise ContractError("approved engine lacks preview conformance")
        if not set(pairs).issubset(
            _effective_rule_engine_job_pairs(preview_rule)
        ):
            raise ContractError(
                "approved engine-job pair lacks preview conformance"
            )
        rules.append(
            ProgramSupportRuleV1(
                program=preview_rule.program,
                support_level=SupportLevel.AVAILABLE,
                allowed_jobtypes=jobtypes,
                allowed_engines=engines,
                allowed_engine_job_pairs=pairs,
                compiler_evidence_sha256=(
                    preview_rule.compiler_evidence_sha256
                ),
                preview_evidence_sha256=preview_rule.preview_evidence_sha256,
                preflight_evidence_sha256=(
                    preview_rule.preflight_evidence_sha256
                ),
                verifier_evidence_sha256=(
                    preview_rule.verifier_evidence_sha256
                ),
                execution_evidence_sha256=execution_evidence_sha256,
                rule_ids=(
                    "agent.support.approved_commissioning_execution",
                    "agent.support.conformance_receipts_passed",
                ),
            )
        )
    if set(approved) - {item.program for item in preview_overlay.rules}:
        raise ContractError("approved program lacks a preview support rule")
    return build_support_overlay(
        overlay_id="approved_commissioning_execution",
        registry=registry,
        rules=rules,
    )


def build_program_candidate_proposal(
    *,
    requested_program: str,
    candidate_program: str,
    candidate_engine: str = "",
    source_ids: Iterable[str] = (),
) -> ProgramCandidateProposalV1:
    """Record an advisory candidate; this never establishes support."""

    body = {
        "schema_version": "chemsmart.program-candidate-proposal.v1",
        "requested_program": require_identifier(
            requested_program, "requested_program"
        ),
        "candidate_program": require_identifier(
            candidate_program, "candidate_program"
        ),
        "candidate_engine": (
            require_identifier(candidate_engine, "candidate_engine")
            if candidate_engine
            else ""
        ),
        "source_ids": tuple(sorted(set(str(item) for item in source_ids))),
    }
    return ProgramCandidateProposalV1(
        **body, proposal_sha256=canonical_sha256(body)
    )


def query_capability(
    query: CapabilityQueryV1,
    *,
    registry: ProgramCapabilityRegistryV1 | None = None,
    overlay: ProgramSupportOverlayV1 | None = None,
    live_schema: Any | None = None,
) -> CapabilityQueryReceiptV1:
    registry = registry or load_program_capabilities()
    if live_schema is None:
        from chemsmart.agent.cli_schema import build_live_click_schema

        live_schema = build_live_click_schema()
    live_cli_schema_sha256 = str(live_schema.schema_sha256)
    joined_capability_sha256 = canonical_sha256(
        {
            "registry_sha256": registry.registry_sha256,
            "live_cli_schema_sha256": live_cli_schema_sha256,
        }
    )
    if overlay is not None and (
        overlay.base_registry_sha256 != registry.registry_sha256
    ):
        raise ContractError("support overlay is stale for this registry")
    capability = registry.get(query.program)
    effective_jobtypes: tuple[str, ...] = ()
    effective_engines: tuple[str, ...] = ()
    effective_engine_job_pairs: tuple[tuple[str, str], ...] = ()
    rule_ids: tuple[str, ...] = ()
    status = CapabilityQueryStatus.UNKNOWN_PROGRAM
    overlay_sha256 = overlay.overlay_sha256 if overlay is not None else ""

    if capability is not None:
        effective_jobtypes = capability.jobtypes
        effective_engines = capability.engines
        effective_engine_job_pairs = capability.preview_engine_job_pairs
        rule = overlay.get(query.program) if overlay is not None else None
        level = SupportLevel.DISABLED
        if rule is not None:
            level = rule.support_level
            rule_ids = rule.rule_ids
            if rule.allowed_jobtypes:
                effective_jobtypes = rule.allowed_jobtypes
            if rule.allowed_engines:
                effective_engines = rule.allowed_engines
            rule_pairs = _effective_rule_engine_job_pairs(rule)
            if (
                rule_pairs
                or rule.allowed_jobtypes
                or rule.allowed_engines
                or rule.allowed_engine_job_pairs
            ):
                effective_engine_job_pairs = rule_pairs
            if rule_pairs:
                effective_jobtypes = tuple(
                    sorted({jobtype for _engine, jobtype in rule_pairs})
                )
                effective_engines = tuple(
                    sorted({engine for engine, _jobtype in rule_pairs})
                )
        if overlay is None or rule is None:
            status = CapabilityQueryStatus.DISABLED
            rule_ids = ("agent.support.overlay_required",)
        elif level is SupportLevel.DISABLED:
            status = CapabilityQueryStatus.DISABLED
        elif query.jobtype and query.jobtype not in effective_jobtypes:
            status = CapabilityQueryStatus.UNSUPPORTED_JOBTYPE
        elif query.engine and query.engine not in effective_engines:
            status = CapabilityQueryStatus.UNSUPPORTED_ENGINE
        elif (
            query.jobtype
            and query.engine
            and (query.engine, query.jobtype) not in effective_engine_job_pairs
        ):
            status = CapabilityQueryStatus.UNSUPPORTED_ENGINE_JOB_COMBINATION
        elif level is SupportLevel.REFERENCE_ONLY:
            status = CapabilityQueryStatus.REFERENCE_ONLY
        elif query.jobtype and not (
            live_schema.has_path(("run", query.program, query.jobtype))
            and live_schema.has_path(("sub", query.program, query.jobtype))
        ):
            status = CapabilityQueryStatus.CLI_SCHEMA_MISMATCH
        elif not query.jobtype and not (
            live_schema.has_path(("run", query.program))
            and live_schema.has_path(("sub", query.program))
        ):
            status = CapabilityQueryStatus.CLI_SCHEMA_MISMATCH
        elif level is SupportLevel.PREVIEW_ONLY:
            status = CapabilityQueryStatus.PREVIEW_ONLY
        else:
            status = CapabilityQueryStatus.SUPPORTED

    body = {
        "schema_version": "chemsmart.capability-query-receipt.v1",
        "query": query,
        "registry_sha256": registry.registry_sha256,
        "live_cli_schema_sha256": live_cli_schema_sha256,
        "joined_capability_sha256": joined_capability_sha256,
        "overlay_sha256": overlay_sha256,
        "status": status,
        "capability": capability,
        "effective_jobtypes": effective_jobtypes,
        "effective_engines": effective_engines,
        "rule_ids": rule_ids,
    }
    if effective_engine_job_pairs:
        body["effective_engine_job_pairs"] = effective_engine_job_pairs
    return CapabilityQueryReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def query_environment(
    capability_receipt: CapabilityQueryReceiptV1,
    *,
    targets: Iterable[EnvironmentTargetV1] = (),
    compute_receipts: Iterable[TrustedComputeEnvironmentReceiptV1] = (),
    which: Callable[[str], str | None] = shutil.which,
) -> EnvironmentCapabilityReceiptV1:
    """Observe module/binary presence without importing or executing it."""

    program = capability_receipt.query.program
    engine = capability_receipt.query.engine
    if not engine and capability_receipt.effective_engines:
        if len(capability_receipt.effective_engines) != 1:
            return _environment_receipt(
                capability_receipt,
                engine="",
                status=EnvironmentStatus.NOT_QUERIED,
                rule_ids=("environment.engine.explicit_required",),
            )
        engine = capability_receipt.effective_engines[0]
    target = next(
        (
            item
            for item in targets
            if item.program == program and item.engine == engine
        ),
        None,
    )
    if target is None:
        return _environment_receipt(
            capability_receipt,
            engine=engine,
            status=EnvironmentStatus.NOT_DECLARED,
            rule_ids=("environment.target.not_declared",),
        )

    trusted_compute = next(
        (
            item
            for item in compute_receipts
            if item.program == program and item.engine == engine
        ),
        None,
    )
    if target.target_kind in {"python_module", "compute_receipt"}:
        if trusted_compute is None:
            return _environment_receipt(
                capability_receipt,
                engine=engine,
                status=EnvironmentStatus.NOT_QUERIED,
                target=target,
                observation_method="target_compute_receipt_required",
                rule_ids=("environment.compute_receipt.required",),
            )
        gpu_facts = dict(trusted_compute.gpu_evidence)
        dependencies = dict(trusted_compute.dependency_versions)
        missing_rules = []
        for dependency in target.required_dependencies:
            if not dependencies.get(dependency):
                missing_rules.append(
                    f"environment.dependency.{dependency}_missing"
                )
        for (
            dependency,
            expected_version,
        ) in target.required_dependency_versions:
            if dependencies.get(dependency) != expected_version:
                missing_rules.append(
                    f"environment.dependency.{dependency}_version_mismatch"
                )
        for fact in target.required_gpu_facts:
            if gpu_facts.get(fact) is not True:
                missing_rules.append(f"environment.gpu.{fact}_missing")
        for fact in target.required_gpu_values:
            if gpu_facts.get(fact) in {None, "", 0, False}:
                missing_rules.append(f"environment.gpu.{fact}_missing")
        available = not missing_rules
        return _environment_receipt(
            capability_receipt,
            engine=engine,
            status=(
                EnvironmentStatus.AVAILABLE
                if available
                else EnvironmentStatus.MISSING
            ),
            target=target,
            compute=trusted_compute,
            observation_method="trusted_target_compute_receipt",
            rule_ids=(
                ("environment.compute_receipt.complete",)
                if available
                else tuple(sorted(missing_rules))
            ),
        )

    location = ""
    version = ""
    available = False
    method = "shutil.which"
    resolved = which(target.locator)
    available = bool(resolved)
    location = str(resolved or "")

    return _environment_receipt(
        capability_receipt,
        engine=engine,
        status=(
            EnvironmentStatus.AVAILABLE
            if available
            else EnvironmentStatus.MISSING
        ),
        target=target,
        observed_version=version,
        observed_location_sha256=(
            canonical_sha256({"location": location}) if location else ""
        ),
        observation_method=method,
        rule_ids=(
            (
                "environment.target.available"
                if available
                else "environment.target.missing"
            ),
        ),
    )


def build_trusted_compute_environment_receipt(
    *,
    program: str,
    engine: str,
    compute_interpreter_sha256: str,
    dependency_versions: Mapping[str, str],
    solver_evidence: Mapping[str, bool] | None = None,
    gpu_evidence: Mapping[str, Any] | None = None,
    source_probe_id: str,
) -> TrustedComputeEnvironmentReceiptV1:
    """Normalize a program-owned target-interpreter probe into M2 evidence."""

    body = {
        "schema_version": "chemsmart.compute-environment-receipt.v1",
        "program": require_identifier(program, "program"),
        "engine": require_identifier(engine, "engine"),
        "compute_interpreter_sha256": str(compute_interpreter_sha256),
        "dependency_versions": tuple(
            sorted(
                (require_identifier(name, "dependency"), str(version))
                for name, version in dependency_versions.items()
            )
        ),
        "solver_evidence": tuple(
            sorted(
                (require_identifier(name, "solver"), bool(available))
                for name, available in (solver_evidence or {}).items()
            )
        ),
        "gpu_evidence": tuple(
            sorted(
                (
                    require_identifier(name, "gpu evidence"),
                    (
                        canonical_sha256(value)
                        if isinstance(value, (dict, list, tuple))
                        else value
                    ),
                )
                for name, value in (gpu_evidence or {}).items()
            )
        ),
        "source_probe_id": str(source_probe_id),
    }
    return TrustedComputeEnvironmentReceiptV1(
        **body, evidence_sha256=canonical_sha256(body)
    )


def consume_pyscf_compute_environment_receipt(
    raw_receipt: Mapping[str, Any], *, engine: str
) -> TrustedComputeEnvironmentReceiptV1:
    """Adapt the hash-bound PySCF target-interpreter probe without re-probing."""

    raw = dict(raw_receipt)
    observed_digest = str(raw.pop("receipt_sha256", ""))
    expected_digest = hashlib.sha256(
        json.dumps(raw, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    if not observed_digest or observed_digest != expected_digest:
        raise ContractError("PySCF environment receipt digest mismatch")
    if raw.get("schema_version") != "chemsmart.pyscf-environment.v1":
        raise ContractError("unsupported PySCF environment receipt schema")
    if raw.get("status") != "available":
        raise ContractError("PySCF environment probe is not available")
    if int(raw.get("probe_returncode", -1)) != 0:
        raise ContractError("PySCF environment probe did not exit cleanly")
    interpreter_sha256 = require_sha256(
        str(raw.get("interpreter_sha256") or ""),
        "interpreter_sha256",
    )
    if not raw.get("interpreter") or not raw.get("interpreter_observed"):
        raise ContractError("PySCF environment lacks interpreter identities")
    try:
        interpreter_matches = (
            Path(str(raw["interpreter"])).resolve()
            == Path(str(raw["interpreter_observed"])).resolve()
        )
    except (OSError, TypeError, ValueError) as exc:
        raise ContractError("PySCF interpreter identity is invalid") from exc
    if not interpreter_matches:
        raise ContractError("PySCF compute interpreter differs from probe")
    dependencies_raw = raw.get("dependencies")
    dependencies_raw = (
        dependencies_raw if isinstance(dependencies_raw, Mapping) else {}
    )
    packages_raw = raw.get("packages")
    packages_raw = packages_raw if isinstance(packages_raw, Mapping) else {}
    dependency_versions: dict[str, str] = {
        (
            "libxc_distribution"
            if str(name).lower() == "libxc"
            else str(name).lower()
        ): str(version)
        for name, version in packages_raw.items()
        if version
    }
    for name, detail in dependencies_raw.items():
        if (
            not isinstance(detail, Mapping)
            or detail.get("available") is not True
        ):
            continue
        dependency_name = str(name).lower()
        if dependency_name == "libxc":
            dependency_name = "libxc_distribution"
        dependency_versions.setdefault(
            dependency_name, str(detail.get("version") or "available")
        )
    observed_pyscf_version = str(raw.get("pyscf_version") or "")
    declared_pyscf_version = str(dependency_versions.get("pyscf") or "")
    if (
        observed_pyscf_version
        and declared_pyscf_version
        and observed_pyscf_version != declared_pyscf_version
    ):
        raise ContractError(
            "PySCF runtime version differs from dependency evidence"
        )
    if observed_pyscf_version:
        dependency_versions["pyscf"] = observed_pyscf_version
    observed_libxc_version = str(raw.get("libxc_version") or "")
    if observed_libxc_version:
        # ``pyscf.dft.libxc.libxc_version()`` describes the runtime parser
        # implementation.  A Python distribution named ``libxc`` is separate
        # package evidence and must never be overwritten or interpreted as the
        # runtime functional authority.
        dependency_versions["libxc_runtime"] = observed_libxc_version
    for logical_name, field in (
        ("gpu4pyscf", "gpu4pyscf_distribution"),
        ("cupy", "cupy_distribution"),
    ):
        detail = raw.get(field)
        if isinstance(detail, Mapping) and detail.get("version"):
            dependency_versions[logical_name] = str(detail["version"])
    solvers_raw = raw.get("solver_callables")
    solvers_raw = solvers_raw if isinstance(solvers_raw, Mapping) else {}
    solver_evidence = {
        str(name).lower(): bool(
            detail.get("callable") if isinstance(detail, Mapping) else detail
        )
        for name, detail in solvers_raw.items()
    }
    cupy_detail = raw.get("cupy_distribution")
    cupy_detail = cupy_detail if isinstance(cupy_detail, Mapping) else {}
    cutensor_detail = raw.get("cutensor_distribution")
    cutensor_detail = (
        cutensor_detail if isinstance(cutensor_detail, Mapping) else {}
    )
    gpu4pyscf_detail = raw.get("gpu4pyscf_distribution")
    gpu4pyscf_detail = (
        gpu4pyscf_detail if isinstance(gpu4pyscf_detail, Mapping) else {}
    )
    gpu_evidence = {
        "cuda_available": int(raw.get("cuda_available", 0) or 0) > 0,
        "cuda_driver_version": int(raw.get("cuda_driver_version", 0) or 0),
        "cuda_runtime_version": int(raw.get("cuda_runtime_version", 0) or 0),
        "cupy_version": str(cupy_detail.get("version") or ""),
        "cupy_distribution_name": str(cupy_detail.get("name") or ""),
        "cutensor_version": str(cutensor_detail.get("version") or ""),
        "cutensor_runtime_version": int(raw.get("cutensor_runtime", 0) or 0),
        "cutensor_compatible": raw.get("cutensor_compatible") is True,
        "cutensor_runtime_available": (
            raw.get("cutensor_runtime_available") is True
        ),
        "device_available": int(raw.get("device_count", 0) or 0) > 0,
        "device_count": int(raw.get("device_count", 0) or 0),
        "gpu_model": str(raw.get("device_name") or ""),
        "gpu_uuid": str(raw.get("device_uuid") or ""),
        "gpu4pyscf_distribution": bool(gpu4pyscf_detail),
        "gpu4pyscf_distribution_name": str(gpu4pyscf_detail.get("name") or ""),
        "gpu4pyscf_distribution_suffix_supported": bool(
            str(gpu4pyscf_detail.get("name") or "") == "gpu4pyscf"
            or str(gpu4pyscf_detail.get("name") or "").startswith(
                "gpu4pyscf-cuda"
            )
        ),
        "gpu4pyscf_version": str(gpu4pyscf_detail.get("version") or ""),
        "pyscf_version": str(dependency_versions.get("pyscf") or ""),
    }
    return build_trusted_compute_environment_receipt(
        program="pyscf",
        engine=engine,
        compute_interpreter_sha256=interpreter_sha256,
        dependency_versions=dependency_versions,
        solver_evidence=solver_evidence,
        gpu_evidence=gpu_evidence,
        source_probe_id=(
            f"{raw.get('schema_version', 'unknown')}:{observed_digest}"
        ),
    )


def pyscf_compute_environment_target(engine: str) -> EnvironmentTargetV1:
    """Declare which facts a PySCF CPU/GPU target receipt must contain."""

    normalized = require_identifier(engine, "engine")
    if normalized not in {"cpu", "gpu"}:
        raise ContractError("PySCF engine must be cpu or gpu")
    dependencies = ("h5py", "numpy", "pyscf")
    versions = (("pyscf", "2.14.0"),)
    gpu_facts: tuple[str, ...] = ()
    gpu_values: tuple[str, ...] = ()
    if normalized == "gpu":
        dependencies = tuple(sorted(dependencies + ("cupy", "gpu4pyscf")))
        versions = (("gpu4pyscf", "1.8.0"), ("pyscf", "2.14.0"))
        gpu_facts = (
            "cuda_available",
            "cutensor_compatible",
            "cutensor_runtime_available",
            "device_available",
            "gpu4pyscf_distribution",
            "gpu4pyscf_distribution_suffix_supported",
        )
        gpu_values = (
            "cuda_driver_version",
            "cuda_runtime_version",
            "cupy_distribution_name",
            "cupy_version",
            "cutensor_runtime_version",
            "cutensor_version",
            "gpu4pyscf_distribution_name",
            "gpu4pyscf_version",
            "gpu_model",
            "gpu_uuid",
            "pyscf_version",
        )
    return EnvironmentTargetV1(
        program="pyscf",
        engine=normalized,
        target_kind="compute_receipt",
        locator="chemsmart.jobs.pyscf.environment.probe_compute_environment",
        distribution="pyscf",
        required_dependencies=dependencies,
        required_dependency_versions=versions,
        required_gpu_facts=gpu_facts,
        required_gpu_values=gpu_values,
    )


def resolve_program_engine_binding(
    capability: CapabilityQueryReceiptV1,
    environment: EnvironmentCapabilityReceiptV1,
) -> ResolvedEngineBindingV1:
    """Compatibility wrapper returning the canonical engine binding."""

    return resolve_engine_binding(
        resolve_program_binding(capability), environment
    )


def resolve_program_binding(
    capability: CapabilityQueryReceiptV1,
    *,
    requested_program: str | None = None,
    substitution_receipt_sha256: str = "",
) -> ResolvedProgramBindingV1:
    if capability.status is CapabilityQueryStatus.SUPPORTED:
        state = "resolved"
        rules = ("binding.program.resolved",)
    elif capability.status is CapabilityQueryStatus.PREVIEW_ONLY:
        state = "preview_only"
        rules = ("binding.program.preview_only",)
    else:
        state = "blocked"
        rules = ("binding.program.capability_red",)
    body = {
        "schema_version": "chemsmart.resolved-program-binding.v1",
        "requested_program": require_identifier(
            requested_program or capability.query.program,
            "requested_program",
        ),
        "selected_program": capability.query.program,
        "capability_receipt_sha256": capability.receipt_sha256,
        "substitution_receipt_sha256": substitution_receipt_sha256,
        "joined_capability_sha256": capability.joined_capability_sha256,
        "state": state,
        "rule_ids": rules,
    }
    return ResolvedProgramBindingV1(
        **body, binding_sha256=canonical_sha256(body)
    )


def resolve_engine_binding(
    program_binding: ResolvedProgramBindingV1,
    environment: EnvironmentCapabilityReceiptV1,
) -> ResolvedEngineBindingV1:
    capability_receipt_sha256 = environment.capability_receipt_sha256
    if not capability_receipt_sha256:
        raise ContractError("environment receipt lacks a capability binding")
    if environment.program != program_binding.program:
        raise ContractError("environment program differs from program binding")
    if capability_receipt_sha256 != program_binding.capability_receipt_sha256:
        raise ContractError(
            "environment capability differs from program binding"
        )
    if (
        environment.query.capability_receipt_sha256
        != capability_receipt_sha256
    ):
        raise ContractError("environment query and receipt bindings differ")
    engine = environment.engine
    if program_binding.state == "blocked":
        state = "blocked"
        ready = False
        rules = ("binding.engine.program_blocked",)
    elif environment.status is EnvironmentStatus.AVAILABLE:
        state = "resolved"
        ready = program_binding.state == "resolved"
        rules = ("binding.engine.environment_observed",)
    elif environment.status is EnvironmentStatus.NOT_DECLARED:
        state = "preview_only"
        ready = False
        rules = ("binding.engine.environment_not_declared",)
    else:
        state = "blocked"
        ready = False
        rules = ("binding.engine.environment_unavailable",)
    body = {
        "schema_version": "chemsmart.resolved-engine-binding.v1",
        "program": program_binding.program,
        "requested_engine": environment.query.engine,
        "selected_engine": engine,
        "program_binding_sha256": program_binding.binding_sha256,
        "capability_receipt_sha256": capability_receipt_sha256,
        "environment_receipt_sha256": environment.receipt_sha256,
        "state": state,
        "execution_ready": ready,
        "rule_ids": rules,
    }
    return ResolvedEngineBindingV1(
        **body, binding_sha256=canonical_sha256(body)
    )


def _environment_receipt(
    capability: CapabilityQueryReceiptV1,
    *,
    engine: str,
    status: EnvironmentStatus,
    target: EnvironmentTargetV1 | None = None,
    observed_version: str = "",
    observed_location_sha256: str = "",
    compute: TrustedComputeEnvironmentReceiptV1 | None = None,
    observation_method: str = "none",
    rule_ids: tuple[str, ...] = (),
) -> EnvironmentCapabilityReceiptV1:
    query = ProgramEnvironmentQueryV1(
        capability_receipt_sha256=capability.receipt_sha256,
        program=capability.query.program,
        engine=engine,
    )
    body = {
        "schema_version": "chemsmart.environment-capability-receipt.v1",
        "query": query,
        "capability_receipt_sha256": capability.receipt_sha256,
        "program": capability.query.program,
        "engine": engine,
        "status": status,
        "target_kind": target.target_kind if target is not None else "",
        "locator": target.locator if target is not None else "",
        "observed_version": observed_version,
        "observed_location_sha256": observed_location_sha256,
        "compute_interpreter_sha256": (
            compute.compute_interpreter_sha256 if compute is not None else ""
        ),
        "compute_evidence_sha256": (
            compute.evidence_sha256 if compute is not None else ""
        ),
        "dependency_versions": (
            compute.dependency_versions if compute is not None else ()
        ),
        "solver_evidence": (
            compute.solver_evidence if compute is not None else ()
        ),
        "gpu_evidence": compute.gpu_evidence if compute is not None else (),
        "observation_method": observation_method,
        "rule_ids": rule_ids,
    }
    return EnvironmentCapabilityReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def _require_engine_job_capabilities(
    values: tuple[EngineJobCapabilityV1, ...],
    *,
    engines: tuple[str, ...],
    jobtypes: tuple[str, ...],
    allow_empty: bool,
) -> None:
    if not values and not allow_empty:
        raise ContractError("engine-job capability matrix must not be empty")
    if values != tuple(sorted(set(values))):
        raise ContractError(
            "engine-job capabilities must be sorted and unique"
        )
    for item in values:
        if not isinstance(item, EngineJobCapabilityV1):
            raise ContractError(
                "engine-job capability matrix contains an invalid row"
            )
        if item.engine not in engines:
            raise ContractError(
                "engine-job capability uses an undeclared engine"
            )
        if item.jobtype not in jobtypes:
            raise ContractError(
                "engine-job capability uses an undeclared jobtype"
            )


def _effective_rule_engine_job_pairs(
    rule: ProgramSupportRuleV1,
) -> tuple[tuple[str, str], ...]:
    if rule.allowed_engine_job_pairs:
        return rule.allowed_engine_job_pairs
    return tuple(
        (engine, jobtype)
        for engine in rule.allowed_engines
        for jobtype in rule.allowed_jobtypes
    )


def _require_engine_job_pairs(
    values: tuple[tuple[str, str], ...],
    *,
    engines: tuple[str, ...],
    jobtypes: tuple[str, ...],
    allow_empty: bool,
) -> None:
    if not values and not allow_empty:
        raise ContractError("engine-job pairs must not be empty")
    normalized = tuple(
        sorted(
            {
                (
                    require_identifier(engine, "engine-job engine"),
                    require_identifier(jobtype, "engine-job jobtype"),
                )
                for engine, jobtype in values
            }
        )
    )
    if values != normalized:
        raise ContractError("engine-job pairs must be sorted and unique")
    for engine, jobtype in values:
        if engine not in engines:
            raise ContractError("engine-job pair uses an undeclared engine")
        if jobtype not in jobtypes:
            raise ContractError("engine-job pair uses an undeclared jobtype")


def _normalized_tuple(values: Any, field_name: str) -> tuple[str, ...]:
    if values is None:
        values = ()
    try:
        normalized = tuple(
            sorted(
                {
                    require_identifier(str(value), field_name)
                    for value in values
                }
            )
        )
    except TypeError as exc:
        raise ContractError(f"{field_name} must be iterable") from exc
    return normalized


def _require_sorted_unique(values: tuple[str, ...], field_name: str) -> None:
    normalized = _normalized_tuple(values, field_name)
    if values != normalized:
        raise ContractError(f"{field_name} must be sorted and unique")


def _sha256_existing_file(path: Path) -> str:
    if not path.is_file():
        return ""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


__all__ = [
    "AgentProgramSupportOverlayV1",
    "CapabilityQueryReceiptV1",
    "CapabilityQueryStatus",
    "CapabilityQueryV1",
    "EnvironmentCapabilityReceiptV1",
    "EnvironmentStatus",
    "EnvironmentTargetV1",
    "EngineJobCapabilityV1",
    "ProgramCandidateProposalV1",
    "ProgramComponentConformanceReceiptV1",
    "ProgramCapabilityQueryV1",
    "ProgramCapabilityReceiptV1",
    "ProgramEnvironmentQueryV1",
    "ProgramEnvironmentReceiptV1",
    "ProgramCapabilityRegistryV1",
    "ProgramCapabilityV1",
    "ProgramSupportOverlayV1",
    "ProgramSupportRuleV1",
    "TrustedComputeEnvironmentReceiptV1",
    "ResolvedEngineBindingV1",
    "ResolvedProgramEngineBindingV1",
    "ResolvedProgramBindingV1",
    "SupportLevel",
    "build_support_overlay",
    "build_approved_execution_overlay",
    "build_command_compiled_preview_overlay",
    "build_program_component_conformance_receipt",
    "build_program_candidate_proposal",
    "build_trusted_compute_environment_receipt",
    "consume_pyscf_compute_environment_receipt",
    "load_program_capabilities",
    "query_capability",
    "query_environment",
    "pyscf_compute_environment_target",
    "resolve_engine_binding",
    "resolve_program_binding",
    "resolve_program_engine_binding",
]


# Approved M2 canonical public names.  Earlier internal names remain aliases
# so a staged merge does not require a parallel runtime migration.
AgentProgramSupportOverlayV1 = ProgramSupportOverlayV1
ProgramCapabilityQueryV1 = CapabilityQueryV1
ProgramCapabilityReceiptV1 = CapabilityQueryReceiptV1
ProgramEnvironmentReceiptV1 = EnvironmentCapabilityReceiptV1
ResolvedProgramEngineBindingV1 = ResolvedEngineBindingV1
