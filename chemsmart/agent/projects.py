"""Read, render, and validate project-YAML candidates without executing jobs."""

from __future__ import annotations

import importlib
import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import yaml

from chemsmart.agent._contracts import (
    ContractError,
    TrustedArtifactRefV1,
    canonical_data,
    canonical_sha256,
    file_sha256,
    require_identifier,
)
from chemsmart.agent.capabilities import (
    CapabilityQueryReceiptV1,
    CapabilityQueryStatus,
    ProgramCapabilityRegistryV1,
    load_program_capabilities,
)


@dataclass(frozen=True)
class ProjectSectionV1:
    name: str
    settings: tuple[tuple[str, Any], ...]

    def __post_init__(self) -> None:
        require_identifier(self.name, "project section")
        names = tuple(item[0] for item in self.settings)
        if names != tuple(sorted(set(names))):
            raise ContractError("project setting names must be sorted and unique")
        canonical_data(dict(self.settings))


@dataclass(frozen=True)
class ProjectDocumentV1:
    schema_version: str
    program: str
    sections: tuple[ProjectSectionV1, ...]
    content_sha256: str
    document_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.project-document.v1":
            raise ContractError("unsupported project document schema")
        require_identifier(self.program, "program")
        names = tuple(item.name for item in self.sections)
        if names != tuple(sorted(set(names))):
            raise ContractError("project sections must be sorted and unique")
        expected = canonical_sha256(
            {
                "schema_version": self.schema_version,
                "program": self.program,
                "sections": self.sections,
                "content_sha256": self.content_sha256,
            }
        )
        if self.document_sha256 != expected:
            raise ContractError("project document digest mismatch")


@dataclass(frozen=True)
class ProjectRenderReceiptV1:
    schema_version: str
    program: str
    registry_sha256: str
    document_sha256: str
    rendered_yaml: str
    rendered_sha256: str
    status: str
    rule_ids: tuple[str, ...]
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.project-render-receipt.v1":
            raise ContractError("unsupported project render receipt schema")
        if self.status != "candidate_rendered":
            raise ContractError("rendering alone can only create a candidate")
        expected = canonical_sha256(
            {
                "schema_version": self.schema_version,
                "program": self.program,
                "registry_sha256": self.registry_sha256,
                "document_sha256": self.document_sha256,
                "rendered_yaml": self.rendered_yaml,
                "rendered_sha256": self.rendered_sha256,
                "status": self.status,
                "rule_ids": self.rule_ids,
            }
        )
        if self.receipt_sha256 != expected:
            raise ContractError("project render receipt digest mismatch")


@dataclass(frozen=True)
class ProjectValidationReceiptV1:
    schema_version: str
    project_artifact_id: str
    project_sha256: str
    capability_receipt_sha256: str
    program: str
    jobtype: str
    loader_id: str
    settings_sha256: str
    settings: tuple[tuple[str, Any], ...]
    status: str
    error_class: str
    diagnostic: str
    rule_ids: tuple[str, ...]
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.project-validation-receipt.v1":
            raise ContractError("unsupported project validation receipt schema")
        if self.status not in {"valid", "invalid", "loader_unavailable"}:
            raise ContractError("invalid project validation status")
        expected = canonical_sha256(
            {
                "schema_version": self.schema_version,
                "project_artifact_id": self.project_artifact_id,
                "project_sha256": self.project_sha256,
                "capability_receipt_sha256": self.capability_receipt_sha256,
                "program": self.program,
                "jobtype": self.jobtype,
                "loader_id": self.loader_id,
                "settings_sha256": self.settings_sha256,
                "settings": self.settings,
                "status": self.status,
                "error_class": self.error_class,
                "diagnostic": self.diagnostic,
                "rule_ids": self.rule_ids,
            }
        )
        if self.receipt_sha256 != expected:
            raise ContractError("project validation receipt digest mismatch")


@dataclass(frozen=True)
class PySCFFunctionalResolutionReceiptV1:
    """Host-side XC alias resolution, distinct from target LibXC parsing."""

    schema_version: str
    project_validation_receipt_sha256: str
    project_sha256: str
    jobtype: str
    setting_path: str
    requested_method_kind: str
    requested_literal: str | None
    normalized_requested_literal: str
    applied_xc: str | None
    normalized_applied_xc: str
    status: str
    functional_family: str
    correlation_convention: str
    source: str
    rule_id: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.pyscf-functional-resolution.v1":
            raise ContractError("unsupported PySCF functional resolution schema")
        if self.requested_method_kind not in {"hf", "dft"}:
            raise ContractError("invalid PySCF functional method kind")
        if self.status not in {
            "not_applicable",
            "registered_alias",
            "explicit_variant",
            "literal_preserved",
        }:
            raise ContractError("invalid PySCF functional resolution status")
        body = {
            key: value
            for key, value in self.__dict__.items()
            if key != "receipt_sha256"
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("PySCF functional resolution digest mismatch")

    @property
    def evidence_ref(self) -> str:
        return f"functional_resolution:{self.receipt_sha256}"

    def public_record(self) -> dict[str, Any]:
        return {
            **canonical_data(self),
            "evidence_ref": self.evidence_ref,
            "evidence_scope": "host_xc_resolution_only",
        }


def project_document(
    *, program: str, sections: Mapping[str, Mapping[str, Any]]
) -> ProjectDocumentV1:
    """Create a typed project candidate; this is not loader validation."""

    normalized_sections = tuple(
        ProjectSectionV1(
            name=require_identifier(str(section), "project section"),
            settings=tuple(
                sorted(
                    (
                        require_identifier(str(key), "project setting"),
                        canonical_data(value),
                    )
                    for key, value in settings.items()
                )
            ),
        )
        for section, settings in sorted(sections.items())
    )
    semantic_body = {
        section.name: dict(section.settings) for section in normalized_sections
    }
    content_sha256 = canonical_sha256(semantic_body)
    body = {
        "schema_version": "chemsmart.project-document.v1",
        "program": require_identifier(program, "program"),
        "sections": normalized_sections,
        "content_sha256": content_sha256,
    }
    return ProjectDocumentV1(
        **body, document_sha256=canonical_sha256(body)
    )


def read_project_yaml(
    binding: TrustedArtifactRefV1, *, program: str
) -> ProjectDocumentV1:
    """Read a trusted project artifact and verify its exact-byte identity."""

    path = Path(binding.path)
    if path.stat().st_size != binding.size_bytes:
        raise ContractError("project artifact size differs from its binding")
    if file_sha256(path) != binding.sha256:
        raise ContractError("project artifact digest differs from its binding")
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, Mapping) or not payload:
        raise ContractError("project YAML root must be a non-empty mapping")
    sections = {}
    for section, settings in payload.items():
        if not isinstance(settings, Mapping):
            raise ContractError("every project YAML section must be a mapping")
        sections[str(section)] = dict(settings)
    return project_document(program=program, sections=sections)


def _require_declared_section_shape(
    program: str, payload: Mapping[str, Any]
) -> None:
    """Refuse a project whose top-level sections the program cannot read.

    The section vocabulary differs by program: the route-building programs
    group settings by phase (``gas``/``solv``) while PySCF and xTB key them by
    job type.  A document with the wrong shape used to render and promote
    cleanly, then fail deep inside the loader as ``'NoneType' object has no
    attribute 'items'`` -- a message that names neither the file nor the
    mistake.  Refusing here states the expected sections while the author can
    still act on it.
    """

    from chemsmart.settings.capabilities import PROGRAM_CAPABILITIES

    capability = PROGRAM_CAPABILITIES.get(program)
    declared = set(capability.project_section_names if capability else ())
    if not declared:
        return
    unknown = sorted(set(payload) - declared)
    if unknown:
        raise ContractError(
            f"{program} project YAML has no section(s) "
            f"{unknown}; it accepts {sorted(declared)}"
        )
    if not payload:
        raise ContractError(
            f"{program} project YAML declares no section; it accepts "
            f"{sorted(declared)}"
        )


def render_project_yaml(
    document: ProjectDocumentV1,
    *,
    registry: ProgramCapabilityRegistryV1 | None = None,
) -> ProjectRenderReceiptV1:
    """Render a deterministic candidate view; the program loader remains final."""

    registry = registry or load_program_capabilities()
    capability = registry.get(document.program)
    if capability is None or not capability.supports_project_configuration:
        raise ContractError("program does not declare project-YAML support")
    payload = {
        section.name: dict(section.settings) for section in document.sections
    }
    _require_declared_section_shape(document.program, payload)
    rule_ids = ["project.render.candidate_only"]
    if document.program == "pyscf":
        settings_module = importlib.import_module("chemsmart.settings.pyscf")
        materializer = getattr(
            settings_module, "materialize_canonical_project_sections"
        )
        materialized = materializer(payload)
        if not isinstance(materialized, Mapping) or not materialized:
            raise ContractError(
                "project materializer must return a non-empty mapping"
            )
        payload = {
            str(section): dict(settings)
            for section, settings in materialized.items()
        }
        rule_ids.append("project.render.host_effective_settings")
    rendered = yaml.safe_dump(
        payload,
        sort_keys=True,
        default_flow_style=False,
        allow_unicode=True,
    )
    rendered_sha256 = hashlib.sha256(rendered.encode("utf-8")).hexdigest()
    body = {
        "schema_version": "chemsmart.project-render-receipt.v1",
        "program": document.program,
        "registry_sha256": registry.registry_sha256,
        "document_sha256": document.document_sha256,
        "rendered_yaml": rendered,
        "rendered_sha256": rendered_sha256,
        "status": "candidate_rendered",
        "rule_ids": tuple(sorted(rule_ids)),
    }
    return ProjectRenderReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def validate_project_yaml(
    binding: TrustedArtifactRefV1,
    *,
    capability: CapabilityQueryReceiptV1,
) -> ProjectValidationReceiptV1:
    """Validate through the checked-out program loader without running a job."""

    program = capability.query.program
    jobtype = capability.query.jobtype
    loader_id = ""
    settings_rows: tuple[tuple[str, Any], ...] = ()
    settings_sha256 = ""
    status = "invalid"
    error_class = ""
    diagnostic = ""
    rules: tuple[str, ...]
    if capability.status not in {
        CapabilityQueryStatus.SUPPORTED,
        CapabilityQueryStatus.PREVIEW_ONLY,
    }:
        rules = ("project.capability.not_supported",)
    elif not jobtype:
        rules = ("project.jobtype.required",)
    elif Path(binding.path).stat().st_size != binding.size_bytes:
        rules = ("project.artifact.size_mismatch",)
    elif file_sha256(binding.path) != binding.sha256:
        rules = ("project.artifact.digest_mismatch",)
    else:
        try:
            module = importlib.import_module(f"chemsmart.settings.{program}")
            loader = _yaml_project_loader(module)
            loader_id = f"{loader.__module__}.{loader.__qualname__}"
            project = loader.from_yaml(binding.path)
            method = getattr(project, f"{jobtype}_settings", None)
            if not callable(method):
                raise ContractError(
                    "project loader has no settings method for this jobtype"
                )
            settings = method()
            validator = getattr(settings, "validate", None)
            if callable(validator):
                validator()
            allowed = set(
                capability.capability.project_owned_parameters
                if capability.capability is not None
                else ()
            )
            allowed.update({"jobtype", "engine", "freq"})
            values = {
                key: canonical_data(getattr(settings, key))
                for key in sorted(allowed)
                if hasattr(settings, key)
            }
            settings_rows = tuple(values.items())
            settings_sha256 = canonical_sha256(values)
            status = "valid"
            rules = ("project.loader.valid",)
        except (ImportError, AttributeError) as exc:
            status = "loader_unavailable"
            error_class = type(exc).__name__
            diagnostic = _public_loader_diagnostic(exc, binding.path)
            rules = ("project.loader.unavailable",)
        except Exception as exc:
            error_class = type(exc).__name__
            diagnostic = _public_loader_diagnostic(exc, binding.path)
            rules = ("project.loader.rejected",)
    body = {
        "schema_version": "chemsmart.project-validation-receipt.v1",
        "project_artifact_id": binding.artifact_id,
        "project_sha256": binding.sha256,
        "capability_receipt_sha256": capability.receipt_sha256,
        "program": program,
        "jobtype": jobtype,
        "loader_id": loader_id,
        "settings_sha256": settings_sha256,
        "settings": settings_rows,
        "status": status,
        "error_class": error_class,
        "diagnostic": diagnostic,
        "rule_ids": rules,
    }
    return ProjectValidationReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def project_scientific_materializations(
    receipt: ProjectValidationReceiptV1,
) -> tuple[PySCFFunctionalResolutionReceiptV1, ...]:
    """Derive preview-safe scientific resolution evidence from one loader receipt."""

    if receipt.status != "valid" or receipt.program != "pyscf":
        return ()
    from chemsmart.jobs.pyscf.settings import describe_functional_resolution

    values = dict(receipt.settings)
    resolution = describe_functional_resolution(
        values.get("functional"), ab_initio=values.get("ab_initio")
    )
    if resolution["status"] == "missing":
        return ()
    body = {
        "schema_version": resolution["schema_version"],
        "project_validation_receipt_sha256": receipt.receipt_sha256,
        "project_sha256": receipt.project_sha256,
        "jobtype": receipt.jobtype,
        "setting_path": f"{receipt.jobtype}.functional",
        "requested_method_kind": resolution["requested_method_kind"],
        "requested_literal": resolution["requested_literal"],
        "normalized_requested_literal": resolution[
            "normalized_requested_literal"
        ],
        "applied_xc": resolution["applied_xc"],
        "normalized_applied_xc": resolution["normalized_applied_xc"],
        "status": resolution["status"],
        "functional_family": resolution["functional_family"],
        "correlation_convention": resolution["correlation_convention"],
        "source": resolution["source"],
        "rule_id": resolution["rule_id"],
    }
    return (
        PySCFFunctionalResolutionReceiptV1(
            **body, receipt_sha256=canonical_sha256(body)
        ),
    )


def _yaml_project_loader(module: Any) -> type:
    candidates = []
    for name, value in vars(module).items():
        if (
            isinstance(value, type)
            and name.startswith("Yaml")
            and name.endswith("ProjectSettings")
            and callable(getattr(value, "from_yaml", None))
        ):
            candidates.append(value)
    if len(candidates) != 1:
        raise ContractError("program settings module has no unique YAML loader")
    return candidates[0]


def _public_loader_diagnostic(exc: Exception, project_path: str) -> str:
    """Return a bounded, path-free counterexample for model repair."""

    message = " ".join(str(exc).split())
    if project_path:
        message = message.replace(str(project_path), "<project>")
    return message[:500]


__all__ = [
    "ProjectDocumentV1",
    "ProjectRenderReceiptV1",
    "ProjectSectionV1",
    "ProjectValidationReceiptV1",
    "PySCFFunctionalResolutionReceiptV1",
    "project_document",
    "project_scientific_materializations",
    "read_project_yaml",
    "render_project_yaml",
    "validate_project_yaml",
]
