"""Program-specific generated-artifact inspection through typed verifiers."""

from __future__ import annotations

import importlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Mapping

from chemsmart.agent._contracts import (
    ContractError,
    TrustedArtifactRefV1,
    canonical_data,
    canonical_sha256,
    file_sha256,
    require_identifier,
)


@dataclass(frozen=True)
class ValidationFindingV1:
    rule_id: str
    field: str
    expected: Any
    observed: Any
    evidence_ref: str


@dataclass(frozen=True)
class GeneratedArtifactInspectionReceiptV1:
    schema_version: str
    program: str
    artifact_id: str
    artifact_sha256: str
    project_artifact_id: str
    project_sha256: str
    expected_receipt_sha256: str
    verifier_id: str
    findings: tuple[ValidationFindingV1, ...]
    status: str
    rule_ids: tuple[str, ...]
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.generated-inspection-receipt.v1":
            raise ContractError("unsupported generated inspection schema")
        if self.status not in {"valid", "invalid", "verifier_unavailable"}:
            raise ContractError("invalid generated inspection state")
        expected = canonical_sha256(
            {
                "schema_version": self.schema_version,
                "program": self.program,
                "artifact_id": self.artifact_id,
                "artifact_sha256": self.artifact_sha256,
                "project_artifact_id": self.project_artifact_id,
                "project_sha256": self.project_sha256,
                "expected_receipt_sha256": self.expected_receipt_sha256,
                "verifier_id": self.verifier_id,
                "findings": self.findings,
                "status": self.status,
                "rule_ids": self.rule_ids,
            }
        )
        if self.receipt_sha256 != expected:
            raise ContractError("generated inspection receipt digest mismatch")


def inspect_generated_artifact(
    *,
    program: str,
    settings: Any,
    artifact: TrustedArtifactRefV1,
    project_artifact: TrustedArtifactRefV1 | None,
    expected_receipt: Mapping[str, Any] | None,
    verifier: Callable[[Any, str | Path], list[Any]] | None = None,
) -> GeneratedArtifactInspectionReceiptV1:
    """Inspect a native/structured artifact without engine-specific guessing.

    For PySCF, the conventional resolver selects
    ``chemsmart.jobs.pyscf.validation.verify_provenance`` and therefore reads
    the structured HDF5 result rather than the human log.  Other programs may
    expose the same verifier name without changing this agent module.
    """

    normalized_program = require_identifier(program, "program")
    _verify_binding(artifact)
    if project_artifact is not None:
        _verify_binding(project_artifact)
    artifact_sha256 = artifact.sha256
    expected = dict(expected_receipt or {})
    expected_receipt_sha256 = canonical_sha256(expected) if expected else ""
    if normalized_program == "pyscf":
        _validate_pyscf_expected_receipt(expected, project_artifact)
    verifier_id = ""
    findings: tuple[ValidationFindingV1, ...] = ()
    status = "verifier_unavailable"
    rules = ("inspection.verifier.unavailable",)
    try:
        selected = verifier or _resolve_verifier(normalized_program)
        verifier_id = f"{selected.__module__}.{selected.__qualname__}"
        if normalized_program == "pyscf":
            raw_findings = selected(
                settings,
                artifact.path,
                expected_receipt=expected,
            )
        else:
            raw_findings = selected(settings, artifact.path)
        findings = tuple(
            ValidationFindingV1(
                rule_id=str(item.rule_id),
                field=str(item.field),
                expected=canonical_data(item.expected),
                observed=canonical_data(item.observed),
                evidence_ref=str(item.evidence_ref),
            )
            for item in raw_findings
        )
        status = "valid" if not findings else "invalid"
        rules = (
            (
                "inspection.provenance.valid"
                if not findings
                else "inspection.provenance.findings"
            ),
        )
    except (ImportError, AttributeError):
        pass
    body = {
        "schema_version": "chemsmart.generated-inspection-receipt.v1",
        "program": normalized_program,
        "artifact_id": artifact.artifact_id,
        "artifact_sha256": artifact_sha256,
        "project_artifact_id": (
            project_artifact.artifact_id
            if project_artifact is not None
            else ""
        ),
        "project_sha256": (
            project_artifact.sha256 if project_artifact is not None else ""
        ),
        "expected_receipt_sha256": expected_receipt_sha256,
        "verifier_id": verifier_id,
        "findings": findings,
        "status": status,
        "rule_ids": rules,
    }
    return GeneratedArtifactInspectionReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def _resolve_verifier(program: str) -> Callable[[Any, str | Path], list[Any]]:
    module = importlib.import_module(f"chemsmart.jobs.{program}.validation")
    verifier = getattr(module, "verify_provenance")
    if not callable(verifier):
        raise AttributeError("verify_provenance is not callable")
    return verifier


def _verify_binding(binding: TrustedArtifactRefV1) -> None:
    path = Path(binding.path)
    if not path.is_file() or path.stat().st_size != binding.size_bytes:
        raise ContractError("generated artifact size differs from binding")
    if file_sha256(path) != binding.sha256:
        raise ContractError("generated artifact digest differs from binding")


def _validate_pyscf_expected_receipt(
    expected: Mapping[str, Any],
    project_artifact: TrustedArtifactRefV1 | None,
) -> None:
    required = {
        "script_sha256",
        "input_receipt_sha256",
        "environment_receipt_sha256",
        "input_geometry_sha256",
        "requested_settings_sha256",
        "project_yaml_digest",
    }
    missing = sorted(required.difference(expected))
    if missing:
        raise ContractError(
            "PySCF result verification lacks bound run receipt fields: "
            + ", ".join(missing)
        )
    if project_artifact is None:
        raise ContractError(
            "PySCF result verification requires project binding"
        )
    if expected.get("project_yaml_digest") != project_artifact.sha256:
        raise ContractError("PySCF expected receipt uses another project")


__all__ = [
    "GeneratedArtifactInspectionReceiptV1",
    "ValidationFindingV1",
    "inspect_generated_artifact",
]
