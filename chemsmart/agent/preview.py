"""Isolated, exact-argv ChemSmart safe preview execution.

This module invokes the checked-out Click tree only with compiler-proven
``--fake``/``--no-scratch`` and, for submissions, ``--test``.  It hashes the
files actually emitted inside an isolated directory and re-hashes every bound
external input before and after invocation.  A caller cannot supply an exit
status or claim an artifact hash.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import click
from click.testing import CliRunner

from chemsmart.agent._contracts import (
    AuxiliaryArtifactBindingV1,
    ContractError,
    TrustedArtifactRefV1,
    canonical_sha256,
    file_sha256,
    require_auxiliary_artifact_bindings,
)
from chemsmart.agent.commands import CanonicalCommandInvocationV1
from chemsmart.agent.program_verifiers import (
    ProgramPreviewExpectationV1,
    validate_preview_workspace,
)


@dataclass(frozen=True)
class PreviewArtifactV1:
    relative_path: str
    size_bytes: int
    sha256: str


@dataclass(frozen=True)
class SafePreviewReceiptV1:
    schema_version: str
    invocation_sha256: str
    observed_argv_sha256: str
    input_sha256: str
    project_sha256: str
    exit_status: int
    fake_mode: bool
    no_scratch_mode: bool
    scheduler_test_mode: bool
    artifacts: tuple[PreviewArtifactV1, ...]
    artifact_set_sha256: str
    output_sha256: str
    exception_class: str
    program_validation_receipt_sha256: str
    program_validation_status: str
    critical_finding_sha256s: tuple[str, ...]
    status: str
    rule_ids: tuple[str, ...]
    receipt_sha256: str
    #: The findings those digests identify.  A digest cannot be acted on: a
    #: live session was told six times that its command was invalid, given
    #: only hashes, and recompiled blindly until it abandoned two programs.
    #: The validator computes these bodies in order to hash them, so keeping
    #: them costs nothing and is the difference between a refusal and a
    #: repair.  Excluded from the digest when empty so receipts recorded
    #: before this field keep their identity.
    critical_findings: tuple[Any, ...] = ()
    auxiliary_input_bindings: tuple[AuxiliaryArtifactBindingV1, ...] = ()

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.safe-preview-receipt.v1":
            raise ContractError("unsupported safe preview receipt schema")
        if self.status not in {"previewed", "failed"}:
            raise ContractError("invalid safe preview state")
        if self.status == "previewed" and (
            self.exit_status != 0
            or not self.fake_mode
            or not self.no_scratch_mode
            or not self.artifacts
            or self.program_validation_status != "valid"
            or not self.program_validation_receipt_sha256
            or self.critical_finding_sha256s
        ):
            raise ContractError(
                "previewed requires fake, no-scratch, emitted-artifact evidence"
            )
        require_auxiliary_artifact_bindings(self.auxiliary_input_bindings)
        body = {
            "schema_version": self.schema_version,
            "invocation_sha256": self.invocation_sha256,
            "observed_argv_sha256": self.observed_argv_sha256,
            "input_sha256": self.input_sha256,
            "project_sha256": self.project_sha256,
            "exit_status": self.exit_status,
            "fake_mode": self.fake_mode,
            "no_scratch_mode": self.no_scratch_mode,
            "scheduler_test_mode": self.scheduler_test_mode,
            "artifacts": self.artifacts,
            "artifact_set_sha256": self.artifact_set_sha256,
            "output_sha256": self.output_sha256,
            "exception_class": self.exception_class,
            "program_validation_receipt_sha256": (
                self.program_validation_receipt_sha256
            ),
            "program_validation_status": self.program_validation_status,
            "critical_finding_sha256s": self.critical_finding_sha256s,
            "status": self.status,
            "rule_ids": self.rule_ids,
        }
        if self.auxiliary_input_bindings:
            body["auxiliary_input_bindings"] = self.auxiliary_input_bindings
        expected = canonical_sha256(body)
        if self.receipt_sha256 != expected:
            raise ContractError("safe preview receipt digest mismatch")


def execute_safe_preview(
    invocation: CanonicalCommandInvocationV1,
    *,
    input_artifact: TrustedArtifactRefV1,
    project_artifact: TrustedArtifactRefV1 | None,
    expectation: ProgramPreviewExpectationV1,
    auxiliary_input_artifacts: (
        Mapping[str, TrustedArtifactRefV1] | None
    ) = None,
    root: click.Command | None = None,
    runner: CliRunner | None = None,
) -> SafePreviewReceiptV1:
    """Execute one compiler-owned safe preview and observe actual artifacts."""

    auxiliary_input_artifacts = dict(auxiliary_input_artifacts or {})
    auxiliary_bindings = _auxiliary_input_bindings(auxiliary_input_artifacts)
    _validate_preview_bindings(
        invocation,
        input_artifact=input_artifact,
        project_artifact=project_artifact,
        expectation=expectation,
        auxiliary_input_bindings=auxiliary_bindings,
    )
    fake = _option_enabled(invocation, "fake", "--fake")
    no_scratch = _option_enabled(invocation, "scratch", "--no-scratch")
    scheduler_test = _option_enabled(invocation, "test", "--test")
    rules = []
    if not fake:
        rules.append("preview.fake_flag_missing")
    if not no_scratch:
        rules.append("preview.no_scratch_flag_missing")
    if invocation.command_path[0] == "sub" and not scheduler_test:
        rules.append("preview.scheduler_test_flag_missing")
    if invocation.command_path[0] == "run" and scheduler_test:
        rules.append("preview.scheduler_test_flag_unexpected")
    if rules:
        raise ContractError(", ".join(sorted(rules)))

    if root is None:
        from chemsmart.cli.main import entry_point

        root = entry_point
    runner = runner or CliRunner()
    before_input = _verified_artifact_hash(input_artifact)
    before_project = (
        _verified_artifact_hash(project_artifact)
        if project_artifact is not None
        else ""
    )
    before_auxiliary = tuple(
        _verified_artifact_hash(artifact)
        for _name, artifact in sorted(auxiliary_input_artifacts.items())
    )
    artifacts: tuple[PreviewArtifactV1, ...] = ()
    exit_status = -1
    output_sha256 = canonical_sha256({"output": ""})
    exception_class = ""
    program_validation_receipt_sha256 = ""
    program_validation_status = "invalid"
    critical_finding_sha256s: tuple[str, ...] = ()
    critical_findings: tuple[Any, ...] = ()
    with runner.isolated_filesystem() as workspace:
        result = runner.invoke(
            root,
            list(invocation.argv[1:]),
            catch_exceptions=True,
        )
        exit_status = int(result.exit_code)
        output_sha256 = canonical_sha256({"output": result.output})
        if result.exception is not None:
            exception_class = type(result.exception).__name__
        artifacts = _collect_preview_artifacts(Path(workspace))
        program_validation = validate_preview_workspace(
            expectation, Path(workspace)
        )
        program_validation_receipt_sha256 = program_validation.receipt_sha256
        program_validation_status = program_validation.status
        critical_finding_sha256s = tuple(
            sorted(
                canonical_sha256(item) for item in program_validation.findings
            )
        )
        # Keep what the digests identify.  Hashing a finding and discarding
        # its body is what left a live session recompiling against opaque
        # hashes; the digests above still bind these bodies for integrity.
        critical_findings = tuple(program_validation.findings)
        rules.extend(
            _validate_staged_auxiliary_inputs(
                Path(workspace), auxiliary_input_artifacts
            )
        )

    after_input = _verified_artifact_hash(input_artifact)
    after_project = (
        _verified_artifact_hash(project_artifact)
        if project_artifact is not None
        else ""
    )
    after_auxiliary = tuple(
        _verified_artifact_hash(artifact)
        for _name, artifact in sorted(auxiliary_input_artifacts.items())
    )
    if before_input != after_input:
        rules.append("preview.input_artifact_mutated")
    if before_project != after_project:
        rules.append("preview.project_artifact_mutated")
    if before_auxiliary != after_auxiliary:
        rules.append("preview.auxiliary_input_artifact_mutated")
    if exit_status != 0:
        rules.append("preview.click_invocation_failed")
    if not artifacts:
        rules.append("preview.generated_artifact_missing")
    if program_validation_status != "valid":
        rules.append("preview.program_validator.red")
    status = "previewed" if not rules else "failed"
    if status == "previewed":
        rules.append("preview.click_exact_argv_observed")
    artifact_set_sha256 = canonical_sha256(artifacts)
    body = {
        "schema_version": "chemsmart.safe-preview-receipt.v1",
        "invocation_sha256": invocation.invocation_sha256,
        "observed_argv_sha256": canonical_sha256(invocation.argv),
        "input_sha256": after_input,
        "project_sha256": after_project,
        "exit_status": exit_status,
        "fake_mode": fake,
        "no_scratch_mode": no_scratch,
        "scheduler_test_mode": scheduler_test,
        "artifacts": artifacts,
        "artifact_set_sha256": artifact_set_sha256,
        "output_sha256": output_sha256,
        "exception_class": exception_class,
        "program_validation_receipt_sha256": (
            program_validation_receipt_sha256
        ),
        "program_validation_status": program_validation_status,
        "critical_finding_sha256s": critical_finding_sha256s,
        "status": status,
        "rule_ids": tuple(sorted(set(rules))),
    }
    if auxiliary_bindings:
        body["auxiliary_input_bindings"] = auxiliary_bindings
    return SafePreviewReceiptV1(
        **body,
        receipt_sha256=canonical_sha256(body),
        critical_findings=critical_findings,
    )


def _validate_preview_bindings(
    invocation: CanonicalCommandInvocationV1,
    *,
    input_artifact: TrustedArtifactRefV1,
    project_artifact: TrustedArtifactRefV1 | None,
    expectation: ProgramPreviewExpectationV1,
    auxiliary_input_bindings: tuple[AuxiliaryArtifactBindingV1, ...],
) -> None:
    if invocation.status != "compiled":
        raise ContractError("only a compiled invocation can be previewed")
    if input_artifact.sha256 != invocation.input_sha256:
        raise ContractError("preview input binding differs from invocation")
    observed_project = (
        project_artifact.sha256 if project_artifact is not None else ""
    )
    if observed_project != invocation.project_sha256:
        raise ContractError("preview project binding differs from invocation")
    if not invocation.argv or invocation.argv[0] != "chemsmart":
        raise ContractError(
            "preview argv must be compiler-owned ChemSmart argv"
        )
    if (expectation.program, expectation.jobtype) != (
        invocation.command_path[1],
        invocation.command_path[2],
    ):
        raise ContractError("preview expectation differs from command path")
    if expectation.input_artifact.sha256 != invocation.input_sha256:
        raise ContractError("preview expectation differs from input binding")
    if expectation.project_sha256 != invocation.project_sha256:
        raise ContractError("preview expectation differs from project binding")
    if auxiliary_input_bindings != invocation.auxiliary_input_bindings:
        raise ContractError(
            "preview auxiliary inputs differ from compiled invocation"
        )


def _auxiliary_input_bindings(
    artifacts: Mapping[str, TrustedArtifactRefV1],
) -> tuple[AuxiliaryArtifactBindingV1, ...]:
    return tuple(
        AuxiliaryArtifactBindingV1(
            parameter_name=parameter_name,
            artifact_id=artifact.artifact_id,
            artifact_sha256=artifact.sha256,
        )
        for parameter_name, artifact in sorted(artifacts.items())
    )


def _validate_staged_auxiliary_inputs(
    workspace: Path,
    artifacts: Mapping[str, TrustedArtifactRefV1],
) -> tuple[str, ...]:
    """Confirm the generated job staged the exact bound auxiliary files."""

    findings = []
    for parameter_name, artifact in sorted(artifacts.items()):
        basename = Path(artifact.path).name
        candidates = tuple(
            path
            for path in workspace.rglob(basename)
            if path.is_file() and not path.is_symlink()
        )
        if not any(
            path.stat().st_size == artifact.size_bytes
            and file_sha256(path) == artifact.sha256
            for path in candidates
        ):
            findings.append(
                f"preview.auxiliary_input_not_staged.{parameter_name}"
            )
    return tuple(findings)


def _verified_artifact_hash(binding: TrustedArtifactRefV1) -> str:
    path = Path(binding.path)
    if not path.is_file() or path.stat().st_size != binding.size_bytes:
        raise ContractError("trusted preview artifact size changed")
    digest = file_sha256(path)
    if digest != binding.sha256:
        raise ContractError("trusted preview artifact digest changed")
    return digest


def _option_enabled(
    invocation: CanonicalCommandInvocationV1,
    parameter_name: str,
    expected_flag: str,
) -> bool:
    return any(
        item.parameter_name == parameter_name
        and item.flag == expected_flag
        and not item.values
        for item in invocation.scoped_options
    )


def _collect_preview_artifacts(
    workspace: Path,
) -> tuple[PreviewArtifactV1, ...]:
    rows = []
    for path in sorted(workspace.rglob("*")):
        if path.is_symlink():
            raise ContractError("safe preview emitted a symbolic link")
        if not path.is_file():
            continue
        relative = path.relative_to(workspace).as_posix()
        rows.append(
            PreviewArtifactV1(
                relative_path=relative,
                size_bytes=path.stat().st_size,
                sha256=file_sha256(path),
            )
        )
    return tuple(rows)


__all__ = [
    "PreviewArtifactV1",
    "SafePreviewReceiptV1",
    "execute_safe_preview",
]
