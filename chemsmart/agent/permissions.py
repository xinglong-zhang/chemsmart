"""Hash-bound permission decisions for agent actions."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from threading import Lock

from chemsmart.agent._contracts import (
    ContractError,
    canonical_sha256,
    require_sha256,
)


class PermissionDecision(str, Enum):
    AUTO_ALLOW = "auto_allow"
    NEEDS_USER = "needs_user"
    ALLOW_ONCE = "allow_once"
    DENY = "deny"


_READ_ONLY_ACTIONS = frozenset(
    {
        "capability_query",
        "environment_query",
        "project_read",
        "project_render",
        "project_validate",
        "command_compile",
        "command_inspect",
        "artifact_inspect",
        "result_extract",
        "thermochemistry_derive",
        "quantity_evaluate",
        "fixture_safe_preview",
    }
)

_MATERIAL_ACTIONS = frozenset(
    {
        "project_write",
        "execute_local",
        "submit_hpc",
        "cancel",
        "retry",
        "paid_network",
        "publish",
        "overwrite",
    }
)


@dataclass(frozen=True)
class PermissionRequestV1:
    schema_version: str
    request_id: str
    action: str
    scope_sha256: str
    command_sha256: str = ""
    input_sha256s: tuple[str, ...] = ()
    project_sha256: str = ""
    environment_sha256: str = ""
    provider: str = ""
    quota_scope: str = ""
    request_sha256: str = ""

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.permission-request.v1":
            raise ContractError("unsupported permission request schema")
        if self.action not in _READ_ONLY_ACTIONS | _MATERIAL_ACTIONS:
            raise ContractError("unknown permission action")
        if not self.scope_sha256:
            raise ContractError("permission scope binding must not be empty")
        require_sha256(self.scope_sha256, "permission scope")
        for field_name, digest in (
            ("command_sha256", self.command_sha256),
            ("project_sha256", self.project_sha256),
            ("environment_sha256", self.environment_sha256),
        ):
            if digest:
                require_sha256(digest, field_name)
        for digest in self.input_sha256s:
            require_sha256(digest, "input_sha256")
        if self.action in {"execute_local", "submit_hpc", "retry"}:
            required = (
                self.command_sha256,
                self.input_sha256s,
                self.project_sha256,
                self.environment_sha256,
            )
            if not all(required):
                raise ContractError(
                    "calculation permission requires command, input, project, "
                    "and environment bindings"
                )
        if self.action == "paid_network" and not (
            self.provider and self.quota_scope
        ):
            raise ContractError(
                "paid-network permission requires provider and quota scope"
            )
        body = {
            "schema_version": self.schema_version,
            "request_id": self.request_id,
            "action": self.action,
            "scope_sha256": self.scope_sha256,
            "command_sha256": self.command_sha256,
            "input_sha256s": self.input_sha256s,
            "project_sha256": self.project_sha256,
            "environment_sha256": self.environment_sha256,
            "provider": self.provider,
            "quota_scope": self.quota_scope,
        }
        if self.request_sha256 != canonical_sha256(body):
            raise ContractError("permission request digest mismatch")


@dataclass(frozen=True)
class ApprovalResolutionV1:
    schema_version: str
    approval_id: str
    permission_request_sha256: str
    decision: PermissionDecision
    actor: str
    one_shot: bool
    consumed: bool
    resolution_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.approval-resolution.v1":
            raise ContractError("unsupported approval resolution schema")
        if self.decision is PermissionDecision.ALLOW_ONCE and not self.one_shot:
            raise ContractError("material approval must be one-shot")
        body = {
            "schema_version": self.schema_version,
            "approval_id": self.approval_id,
            "permission_request_sha256": self.permission_request_sha256,
            "decision": self.decision,
            "actor": self.actor,
            "one_shot": self.one_shot,
            "consumed": self.consumed,
        }
        if self.resolution_sha256 != canonical_sha256(body):
            raise ContractError("approval resolution digest mismatch")


@dataclass(frozen=True)
class PermissionReceiptV1:
    schema_version: str
    permission_request_sha256: str
    approval_resolution_sha256: str
    decision: PermissionDecision
    reason: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.permission-receipt.v1":
            raise ContractError("unsupported permission receipt schema")
        body = {
            "schema_version": self.schema_version,
            "permission_request_sha256": self.permission_request_sha256,
            "approval_resolution_sha256": self.approval_resolution_sha256,
            "decision": self.decision,
            "reason": self.reason,
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("permission receipt digest mismatch")


def build_permission_request(
    *, request_id: str, action: str, scope_sha256: str, **bindings: object
) -> PermissionRequestV1:
    body = {
        "schema_version": "chemsmart.permission-request.v1",
        "request_id": request_id,
        "action": action,
        "scope_sha256": scope_sha256,
        "command_sha256": str(bindings.get("command_sha256", "")),
        "input_sha256s": tuple(bindings.get("input_sha256s", ())),
        "project_sha256": str(bindings.get("project_sha256", "")),
        "environment_sha256": str(bindings.get("environment_sha256", "")),
        "provider": str(bindings.get("provider", "")),
        "quota_scope": str(bindings.get("quota_scope", "")),
    }
    return PermissionRequestV1(
        **body, request_sha256=canonical_sha256(body)
    )


def resolve_permission(
    request: PermissionRequestV1,
    *,
    approval: ApprovalResolutionV1 | None = None,
) -> PermissionReceiptV1:
    if request.action in _MATERIAL_ACTIONS and approval is not None:
        raise ContractError(
            "material approvals must be consumed by the persistent event store"
        )
    return _evaluate_permission(request, approval=approval)


def _evaluate_permission(
    request: PermissionRequestV1,
    *,
    approval: ApprovalResolutionV1 | None = None,
) -> PermissionReceiptV1:
    if request.action in _READ_ONLY_ACTIONS:
        decision = PermissionDecision.AUTO_ALLOW
        reason = "policy.read_only_or_fixture"
        approval_sha = ""
    elif approval is None:
        decision = PermissionDecision.NEEDS_USER
        reason = "policy.explicit_approval_required"
        approval_sha = ""
    elif approval.permission_request_sha256 != request.request_sha256:
        decision = PermissionDecision.DENY
        reason = "policy.approval_binding_mismatch"
        approval_sha = approval.resolution_sha256
    elif approval.consumed:
        decision = PermissionDecision.DENY
        reason = "policy.approval_already_consumed"
        approval_sha = approval.resolution_sha256
    elif approval.decision is PermissionDecision.ALLOW_ONCE:
        decision = PermissionDecision.ALLOW_ONCE
        reason = "policy.exact_one_shot_approval"
        approval_sha = approval.resolution_sha256
    else:
        decision = PermissionDecision.DENY
        reason = "policy.user_denied"
        approval_sha = approval.resolution_sha256
    body = {
        "schema_version": "chemsmart.permission-receipt.v1",
        "permission_request_sha256": request.request_sha256,
        "approval_resolution_sha256": approval_sha,
        "decision": decision,
        "reason": reason,
    }
    return PermissionReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def build_approval_resolution(
    *,
    approval_id: str,
    request: PermissionRequestV1,
    allow: bool,
    actor: str,
) -> ApprovalResolutionV1:
    body = {
        "schema_version": "chemsmart.approval-resolution.v1",
        "approval_id": approval_id,
        "permission_request_sha256": request.request_sha256,
        "decision": (
            PermissionDecision.ALLOW_ONCE if allow else PermissionDecision.DENY
        ),
        "actor": actor,
        "one_shot": bool(allow),
        "consumed": False,
    }
    return ApprovalResolutionV1(
        **body, resolution_sha256=canonical_sha256(body)
    )


class PermissionLedgerV1:
    """Compatibility ledger; material consumption moved to RuntimeEventStore."""

    def __init__(self) -> None:
        self._consumed_approval_ids: set[str] = set()
        self._lock = Lock()

    def resolve(
        self,
        request: PermissionRequestV1,
        *,
        approval: ApprovalResolutionV1 | None = None,
    ) -> PermissionReceiptV1:
        if request.action in _READ_ONLY_ACTIONS or approval is None:
            return resolve_permission(request, approval=approval)
        raise ContractError(
            "material approvals require persistent, crash-stable "
            "RuntimeEventStore consumption"
        )

    def _unsafe_legacy_resolve(
        self,
        request: PermissionRequestV1,
        *,
        approval: ApprovalResolutionV1,
    ) -> PermissionReceiptV1:
        """Unexported historical behavior for migration-only replay tests."""

        with self._lock:
            if approval.approval_id in self._consumed_approval_ids:
                consumed_body = {
                    "schema_version": "chemsmart.approval-resolution.v1",
                    "approval_id": approval.approval_id,
                    "permission_request_sha256": (
                        approval.permission_request_sha256
                    ),
                    "decision": approval.decision,
                    "actor": approval.actor,
                    "one_shot": approval.one_shot,
                    "consumed": True,
                }
                consumed = ApprovalResolutionV1(
                    **consumed_body,
                    resolution_sha256=canonical_sha256(consumed_body),
                )
                return _evaluate_permission(request, approval=consumed)
            receipt = _evaluate_permission(request, approval=approval)
            if receipt.decision is PermissionDecision.ALLOW_ONCE:
                self._consumed_approval_ids.add(approval.approval_id)
            return receipt


__all__ = [
    "ApprovalResolutionV1",
    "PermissionDecision",
    "PermissionReceiptV1",
    "PermissionRequestV1",
    "PermissionLedgerV1",
    "build_approval_resolution",
    "build_permission_request",
    "resolve_permission",
]
