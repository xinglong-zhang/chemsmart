"""Shared result model for deterministic command contracts."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Literal


ContractSeverity = Literal["warn", "reject"]


@dataclass(frozen=True)
class CommandContractIssue:
    rule_id: str
    severity: ContractSeverity
    message: str
    evidence: dict[str, Any] = field(default_factory=dict)
    missing_info: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def reject(
    rule_id: str,
    message: str,
    evidence: dict[str, Any],
    missing_info: tuple[str, ...] = (),
) -> CommandContractIssue:
    return CommandContractIssue(
        rule_id=rule_id,
        severity="reject",
        message=message,
        evidence=evidence,
        missing_info=missing_info,
    )


__all__ = ["CommandContractIssue", "ContractSeverity", "reject"]
