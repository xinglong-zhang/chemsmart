from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Literal

HarnessVerdict = Literal["ok", "warn", "reject"]
IssueSeverity = Literal["warn", "reject"]


@dataclass(frozen=True)
class InvariantIssue:
    rule_id: str
    severity: IssueSeverity
    message: str
    evidence: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class InvariantResult:
    rule_id: str
    verdict: HarnessVerdict = "ok"
    issues: tuple[InvariantIssue, ...] = ()
    context: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class HarnessResult:
    verdict: HarnessVerdict = "ok"
    issues: tuple[InvariantIssue, ...] = ()
    rule_results: tuple[InvariantResult, ...] = ()

    @classmethod
    def from_rule_results(
        cls, rule_results: list[InvariantResult] | tuple[InvariantResult, ...]
    ) -> "HarnessResult":
        results = tuple(rule_results)
        issues = tuple(issue for result in results for issue in result.issues)
        verdict: HarnessVerdict = "ok"
        if any(issue.severity == "reject" for issue in issues):
            verdict = "reject"
        elif any(issue.severity == "warn" for issue in issues):
            verdict = "warn"
        return cls(verdict=verdict, issues=issues, rule_results=results)

    @property
    def failed_rule_ids(self) -> list[str]:
        seen: set[str] = set()
        failed: list[str] = []
        for issue in self.issues:
            if issue.rule_id not in seen:
                seen.add(issue.rule_id)
                failed.append(issue.rule_id)
        return failed

    def to_dict(self) -> dict[str, Any]:
        return {
            "verdict": self.verdict,
            "issues": [issue.to_dict() for issue in self.issues],
            "rule_results": [result.to_dict() for result in self.rule_results],
            "failed_rule_ids": self.failed_rule_ids,
        }
