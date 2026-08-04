"""Deterministic convention rules declared by domain-knowledge skills.

A convention rule states how a quantity is *expressed* -- its unit, its sign
convention, its standard state, the companion quantity it requires.  It never
states how *accurate* a result must be, and it never establishes readiness.
Both refusals are enforced in ``__post_init__`` so an overlay file cannot widen
the authority of a skill by editing text.
"""

from __future__ import annotations

from dataclasses import dataclass

from chemsmart.agent._contracts import (
    ContractError,
    canonical_sha256,
    require_identifier,
)

#: Scopes a convention rule may bind to.  ``thermochemistry`` is consumed by
#: :mod:`chemsmart.analysis.thermochemistry`; ``result_validation`` by
#: :func:`chemsmart.jobs.pyscf.validation.validate_pyscf_result`.
CONVENTION_SCOPES = ("result_validation", "thermochemistry")


@dataclass(frozen=True)
class DeterministicConventionRuleV1:
    """One host-consumable expression convention contributed by a skill."""

    schema_version: str
    rule_id: str
    skill_id: str
    scope: str
    subject: str
    convention: str
    rule_sha256: str
    readiness_authority: bool = False
    accuracy_authority: bool = False

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.deterministic-convention-rule.v1":
            raise ContractError("unsupported convention rule schema")
        if self.readiness_authority:
            raise ContractError("convention rules cannot establish readiness")
        if self.accuracy_authority:
            raise ContractError(
                "convention rules cannot assert accuracy requirements"
            )
        require_identifier(self.rule_id, "rule_id")
        require_identifier(self.skill_id, "skill_id")
        if self.scope not in CONVENTION_SCOPES:
            raise ContractError(
                f"convention rule scope must be one of {CONVENTION_SCOPES}"
            )
        for value, name in (
            (self.subject, "subject"),
            (self.convention, "convention"),
        ):
            if not str(value).strip():
                raise ContractError(
                    f"convention rule {name} must not be empty"
                )
        if self.rule_sha256 != canonical_sha256(self._body()):
            raise ContractError("convention rule digest mismatch")

    def _body(self) -> dict[str, object]:
        return {
            "schema_version": self.schema_version,
            "rule_id": self.rule_id,
            "skill_id": self.skill_id,
            "scope": self.scope,
            "subject": self.subject,
            "convention": self.convention,
            "readiness_authority": False,
            "accuracy_authority": False,
        }

    def as_dict(self) -> dict[str, object]:
        """Return a detached canonical mapping including the digest."""

        body = self._body()
        body["rule_sha256"] = self.rule_sha256
        return body


def build_convention_rule(
    *,
    rule_id: str,
    skill_id: str,
    scope: str,
    subject: str,
    convention: str,
) -> DeterministicConventionRuleV1:
    """Construct a digest-bound convention rule."""

    body = {
        "schema_version": "chemsmart.deterministic-convention-rule.v1",
        "rule_id": require_identifier(rule_id, "rule_id"),
        "skill_id": require_identifier(skill_id, "skill_id"),
        "scope": str(scope),
        "subject": str(subject),
        "convention": str(convention),
        "readiness_authority": False,
        "accuracy_authority": False,
    }
    return DeterministicConventionRuleV1(
        **body, rule_sha256=canonical_sha256(body)
    )


def rules_for_scope(
    rules: tuple[DeterministicConventionRuleV1, ...], scope: str
) -> tuple[DeterministicConventionRuleV1, ...]:
    """Return the subset of ``rules`` bound to ``scope``, digest-ordered."""

    selected = [item for item in rules if item.scope == scope]
    return tuple(sorted(selected, key=lambda item: item.rule_sha256))


__all__ = [
    "CONVENTION_SCOPES",
    "DeterministicConventionRuleV1",
    "build_convention_rule",
    "rules_for_scope",
]
