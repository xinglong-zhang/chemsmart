"""Version-controlled convention rules contributed by builtin skills.

These are declared in Python, never in a user overlay, so editing a
``SKILL.md`` cannot introduce or alter one.

Scope note for this round: consumers **record** the applicable conventions
alongside a reported quantity so the convention in force is auditable. None of
them changes a scientific verdict. Rules that would change a verdict require
execution evidence and are deliberately absent.
"""

from __future__ import annotations

from chemsmart.agent.skills.rules import (
    DeterministicConventionRuleV1,
    build_convention_rule,
    rules_for_scope,
)

_SKILL_ID = "scientific-conventions"


BUILTIN_CONVENTION_RULES: tuple[DeterministicConventionRuleV1, ...] = (
    build_convention_rule(
        rule_id="standard-state-pressure",
        skill_id=_SKILL_ID,
        scope="thermochemistry",
        subject="standard state pressure",
        convention=(
            "The modern standard state is 1 bar. A value computed at 1 atm "
            "differs by R*ln(1.01325) in the standard entropy; report which "
            "pressure the calculation used."
        ),
    ),
    build_convention_rule(
        rule_id="rotational-symmetry-number",
        skill_id=_SKILL_ID,
        scope="thermochemistry",
        subject="rotational symmetry number",
        convention=(
            "Rotational entropy uses the symmetry number sigma taken from the "
            "point group; omitting it inflates the entropy by R*ln(sigma)."
        ),
    ),
    build_convention_rule(
        rule_id="term-value-orientation",
        skill_id=_SKILL_ID,
        scope="result_validation",
        subject="term value orientation",
        convention=(
            "A term value is E(upper state) minus E(ground state) and is "
            "non-negative; a negative value means the states were ordered the "
            "wrong way round."
        ),
    ),
)


BUILTIN_CONVENTION_RULE_SHA256S: tuple[str, ...] = tuple(
    sorted(item.rule_sha256 for item in BUILTIN_CONVENTION_RULES)
)


def conventions_for_scope(
    scope: str,
    rules: tuple[
        DeterministicConventionRuleV1, ...
    ] = BUILTIN_CONVENTION_RULES,
) -> tuple[DeterministicConventionRuleV1, ...]:
    """Return the builtin convention rules bound to ``scope``."""

    return rules_for_scope(rules, scope)


__all__ = [
    "BUILTIN_CONVENTION_RULES",
    "BUILTIN_CONVENTION_RULE_SHA256S",
    "conventions_for_scope",
]
