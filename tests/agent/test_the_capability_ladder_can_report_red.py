"""The ladder is the instrument; an instrument that cannot fail lies.

``chemsmart agent capabilities`` exists to say what the agent can do and
how far each capability has climbed -- declared, wired, advertised,
tested, qualified. It was built to prevent exactly the defect class this
campaign keeps finding: a surface advertised to the model that is not
executable end to end.

It could not detect any of it. Measured on 2026-09-11, before this test:

- ``setting`` was not a capability kind at all, so project-owned
  parameters -- the load-bearing surface of the hub thesis, where the
  model's whole method vocabulary travels -- were absent from the
  ladder. Two tests already carried ``setting:`` markers that matched
  nothing.
- the marker regex was ``capability\\(([^)]*)\\)``, which also matched
  ``program_capability(...)`` and ``engine_capability(...)``, so the
  instrument's own input carried sixteen tokens no capability could ever
  have: ``cpu``, ``opt``, ``sp``, ``ts``, ``xtb``, ``2``, ``3``.
- twenty well-formed markers matched no capability, including **all six
  anomaly-sensor predicates** -- the charter's "the anomaly has standing"
  mechanism, pinned against the expression-predicate vocabulary it does
  not belong to -- and seven code gates from ``CONDUCT.md`` section 2's
  own list, which had no registry.
- ``wired_by`` was the constant string ``"reader accessor"`` for all 393
  selectors, so the rung would have reported ``wired`` for a selector
  whose accessor had been deleted.

A rung that cannot fail measures nothing, and a marker that matches
nothing certifies nothing. This test holds the instrument to the
property it is for.
"""

from pathlib import Path

import pytest

from chemsmart.agent.capability_registry import (
    CAPABILITY_KINDS,
    build_capability_registry,
)
from chemsmart.agent.capability_registry import (
    test_markers as read_markers,  # not a test: pytest would collect it
)

#: Kinds whose members genuinely share one wiring mechanism, so a single
#: `wired_by` string is the truth rather than a decoration: every
#: expression operation is joined by `evaluate_quantity_expression`,
#: every guide by `open_guide`, and so on. A kind absent from this set
#: must compute the rung per capability. Adding a kind therefore forces
#: the decision at declaration time, which is the whole point.
UNIFORM_WIRING = frozenset(
    {
        "operation",
        "predicate",
        "constant",
        "skill",
        "guide",
        "rule",
        "tool",
        "program_jobtype",
    }
)

TESTS_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def ladder():
    return build_capability_registry(tests_root=TESTS_ROOT, host_store=None)


@pytest.fixture(scope="module")
def markers():
    return read_markers(TESTS_ROOT)


@pytest.mark.capability("gate:capability.marker_names_a_capability")
def test_every_marker_names_a_capability(ladder, markers):
    """A marker that matches nothing records coverage against nothing."""

    keys = {item.key for item in ladder}
    orphans = sorted(
        pattern
        for pattern in markers
        if not pattern.endswith(":*") and pattern not in keys
    )
    assert not orphans, (
        "these capability markers match no capability, so the tests "
        "carrying them certify nothing:\n  " + "\n  ".join(orphans)
    )


@pytest.mark.capability("gate:capability.marker_names_a_capability")
def test_no_marker_is_a_bare_token(markers):
    """The marker regex must read markers, not every call it resembles."""

    bare = sorted(pattern for pattern in markers if ":" not in pattern)
    assert not bare, (
        "these are not capability ids; the marker pattern is picking up "
        f"unrelated call sites: {bare}"
    )


@pytest.mark.capability("gate:capability.marker_names_a_capability")
def test_every_marker_kind_is_a_declared_kind(markers):
    kinds = {pattern.split(":", 1)[0] for pattern in markers if ":" in pattern}
    unknown = sorted(kinds - set(CAPABILITY_KINDS))
    assert not unknown, (
        f"markers name kinds the registry does not carry: {unknown}. "
        f"Declared kinds are {CAPABILITY_KINDS}"
    )


@pytest.mark.capability("gate:capability.wildcard_is_not_sole_coverage")
def test_a_wildcard_is_never_a_capabilitys_only_coverage(ladder, markers):
    """``selector:*`` certified all 393 selectors from two files.

    A blanket marker is legitimate for a family test, and it must not be
    the *only* thing standing behind a capability -- otherwise one marker
    reports coverage for hundreds of surfaces nobody exercised. This is
    the saturation that let every defect in this campaign's ladder pass.
    """

    wildcards = {p for p in markers if p.endswith(":*")}
    assert wildcards, "expected family wildcards to exist"
    exact = {p for p in markers if not p.endswith(":*")}

    # The property is about the *rung*, not about writing 557 markers: a
    # family test is real evidence about the family and no evidence
    # about any one member, so it must not lift a capability to
    # `tested`. Before this split, nine blanket markers reported
    # `tested` for 557 of 712 capabilities.
    lifted = sorted(
        item.key
        for item in ladder
        if item.status == "tested"
        and item.key not in exact
        and not item.qualified_by
    )
    assert not lifted, (
        f"{len(lifted)} capabilities reach status 'tested' without any "
        "test naming them, so a wildcard is lifting the rung. First "
        f"ten: {lifted[:10]}"
    )
    covered = [item for item in ladder if item.family_tested_by]
    assert covered, (
        "no capability records family coverage, so the wildcard split "
        "is not wired: family_tested_by should carry what the blanket "
        "markers actually cover"
    )


@pytest.mark.capability("gate:capability.rung_is_computed")
def test_the_ladder_is_not_saturated(ladder):
    """The instrument must be able to report every rung, not just the top.

    Before the wildcard split, 557 of 712 capabilities reported
    ``tested`` and **zero** could report ``declared``, ``wired`` or
    ``advertised`` -- nine blanket markers lifted almost everything to
    the top rung, so the ladder had one value and no discriminating
    power. A distribution with one bucket is not a measurement.

    Deliberately not asserted here: that ``wired_by`` *differs* between
    members of a kind. A computed rung may legitimately return the same
    answer for every member -- all seven anomaly signals are raised from
    ``tool_runtime.py``, so uniformity there is the truth rather than a
    decoration. Variance is not the property; being computed is, and
    that is held by the kinds' own construction sites plus
    ``UNIFORM_WIRING`` above, which forces the decision when a kind is
    added.
    """

    counts: dict[str, int] = {}
    for item in ladder:
        counts[item.status] = counts.get(item.status, 0) + 1
    assert len(counts) > 1, (
        "every capability reports the same status, so the ladder "
        f"discriminates nothing: {counts}"
    )
    below = sum(
        value
        for status, value in counts.items()
        if status in {"declared", "wired", "advertised"}
    )
    assert below, (
        "no capability sits below 'tested', which is what a saturated "
        f"instrument looks like: {counts}"
    )


@pytest.mark.capability("gate:capability.rung_is_computed")
def test_an_unwired_capability_can_say_so(ladder):
    """A declared surface nothing raises must be visible as such.

    Concretely: ``CONDUCT.md`` section 2 lists the code gates in prose,
    and four of the seven now declared in ``CODE_GATES`` are raised
    nowhere in the agent package under their own id. Whether each is
    genuinely absent or merely implemented under another name is a
    question for a reader -- the point is that the ladder now asks it
    instead of reporting silence.
    """

    unwired = [item for item in ladder if not item.wired_by]
    assert unwired, (
        "no capability can report itself unwired, so the wired rung "
        "cannot fail"
    )
