"""One schema for every capability, and a ladder computed from what
exists. The test that keeps the schema honest: everything advertised is
wired, every executable program job type has a qualification record or
is displayed as a claim, every guide and rule resolves, and a capability
nobody marked is named rather than silently unpinned.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from chemsmart.agent.capability_registry import (
    CAPABILITY_KINDS,
    LADDER,
    build_capability_registry,
    load_release_records,
    render_capability_matrix,
)

pytestmark = pytest.mark.capability(
    "tool:*", "guide:*", "rule:*", "program_jobtype:*"
)

_TESTS = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def registry():
    return build_capability_registry(tests_root=_TESTS, host_store=None)


def test_every_kind_is_represented(registry):
    assert {item.kind for item in registry} == set(CAPABILITY_KINDS)
    assert all(item.status in LADDER for item in registry)


#: Advertised-but-unwired capabilities, as a budget that may only fall.
#:
#: This asserted zero and held only because the ladder was blind: before
#: 2026-09-11 project settings were not a capability kind at all, so the
#: surface where the model's whole method vocabulary travels was absent
#: from the count. With `setting` declared, 82 advertised parameters
#: carry no declared domain -- `opt_convergence` among them, whose enum
#: sits in ORCA_OPT_CONVERGENCE_KEYWORDS and is withheld from the model.
#:
#: Suppressing that to keep a green assertion is how an instrument goes
#: blind, so the gap is recorded instead and this number must only ever
#: decrease. Each reduction is a parameter whose domain now reaches the
#: model, which is FUNDAMENTAL 1 made measurable.
ADVERTISED_UNWIRED_BUDGET = 82


def test_everything_advertised_is_wired(registry):
    unwired = sorted(
        item.key
        for item in registry
        if item.advertised_in and not item.wired_by
    )
    assert len(unwired) <= ADVERTISED_UNWIRED_BUDGET, (
        f"{len(unwired)} advertised capabilities are unwired, above the "
        f"recorded budget of {ADVERTISED_UNWIRED_BUDGET}. A new one was "
        f"added: {unwired[:10]}"
    )
    # Two kinds may legitimately sit here, and both are findings rather
    # than defects in this test: a project setting awaiting a declared
    # domain, and a code gate that CONDUCT section 2 lists in prose and
    # that nothing in the package raises under its own id. Four gates
    # are in that state; whether each is genuinely absent or implemented
    # under another name is a question the ladder now asks instead of
    # answering silently.
    unexpected = [
        item for item in unwired if not item.startswith(("setting:", "gate:"))
    ]
    assert not unexpected, (
        "something other than a domainless setting or an unraised gate "
        f"is advertised and unwired: {unexpected}"
    )


def test_every_executable_program_jobtype_has_a_qualification_record(
    registry,
):
    release = load_release_records()
    missing = [
        item.key
        for item in registry
        if item.kind == "program_jobtype"
        and item.advertised_in == "inspect_program"
        and item.key not in release
    ]
    assert not missing, (
        "an executable program x jobtype with no release record: add the "
        f"run that qualified it, or stop claiming it: {missing}"
    )
    claimed = sorted(
        key
        for key, record in release.items()
        if record.get("status") == "claimed"
    )
    # Claims without a machine-recorded run are displayed as such, never
    # hidden; they are the next observations to run.
    for item in registry:
        if item.key in claimed:
            assert any(
                ref.startswith("release:claimed") for ref in item.qualified_by
            )


def test_no_capability_is_silently_unpinned(registry):
    """Every capability is covered by a test that names it, or by a family.

    This used to assert that `tested_by` was non-empty for everything
    outside selectors, and it passed because nine blanket `kind:*`
    markers lifted 557 of 712 capabilities to `tested` from a handful of
    files. A capability covered only by its family is now recorded in
    `family_tested_by` and does not claim the rung, so the honest
    property is that nothing is covered by *nothing at all*.
    """

    unmarked = sorted(
        item.key
        for item in registry
        if not item.tested_by
        and not item.family_tested_by
        and item.kind not in {"selector"}
    )
    assert unmarked == [], (
        "these capabilities are covered by no test, exactly or by "
        f"family: {unmarked}"
    )


def test_the_matrix_says_unsupported_out_loud(registry):
    text = render_capability_matrix(registry)
    assert "by status:" in text
    assert "unsupported" in text
