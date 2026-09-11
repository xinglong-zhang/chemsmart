"""A model-proposed precision became a user-authorised one, unseen.

The owner's ruling is that the session restates the tolerance and the
host freezes it. That only holds if the human sees the restatement: the
declared-observable panel rendered the meaning, the unit and the
expectation band, and never the required precision or where the session
got it. So the one decision covered a contract nobody was shown.

The absence matters as much as the value. "Declare no tolerance at all"
is the cheapest escape from an obligation and the least visible one, so
a requested observable with no required precision says so in the panel,
on a page a reviewer reads before granting anything.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from chemsmart.agent.tui.review import _declared_observable_panel

pytestmark = pytest.mark.capability("rule:wake.restate_observable")


def _rendered(records):
    from rich.console import Console

    panel = _declared_observable_panel(
        SimpleNamespace(requested_observable_declarations=records)
    )
    console = Console(width=200, record=True)
    console.print(panel)
    return console.export_text()


def test_a_declared_tolerance_and_its_source_are_shown():
    text = _rendered(
        (
            {
                "observable_id": "e-couple",
                "unit": "V",
                "meaning": "the +1/0 couple versus ferrocene",
                "required_tolerance": 0.2,
                "tolerance_basis": "Task: 'good to about plus or minus 0.2 V'",
            },
        )
    )
    assert "required precision: +/-0.2 V" in text
    assert "good to about plus or minus 0.2 V" in text


def test_an_absent_tolerance_is_shown_as_an_absence():
    text = _rendered(
        (
            {
                "observable_id": "e-couple",
                "unit": "V",
                "meaning": "the +1/0 couple versus ferrocene",
            },
        )
    )
    assert "none declared" in text
    assert "if the task states one, it is not on this contract" in text


def test_a_diagnostic_is_owed_no_precision_and_is_not_marked_for_one():
    """A diagnostic is the session's own prediction and never owed, so
    an absent tolerance on one is not a gap the reviewer should read as
    a missing obligation."""

    text = _rendered(
        (
            {
                "observable_id": "spinsq",
                "unit": "1",
                "meaning": "the spin-contamination diagnostic",
                "role": "diagnostic",
            },
        )
    )
    assert "none declared" not in text
