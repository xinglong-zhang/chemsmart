"""A run that simply stops has a reason, and it is in the last lines.

The engine-failure reader quotes a window around an abort marker. That is the
right shape for a program that announces its failure, and it collects nothing
at all for a program whose output just ends: with no marker, the window never
opens, and the session is handed `error_class: incomplete_output` with
`engine_lines: []` -- a classification with no reason attached, which no
session can act on.

The fixture is the tail of a real approved ORCA scan launched through the
Agent's own path. ORCA had refused the very first scan point because the
requested starting dihedral was about 120 degrees away from the geometry it was
given, and it said so plainly on its last line. None of that reached the typed
evidence.

The loop already keeps the recent substantive lines, for the backwards window.
This is the case where those lines are the whole of what the engine said.
"""

from __future__ import annotations

from chemsmart.io.native_failure import summarize_orca_native_failure

_REFUSED_SCAN = (
    "tests/data/ORCATests/outputs/hooh_scan_constraint_refusal_excerpt.out"
)
_COMPLETED = "tests/data/ORCATests/outputs/CO2.out"


def _summary(path):
    with open(path, encoding="utf-8", errors="replace") as handle:
        return summarize_orca_native_failure(handle.read().splitlines())


def test_the_failure_is_still_classified_as_an_incomplete_output():
    """The class is unchanged; only the evidence attached to it grows."""

    summary = _summary(_REFUSED_SCAN)

    assert summary.error_class == "incomplete_output"
    assert summary.termination_state == "incomplete"


def test_the_engine_gets_to_say_why_it_stopped():
    summary = _summary(_REFUSED_SCAN)

    assert summary.engine_lines
    quoted = " ".join(summary.engine_lines)
    assert "could not impose initial constraints" in quoted


def test_the_quoted_lines_are_still_redacted():
    """The reason travels; the build path it names does not."""

    summary = _summary(_REFUSED_SCAN)

    quoted = " ".join(summary.engine_lines)
    assert "qcgstep.cpp" not in quoted
    assert "<redacted>" in quoted


def test_a_completed_run_is_not_turned_into_a_failure():
    """The control: quoting a tail must not invent a failure that never was."""

    assert _summary(_COMPLETED) is None
