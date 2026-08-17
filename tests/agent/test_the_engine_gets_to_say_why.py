"""An unanticipated engine failure must still reach whoever has to re-plan.

The classifier's vocabulary of failure classes is closed.  A failure nobody
wrote a rule for lands in ``native_runtime``, whose canonical template says only
that an error occurred -- exactly as useful as the empty string it replaced.

That is not hypothetical.  Four coupled-cluster nodes once failed because
canonical CCSD under the RIJK approximation is not size-consistent, which is
qualitatively fatal for the interaction energy being computed.  ORCA said so in
plain words and named the remedy.  No layer above it carried a single one of
those words, so the reason had to be recovered by opening the file by hand.

These tests execute that path end to end rather than reading the code: an
unmatched failure keeps its stable class *and* carries the engine's own lines,
those lines survive into the failed node's summary, host locators are redacted
out of them, and xTB -- one of the three programs this release executes -- is
covered too.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from chemsmart.agent.executor import _execution_failure_summary
from chemsmart.io.native_failure import (
    NativeFailureSummaryV1,
    summarize_orca_native_failure,
    summarize_xtb_native_failure,
)

_FIXTURES = Path("tests/data/agent/native_failures")


def _lines(name: str) -> tuple[str, ...]:
    return tuple((_FIXTURES / name).read_text(encoding="utf-8").splitlines())


class _Receipt:
    """The fields the failure summary actually reads off an execution receipt."""

    def __init__(self, **fields):
        self.execution_state = "failed"
        self.wrapper_exit_status = 0
        self.child_exit_status = 1
        self.findings = ("orca.native_failure.native_runtime",)
        self.result_validation_receipt_sha256 = "v" * 64
        for name, value in fields.items():
            setattr(self, name, value)


class _Validation:
    def __init__(self, observations):
        self.observations = observations


class _Host:
    def __init__(self, summary):
        self.result_validation_receipts = {
            "v"
            * 64: _Validation(
                {"orca": {"native_failure": dict(summary.as_dict())}}
            )
        }


def test_an_unclassified_failure_still_carries_the_engine_s_own_words():
    summary = summarize_orca_native_failure(
        _lines("orca_unclassified_size_consistency.txt")
    )

    assert summary is not None
    assert summary.termination_state == "error_termination"
    # No rule matches this one -- that is precisely the case under test.
    assert summary.error_class == "native_runtime"
    quoted = " ".join(summary.engine_lines).lower()
    assert "size-consistent" in quoted, "the rule that was broken"
    assert "non-covalent binding" in quoted, "why it matters for this quantity"
    assert "turn off rijk" in quoted, "the remedy the engine named"


def test_the_diagnosis_before_an_abort_marker_is_not_lost():
    """Engines diagnose, then abort; the abort line is the least informative.

    Observed on a real run: the terminating line was quoted correctly and said
    only "Error (ORCA_MAIN): ... aborting the run", with its own reason elided,
    while the four lines naming the rule and its remedy sat above it.
    """

    summary = summarize_orca_native_failure(
        _lines("orca_unclassified_size_consistency.txt")
    )

    assert summary.engine_lines[0].startswith("WARNING: Canonical CCSD")
    assert "aborting the run" in summary.engine_lines[-1]


def test_host_locators_are_redacted_from_quoted_lines():
    summary = summarize_orca_native_failure(
        _lines("orca_unclassified_size_consistency.txt")
    )

    quoted = " ".join(summary.engine_lines)
    assert "/opt/chem" not in quoted
    assert "<redacted>" in quoted


def test_the_existing_credential_and_url_scrubbing_still_holds():
    """The auxiliary-basis fixture carries a bearer token and a URL."""

    summary = summarize_orca_native_failure(_lines("orca_auxiliary_basis.txt"))

    quoted = " ".join(summary.engine_lines)
    assert "Bearer" not in quoted
    assert "https://" not in quoted
    assert "/secret/path" not in quoted


def test_the_reason_reaches_the_failed_node_summary():
    """The join that was missing: receipt observations -> node failure text."""

    summary = summarize_orca_native_failure(
        _lines("orca_unclassified_size_consistency.txt")
    )

    text = _execution_failure_summary(_Receipt(), _Host(summary))

    assert "execution_state=failed" in text
    assert "child_exit_status=1" in text
    assert "engine reported (verbatim):" in text
    assert "turn off rijk" in text.lower(), "the remedy must reach the node"


def test_a_summary_without_a_host_still_reports_what_it_can():
    """Losing the engine's words must never cost the rest of the summary."""

    text = _execution_failure_summary(_Receipt())

    assert "execution_state=failed" in text
    assert "engine reported" not in text


def test_a_successful_node_reports_no_failure():
    assert _execution_failure_summary(_Receipt(execution_state="validated")) == ""


def test_xtb_failures_are_classified_and_quoted():
    """xTB is executable in this release; it was not covered at all."""

    summary = summarize_xtb_native_failure(_lines("xtb_scc_not_converged.txt"))

    assert summary is not None
    assert summary.program == "xtb"
    assert summary.error_class == "scf_convergence"
    quoted = " ".join(summary.engine_lines).lower()
    assert "not converged" in quoted
    assert "iterations" in quoted, "the remedy the engine named must survive"


def test_a_normally_terminated_run_produces_no_summary():
    assert (
        summarize_xtb_native_failure(
            ("normal xtb output", " * finished run on 2026/08/18 at 02:00:00")
        )
        is None
    )


def test_an_unredacted_engine_line_cannot_be_constructed():
    """The contract refuses to carry a locator, not merely avoid making one."""

    with pytest.raises(ValueError, match="not redacted"):
        NativeFailureSummaryV1(
            schema_version="chemsmart.native-failure-summary.v1",
            program="orca",
            termination_state="error_termination",
            error_class="native_runtime",
            diagnostic_lines=(),
            engine_lines=("failed reading /home/someone/private/run.inp",),
        )
