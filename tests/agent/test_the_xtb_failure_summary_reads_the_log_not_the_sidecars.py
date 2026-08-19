"""Every healthy live xTB run failed validation: the summary read sidecars.

The first real xTB execution through `agent execute` -- three normally
terminated single points -- came back with
`xtb.native_failure.incomplete_output` on every node. The failure
summarizer selected artifacts of kind "program_output", which for xTB are
the charges/wbo/xtbrestart/topology sidecars; the actual log, kind
"xtb_output", was excluded, so the classifier never saw the
"finished run" sentinel and called the run incomplete. The summary now
reads the log.
"""

from __future__ import annotations

from chemsmart.io.native_failure import summarize_xtb_native_failure

_HEALTHY_LOG = (
    "          | HOMO-LUMO GAP              14.629673576893 eV   |",
    "------------------------------------------------------------",
    " * finished run on 2026/08/20 at 02:14:31.690",
    " total:",
    " * wall-time:     0 d,  0 h,  0 min,  0.005 sec",
)
_SIDECAR = ("0.1234", "-0.5678", "0.4444")


def test_a_finished_log_is_not_a_failure():
    assert summarize_xtb_native_failure(iter(_HEALTHY_LOG)) is None


def test_a_sidecar_alone_would_have_read_as_incomplete():
    summary = summarize_xtb_native_failure(iter(_SIDECAR))
    assert summary is not None
    assert summary.error_class == "incomplete_output"
