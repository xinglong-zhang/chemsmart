"""The terminal shows the report's science; the file keeps its receipts.

The parser is keyed to the host renderer's exact headings; a source-pin
assertion fails loudly if the renderer drifts, and unknown sections pass
through instead of being dropped.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

pytest.importorskip("textual")

from rich.console import Console  # noqa: E402

from chemsmart.agent.tui.report import (  # noqa: E402
    looks_like_host_report,
    render_report_for_humans,
)

_TOOLCHAIN_MD = """# Host-validated structured analysis

Completion receipt: `{h1}`
Toolchain plan: `{h2}`

## Thermochemical conditions

| Stage | Temperature (K) | Pressure (atm) |
|---|---:|---:|
| thermo | `298.15` | `1.0` |

Claim record: `{h3}`

## Host-rendered numerical claims

| Claim | Value | Unit | Source receipt |
|---|---:|---|---|
| gibbs | `-5.070370761845` | hartree | `{h4}` |

Scientific decision: not recorded -- interpretation is a session act; this
run executed extraction, thermochemistry, expressions, validation verdicts,
and claim rendering only.
""".format(h1="1" * 64, h2="2" * 64, h3="3" * 64, h4="4" * 64)

_POLICY_MD = """# Host-validated structured analysis

Completion receipt: `{h1}`
Claim record: `{h2}`

## Task-owned conditions

| Condition | Value | Unit | Origin | Evidence |
|---|---:|---|---|---|
| temperature | `298.15` | K | task | `{h3}` |

## Host-rendered numerical claims

| Claim | Value | Unit | Source receipt |
|---|---:|---|---|
| pka | `4.76` | 1 | `{h4}` |

## Method rationale

- The anchor ladder cancels the proton term.

## Uncertainties

- Absolute pKa carries the documented +2 systematic.
""".format(h1="5" * 64, h2="6" * 64, h3="7" * 64, h4="8" * 64)


def _flatten(markdown: str) -> str:
    console = Console(record=True, width=160)
    console.print(render_report_for_humans(markdown))
    return console.export_text()


def test_the_toolchain_report_reads_as_science():
    text = _flatten(_TOOLCHAIN_MD)

    assert "Analysis results" in text
    assert "gibbs" in text and "-5.070370761845" in text and "hartree" in text
    assert "298.15" in text and "1.0" in text
    assert "Source receipt" not in text
    assert "Completion receipt" not in text
    assert "No interpretation was recorded in this run" in text
    assert re.search(r"[0-9a-f]{64}", text) is None


def test_the_policy_report_keeps_decision_prose_verbatim():
    text = _flatten(_POLICY_MD)

    assert "The anchor ladder cancels the proton term." in text
    assert "documented +2 systematic" in text
    assert "pka" in text and "4.76" in text
    # Every Evidence cell was a digest, so the column is dropped whole.
    assert "Evidence" not in text
    assert re.search(r"[0-9a-f]{64}", text) is None


def test_a_human_evidence_column_survives():
    md = _POLICY_MD.replace("`" + "7" * 64 + "`", "measured on this host")
    text = _flatten(md)

    assert "Evidence" in text
    assert "measured on this host" in text


def test_unknown_sections_pass_through_with_digests_elided():
    md = _TOOLCHAIN_MD + "\n## A future section\n\nkeep me `" + "9" * 64 + "`\n"
    text = _flatten(md)

    assert "A future section" in text
    assert "keep me" in text
    assert re.search(r"[0-9a-f]{64}", text) is None


def test_detection_only_claims_host_reports():
    assert looks_like_host_report(_TOOLCHAIN_MD)
    assert not looks_like_host_report("# My notes\n\nhello")


def test_the_parser_is_pinned_to_the_live_renderer():
    source = Path("chemsmart/agent/tool_runtime.py").read_text()
    for heading in (
        '"# Host-validated structured analysis"',
        '"## Host-rendered numerical claims"',
        '"## Thermochemical conditions"',
        '"## Task-owned conditions"',
        '"| Claim | Value | Unit | Source receipt |"',
        '"Scientific decision: not recorded -- interpretation is "',
    ):
        assert heading in source, f"host renderer drifted: {heading}"
