"""The brand: two-tone wordmark, the spectrum rule, honest degradation."""

from __future__ import annotations

import pytest

pytest.importorskip("textual")

from rich.console import Console  # noqa: E402

from chemsmart.agent.tui.panels import (  # noqa: E402
    SPECTRUM_COLORS,
    spectrum_rule,
    wordmark,
)


def test_the_wordmark_is_the_logo_pair():
    mark = wordmark()
    spans = {(span.style) for span in mark.spans}
    assert "bold magenta" in spans
    assert "bold green" in spans
    assert mark.plain.startswith("CHEMSMART Agent")


def test_the_spectrum_rule_carries_four_segments():
    rule = spectrum_rule(24)
    styles = [span.style for span in rule.spans]
    assert styles == list(SPECTRUM_COLORS)
    assert set(rule.plain) <= set("▂▄▆█▁")
    assert len(rule.plain) == 24


def test_plain_mode_degrades_to_a_dim_rule():
    rule = spectrum_rule(24, plain=True)
    assert rule.plain == "─" * 24
    assert all(span.style == "dim" for span in rule.spans)


def test_both_render_without_color_support():
    console = Console(record=True, width=60, no_color=True)
    console.print(wordmark())
    console.print(spectrum_rule(44))
    console.print(spectrum_rule(44, plain=True))
    text = console.export_text()
    assert "CHEMSMART Agent" in text
    assert "▂" in text and "─" in text
