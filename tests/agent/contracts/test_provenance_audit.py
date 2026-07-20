"""Focused tests for the source-provenance comparison primitives."""

from __future__ import annotations

import json
from pathlib import Path

from scripts.review.audit_agent_provenance import (
    TextDocument,
    exact_line_matches,
    exact_token_matches,
    near_token_matches,
)

REPO_ROOT = Path(__file__).resolve().parents[3]


def test_exact_line_match_requires_six_substantial_nonblank_lines() -> None:
    lines = [
        (
            f"result_{index} = transform(value_{index}, option_{index}, "
            f"context_{index}, registry_{index}, provenance_{index})"
        )
        for index in range(8)
    ]
    candidate = TextDocument("chemsmart", "agent.py", "\n".join(lines))
    source = TextDocument(
        "reference",
        "source.py",
        "header = True\n\n" + "\n".join(lines) + "\nfooter = True\n",
    )

    matches = exact_line_matches([candidate], [source])

    assert len(matches) == 1
    assert matches[0].units == 8
    assert matches[0].tokens >= 50


def test_exact_token_match_finds_long_single_line_sequence() -> None:
    shared = " ".join(f"token_{index}" for index in range(70))
    candidate = TextDocument("chemsmart", "agent.md", shared)
    source = TextDocument("reference", "guide.md", f"prefix {shared} suffix")

    matches = exact_token_matches([candidate], [source])

    assert len(matches) == 1
    assert matches[0].tokens == 70


def test_near_token_match_flags_large_high_similarity_unit() -> None:
    candidate_tokens = [f"workflow_{index}" for index in range(120)]
    source_tokens = list(candidate_tokens)
    for index in range(0, 120, 20):
        source_tokens[index] = f"replacement_{index}"
    candidate = TextDocument(
        "chemsmart",
        "agent.txt",
        " ".join(candidate_tokens),
    )
    source = TextDocument(
        "reference",
        "source.txt",
        " ".join(source_tokens),
    )

    matches = near_token_matches([candidate], [source])

    assert len(matches) == 1
    assert matches[0].tokens == 120
    assert matches[0].similarity >= 0.85


def test_punctuation_separator_is_not_a_lexical_token_match() -> None:
    separator = "# " + "- " * 100
    candidate = TextDocument("chemsmart", "agent.py", separator)
    source = TextDocument("reference", "source.py", separator)

    assert exact_token_matches([candidate], [source]) == []


def test_frozen_provenance_report_covers_all_production_files() -> None:
    path = REPO_ROOT / "docs/review/agent-provenance-audit.json"
    report = json.loads(path.read_text(encoding="utf-8"))

    assert report["audit"]["production_python_files"] == 151
    assert report["audit"]["candidate_documents"] == 176
    assert sum(row["lines"] for row in report["lineage"]) == 44691
    assert all(row["introduction"] for row in report["lineage"])
    assert all(row["introduced_after_boundary"] for row in report["lineage"])
    assert report["matches"] == {
        "exact_lines": [],
        "exact_tokens": [],
        "near_tokens": [],
    }
