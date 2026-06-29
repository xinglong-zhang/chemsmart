from __future__ import annotations

import pytest

from chemsmart.agent.error_summary import summarize_log


@pytest.mark.parametrize(
    ("kind", "line"),
    [
        (
            "oom_killed",
            "slurmstepd: error: Detected 1 oom-kill event(s) in StepId=7.0",
        ),
        (
            "walltime_exceeded",
            "slurmstepd: error: *** JOB 123 DUE TO TIME LIMIT ***",
        ),
        (
            "missing_module",
            "module: command not found",
        ),
        (
            "node_failure",
            "Job requeued because of node failure on node042",
        ),
        (
            "scheduler_reject",
            "sbatch: error: Batch job submission failed: Invalid account",
        ),
        (
            "segfault",
            "Segmentation fault (core dumped)",
        ),
    ],
)
def test_summarize_log_detects_each_supported_kind(kind, line):
    result = summarize_log(line)

    assert len(result) == 1
    assert result[0].kind == kind
    assert result[0].line == line
    assert result[0].line_no == 1
    assert result[0].severity == "error"


def test_summarize_log_returns_empty_for_empty_text():
    assert summarize_log("") == []


def test_summarize_log_handles_multi_kind_text_and_dedupes_nearby_lines():
    text = "\n".join(
        [
            "module: command not found",
            "ERROR: Unable to locate a modulefile for 'orca'",
            "launch failed requeued held",
            "Segmentation fault (core dumped)",
        ]
    )

    result = summarize_log(text)

    assert [item.kind for item in result] == [
        "missing_module",
        "node_failure",
        "segfault",
    ]
    assert [item.line_no for item in result] == [1, 3, 4]


def test_summarize_log_honors_max_signatures_cap():
    text = "\n".join(
        [
            "Segmentation fault (core dumped)",
            "line 2",
            "line 3",
            "line 4",
            "slurmstepd: error: Detected 1 oom-kill event(s)",
        ]
    )

    result = summarize_log(text, max_signatures=1)

    assert len(result) == 1
    assert result[0].kind == "segfault"
