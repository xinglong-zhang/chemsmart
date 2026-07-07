from __future__ import annotations

from chemsmart.agent.tui.tool_meta import format_assumptions_banner


def test_format_assumptions_banner_returns_none_without_entities():
    assert format_assumptions_banner(None, "ok") is None


def test_format_assumptions_banner_renders_scheduler_server_and_high_confidence():
    banner = format_assumptions_banner(
        {"last_scheduler": "pbs", "last_server": "chemnode1"},
        "ok",
    )

    assert banner is not None
    assert "PBS" in banner
    assert "chemnode1" in banner
    assert "conf=high" in banner


def test_format_assumptions_banner_renders_job_log_and_low_confidence():
    banner = format_assumptions_banner(
        {
            "last_job_id": "4.chemnode1",
            "last_log_path": "/scratch/run.log",
        },
        "error",
    )

    assert banner is not None
    assert "job 4.chemnode1" in banner
    assert "log run.log" in banner
    assert "conf=low" in banner


def test_format_assumptions_banner_maps_status_to_confidence():
    expectations = {
        "ok": "high",
        "partial": "med",
        "error": "low",
        "denied": "low",
        "skipped": "low",
        "ask_user": "low",
        None: "unknown",
    }

    for status, confidence in expectations.items():
        banner = format_assumptions_banner(
            {"last_server": "chemnode1"},
            status,
        )

        assert banner is not None
        assert f"conf={confidence}" in banner
