from __future__ import annotations

from click.testing import CliRunner

from chemsmart.agent.cli import agent, sanitize_inline_cli_output


def test_agent_wizard_refresh_prints_rich_summary(monkeypatch):
    monkeypatch.setattr(
        "chemsmart.agent.cli.run_wizard_refresh",
        lambda name, force=False: {
            "server_name": name,
            "cache_path": "/tmp/perlmutter.cache.json",
            "status": "fresh",
            "host": "login.cluster",
            "mode": "ssh",
            "scheduler": "SLURM",
            "probed_at": "2026-05-10T00:00:00Z",
            "node_summary": {
                "selected_queue": "debug",
                "cpu": 16,
                "mem_gb": 64,
                "gpu": 0,
                "queue_count": 1,
                "enabled_queue_count": 1,
                "started_queue_count": 1,
                "gpu_queue_count": 0,
                "project": "chem-123",
                "scratch_dir": "/scratch/user",
                "scratch_writable": True,
            },
            "program_candidates": {
                "gaussian": {
                    "source": "path",
                    "exefolder": "/apps/gaussian",
                    "module_candidates": [],
                }
            },
            "last_error": None,
        },
    )

    result = CliRunner().invoke(
        agent,
        ["wizard-refresh", "perlmutter", "--force"],
        catch_exceptions=False,
    )

    output = sanitize_inline_cli_output(result.output)
    assert result.exit_code == 0
    assert "Wizard refresh: perlmutter" in output
    assert "fresh" in output
    assert "login.cluster" in output
    assert "debug" in output
    assert "/apps/gaussian" in output
