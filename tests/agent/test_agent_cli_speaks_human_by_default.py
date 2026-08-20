"""The agent CLI speaks human by default; --json is the machine channel.

Digests live on in the durable records and the --json payloads; the default
stdout of plan/review/run never shows a bare 64-hex value.
"""

from __future__ import annotations

import re
from types import SimpleNamespace

from click.testing import CliRunner

from chemsmart.cli.agent import agent

_HEX = re.compile(r"[0-9a-f]{64}")


class _Result:
    terminal_state = "waiting_for_approval"
    successful_tool_calls = 27
    failed_tool_calls = 1
    final_text = "The workflow is ready. digest " + "a" * 64
    prepared_execution = object()

    def public_summary_json(self):
        return '{"machine": "payload", "sha": "' + "b" * 64 + '"}'


def _invoke_plan(tmp_path, monkeypatch, extra=()):
    workspace = tmp_path / "workspace"
    workspace.mkdir(exist_ok=True)
    secret = tmp_path / "secret.env"
    secret.write_text("KEY=value\n", encoding="utf-8")

    import chemsmart.agent.live_session as live_session

    monkeypatch.setattr(
        live_session, "run_live_agent_session", lambda **_kw: _Result()
    )
    return CliRunner().invoke(
        agent,
        [
            "plan",
            "--task",
            "preview this calculation",
            "--secret-file",
            str(secret),
            "--workspace",
            str(workspace),
            *extra,
        ],
    )


def test_plan_json_is_the_exact_machine_payload(tmp_path, monkeypatch):
    result = _invoke_plan(tmp_path, monkeypatch, extra=("--json",))

    assert result.exit_code == 0, result.output
    assert result.output.strip() == _Result().public_summary_json()


def test_plan_default_speaks_words_without_digests(tmp_path, monkeypatch):
    result = _invoke_plan(tmp_path, monkeypatch)

    assert result.exit_code == 0, result.output
    assert "outcome: ready for your review" in result.output
    assert "steps: 27 steps, 1 refused" in result.output
    assert "The workflow is ready." in result.output
    assert _HEX.search(result.output) is None


def test_review_summary_names_artifacts_and_drops_digests(
    tmp_path, monkeypatch
):
    review_file = tmp_path / "review.json"
    review_file.write_text("{}", encoding="utf-8")
    workspace = tmp_path / "workspace"
    workspace.mkdir()

    import chemsmart.agent.live_session as live_session

    monkeypatch.setattr(
        live_session,
        "inspect_workflow_execution_replay",
        lambda **_kw: {
            "review_sha256": "a" * 64,
            "task_spec_sha256": "b" * 64,
            "workflow_id": "water-sp",
            "node_count": 1,
            "non_executable_node_ids": [],
            "approved_artifacts_present": [
                "sp-water:project:" + "c" * 16,
                "sp-water:input:" + "d" * 16,
            ],
            "missing_approved_artifacts": [],
            "previously_consumed_approval_ids": [],
            "canonical_review": "{}",
        },
    )
    result = CliRunner().invoke(
        agent,
        [
            "review",
            "--review-file",
            str(review_file),
            "--workspace",
            str(workspace),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "sp-water · project settings" in result.output
    assert "sp-water · input geometry" in result.output
    assert "review_sha256" not in result.output
    assert "task_spec_sha256" not in result.output
    assert _HEX.search(result.output) is None
