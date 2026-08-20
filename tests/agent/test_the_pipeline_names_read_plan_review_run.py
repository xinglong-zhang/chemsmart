"""The public agent surface is one pipeline: plan, review, run, tui.

The production close renamed the file-based decision command (`replay` ->
`review`) and gave the provider-free executor the name a scientist reaches
for (`execute`, hidden -> `run`, public), while the old plan-shaped `run`
and the hidden legacy `approve` were removed. These tests pin the mapping of
each name to its authority and that the retired names are gone rather than
aliased.
"""

from types import SimpleNamespace

from click.testing import CliRunner

from chemsmart.cli.agent import agent


def test_review_records_the_decision_and_prints_the_run_step(
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
        lambda **_kwargs: {
            "review_sha256": "a" * 64,
            "workflow_id": "water-sp",
            "missing_approved_artifacts": [],
            "canonical_review": "{}",
        },
    )
    recorded = {}

    def fake_resolve(**kwargs):
        recorded.update(kwargs)
        return (
            SimpleNamespace(
                decision="approve",
                approval_id=kwargs["approval_id"],
                resolution_sha256="b" * 64,
                review_sha256="a" * 64,
            ),
            SimpleNamespace(bundle_sha256="c" * 64),
        )

    monkeypatch.setattr(
        live_session, "resolve_workflow_execution_review", fake_resolve
    )

    result = CliRunner().invoke(
        agent,
        [
            "review",
            "--review-file",
            str(review_file),
            "--workspace",
            str(workspace),
            "--decision",
            "approve",
            "--actor",
            "tester",
        ],
    )

    assert result.exit_code == 0, result.output
    assert recorded["decision"] == "approve"
    assert "chemsmart agent run --approval-file" in result.output
    assert "--run-directory" in result.output


def test_run_is_the_provider_free_executor(tmp_path, monkeypatch):
    approval = tmp_path / "bundle.json"
    approval.write_text("{}", encoding="utf-8")
    workspace = tmp_path / "workspace"
    workspace.mkdir()

    import chemsmart.agent.executor as executor

    calls = {}

    def fake_execute(**kwargs):
        calls.update(kwargs)
        return SimpleNamespace(
            workflow_id="water-sp",
            plan_sha256="a" * 64,
            approval_sha256="b" * 64,
            status="completed",
            provider_calls=0,
            non_executable_node_ids=(),
            run_directory=str(tmp_path / "run"),
            analysis_status="",
            analysis_report_path="",
            analysis_completion_receipt_sha256s=(),
            analysis_nodes=(),
            nodes=(),
        )

    monkeypatch.setattr(executor, "execute_approved_workflow", fake_execute)

    result = CliRunner().invoke(
        agent,
        [
            "run",
            "--approval-file",
            str(approval),
            "--workspace",
            str(workspace),
            "--run-directory",
            str(tmp_path / "run"),
            "--json",
        ],
    )

    assert result.exit_code == 0, result.output
    assert calls["approval_file"] == approval
    assert '"provider_calls": 0' in result.output

    human = CliRunner().invoke(
        agent,
        [
            "run",
            "--approval-file",
            str(approval),
            "--workspace",
            str(workspace),
            "--run-directory",
            str(tmp_path / "run"),
        ],
    )

    assert human.exit_code == 0, human.output
    assert "workflow: water-sp" in human.output
    assert "status: completed" in human.output
    assert "plan_sha256" not in human.output
    import re

    assert re.search(r"[0-9a-f]{64}", human.output) is None


def test_the_retired_command_names_are_gone_not_aliased():
    for retired in ("replay", "approve", "execute"):
        result = CliRunner().invoke(agent, [retired, "--help"])
        assert result.exit_code == 2, (retired, result.output)
