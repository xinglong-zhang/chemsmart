from __future__ import annotations

from chemsmart.agent.core import AgentSession, Step


class DummyJob:
    def __repr__(self) -> str:
        return "DummyJob()"


class RecordingRegistry:
    def __init__(self, results):
        self.results = dict(results)
        self.calls = []

    def call(self, name, args=None):
        args = dict(args or {})
        self.calls.append((name, args))
        return self.results[name]


def test_execute_step_puts_results_into_handle_store_and_resolves_handles(
    tmp_path,
):
    job = DummyJob()
    registry = RecordingRegistry(
        {
            "build_job": job,
            "dry_run_input": {
                "inputfile": str(tmp_path / "water.com"),
                "content": "#p opt\n",
            },
            "consume_legacy_ref": {"ok": True},
        }
    )
    session = AgentSession(registry=registry, session_root=tmp_path)
    session._start_new_session("prepare a job")

    result = session._execute_step(
        0,
        Step(tool="build_job", args={}, rationale="Create the job."),
        prior_results=[],
    )

    assert result is job
    entries = session.decision_log.read_all()
    tool_result = next(
        entry for entry in entries if entry["kind"] == "tool_result"
    )
    handle_id = tool_result["payload"]["handle_id"]
    assert handle_id.startswith("job_")
    assert session.handle_store.get(handle_id) is job
    assert session.handle_store.get_summary(handle_id)["type"] == "DummyJob"

    session._execute_step(
        1,
        Step(
            tool="dry_run_input",
            args={"job": handle_id},
            rationale="Render the input file.",
        ),
        prior_results=[],
    )

    assert registry.calls[1][0] == "dry_run_input"
    assert registry.calls[1][1]["job"] is job


def test_legacy_step_refs_still_resolve_with_handle_support(tmp_path):
    job = DummyJob()
    registry = RecordingRegistry({"consume_legacy_ref": {"ok": True}})
    session = AgentSession(registry=registry, session_root=tmp_path)
    session._start_new_session("reuse prior step output")

    session._execute_step(
        0,
        Step(
            tool="consume_legacy_ref",
            args={"job": "$step1"},
            rationale="Keep legacy step references working.",
        ),
        prior_results=[job],
    )

    assert registry.calls[0][1]["job"] is job
