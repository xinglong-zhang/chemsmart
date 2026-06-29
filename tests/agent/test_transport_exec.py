from __future__ import annotations

from chemsmart.agent.transport import MockExecTransport, SubmitTransport


def test_mock_exec_transport_records_exec_calls():
    transport = MockExecTransport()

    result = transport.exec(
        command="sinfo --json",
        server="cluster-a",
        timeout_s=9,
    )

    assert transport.calls == [
        {
            "command": "sinfo --json",
            "server": "cluster-a",
            "timeout_s": 9,
        }
    ]
    assert result.command == "sinfo --json"
    assert result.server == "cluster-a"
    assert result.returncode == 0


def test_submit_transport_subclasses_remain_instantiable_without_exec():
    class DummySubmitTransport(SubmitTransport):
        def submit(self, script_path, working_dir, server):
            return {
                "submitted_at_path": script_path,
                "command_executed": "qsub job.sh",
                "returncode": 0,
                "job_id": None,
            }

    transport = DummySubmitTransport()

    result = transport.submit("job.sh", ".", None)

    assert result["command_executed"] == "qsub job.sh"
