from __future__ import annotations

from chemsmart.agent.permissions import RuntimePermissionMode
from chemsmart.agent.registry import ToolRegistry
from chemsmart.agent.tool_protocol import RuntimeToolMetadata
from chemsmart.agent.tools_hpc import scheduler_query
from chemsmart.agent.transport import ExecResult, MockExecTransport


def test_scheduler_query_normalizes_slurm_job(monkeypatch):
    transport = MockExecTransport(
        scripted_results=[
            ExecResult(
                returncode=0,
                stdout=(
                    "12345|calc-opt|RUNNING|alice|debug|01:02:03|02:00:00|"
                    "node001\n"
                ),
                stderr="",
                command="",
                server="cluster-a",
                duration_s=0.1,
            )
        ]
    )
    monkeypatch.setattr(
        "chemsmart.agent.tools_hpc._TRANSPORT_FACTORY",
        lambda: transport,
    )

    result = scheduler_query(
        server="cluster-a",
        scheduler="slurm",
        job_id="12345",
    )

    assert result == {
        "job_id": "12345",
        "state": "RUNNING",
        "queue": "debug",
        "user": "alice",
        "name": "calc-opt",
        "runtime_s": 3723,
        "node": "node001",
        "scheduler": "slurm",
        "raw": {
            "job_id": "12345",
            "name": "calc-opt",
            "state": "RUNNING",
            "user": "alice",
            "queue": "debug",
            "runtime": "01:02:03",
            "limit": "02:00:00",
            "node_or_reason": "node001",
        },
    }
    assert transport.calls == [
        {
            "command": ("squeue -h -j 12345 -o " "'%A|%j|%T|%u|%P|%M|%l|%R'"),
            "server": "cluster-a",
            "timeout_s": 15,
        }
    ]


def test_scheduler_query_normalizes_pbs_job(monkeypatch):
    transport = MockExecTransport(
        scripted_results=[
            ExecResult(
                returncode=0,
                stdout=(
                    '{"Jobs":{"221.server":{"job_state":"R","queue":"workq",'
                    '"Job_Owner":"alice@login","Job_Name":"calc",'
                    '"resources_used.walltime":"00:10:05",'
                    '"exec_host":"node002/0*4"}}}'
                ),
                stderr="",
                command="",
                server="cluster-a",
                duration_s=0.1,
            )
        ]
    )
    monkeypatch.setattr(
        "chemsmart.agent.tools_hpc._TRANSPORT_FACTORY",
        lambda: transport,
    )

    result = scheduler_query(
        server="cluster-a",
        scheduler="pbs",
        job_id="221.server",
    )

    assert result["state"] == "RUNNING"
    assert result["queue"] == "workq"
    assert result["user"] == "alice"
    assert result["name"] == "calc"
    assert result["runtime_s"] == 605
    assert result["node"] == "node002"
    assert result["scheduler"] == "pbs"


def test_scheduler_query_normalizes_sge_job(monkeypatch):
    transport = MockExecTransport(
        scripted_results=[
            ExecResult(
                returncode=0,
                stdout=(
                    "job_number: 555\n"
                    "job_name: calc\n"
                    "owner: alice\n"
                    "state: r\n"
                    "master_queue: all.q@node003\n"
                    "running_time: 00:05:00\n"
                ),
                stderr="",
                command="",
                server="cluster-a",
                duration_s=0.1,
            )
        ]
    )
    monkeypatch.setattr(
        "chemsmart.agent.tools_hpc._TRANSPORT_FACTORY",
        lambda: transport,
    )

    result = scheduler_query(
        server="cluster-a",
        scheduler="sge",
        job_id="555",
    )

    assert result == {
        "job_id": "555",
        "state": "RUNNING",
        "queue": "all.q",
        "user": "alice",
        "name": "calc",
        "runtime_s": 300,
        "node": "node003",
        "scheduler": "sge",
        "raw": {
            "job_number": "555",
            "job_name": "calc",
            "owner": "alice",
            "state": "r",
            "master_queue": "all.q@node003",
            "running_time": "00:05:00",
        },
    }


def test_scheduler_query_normalizes_lsf_job(monkeypatch):
    transport = MockExecTransport(
        scripted_results=[
            ExecResult(
                returncode=0,
                stdout=(
                    '{"RECORDS":[{"JOBID":"777","STAT":"RUN",'
                    '"QUEUE":"normal","USER":"alice","JOB_NAME":"calc",'
                    '"RUN_TIME":"75","EXEC_HOST":"node004*8"}]}'
                ),
                stderr="",
                command="",
                server="cluster-a",
                duration_s=0.1,
            )
        ]
    )
    monkeypatch.setattr(
        "chemsmart.agent.tools_hpc._TRANSPORT_FACTORY",
        lambda: transport,
    )

    result = scheduler_query(
        server="cluster-a",
        scheduler="lsf",
        job_id="777",
    )

    assert result["job_id"] == "777"
    assert result["state"] == "RUNNING"
    assert result["queue"] == "normal"
    assert result["user"] == "alice"
    assert result["name"] == "calc"
    assert result["runtime_s"] == 75
    assert result["node"] == "node004"
    assert result["scheduler"] == "lsf"


def test_scheduler_query_rejects_unknown_scheduler():
    result = scheduler_query(server="cluster-a", scheduler="flux")

    assert result == {
        "error": "unknown_scheduler",
        "supported": ["slurm", "pbs", "sge", "lsf"],
    }


def test_scheduler_query_rejects_invalid_job_id(monkeypatch):
    transport = MockExecTransport()
    monkeypatch.setattr(
        "chemsmart.agent.tools_hpc._TRANSPORT_FACTORY",
        lambda: transport,
    )

    result = scheduler_query(
        server="cluster-a",
        scheduler="slurm",
        job_id="12345;rm",
    )

    assert result["error"] == "invalid_job_id"
    assert "probe slot" in result["message"].lower() or "forbidden" in (
        result["message"].lower()
    )
    assert transport.calls == []


def test_scheduler_query_normalizes_slurm_queue(monkeypatch):
    transport = MockExecTransport(
        scripted_results=[
            ExecResult(
                returncode=0,
                stdout=(
                    "PartitionName=debug Default=YES State=UP "
                    "TotalNodes=4 TotalCPUs=128 DefaultTime=01:00:00 "
                    "MaxTime=02:00:00\n"
                ),
                stderr="",
                command="",
                server="cluster-a",
                duration_s=0.1,
            )
        ]
    )
    monkeypatch.setattr(
        "chemsmart.agent.tools_hpc._TRANSPORT_FACTORY",
        lambda: transport,
    )

    result = scheduler_query(server="cluster-a", scheduler="slurm")

    assert result["scheduler"] == "slurm"
    assert result["partition_or_queue"] == "debug"
    assert result["nodes"] == 4
    assert result["total_cpus"] == 128
    assert result["raw"]["queues"][0]["name"] == "debug"


def test_scheduler_query_registry_metadata_and_runtime_modes():
    registry = ToolRegistry.default()
    tool = registry.get_tool("scheduler_query")

    assert tool is not None
    assert tool.metadata == RuntimeToolMetadata(
        read_only=True,
        ui_summary_template="Scheduler query {scheduler} on {server}",
    )

    for mode in (
        RuntimePermissionMode.READ_ONLY,
        RuntimePermissionMode.ACCEPT_EDITS,
        RuntimePermissionMode.BYPASS,
    ):
        names = {item.name for item in registry.assemble_tool_pool(mode)}
        assert "scheduler_query" in names

    plan_names = {
        item.name
        for item in registry.assemble_tool_pool(RuntimePermissionMode.PLAN)
    }
    assert "scheduler_query" not in plan_names
