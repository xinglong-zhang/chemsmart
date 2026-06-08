from __future__ import annotations

from chemsmart.agent.transport import _server_host, build_submit_invocation
from chemsmart.settings.server import Server


def _server(name: str, **server_kwargs) -> Server:
    return Server.from_dict(
        {
            "name": name,
            "SCHEDULER": "PBS",
            "SUBMIT_COMMAND": "qsub",
            **server_kwargs,
        }
    )


def test_server_host_prefers_explicit_host_field():
    server = _server(
        "/tmp/server/SLURM.yaml",
        HOST="login.cluster.example.edu",
    )
    command, executed = build_submit_invocation(
        script_path="/tmp/job.sh",
        working_dir="/tmp",
        server=server,
    )

    assert _server_host(server) == "login.cluster.example.edu"
    assert command == [
        "ssh",
        "login.cluster.example.edu",
        "cd /tmp && qsub /tmp/job.sh",
    ]
    assert (
        executed
        == 'ssh login.cluster.example.edu "cd /tmp && qsub /tmp/job.sh"'
    )


def test_server_host_falls_back_to_yaml_stem_when_host_missing():
    server = _server("/tmp/server/remote-hpc.yaml")

    assert _server_host(server) == "remote-hpc"


def test_build_submit_invocation_keeps_local_special_case():
    command, executed = build_submit_invocation(
        script_path="/tmp/job.sh",
        working_dir="/tmp",
        server=_server("/tmp/server/localhost.yaml"),
    )

    assert command == ["qsub", "/tmp/job.sh"]
    assert executed == "qsub /tmp/job.sh"
