
from chemsmart.agent.wizard import NoTargetError, Topology, detect_topology
from chemsmart.agent.wizard.probe import ProbeResult


class StubRunner:
    def __init__(self, local_results=None, ssh_results=None):
        self.local_results = local_results or {}
        self.ssh_results = ssh_results or {}

    def run_local(self, command, timeout_s=15):
        return self.local_results.get(
            tuple(command),
            ProbeResult(
                command=" ".join(command),
                mode="local",
                host=None,
                returncode=1,
                stdout="",
                stderr="",
                duration_s=0.0,
                truncated=False,
            ),
        )

    def run_ssh(self, host, command, timeout_s=15):
        return self.ssh_results.get(
            (host, command),
            ProbeResult(
                command=command,
                mode="ssh",
                host=host,
                returncode=1,
                stdout="",
                stderr="",
                duration_s=0.0,
                truncated=False,
            ),
        )


def test_detect_topology_mode_a_on_local_scheduler_probe():
    runner = StubRunner(
        local_results={
            ("env",): ProbeResult(
                command="env",
                mode="local",
                host=None,
                returncode=0,
                stdout="PATH=/bin\n",
                stderr="",
                duration_s=0.0,
                truncated=False,
            ),
            ("sinfo",): ProbeResult(
                command="sinfo",
                mode="local",
                host=None,
                returncode=0,
                stdout="ok",
                stderr="",
                duration_s=0.0,
                truncated=False,
            ),
        }
    )

    assert detect_topology(runner) == Topology(
        mode="A", host="localhost", evidence=["local:sinfo"]
    )


def test_detect_topology_mode_b_on_hint_when_local_probes_fail():
    runner = StubRunner(
        local_results={
            ("env",): ProbeResult(
                command="env",
                mode="local",
                host=None,
                returncode=0,
                stdout="PATH=/bin\n",
                stderr="",
                duration_s=0.0,
                truncated=False,
            )
        }
    )

    assert detect_topology(runner, ssh_host_hint="cluster") == Topology(
        mode="B", host="cluster", evidence=["ssh_host_hint:cluster"]
    )


def test_detect_topology_raises_when_no_target_is_available():
    runner = StubRunner(
        local_results={
            ("env",): ProbeResult(
                command="env",
                mode="local",
                host=None,
                returncode=0,
                stdout="PATH=/bin\n",
                stderr="",
                duration_s=0.0,
                truncated=False,
            )
        }
    )

    try:
        detect_topology(runner)
    except NoTargetError:
        assert True
    else:
        assert False
