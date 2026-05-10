from chemsmart.agent.wizard import ScratchFinding, Topology, discover_scratch
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


def _result(command, stdout="", returncode=0):
    return ProbeResult(
        command=command,
        mode="local",
        host=None,
        returncode=returncode,
        stdout=stdout,
        stderr="",
        duration_s=0.0,
        truncated=False,
    )


def test_discover_scratch_prefers_scratch_env_var():
    runner = StubRunner(
        local_results={
            ("printf", "%s\\n", "$SCRATCH", "$WORK", "$TMPDIR"): _result(
                "printf %s\\n $SCRATCH $WORK $TMPDIR",
                "/scratch/user\n/work/user\n/tmp/user\n",
            ),
            ("test", "-d", "~/scratch", "-a", "-w", "~/scratch"): _result(
                "test -d ~/scratch -a -w ~/scratch",
                returncode=1,
            ),
            ("test", "-w", "/scratch/user"): _result("test -w /scratch/user"),
        }
    )

    assert discover_scratch(runner, Topology("A", "localhost", [])) == (
        ScratchFinding(
            path="/scratch/user",
            source="env:SCRATCH",
            writable=True,
            candidates=[
                ("env:SCRATCH", "/scratch/user"),
                ("env:WORK", "/work/user"),
                ("env:TMPDIR", "/tmp/user"),
            ],
        )
    )


def test_discover_scratch_falls_back_to_home_scratch():
    runner = StubRunner(
        local_results={
            ("printf", "%s\\n", "$SCRATCH", "$WORK", "$TMPDIR"): _result(
                "printf %s\\n $SCRATCH $WORK $TMPDIR",
                "\n\n\n",
            ),
            ("test", "-d", "~/scratch", "-a", "-w", "~/scratch"): _result(
                "test -d ~/scratch -a -w ~/scratch"
            ),
            ("test", "-w", "~/scratch"): _result("test -w ~/scratch"),
        }
    )

    assert discover_scratch(runner, Topology("A", "localhost", [])) == (
        ScratchFinding(
            path="~/scratch",
            source="home:~/scratch",
            writable=True,
            candidates=[("home:~/scratch", "~/scratch")],
        )
    )


def test_discover_scratch_returns_none_when_no_candidate_exists():
    runner = StubRunner(
        local_results={
            ("printf", "%s\\n", "$SCRATCH", "$WORK", "$TMPDIR"): _result(
                "printf %s\\n $SCRATCH $WORK $TMPDIR",
                "\n\n\n",
            ),
            ("test", "-d", "~/scratch", "-a", "-w", "~/scratch"): _result(
                "test -d ~/scratch -a -w ~/scratch",
                returncode=1,
            ),
        }
    )

    assert discover_scratch(runner, Topology("A", "localhost", [])) == (
        ScratchFinding(
            path=None,
            source="none",
            writable=False,
            candidates=[],
        )
    )
