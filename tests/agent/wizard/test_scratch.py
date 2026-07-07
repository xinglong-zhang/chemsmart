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


def test_discover_scratch_prefers_scratch_env_var(monkeypatch):
    monkeypatch.delenv("TMPDIR", raising=False)
    monkeypatch.delenv("WORK", raising=False)
    monkeypatch.delenv("SCRATCH", raising=False)
    monkeypatch.delenv("HOME", raising=False)
    runner = StubRunner(
        local_results={
            ("printenv", "SCRATCH"): _result(
                "printenv SCRATCH",
                "/scratch/user\n",
            ),
            ("printenv", "WORK"): _result("printenv WORK", "/work/user\n"),
            ("printenv", "TMPDIR"): _result("printenv TMPDIR", "/tmp/user\n"),
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


def test_discover_scratch_falls_back_to_home_scratch(monkeypatch):
    monkeypatch.delenv("TMPDIR", raising=False)
    monkeypatch.delenv("WORK", raising=False)
    monkeypatch.delenv("SCRATCH", raising=False)
    monkeypatch.delenv("HOME", raising=False)
    runner = StubRunner(
        local_results={
            ("printenv", "SCRATCH"): _result("printenv SCRATCH", "\n", 1),
            ("printenv", "WORK"): _result("printenv WORK", "\n", 1),
            ("printenv", "TMPDIR"): _result("printenv TMPDIR", "\n", 1),
            ("printenv", "HOME"): _result("printenv HOME", "/home/tester\n"),
            (
                "test",
                "-d",
                "/home/tester/scratch",
                "-a",
                "-w",
                "/home/tester/scratch",
            ): _result(
                "test -d /home/tester/scratch -a -w /home/tester/scratch"
            ),
            ("test", "-w", "/home/tester/scratch"): _result(
                "test -w /home/tester/scratch"
            ),
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


def test_discover_scratch_returns_none_when_no_candidate_exists(monkeypatch):
    monkeypatch.delenv("TMPDIR", raising=False)
    monkeypatch.delenv("WORK", raising=False)
    monkeypatch.delenv("SCRATCH", raising=False)
    monkeypatch.delenv("HOME", raising=False)
    runner = StubRunner(
        local_results={
            ("printenv", "SCRATCH"): _result("printenv SCRATCH", "\n", 1),
            ("printenv", "WORK"): _result("printenv WORK", "\n", 1),
            ("printenv", "TMPDIR"): _result("printenv TMPDIR", "\n", 1),
            ("printenv", "HOME"): _result("printenv HOME", "/home/tester\n"),
            (
                "test",
                "-d",
                "/home/tester/scratch",
                "-a",
                "-w",
                "/home/tester/scratch",
            ): _result(
                "test -d /home/tester/scratch -a -w /home/tester/scratch",
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
