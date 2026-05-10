from chemsmart.agent.wizard import ProjectFinding, Topology, discover_project
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


def test_discover_project_prefers_sbatch_account_for_slurm():
    runner = StubRunner(
        local_results={
            (
                "printf",
                "%s\\n",
                "$SBATCH_ACCOUNT",
                "$SLURM_ACCOUNT",
            ): _result(
                "printf %s\\n $SBATCH_ACCOUNT $SLURM_ACCOUNT",
                "chem123\nbackup456\n",
            )
        }
    )

    assert discover_project(
        runner,
        Topology("A", "localhost", []),
        "SLURM",
    ) == ProjectFinding(
        project="chem123",
        source="env:SBATCH_ACCOUNT",
        candidates=["chem123", "backup456"],
    )


def test_discover_project_uses_sacctmgr_default_account_for_slurm():
    runner = StubRunner(
        local_results={
            (
                "printf",
                "%s\\n",
                "$SBATCH_ACCOUNT",
                "$SLURM_ACCOUNT",
            ): _result(
                "printf %s\\n $SBATCH_ACCOUNT $SLURM_ACCOUNT",
                "\n\n",
            ),
            (
                "sacctmgr",
                "-n",
                "-p",
                "show",
                "user",
                "$USER",
                "format=DefaultAccount,Account",
            ): _result(
                "sacctmgr -n -p show user $USER format=DefaultAccount,Account",
                "chem-main|chem-main,chem-alt\n",
            ),
        }
    )

    assert discover_project(
        runner,
        Topology("A", "localhost", []),
        "SLURM",
    ) == ProjectFinding(
        project="chem-main",
        source="sacctmgr",
        candidates=["chem-main"],
    )


def test_discover_project_prefers_pbs_account_for_pbs():
    runner = StubRunner(
        local_results={
            ("printf", "%s\\n", "$PBS_ACCOUNT"): _result(
                "printf %s\\n $PBS_ACCOUNT",
                "pbs-chem\n",
            )
        }
    )

    assert discover_project(
        runner,
        Topology("A", "localhost", []),
        "PBS",
    ) == ProjectFinding(
        project="pbs-chem",
        source="env:PBS_ACCOUNT",
        candidates=["pbs-chem"],
    )


def test_discover_project_returns_none_when_nothing_is_available():
    runner = StubRunner(
        local_results={
            (
                "printf",
                "%s\\n",
                "$SBATCH_ACCOUNT",
                "$SLURM_ACCOUNT",
            ): _result(
                "printf %s\\n $SBATCH_ACCOUNT $SLURM_ACCOUNT",
                "\n\n",
            ),
            (
                "sacctmgr",
                "-n",
                "-p",
                "show",
                "user",
                "$USER",
                "format=DefaultAccount,Account",
            ): _result(
                "sacctmgr -n -p show user $USER format=DefaultAccount,Account",
                returncode=1,
            ),
            (
                "sacctmgr",
                "-n",
                "-p",
                "show",
                "user",
                "hongjiseung",
                "format=DefaultAccount,Account",
            ): _result(
                "sacctmgr -n -p show user hongjiseung format=DefaultAccount,Account",
                returncode=1,
            ),
            ("groups",): _result("groups", "", returncode=1),
        }
    )

    assert discover_project(
        runner,
        Topology("A", "localhost", []),
        "SLURM",
    ) == ProjectFinding(
        project=None,
        source="none",
        candidates=[],
    )
