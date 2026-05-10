from __future__ import annotations

from chemsmart.agent.transport import build_submit_invocation
from chemsmart.agent.wizard.orchestrator import run_wizard
from chemsmart.agent.wizard.probe import ProbeResult
from chemsmart.agent.wizard.write import write_server_yaml
from chemsmart.settings.server import Server


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


def _result(
    command, *, mode="local", host=None, stdout="", stderr="", returncode=0
):
    return ProbeResult(
        command=command,
        mode=mode,
        host=host,
        returncode=returncode,
        stdout=stdout,
        stderr=stderr,
        duration_s=0.0,
        truncated=False,
    )


def test_mode_a_round_trip(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))
    runner = StubRunner(
        local_results={
            ("printenv",): _result(
                "printenv",
                stdout="SLURM_CLUSTER_NAME=perlmutter\nPATH=/usr/bin\n",
            ),
            ("sinfo", "--json"): _result(
                "sinfo --json", stdout='{"partitions": []}'
            ),
            (
                "scontrol",
                "show",
                "partition",
                "--oneliner",
            ): _result(
                "scontrol show partition --oneliner",
                stdout=(
                    "PartitionName=debug Default=YES MaxTime=2-00:00:00 "
                    "DefaultTime=08:00:00 DefMemPerNode=65536 "
                    "MaxCPUsPerNode=32 State=UP TRES=cpu=32,gres/gpu=0\n"
                ),
            ),
            ("type", "module"): _result("type module", returncode=1),
            ("which", "module"): _result("which module", returncode=1),
            ("printenv", "CONDA_PREFIX"): _result(
                "printenv CONDA_PREFIX",
                stdout="\n",
                returncode=1,
            ),
            ("conda", "info", "--base"): _result(
                "conda info --base",
                returncode=1,
            ),
            ("command", "-v", "g16"): _result(
                "command -v g16",
                stdout="/apps/bin/g16\n",
            ),
            ("readlink", "-f", "/apps/bin/g16"): _result(
                "readlink -f /apps/bin/g16",
                stdout="/apps/gaussian/g16\n",
            ),
            ("command", "-v", "orca"): _result(
                "command -v orca",
                returncode=1,
            ),
            ("which", "orca"): _result("which orca", returncode=1),
            ("command", "-v", "nciplot"): _result(
                "command -v nciplot",
                returncode=1,
            ),
            ("which", "nciplot"): _result("which nciplot", returncode=1),
            ("module", "-t", "avail"): _result(
                "module -t avail",
                returncode=1,
            ),
            ("printenv", "SCRATCH"): _result(
                "printenv SCRATCH",
                stdout="/scratch/user\n",
            ),
            ("printenv", "WORK"): _result(
                "printenv WORK",
                stdout="\n",
                returncode=1,
            ),
            ("printenv", "TMPDIR"): _result(
                "printenv TMPDIR",
                stdout="\n",
                returncode=1,
            ),
            ("test", "-w", "/scratch/user"): _result("test -w /scratch/user"),
            ("printenv", "SBATCH_ACCOUNT"): _result(
                "printenv SBATCH_ACCOUNT",
                stdout="chem123\n",
            ),
            ("printenv", "SLURM_ACCOUNT"): _result(
                "printenv SLURM_ACCOUNT",
                stdout="\n",
                returncode=1,
            ),
            ("printenv", "CHEMSMART_BIN"): _result(
                "printenv CHEMSMART_BIN",
                stdout="\n",
                returncode=1,
            ),
        }
    )

    outcome = run_wizard(runner, server_name="perlmutter", write=False)

    written_path = write_server_yaml("perlmutter", outcome.plan.text)
    server = Server.from_yaml(written_path)
    invocation, _ = build_submit_invocation(
        script_path="/tmp/job.sh",
        working_dir="/tmp",
        server=server,
    )

    assert outcome.validation.ok is True
    assert server.kwargs["HOST"] == "localhost"
    assert "ssh" not in invocation
    assert invocation == ["sbatch", "/tmp/job.sh"]


def test_mode_b_round_trip(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))
    host = "cluster.example.edu"
    runner = StubRunner(
        local_results={
            ("printenv",): _result("printenv", stdout="PATH=/usr/bin\n"),
            ("sinfo",): _result("sinfo", returncode=1),
            ("qstat",): _result("qstat", returncode=1),
            ("bqueues",): _result("bqueues", returncode=1),
            ("qconf",): _result("qconf", returncode=1),
        },
        ssh_results={
            (host, "sinfo --json"): _result(
                "sinfo --json",
                mode="ssh",
                host=host,
                stdout='{"partitions": []}',
            ),
            (host, "scontrol show partition --oneliner"): _result(
                "scontrol show partition --oneliner",
                mode="ssh",
                host=host,
                stdout=(
                    "PartitionName=normal Default=YES MaxTime=2-00:00:00 "
                    "DefaultTime=08:00:00 DefMemPerNode=65536 "
                    "MaxCPUsPerNode=32 State=UP TRES=cpu=32,gres/gpu=0\n"
                ),
            ),
            (host, "type module"): _result(
                "type module",
                mode="ssh",
                host=host,
                returncode=1,
            ),
            (host, "which module"): _result(
                "which module",
                mode="ssh",
                host=host,
                returncode=1,
            ),
            (host, "printenv CONDA_PREFIX"): _result(
                "printenv CONDA_PREFIX",
                mode="ssh",
                host=host,
                stdout="\n",
                returncode=1,
            ),
            (host, "conda info --base"): _result(
                "conda info --base",
                mode="ssh",
                host=host,
                returncode=1,
            ),
            (host, "command -v g16"): _result(
                "command -v g16",
                mode="ssh",
                host=host,
                stdout="/apps/bin/g16\n",
            ),
            (host, "readlink -f /apps/bin/g16"): _result(
                "readlink -f /apps/bin/g16",
                mode="ssh",
                host=host,
                stdout="/apps/gaussian/g16\n",
            ),
            (host, "command -v orca"): _result(
                "command -v orca",
                mode="ssh",
                host=host,
                returncode=1,
            ),
            (host, "command -v nciplot"): _result(
                "command -v nciplot",
                mode="ssh",
                host=host,
                returncode=1,
            ),
            (host, "module -t avail"): _result(
                "module -t avail",
                mode="ssh",
                host=host,
                returncode=1,
            ),
            (host, "printenv SCRATCH"): _result(
                "printenv SCRATCH",
                mode="ssh",
                host=host,
                stdout="/scratch/user\n",
            ),
            (host, "printenv WORK"): _result(
                "printenv WORK",
                mode="ssh",
                host=host,
                stdout="\n",
                returncode=1,
            ),
            (host, "printenv TMPDIR"): _result(
                "printenv TMPDIR",
                mode="ssh",
                host=host,
                stdout="\n",
                returncode=1,
            ),
            (host, "test -w /scratch/user"): _result(
                "test -w /scratch/user",
                mode="ssh",
                host=host,
            ),
            (host, "printenv SBATCH_ACCOUNT"): _result(
                "printenv SBATCH_ACCOUNT",
                mode="ssh",
                host=host,
                stdout="chem123\n",
            ),
            (host, "printenv SLURM_ACCOUNT"): _result(
                "printenv SLURM_ACCOUNT",
                mode="ssh",
                host=host,
                stdout="\n",
                returncode=1,
            ),
        },
    )

    outcome = run_wizard(
        runner,
        server_name="remote-perlmutter",
        ssh_host_hint=host,
        write=False,
    )

    written_path = write_server_yaml("remote-perlmutter", outcome.plan.text)
    server = Server.from_yaml(written_path)
    invocation, _ = build_submit_invocation(
        script_path="/tmp/job.sh",
        working_dir="/tmp",
        server=server,
    )

    assert outcome.validation.ok is True
    assert server.kwargs["HOST"] == host
    assert invocation[0] == "ssh"
    assert invocation[1] == host


def test_legacy_yaml_no_host_falls_to_filename(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))
    target = tmp_path / ".chemsmart" / "server" / "legacy-cluster.yaml"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        """SERVER:
  SCHEDULER: SLURM
  QUEUE_NAME: debug
  NUM_HOURS: 8
  MEM_GB: 64
  NUM_CORES: 16
  SUBMIT_COMMAND: sbatch
GAUSSIAN:
  EXEFOLDER: /apps/gaussian
  LOCAL_RUN: true
  SCRATCH: true
  CONDA_ENV: ""
  MODULES: ""
  ENVARS: ""
""",
        encoding="utf-8",
    )

    server = Server.from_yaml(str(target))
    invocation, _ = build_submit_invocation(
        script_path="/tmp/job.sh",
        working_dir="/tmp",
        server=server,
    )

    assert invocation[0] == "ssh"
    assert invocation[1] == "legacy-cluster"
