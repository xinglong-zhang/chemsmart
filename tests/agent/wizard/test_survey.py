import pytest

from chemsmart.agent.wizard import (
    AmbiguousSchedulerError,
    ScheduleSurvey,
    Topology,
    run_schedule_survey,
)
from chemsmart.agent.wizard.parsers import QueueFacts
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


def _result(command, stdout, returncode=0, mode="local", host=None):
    return ProbeResult(
        command=command,
        mode=mode,
        host=host,
        returncode=returncode,
        stdout=stdout,
        stderr="",
        duration_s=0.0,
        truncated=False,
    )


def test_run_schedule_survey_slurm_only_path():
    runner = StubRunner(
        local_results={
            ("sinfo", "--json"): _result(
                "sinfo --json",
                '{"partitions": [{"name": "debug", "default": true, '
                '"max_time": "02:00:00", "state": "UP"}]}',
            ),
            ("scontrol", "show", "partition", "--oneliner"): _result(
                "scontrol show partition --oneliner",
                "PartitionName=debug Default=YES MaxTime=02:00:00 "
                "DefaultTime=01:00:00 DefMemPerNode=65536 "
                "MaxCPUsPerNode=32 State=UP TRES=cpu=32,gres/gpu=0\n",
            ),
        }
    )

    assert run_schedule_survey(
        runner, Topology("A", "localhost", [])
    ) == ScheduleSurvey(
        scheduler="SLURM",
        submit_command="sbatch",
        queues=[
            QueueFacts(
                name="debug",
                default=True,
                max_walltime_hours=2,
                default_walltime_hours=1,
                default_mem_gb=64,
                default_cores=32,
                gpus_per_node=0,
                enabled=True,
                started=True,
            )
        ],
        chosen_queue="debug",
        evidence={
            "sinfo --json": "parsed",
            "scontrol show partition --oneliner": "parsed",
        },
    )


def test_run_schedule_survey_raises_on_ambiguous_scheduler():
    runner = StubRunner(
        local_results={
            ("sinfo", "--json"): _result(
                "sinfo --json",
                '{"partitions": [{"name": "debug", "default": true, '
                '"max_time": "02:00:00", "state": "UP"}]}',
            ),
            ("qstat", "-Q", "-f", "-F", "json"): _result(
                "qstat -Q -f -F json",
                '{"Queue": {"workq": {"default_queue": "True", '
                '"enabled": "True", "started": "True"}}}',
            ),
        }
    )

    with pytest.raises(AmbiguousSchedulerError):
        run_schedule_survey(runner, Topology("A", "localhost", []))
