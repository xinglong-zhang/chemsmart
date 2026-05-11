from pathlib import Path

import pytest

from chemsmart.agent.wizard import (
    AmbiguousSchedulerError,
    ScheduleSurvey,
    Topology,
    run_schedule_survey,
)
from chemsmart.agent.wizard.parsers import QueueFacts
from chemsmart.agent.wizard.probe import ProbeResult

FIXTURE_DIR = Path(__file__).parent / "fixtures"
SLURM_FIXTURE_DIR = FIXTURE_DIR / "slurm"


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
    node_payload = (
        SLURM_FIXTURE_DIR / "scontrol_show_node_768.txt"
    ).read_text()
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
                "MaxCPUsPerNode=32 Nodes=chemslurm1 "
                "State=UP TRES=cpu=32,gres/gpu=0\n",
            ),
            ("scontrol", "show", "node", "chemslurm1"): _result(
                "scontrol show node chemslurm1",
                node_payload,
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
                mem_mb=768,
            )
        ],
        chosen_queue="debug",
        evidence={
            "sinfo --json": "parsed",
            "scontrol show partition --oneliner": "parsed",
            "scontrol show node chemslurm1": "parsed",
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


def test_run_schedule_survey_merges_sge_total_slots():
    runner = StubRunner(
        local_results={
            ("qconf", "-sql"): _result(
                "qconf -sql",
                "all.q\n20core.q\n40core.q\n",
            ),
            ("qstat", "-g", "c"): _result(
                "qstat -g c",
                "CLUSTER QUEUE                   CQLOAD   USED   RES  AVAIL  TOTAL "
                "aoACDS  cdsuE\n"
                "--------------------------------------------------------------------------------\n"
                "all.q                            0.00      0     0      0      0      0      0\n"
                "20core.q                         0.10     70     0     30    100      0      0\n"
                "40core.q                         0.20    320     0    240    560      0      0\n",
            ),
            ("qhost",): _result("qhost", returncode=1, stdout=""),
            ("qconf", "-sq", "all.q"): _result(
                "qconf -sq all.q",
                "qname all.q\nslots 1\nh_rt 86400\nh_vmem INFINITY\n"
                "state enabled\n",
            ),
            ("qconf", "-sq", "20core.q"): _result(
                "qconf -sq 20core.q",
                "qname 20core.q\nslots 20\nh_rt 86400\nh_vmem 62.7G\n"
                "state enabled\n",
            ),
            ("qconf", "-sq", "40core.q"): _result(
                "qconf -sq 40core.q",
                "qname 40core.q\nslots 40\nh_rt 86400\nh_vmem 187.4G\n"
                "state enabled\n",
            ),
        }
    )

    survey = run_schedule_survey(runner, Topology("A", "localhost", []))

    assert survey.scheduler == "SGE"
    assert survey.chosen_queue == "20core.q"
    assert survey.evidence["qstat -g c"] == "parsed"
    assert {queue.name: queue.slots_total for queue in survey.queues} == {
        "all.q": 0,
        "20core.q": 100,
        "40core.q": 560,
    }


def test_run_schedule_survey_slurm_uses_node_cpus_when_sinfo_json_fails():
    node_payload = (
        SLURM_FIXTURE_DIR / "scontrol_show_node_chemnode2.txt"
    ).read_text()
    runner = StubRunner(
        local_results={
            ("sinfo", "--json"): _result(
                "sinfo --json",
                "",
                returncode=1,
            ),
            ("scontrol", "show", "partition", "--oneliner"): _result(
                "scontrol show partition --oneliner",
                "PartitionName=workq Default=YES MaxTime=1-00:00:00 "
                "MaxCPUsPerNode=UNLIMITED TotalCPUs=2 Nodes=chemnode2 "
                "State=UP TRES=cpu=2,mem=3000M,node=1\n",
            ),
            ("scontrol", "show", "node", "chemnode2"): _result(
                "scontrol show node chemnode2",
                node_payload,
            ),
        }
    )

    survey = run_schedule_survey(runner, Topology("A", "localhost", []))

    assert survey.scheduler == "SLURM"
    assert survey.chosen_queue == "workq"
    assert survey.evidence == {
        "scontrol show partition --oneliner": "parsed",
        "scontrol show node chemnode2": "parsed",
    }
    assert survey.queues == [
        QueueFacts(
            name="workq",
            default=True,
            max_walltime_hours=24,
            default_walltime_hours=None,
            default_mem_gb=None,
            default_cores=2,
            gpus_per_node=None,
            enabled=True,
            started=True,
            slots_total=None,
            mem_mb=3000,
        )
    ]


def test_run_schedule_survey_pbs_openpbs_live_shape_uses_node_resources():
    runner = StubRunner(
        local_results={
            ("qstat", "-Q", "-f", "-F", "json"): _result(
                "qstat -Q -f -F json",
                (FIXTURE_DIR / "pbs" / "qstat_qf_json.txt").read_text(),
            ),
            ("qstat", "-Q", "-f"): _result(
                "qstat -Q -f",
                (FIXTURE_DIR / "pbs" / "qstat_qf_text.txt").read_text(),
            ),
            ("qmgr", "-c", "list server"): _result(
                'qmgr -c "list server"',
                (FIXTURE_DIR / "pbs" / "qmgr_list_server.txt").read_text(),
            ),
            ("pbsnodes", "-av"): _result(
                "pbsnodes -av",
                (FIXTURE_DIR / "pbs" / "pbsnodes_av.txt").read_text(),
            ),
        }
    )

    survey = run_schedule_survey(runner, Topology("A", "localhost", []))

    assert survey.scheduler == "PBS"
    assert survey.chosen_queue == "workq"
    assert survey.evidence == {
        "qstat -Q -f -F json": "parsed",
        "qstat -Q -f": "parsed",
        'qmgr -c "list server"': "parsed",
        "pbsnodes -av": "parsed",
    }
    assert survey.queues == [
        QueueFacts(
            name="workq",
            default=True,
            max_walltime_hours=None,
            default_walltime_hours=None,
            default_mem_gb=3,
            default_cores=2,
            gpus_per_node=None,
            enabled=True,
            started=True,
            slots_total=None,
        )
    ]
