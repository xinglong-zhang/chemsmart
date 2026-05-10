from chemsmart.agent.wizard.parsers import (
    QueueFacts,
    parse_lsf_bqueues_l,
    parse_pbs_qstat_qf_json,
    parse_pbs_qstat_qf_text,
    parse_sge_qconf_sq,
    parse_sge_qconf_sql,
    parse_slurm_scontrol_partition_oneliner,
    parse_slurm_sinfo_json,
)


def test_parse_slurm_sinfo_json():
    assert parse_slurm_sinfo_json(
        '{"partitions": [{"name": "debug", "default": true, '
        '"max_time": "02:00:00", "state": "UP"}]}'
    ) == [
        QueueFacts(
            name="debug",
            default=True,
            max_walltime_hours=2,
            default_walltime_hours=None,
            default_mem_gb=None,
            default_cores=None,
            gpus_per_node=None,
            enabled=True,
            started=True,
        )
    ]


def test_parse_slurm_scontrol_partition_oneliner():
    assert parse_slurm_scontrol_partition_oneliner(
        "PartitionName=normal Default=YES MaxTime=2-00:00:00 "
        "DefaultTime=08:00:00 DefMemPerNode=65536 MaxCPUsPerNode=32 "
        "State=UP TRES=cpu=32,gres/gpu=0\n"
    ) == [
        QueueFacts(
            name="normal",
            default=True,
            max_walltime_hours=48,
            default_walltime_hours=8,
            default_mem_gb=64,
            default_cores=32,
            gpus_per_node=0,
            enabled=True,
            started=True,
        )
    ]


def test_parse_pbs_qstat_qf_json():
    assert parse_pbs_qstat_qf_json(
        '{"Queue": {"workq": {"default_queue": "True", '
        '"resources_default.walltime": "12:00:00", '
        '"resources_max.walltime": "48:00:00", '
        '"resources_default.mem": "64gb", '
        '"resources_default.ncpus": "16", '
        '"resources_available.ngpus": "0", '
        '"enabled": "True", "started": "True"}}}'
    ) == [
        QueueFacts(
            name="workq",
            default=True,
            max_walltime_hours=48,
            default_walltime_hours=12,
            default_mem_gb=64,
            default_cores=16,
            gpus_per_node=0,
            enabled=True,
            started=True,
        )
    ]


def test_parse_pbs_qstat_qf_text():
    assert parse_pbs_qstat_qf_text(
        "Queue workq\n"
        "    default_queue = True\n"
        "    resources_default.walltime = 12:00:00\n"
        "    resources_max.walltime = 48:00:00\n"
        "    resources_default.mem = 64gb\n"
        "    resources_default.ncpus = 16\n"
        "    enabled = True\n"
        "    started = True\n"
    ) == [
        QueueFacts(
            name="workq",
            default=True,
            max_walltime_hours=48,
            default_walltime_hours=12,
            default_mem_gb=64,
            default_cores=16,
            gpus_per_node=None,
            enabled=True,
            started=True,
        )
    ]


def test_parse_lsf_bqueues_l():
    assert parse_lsf_bqueues_l(
        "QUEUE: normal\n"
        "Open:Active\n"
        "DEFAULT_QUEUE = Y\n"
        "RUNLIMIT = 24.0 h\n"
        "MEMLIMIT = 64 G\n"
        "PROCLIMIT = 16\n"
        "GPU REQUIREMENT DETAILS: num=0\n"
    ) == [
        QueueFacts(
            name="normal",
            default=True,
            max_walltime_hours=24,
            default_walltime_hours=24,
            default_mem_gb=64,
            default_cores=16,
            gpus_per_node=0,
            enabled=True,
            started=True,
        )
    ]


def test_parse_sge_qconf_sql():
    assert parse_sge_qconf_sql("# comment\nall.q\ngpu.q\n") == [
        "all.q",
        "gpu.q",
    ]


def test_parse_sge_qconf_sq():
    assert parse_sge_qconf_sq(
        "qname all.q\n"
        "slots 8\n"
        "h_rt 86400\n"
        "h_vmem 8G\n"
        "complex_values gpu=2\n"
        "state enabled\n"
    ) == QueueFacts(
        name="all.q",
        default=False,
        max_walltime_hours=24,
        default_walltime_hours=None,
        default_mem_gb=8,
        default_cores=8,
        gpus_per_node=2,
        enabled=True,
        started=True,
    )
