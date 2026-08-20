"""Canned real outputs drive the discovery parsers (the aiida method).

Each blob is a module constant with its expected counts declared beside
it; the parsers are pure functions over (returncode, stdout, stderr), so
no cluster is needed to prove scheduler behavior.
"""

from __future__ import annotations

import pytest

from chemsmart.settings.probe import ProbeUnitError
from chemsmart.settings.probe.pbs import (
    parse_pbsnodes,
    parse_qstat_queues,
)
from chemsmart.settings.probe.slurm import (
    parse_sinfo,
    parse_sinfo_version,
)

# Recorded on this host (slurm 26.05.2), plus a synthetic centre with a
# GPU partition, an UNLIMITED limit, and a heterogeneous memory marker.
_SINFO = """\
compute*|up|2-00:00:00|6|60000|1|(null)
standard|up|1-00:00:00|128|515000+|412|(null)
gpu|up|12:00:00|64|1031000|24|gpu:a100:4
debug|down|30|8|64000|2|(null)
longrun|up|UNLIMITED|128|515000|64|(null)
"""
_SINFO_QUEUES = 5
_SINFO_DEFAULT = "compute"
_SINFO_DOWN = {"debug"}


def test_sinfo_partitions_parse_to_typed_queue_facts():
    queues = parse_sinfo(0, _SINFO, "")

    assert len(queues) == _SINFO_QUEUES
    by_name = {queue.name: queue for queue in queues}
    assert [q.name for q in queues if q.is_default] == [_SINFO_DEFAULT]
    assert {q.name for q in queues if not q.available} == _SINFO_DOWN

    compute = by_name["compute"]
    assert compute.max_time_seconds == 172800
    assert compute.cores_per_node == 6
    assert compute.mem_kb_per_node == 60000 * 1024
    assert compute.node_count == 1

    standard = by_name["standard"]
    assert standard.mem_kb_per_node == 515000 * 1024  # trailing + stripped

    gpu = by_name["gpu"]
    assert gpu.gres == "gpu:a100:4"

    assert by_name["longrun"].max_time_seconds is None  # UNLIMITED
    assert by_name["debug"].max_time_seconds == 30 * 60  # bare minutes


def test_a_failing_sinfo_refuses_instead_of_inventing_queues():
    with pytest.raises(ProbeUnitError, match="sinfo exited 1"):
        parse_sinfo(1, "", "slurm_load_partitions: unable to contact")


def test_sinfo_version_line_is_taken_verbatim():
    assert parse_sinfo_version(0, "slurm 26.05.2\n", "") == "slurm 26.05.2"
    assert parse_sinfo_version(1, "", "") == ""


# Shaped after real qstat -Q -f output; stanza folding adapted from
# aiida-core (MIT).
_QSTAT_QF = """\
Queue: workq
    queue_type = Execution
    resources_max.ncpus = 64
    resources_max.mem = 250gb
    resources_max.walltime = 72:00:00
    enabled = True
    started = True

Queue: gpuq
    queue_type = Execution
    resources_max.walltime = 24:00:00
    enabled = True
    started = False
"""
_QSTAT_QUEUES = 2
_QSTAT_STOPPED = {"gpuq"}


def test_qstat_queue_stanzas_parse_with_limits():
    queues = parse_qstat_queues(0, _QSTAT_QF, "")

    assert len(queues) == _QSTAT_QUEUES
    by_name = {queue.name: queue for queue in queues}
    workq = by_name["workq"]
    assert workq.cores_per_node == 64
    assert workq.mem_kb_per_node == 250 * 1024 * 1024
    assert workq.max_time_seconds == 72 * 3600
    assert workq.available is True
    assert {q.name for q in queues if not q.available} == _QSTAT_STOPPED


_PBSNODES = """\
node001
    Mom = node001.cluster
    resources_available.ncpus = 64
    resources_available.mem = 257698037kb
    state = free

node002
    Mom = node002.cluster
    resources_available.ncpus = 64
    resources_available.mem = 257698037kb
    state = job-busy
"""


def test_pbsnodes_reports_node_capability_and_count():
    cores, mem_kb, count = parse_pbsnodes(0, _PBSNODES, "")

    assert cores == 64
    assert mem_kb == 257698037
    assert count == 2


def test_dead_pbsnodes_degrades_to_unknown_capability():
    assert parse_pbsnodes(1, "", "connection refused") == (None, None, 0)
