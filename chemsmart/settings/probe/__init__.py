"""Read-only host and scheduler probes behind ``chemsmart wizard``.

Every probe is a pure pair -- a command builder and a
``parse(returncode, stdout, stderr)`` function returning typed facts -- so
each scheduler behavior is testable from canned real outputs without a
cluster (the aiida-core testing method).
"""

from .detect import detect_scheduler
from .facts import (
    DetectionV1,
    HostFactsV1,
    QueueFactsV1,
    SchedulerFactsV1,
)
from .localhost import gather_host_facts
from .units import ProbeUnitError, parse_memory_kb, parse_walltime_seconds

__all__ = [
    "DetectionV1",
    "HostFactsV1",
    "ProbeUnitError",
    "QueueFactsV1",
    "SchedulerFactsV1",
    "detect_scheduler",
    "gather_host_facts",
    "parse_memory_kb",
    "parse_walltime_seconds",
]
