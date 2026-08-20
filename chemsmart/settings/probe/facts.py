"""Typed facts a host or scheduler probe reports; nothing else crosses."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class QueueFactsV1:
    """One queue/partition as the scheduler advertises it."""

    name: str
    is_default: bool = False
    available: bool = True
    max_time_seconds: int | None = None
    cores_per_node: int | None = None
    mem_kb_per_node: int | None = None
    node_count: int | None = None
    gres: str = ""


@dataclass(frozen=True)
class SchedulerFactsV1:
    """What was actually observed about the batch system."""

    scheduler: str  # SLURM | PBS | SGE | LSF | HTCondor
    version: str = ""
    submit_path: str = ""
    queues: tuple[QueueFactsV1, ...] = ()
    evidence: tuple[str, ...] = ()


@dataclass(frozen=True)
class HostFactsV1:
    """The login host itself, for the local fallback and the report."""

    hostname: str = ""
    os_name: str = ""
    os_release: str = ""
    machine: str = ""
    cpu_count: int | None = None
    mem_kb: int | None = None
    scratch_candidates: tuple[str, ...] = ()
    has_module_command: bool = False


@dataclass(frozen=True)
class DetectionV1:
    """The detection verdict plus the evidence it rests on."""

    scheduler: str | None
    evidence: tuple[str, ...] = field(default_factory=tuple)


__all__ = [
    "DetectionV1",
    "HostFactsV1",
    "QueueFactsV1",
    "SchedulerFactsV1",
]
