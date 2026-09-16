"""What the scheduler is asked for, and who decides it.

Two objects claimed to own one question and never met. The execution
envelope carries ``cores`` and ``memory_gb``, is digest-sealed as
``resource_sha256``, and is what the human reads in the review's bounds
panel before approving. The operator's server profile carries
``NUM_CORES`` and ``MEM_GB``, and was the only thing an ``#SBATCH`` line
could see, because ``Submitter`` holds no envelope.

The evidence that this mattered: CUHK goal ``sm3-run3-cuhk`` (job
2134343, 2026-09-15) displayed and approved 4 cores / 8 GB on every node
panel and was allocated 32 cores / 160 GB, for a water molecule that ran
31 s at 424 MB peak RSS. The approved digest covered numbers that never
reached the scheduler.

The profile is a **ceiling**, not a request (owner ruling, 2026-09-16).
The envelope says what this work needs; the profile says how large the
machine will let it get, and keeps owning queue, account and the host's
extra commands.

A request above that ceiling is **refused**, not quietly shrunk (owner
ruling, 2026-09-16, revised on evidence the first ruling did not have).
How much memory a correlated calculation gets is not a resource
preference: it selects the algorithm -- integral-direct against
conventional, disk against in-core -- and sometimes decides whether the
calculation is possible at all. Shrinking it is a method decision, and
the host does not make those. A refusal naming both numbers sends the
choice back to the session, which can re-plan the method deliberately;
a clamp sends it to the engine, which discovers it as an OOM kill that
arrives wearing a scientific failure's word.

Wall time is the envelope's too, and here the profile is a cap rather
than a refusal: a shorter wall clock does not change the method, and the
host's own per-node timeout sits far below both. The asymmetry is
deliberate and is stated where it is applied.

For a **sealed** job the ceiling is the server's maximum available CPU
cores and its maximum memory less ``SEALED_MEMORY_HEADROOM_GB``, so a
clamped request can never claim a node's entire RAM. This is a persistent
default for sealed jobs until a later owner decision replaces it. It is a
bound, not a target, and it is a different axis from the concurrency
default: a cohort's footprint is this ceiling *times* the number of
calculations allowed to run at once, and the review displays both so the
product is visible before approval rather than inferred after.

The ceiling reads the server profile and nothing else. Deriving it from a
live scheduler probe would give one question two authorities, which is
the defect class this module exists to remove -- a profile that overstates
its machine is a fact about that profile, visible where it is written.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Optional

from chemsmart.agent._contracts import ContractError

#: Memory a sealed job's ceiling leaves to the node: the operating
#: system, the filesystem cache and the scheduler's own accounting. A
#: ceiling equal to a node's whole RAM is a ceiling that OOMs.
SEALED_MEMORY_HEADROOM_GB = 6.0


def _as_number(value: float) -> Any:
    """Keep an integral value integral.

    ``#SBATCH --mem=46.0G`` is not a memory specification Slurm accepts,
    and an unsealed request must render exactly the bytes the profile
    rendered before this module existed.
    """

    number = float(value)
    return int(number) if number.is_integer() else number


@dataclass(frozen=True)
class SchedulerRequestV1:
    """One resolved allocation request, with what was asked and allowed.

    ``cores``/``memory_gb``/``gpu_count`` are what the scheduler is asked
    for. The ``requested_*`` and ``ceiling_*`` fields are kept beside them
    because a number and the bound that produced it are different facts,
    and a reader of the applied value is entitled to both.
    """

    cores: int
    memory_gb: Any
    gpu_count: int
    hours: int
    requested_cores: int
    requested_memory_gb: Any
    requested_gpu_count: int
    requested_hours: int
    ceiling_cores: int
    ceiling_memory_gb: Any
    ceiling_gpu_count: int
    ceiling_hours: int
    sealed: bool
    observations: tuple[str, ...] = ()

    @property
    def memory_directive(self) -> str:
        """The memory this request asks for, as a scheduler will parse it.

        Every scheduler here takes an integer with a unit suffix, so a
        fractional gigabyte is expressed in megabytes rather than rounded:
        ``#SBATCH --mem=87.5G`` is not a specification Slurm accepts, and
        a node with 93.5 GB minus the sealed headroom is exactly how one
        arises. Flooring to 87G would quietly shrink an approved
        allocation, so nothing is rounded -- the same quantity is simply
        named in a unit that can hold it.

        Integral gigabytes keep the ``NG`` spelling they always had, so a
        profile of whole gigabytes renders byte-identically.
        """

        value = float(self.memory_gb)
        if value.is_integer():
            return f"{int(value)}G"
        megabytes = int(round(value * 1024))
        return f"{megabytes}M"

    @property
    def clamped(self) -> bool:
        """Whether any dimension was reduced to its ceiling.

        Only wall time can be: cores, memory and GPUs are refused above
        the ceiling rather than reduced.
        """

        return bool(self.observations)

    def public_record(self) -> dict:
        return {
            "applied": {
                "cores": self.cores,
                "memory_gb": self.memory_gb,
                "gpu_count": self.gpu_count,
                "hours": self.hours,
            },
            "requested": {
                "cores": self.requested_cores,
                "memory_gb": self.requested_memory_gb,
                "gpu_count": self.requested_gpu_count,
                "hours": self.requested_hours,
            },
            "ceiling": {
                "cores": self.ceiling_cores,
                "memory_gb": self.ceiling_memory_gb,
                "gpu_count": self.ceiling_gpu_count,
                "hours": self.ceiling_hours,
            },
            "sealed": self.sealed,
            "clamped": self.clamped,
            "observations": list(self.observations),
        }


def _ceiling(server: Any, *, sealed: bool) -> tuple[int, float, int]:
    cores = int(getattr(server, "num_cores", 0) or 0)
    memory_gb = float(getattr(server, "mem_gb", 0.0) or 0.0)
    gpus = int(getattr(server, "num_gpus", 0) or 0)
    if cores < 1:
        raise ContractError(
            "the server profile declares no CPU cores, so it can bound "
            "nothing; set NUM_CORES to the partition's own core count"
        )
    if sealed:
        memory_gb = memory_gb - SEALED_MEMORY_HEADROOM_GB
        if memory_gb <= 0:
            raise ContractError(
                "the server profile declares "
                f"{getattr(server, 'mem_gb', 0)} GB, which leaves nothing "
                f"above the {SEALED_MEMORY_HEADROOM_GB} GB sealed-job "
                "headroom, so no sealed request can be satisfied on it. "
                "Raise MEM_GB to describe the machine, or dispatch this "
                "goal with --unsealed, which drops the headroom and holds "
                "the request to the profile alone."
            )
    return cores, memory_gb, gpus


def _requested_hours(resources: Any, envelope: Any) -> int:
    """The wall clock this allocation needs, from the approved envelope.

    The job runs one whole cycle inside the allocation, so the bound is
    the episode window plus the postprocessing reserve the human granted,
    rounded up to the hour Slurm and PBS both want. Where no envelope
    record is available the node timeout is the only bound the host has,
    and it is used rather than falling back to the operator's default --
    which on CUHK is 24 h for a 31-second job, the same defect as cores
    and memory in the dimension a scheduler actually schedules on.
    """

    seconds = 0.0
    for field in ("episode_wall_time_seconds", "postprocess_reserve_seconds"):
        value = getattr(envelope, field, None)
        if value:
            seconds += float(value)
    if seconds <= 0:
        seconds = float(getattr(resources, "node_timeout_seconds", 0) or 0)
    if seconds <= 0:
        return 0
    return max(1, math.ceil(seconds / 3600.0))


def resolve_scheduler_request(
    *,
    resources: Optional[Any],
    server: Any,
    sealed: bool = False,
    envelope: Optional[Any] = None,
) -> SchedulerRequestV1:
    """Resolve one allocation request from an envelope and a profile.

    Args:
        resources: The approved ``ExecutionResourceSpecV1``, or ``None``
            when no envelope governs this submission (the human ``sub``
            path), in which case the profile's own numbers are the
            request and nothing is clamped.
        server: The server profile, which supplies the ceiling.
        sealed: Whether this is a sealed job, whose memory ceiling is the
            profile's maximum less ``SEALED_MEMORY_HEADROOM_GB``.

    Returns:
        SchedulerRequestV1: What to ask for, what was asked, what was
        allowed, and one observation per clamped dimension.

    Raises:
        ContractError: If the profile cannot bound anything, or if its
            memory leaves nothing above the sealed headroom.
    """

    ceiling_cores, ceiling_memory, ceiling_gpus = _ceiling(
        server, sealed=sealed
    )
    ceiling_hours = int(getattr(server, "num_hours", 0) or 0)

    if resources is None:
        requested_cores, requested_memory, requested_gpus = (
            ceiling_cores,
            ceiling_memory,
            ceiling_gpus,
        )
        requested_hours = ceiling_hours
    else:
        requested_cores = int(resources.cores)
        requested_memory = float(resources.memory_gb)
        requested_gpus = int(resources.gpu_count)
        requested_hours = _requested_hours(resources, envelope)

    # Cores, memory and GPUs are refused above the ceiling rather than
    # reduced: each of them can change what calculation is actually run,
    # and choosing that is the session's, not the host's.
    excesses = []
    if requested_cores > ceiling_cores:
        excesses.append(
            f"cores: {requested_cores} requested, ceiling {ceiling_cores}"
        )
    if requested_memory > ceiling_memory:
        excesses.append(
            f"memory: {_as_number(requested_memory)} GB requested, ceiling "
            f"{_as_number(ceiling_memory)} GB"
            + (
                f" (the profile's "
                f"{_as_number(ceiling_memory + SEALED_MEMORY_HEADROOM_GB)} GB "
                f"less the {_as_number(SEALED_MEMORY_HEADROOM_GB)} GB "
                "sealed-job headroom)"
                if sealed
                else ""
            )
        )
    if requested_gpus > ceiling_gpus:
        excesses.append(
            f"GPUs: {requested_gpus} requested, ceiling {ceiling_gpus}"
        )
    if excesses:
        raise ContractError(
            "this allocation asks for more than the server profile allows: "
            + "; ".join(excesses)
            + ". The host does not shrink it: how much memory or how many "
            "cores a calculation gets can decide which algorithm runs and "
            "whether it runs at all, so the method is re-planned by the "
            "session rather than quietly changed here. Lower the "
            "envelope's resources, or name a server profile that "
            "describes a larger machine."
        )

    # Wall time is capped rather than refused: a shorter wall clock does
    # not change the method, and the host's own per-node timeout sits far
    # below both numbers. The cap is recorded like any other observation.
    observations = []
    applied_hours = requested_hours
    if ceiling_hours and requested_hours > ceiling_hours:
        applied_hours = ceiling_hours
        observations.append(
            f"wall time capped by the server profile: requested "
            f"{requested_hours} h, ceiling {ceiling_hours} h, applied "
            f"{applied_hours} h"
        )

    return SchedulerRequestV1(
        cores=int(requested_cores),
        memory_gb=_as_number(requested_memory),
        gpu_count=int(requested_gpus),
        hours=int(applied_hours),
        requested_cores=int(requested_cores),
        requested_memory_gb=_as_number(requested_memory),
        requested_gpu_count=int(requested_gpus),
        requested_hours=int(requested_hours),
        ceiling_cores=int(ceiling_cores),
        ceiling_memory_gb=_as_number(ceiling_memory),
        ceiling_gpu_count=int(ceiling_gpus),
        ceiling_hours=int(ceiling_hours),
        sealed=bool(sealed),
        observations=tuple(observations),
    )
