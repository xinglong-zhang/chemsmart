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
machine will let it get, and keeps owning queue, account, wall time and
the host's extra commands. A clamp is a displayed observation naming both
numbers, never a silent edit.

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
    requested_cores: int
    requested_memory_gb: Any
    requested_gpu_count: int
    ceiling_cores: int
    ceiling_memory_gb: Any
    ceiling_gpu_count: int
    sealed: bool
    observations: tuple[str, ...] = ()

    @property
    def clamped(self) -> bool:
        """Whether any dimension was reduced to its ceiling."""

        return bool(self.observations)

    def public_record(self) -> dict:
        return {
            "applied": {
                "cores": self.cores,
                "memory_gb": self.memory_gb,
                "gpu_count": self.gpu_count,
            },
            "requested": {
                "cores": self.requested_cores,
                "memory_gb": self.requested_memory_gb,
                "gpu_count": self.requested_gpu_count,
            },
            "ceiling": {
                "cores": self.ceiling_cores,
                "memory_gb": self.ceiling_memory_gb,
                "gpu_count": self.ceiling_gpu_count,
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
                "headroom; raise MEM_GB or run this goal unsealed rather "
                "than requesting a negative allocation"
            )
    return cores, memory_gb, gpus


def resolve_scheduler_request(
    *,
    resources: Optional[Any],
    server: Any,
    sealed: bool = False,
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

    if resources is None:
        applied_cores, applied_memory, applied_gpus = (
            ceiling_cores,
            ceiling_memory,
            ceiling_gpus,
        )
        requested_cores, requested_memory, requested_gpus = (
            ceiling_cores,
            ceiling_memory,
            ceiling_gpus,
        )
    else:
        requested_cores = int(resources.cores)
        requested_memory = float(resources.memory_gb)
        requested_gpus = int(resources.gpu_count)
        applied_cores = min(requested_cores, ceiling_cores)
        applied_memory = min(requested_memory, ceiling_memory)
        applied_gpus = min(requested_gpus, ceiling_gpus)

    observations = []
    if applied_cores != requested_cores:
        observations.append(
            f"cores clamped to the server ceiling: requested "
            f"{requested_cores}, ceiling {ceiling_cores}, "
            f"applied {applied_cores}"
        )
    if applied_memory != requested_memory:
        observations.append(
            f"memory clamped to the server ceiling: requested "
            f"{_as_number(requested_memory)} GB, ceiling "
            f"{_as_number(ceiling_memory)} GB, applied "
            f"{_as_number(applied_memory)} GB"
            + (
                f" (the profile's {_as_number(ceiling_memory + SEALED_MEMORY_HEADROOM_GB)} GB "
                f"less the {_as_number(SEALED_MEMORY_HEADROOM_GB)} GB "
                "sealed-job headroom)"
                if sealed
                else ""
            )
        )
    if applied_gpus != requested_gpus:
        observations.append(
            f"GPUs clamped to the server ceiling: requested "
            f"{requested_gpus}, ceiling {ceiling_gpus}, "
            f"applied {applied_gpus}"
        )

    return SchedulerRequestV1(
        cores=int(applied_cores),
        memory_gb=_as_number(applied_memory),
        gpu_count=int(applied_gpus),
        requested_cores=int(requested_cores),
        requested_memory_gb=_as_number(requested_memory),
        requested_gpu_count=int(requested_gpus),
        ceiling_cores=int(ceiling_cores),
        ceiling_memory_gb=_as_number(ceiling_memory),
        ceiling_gpu_count=int(ceiling_gpus),
        sealed=bool(sealed),
        observations=tuple(observations),
    )
