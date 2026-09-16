"""What the scheduler is asked for: the operator's profile, resolved.

The server profile is the authority for how much machine a job gets, and
the agent does not argue with it. ``NUM_CORES``, ``MEM_GB``, ``NUM_GPUS``
and ``NUM_HOURS`` in ``server.yaml`` / ``SLURM.yaml`` are what reaches the
``#SBATCH`` line and what the engine is told it may use, and an operator
who wants a different allocation edits that file rather than arguing with
the agent about it (owner ruling, 2026-09-16).

Two earlier designs are recorded here because the repository record named
them owner rulings and they were mine to propose:

- The envelope's ``cores``/``memory_gb`` were made the scheduler request
  with the profile as a ceiling. That put a scientific object in charge of
  an operator's resource policy.
- A request above that ceiling was then **refused**. A refusal after a
  goal is under way is an agent-tool failure: it teaches a session to
  build workarounds and carry control logic for something the host should
  simply have decided. There is no refusal here now, and no clamp: there
  is one number, from one file, used by the scheduler and the engine
  alike.

What the sealed default *is*: the server's maximum available CPU cores,
and its maximum memory less ``SEALED_MEMORY_HEADROOM_GB``, so a job never
claims a node's entire RAM and leaves the operating system, the
filesystem cache and the scheduler's own accounting somewhere to live.
``--unsealed`` drops the headroom and takes the profile's memory whole.

The envelope's own resources are still recorded beside the applied
numbers, because a reader of a receipt is entitled to see what the
approved episode asked for and what the machine gave it. They are an
observation and never an action.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

#: Memory a sealed job leaves to the node: the operating system, the
#: filesystem cache and the scheduler's own accounting. A request equal to
#: a node's whole RAM is a request that OOMs.
SEALED_MEMORY_HEADROOM_GB = 6.0


def _as_number(value: float) -> Any:
    """Keep an integral value integral.

    A profile of whole gigabytes must render the bytes it always
    rendered; ``46.0`` is not what ``46`` spelled.
    """

    number = float(value)
    return int(number) if number.is_integer() else number


@dataclass(frozen=True)
class SchedulerRequestV1:
    """One resolved allocation: the operator's numbers, and the episode's.

    ``cores``/``memory_gb``/``gpu_count``/``hours`` are what the scheduler
    is asked for and what the engine is told it has -- one set of numbers
    for both, so an allocation and the program running inside it can never
    disagree. ``envelope_*`` is what the approved episode asked for, kept
    beside them as an observation.
    """

    cores: int
    memory_gb: Any
    gpu_count: int
    hours: int
    sealed: bool
    memory_headroom_gb: Any = 0
    envelope_cores: Optional[int] = None
    envelope_memory_gb: Any = None
    observations: tuple[str, ...] = ()

    @property
    def memory_directive(self) -> str:
        """The memory this request asks for, as a scheduler will parse it.

        Every scheduler here takes an integer with a unit suffix, so a
        fractional gigabyte is named in megabytes rather than rounded:
        ``--mem=87.5G`` is not a specification Slurm accepts, and a node
        of 93.5 GB less the sealed headroom is exactly how one arises.
        Rounding would quietly change the allocation, so the same quantity
        is simply named in a unit that can hold it. Integral gigabytes
        keep the ``NG`` spelling they always had.
        """

        value = float(self.memory_gb)
        if value.is_integer():
            return f"{int(value)}G"
        return f"{int(round(value * 1024))}M"

    def public_record(self) -> dict:
        return {
            "applied": {
                "cores": self.cores,
                "memory_gb": self.memory_gb,
                "gpu_count": self.gpu_count,
                "hours": self.hours,
            },
            "source": "server_profile",
            "sealed": self.sealed,
            "memory_headroom_gb": self.memory_headroom_gb,
            "envelope": {
                "cores": self.envelope_cores,
                "memory_gb": self.envelope_memory_gb,
            },
            "observations": list(self.observations),
        }


def resolve_scheduler_request(
    *,
    resources: Optional[Any] = None,
    server: Any,
    sealed: bool = False,
    envelope: Optional[Any] = None,
) -> SchedulerRequestV1:
    """Resolve the allocation this submission asks the scheduler for.

    The numbers are the server profile's. Nothing is refused and nothing
    is clamped: an operator who wants a different allocation edits the
    profile, which is the file that owns the question.

    Args:
        resources: The approved ``ExecutionResourceSpecV1``, recorded as
            an observation beside the applied numbers. It does not change
            them.
        server: The server profile, which supplies every applied number.
        sealed: Whether the sealed-job memory headroom applies.
        envelope: The approved execution envelope, accepted so callers do
            not have to know whether it is consulted. It is not.

    Returns:
        SchedulerRequestV1: What to ask the scheduler for, and what the
        approved episode asked for beside it.
    """

    cores = int(getattr(server, "num_cores", 0) or 0)
    memory_gb = float(getattr(server, "mem_gb", 0.0) or 0.0)
    gpus = int(getattr(server, "num_gpus", 0) or 0)
    hours = int(getattr(server, "num_hours", 0) or 0)

    headroom = 0.0
    if sealed and memory_gb > SEALED_MEMORY_HEADROOM_GB:
        headroom = SEALED_MEMORY_HEADROOM_GB
        memory_gb -= SEALED_MEMORY_HEADROOM_GB

    observations: list[str] = []
    envelope_cores = None
    envelope_memory = None
    if resources is not None:
        envelope_cores = int(resources.cores)
        envelope_memory = _as_number(float(resources.memory_gb))
        if envelope_cores != cores or float(envelope_memory) != memory_gb:
            # An observation, never an action: the operator's profile
            # decides, and a reader of the receipt sees both numbers
            # rather than inferring one from the other.
            observations.append(
                "the approved episode asked for "
                f"{envelope_cores} cores / {envelope_memory} GB; the "
                f"server profile allocates {cores} cores / "
                f"{_as_number(memory_gb)} GB"
                + (
                    f" (its {_as_number(memory_gb + headroom)} GB less the "
                    f"{_as_number(headroom)} GB sealed-job headroom)"
                    if headroom
                    else ""
                )
            )

    return SchedulerRequestV1(
        cores=cores,
        memory_gb=_as_number(memory_gb),
        gpu_count=gpus,
        hours=hours,
        sealed=bool(sealed),
        memory_headroom_gb=_as_number(headroom),
        envelope_cores=envelope_cores,
        envelope_memory_gb=envelope_memory,
        observations=tuple(observations),
    )
