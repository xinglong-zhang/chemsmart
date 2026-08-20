"""User-owned bounds for one bounded local agent execution session.

The envelope authorizes *where and how much* ChemSmart may execute.  It does
not contain a molecule, method, workflow node, project, or DAG.  Those remain
normal project-YAML and scientific-workflow decisions and are admitted only
after the live ChemSmart loader, compiler, preview, and preflight have observed
them.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import yaml

from chemsmart.agent._contracts import ContractError, require_identifier
from chemsmart.agent.execution import (
    ExecutionResourceSpecV1,
    build_execution_resource_spec,
)

_ROOT_FIELDS = frozenset(
    {
        "schema_version",
        "mode",
        "allowed_program_engines",
        "resources",
        "episode_wall_time_seconds",
        "postprocess_reserve_seconds",
        "max_engine_calls",
        "scratch_root",
    }
)
_RESOURCE_FIELDS = frozenset(
    {
        "execution_target",
        "cores",
        "memory_gb",
        "gpu_count",
        "scratch_policy",
        "node_timeout_seconds",
    }
)


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ContractError(f"{label} must be a mapping")
    if any(not isinstance(key, str) for key in value):
        raise ContractError(f"{label} keys must be strings")
    return value


def _exact_fields(
    value: Mapping[str, Any], expected: frozenset[str], label: str
) -> None:
    missing = sorted(expected.difference(value))
    unknown = sorted(set(value).difference(expected))
    if missing:
        raise ContractError(f"{label} is missing fields: {', '.join(missing)}")
    if unknown:
        raise ContractError(
            f"{label} contains unsupported fields: {', '.join(unknown)}"
        )


def _integer(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ContractError(f"{label} must be an integer")
    return value


def _number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ContractError(f"{label} must be a number")
    result = float(value)
    if not math.isfinite(result):
        raise ContractError(f"{label} must be finite")
    return result


@dataclass(frozen=True)
class BoundedExecutionEnvelopeV1:
    """Science-free authorization for a bounded local execution episode."""

    schema_version: str
    mode: str
    allowed_program_engines: tuple[tuple[str, tuple[str, ...]], ...]
    resources: ExecutionResourceSpecV1
    episode_wall_time_seconds: int
    postprocess_reserve_seconds: int
    max_engine_calls: int
    scratch_root: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.bounded-execution-envelope.v1":
            raise ContractError(
                "unsupported bounded execution envelope schema"
            )
        if self.mode == "continuous-local":
            object.__setattr__(self, "mode", "bounded-local")
        if self.mode != "bounded-local":
            raise ContractError("bounded execution mode must be bounded-local")
        if not self.allowed_program_engines:
            raise ContractError(
                "bounded execution requires an engine allowlist"
            )
        normalized = tuple(
            sorted(
                (
                    require_identifier(program, "allowed program"),
                    tuple(
                        sorted(
                            {
                                require_identifier(engine, "allowed engine")
                                for engine in engines
                            }
                        )
                    ),
                )
                for program, engines in self.allowed_program_engines
            )
        )
        if normalized != self.allowed_program_engines:
            raise ContractError(
                "allowed program engines must be sorted, unique, and normalized"
            )
        if any(not engines for _program, engines in normalized):
            raise ContractError("every allowed program requires an engine")
        if len({program for program, _engines in normalized}) != len(
            normalized
        ):
            raise ContractError("allowed programs must be unique")
        if self.resources.execution_target != "run":
            raise ContractError("bounded local execution requires target run")
        if self.resources.scratch_policy != "server":
            raise ContractError(
                "bounded execution requires the server scratch policy"
            )
        permits_gpu = any(
            engine != "cpu"
            for _program, engines in normalized
            for engine in engines
        )
        if self.resources.gpu_count and not permits_gpu:
            raise ContractError(
                "GPU resources require an explicitly allowed non-CPU engine"
            )
        if permits_gpu and not self.resources.gpu_count:
            raise ContractError(
                "a non-CPU engine requires a positive GPU allocation"
            )
        if self.episode_wall_time_seconds < 1:
            raise ContractError("episode wall time must be positive")
        if self.postprocess_reserve_seconds < 0:
            raise ContractError("postprocess reserve cannot be negative")
        if self.postprocess_reserve_seconds >= self.episode_wall_time_seconds:
            raise ContractError(
                "postprocess reserve must be smaller than episode wall time"
            )
        if self.resources.node_timeout_seconds > (
            self.episode_wall_time_seconds - self.postprocess_reserve_seconds
        ):
            raise ContractError(
                "node timeout exceeds the executable part of the episode"
            )
        if self.max_engine_calls < 1:
            raise ContractError("max_engine_calls must be positive")
        scratch = Path(self.scratch_root)
        if not scratch.is_absolute() or scratch == Path("/"):
            raise ContractError("scratch_root must be a narrow absolute path")
        if any(character.isspace() for character in self.scratch_root) or any(
            character in self.scratch_root for character in "'\"`$;&|<>"
        ):
            raise ContractError(
                "scratch_root must be safe for program environment materialization"
            )

    def allows(self, program: str, engine: str) -> bool:
        """Return whether a concrete program/engine pair was authorized."""

        allowed = dict(self.allowed_program_engines)
        return str(engine).strip().lower() in allowed.get(
            str(program).strip().lower(), ()
        )

    def public_record(self) -> dict[str, Any]:
        """Expose operating bounds without leaking an absolute host path."""

        return {
            "schema_version": "chemsmart.public-bounded-execution-envelope.v1",
            "mode": self.mode,
            "allowed_program_engines": {
                program: engines
                for program, engines in self.allowed_program_engines
            },
            "resources": {
                "execution_target": self.resources.execution_target,
                "cores": self.resources.cores,
                "memory_gb": self.resources.memory_gb,
                "gpu_count": self.resources.gpu_count,
                "scratch_policy": self.resources.scratch_policy,
                "node_timeout_seconds": self.resources.node_timeout_seconds,
            },
            "episode_wall_time_seconds": self.episode_wall_time_seconds,
            "postprocess_reserve_seconds": self.postprocess_reserve_seconds,
            "max_engine_calls": self.max_engine_calls,
        }


def load_bounded_execution_envelope(
    path: str | Path,
) -> BoundedExecutionEnvelopeV1:
    """Load an exact YAML envelope without interpreting it as shell text."""

    source = Path(path)
    if not source.is_file() or source.is_symlink():
        raise ContractError(
            "execution envelope must be a current regular file"
        )
    try:
        payload = yaml.safe_load(source.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, yaml.YAMLError) as exc:
        raise ContractError(
            "execution envelope is not valid UTF-8 YAML"
        ) from exc
    root = _mapping(payload, "execution envelope")
    _exact_fields(root, _ROOT_FIELDS, "execution envelope")
    resources = _mapping(root["resources"], "execution resources")
    _exact_fields(resources, _RESOURCE_FIELDS, "execution resources")
    raw_allowlist = _mapping(
        root["allowed_program_engines"], "allowed_program_engines"
    )
    allowlist = []
    for raw_program, raw_engines in raw_allowlist.items():
        program = require_identifier(
            str(raw_program).strip().lower(), "allowed program"
        )
        if isinstance(raw_engines, str) or not isinstance(
            raw_engines, (list, tuple)
        ):
            raise ContractError(
                f"allowed engines for {program} must be a sequence"
            )
        engines = tuple(
            sorted(
                {
                    require_identifier(
                        str(engine).strip().lower(), "allowed engine"
                    )
                    for engine in raw_engines
                }
            )
        )
        allowlist.append((program, engines))
    resource_spec = build_execution_resource_spec(
        execution_target=str(resources["execution_target"]).strip().lower(),
        cores=_integer(resources["cores"], "resources.cores"),
        memory_gb=_number(resources["memory_gb"], "resources.memory_gb"),
        gpu_count=_integer(resources["gpu_count"], "resources.gpu_count"),
        scratch_policy=str(resources["scratch_policy"]).strip().lower(),
        node_timeout_seconds=_integer(
            resources["node_timeout_seconds"],
            "resources.node_timeout_seconds",
        ),
    )
    return BoundedExecutionEnvelopeV1(
        schema_version=str(root["schema_version"]).strip(),
        mode=str(root["mode"]).strip().lower(),
        allowed_program_engines=tuple(sorted(allowlist)),
        resources=resource_spec,
        episode_wall_time_seconds=_integer(
            root["episode_wall_time_seconds"], "episode_wall_time_seconds"
        ),
        postprocess_reserve_seconds=_integer(
            root["postprocess_reserve_seconds"],
            "postprocess_reserve_seconds",
        ),
        max_engine_calls=_integer(
            root["max_engine_calls"], "max_engine_calls"
        ),
        scratch_root=str(root["scratch_root"]).strip(),
    )


__all__ = ["BoundedExecutionEnvelopeV1", "load_bounded_execution_envelope"]
