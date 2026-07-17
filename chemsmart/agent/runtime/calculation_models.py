"""Calculation lifecycle contracts and durable receipt storage."""

from __future__ import annotations

import json
import threading
from contextvars import ContextVar, Token
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any, Callable


UTC = timezone.utc


class CalculationStatus(str, Enum):
    VALIDATING = "validating"
    STARTING = "starting"
    RUNNING = "running"
    COMPLETED = "completed"
    CHEMISTRY_FAILED = "chemistry_failed"
    PROCESS_FAILED = "process_failed"
    CANCELLED = "cancelled"
    TIMEOUT = "timeout"

    @property
    def terminal(self) -> bool:
        return self in {
            self.COMPLETED,
            self.CHEMISTRY_FAILED,
            self.PROCESS_FAILED,
            self.CANCELLED,
            self.TIMEOUT,
        }


@dataclass(slots=True)
class CalculationRun:
    run_id: str
    command: str
    cwd: str
    session_id: str = ""
    turn_id: str = ""
    program: str = ""
    kind: str = ""
    label: str = ""
    project: str = ""
    project_path: str = ""
    project_sha256: str = ""
    input_path: str = ""
    output_path: str = ""
    reused_output: bool = False
    method: str = ""
    basis: str = ""
    execution_mode: str = "local"
    status: str = CalculationStatus.VALIDATING.value
    stage: str = "Validating command"
    pid: int | None = None
    started_at: str = ""
    updated_at: str = ""
    ended_at: str = ""
    elapsed_s: float = 0.0
    chemistry_elapsed_s: float | None = None
    returncode: int | None = None
    stdout_path: str = ""
    stderr_path: str = ""
    normal_termination: bool | None = None
    energy: float | None = None
    scf_cycles: int | None = None
    optimization_cycles: int | None = None
    optimization_converged: bool | None = None
    frequencies: list[float] = field(default_factory=list)
    imag_freqs: list[float] = field(default_factory=list)
    frequency_count: int | None = None
    imaginary_frequency_count: int | None = None
    scan_points_completed: int | None = None
    scan_points_total: int | None = None
    energy_min: float | None = None
    energy_max: float | None = None
    neb_images: int | None = None
    qmmm_qm_atoms: int | None = None
    qmmm_mm_atoms: int | None = None
    qmmm_total_atoms: int | None = None
    path_direction: str = ""
    parser_backend: str = "native"
    parser_warnings: list[str] = field(default_factory=list)
    atom_count: int | None = None
    parsed_charge: int | None = None
    parsed_multiplicity: int | None = None
    semantic_verdict: str = ""
    intent_verdict: str = ""
    error: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "CalculationRun":
        known = cls.__dataclass_fields__
        return cls(
            **{key: value for key, value in payload.items() if key in known}
        )


@dataclass(slots=True, frozen=True)
class CalculationEvent:
    kind: str
    run: CalculationRun
    timestamp: str = field(
        default_factory=lambda: datetime.now(UTC).isoformat()
    )

    def to_dict(self) -> dict[str, Any]:
        return {
            "kind": self.kind,
            "timestamp": self.timestamp,
            "run": self.run.to_dict(),
        }


EventSink = Callable[[CalculationEvent], None]


@dataclass(slots=True, frozen=True)
class CalculationContext:
    session_dir: Path | None = None
    session_id: str = ""
    turn_id: str = ""
    semantic_verdict: str = ""
    intent_verdict: str = ""


_CURRENT_CONTEXT: ContextVar[CalculationContext | None] = ContextVar(
    "chemsmart_calculation_context",
    default=None,
)


def current_calculation_context() -> CalculationContext | None:
    return _CURRENT_CONTEXT.get()


def set_calculation_context(context: CalculationContext) -> Token:
    return _CURRENT_CONTEXT.set(context)


def reset_calculation_context(token: Token) -> None:
    _CURRENT_CONTEXT.reset(token)


class CalculationStore:
    """Durable append-only event store for calculations in one session."""

    def __init__(self, session_dir: str | Path) -> None:
        self.session_dir = Path(session_dir)
        self.root = self.session_dir / "calculations"
        self.root.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()

    def run_dir(self, run_id: str) -> Path:
        path = self.root / run_id
        path.mkdir(parents=True, exist_ok=True)
        return path

    def write_run(self, run: CalculationRun) -> None:
        destination = self.run_dir(run.run_id) / "receipt.json"
        temporary = destination.with_suffix(".json.tmp")
        payload = json.dumps(run.to_dict(), indent=2, sort_keys=True)
        with self._lock:
            temporary.write_text(payload, encoding="utf-8")
            temporary.replace(destination)

    def append_event(self, event: CalculationEvent) -> None:
        path = self.run_dir(event.run.run_id) / "events.jsonl"
        line = json.dumps(event.to_dict(), sort_keys=True)
        with self._lock, path.open("a", encoding="utf-8") as handle:
            handle.write(line + "\n")

    def record(self, event: CalculationEvent) -> None:
        self.write_run(event.run)
        self.append_event(event)

    def load(self, run_id: str) -> CalculationRun | None:
        path = self.root / run_id / "receipt.json"
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return None
        return CalculationRun.from_dict(payload)

    def list_runs(self) -> list[CalculationRun]:
        runs: list[CalculationRun] = []
        for path in self.root.glob("*/receipt.json"):
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            runs.append(CalculationRun.from_dict(payload))
        return sorted(runs, key=lambda run: run.updated_at or run.started_at)


def apply_chemistry_summary(
    run: CalculationRun,
    summary: dict[str, Any],
) -> None:
    if not summary:
        return
    for field_name in (
        "normal_termination",
        "energy",
        "scf_cycles",
        "optimization_cycles",
        "optimization_converged",
        "frequency_count",
        "imaginary_frequency_count",
        "scan_points_completed",
        "scan_points_total",
        "energy_min",
        "energy_max",
        "neb_images",
        "qmmm_qm_atoms",
        "qmmm_mm_atoms",
        "qmmm_total_atoms",
        "chemistry_elapsed_s",
        "atom_count",
        "parsed_charge",
        "parsed_multiplicity",
    ):
        setattr(run, field_name, summary.get(field_name))
    run.frequencies = list(summary.get("frequencies") or [])
    run.imag_freqs = list(summary.get("imag_freqs") or [])
    run.path_direction = str(summary.get("path_direction") or "")
    run.parser_backend = str(summary.get("parser_backend") or "native")
    run.parser_warnings = list(summary.get("parser_warnings") or [])


__all__ = [
    "CalculationContext",
    "CalculationEvent",
    "CalculationRun",
    "CalculationStatus",
    "CalculationStore",
    "EventSink",
    "UTC",
    "apply_chemistry_summary",
    "current_calculation_context",
    "reset_calculation_context",
    "set_calculation_context",
]
