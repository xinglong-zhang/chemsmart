"""Observable lifecycle and chemistry summaries for CLI calculations."""

from __future__ import annotations

import json
import mmap
import os
import re
import signal
import subprocess
import threading
import time
import uuid
from contextvars import ContextVar, Token
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Iterable

import yaml

from chemsmart.agent.model_command_parser import parse_model_command
from chemsmart.agent.runtime.cclib_parser import parse_cclib_output

UTC = timezone.utc
EventSink = Callable[["CalculationEvent"], None]
_CURRENT_CONTEXT: ContextVar["CalculationContext | None"] = ContextVar(
    "chemsmart_calculation_context",
    default=None,
)


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
        return cls(**{key: value for key, value in payload.items() if key in known})


@dataclass(slots=True, frozen=True)
class CalculationEvent:
    kind: str
    run: CalculationRun
    timestamp: str = field(default_factory=lambda: datetime.now(UTC).isoformat())

    def to_dict(self) -> dict[str, Any]:
        return {
            "kind": self.kind,
            "timestamp": self.timestamp,
            "run": self.run.to_dict(),
        }


@dataclass(slots=True, frozen=True)
class CalculationContext:
    session_dir: Path | None = None
    session_id: str = ""
    turn_id: str = ""
    semantic_verdict: str = ""
    intent_verdict: str = ""


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


_ACTIVE_PROCESSES: dict[str, subprocess.Popen[str]] = {}
_ACTIVE_LOCK = threading.Lock()
_CANCELLED_RUNS: set[str] = set()


def execute_observed_process(
    argv: list[str],
    *,
    command: str,
    timeout_s: int,
    context: CalculationContext | None = None,
    event_sink: EventSink | None = None,
    cwd: str | Path | None = None,
) -> dict[str, Any]:
    """Run argv while emitting durable, chemistry-aware lifecycle events."""

    workdir = Path(cwd or Path.cwd()).resolve()
    context = context or _CURRENT_CONTEXT.get() or CalculationContext()
    parsed = parse_model_command(command, cwd=str(workdir))
    run_id = f"calc-{uuid.uuid4().hex[:10]}"
    store = CalculationStore(context.session_dir) if context.session_dir else None
    run_dir = (
        store.run_dir(run_id)
        if store is not None
        else Path.home() / ".chemsmart" / "agent" / "calculations" / run_id
    )
    run_dir.mkdir(parents=True, exist_ok=True)
    project_path, project_sha256, method, basis = _project_metadata(
        workdir,
        parsed.program or "",
        parsed.project or "",
    )
    run = CalculationRun(
        run_id=run_id,
        session_id=context.session_id,
        turn_id=context.turn_id,
        command=command,
        cwd=str(workdir),
        program=parsed.program or "",
        kind=parsed.job or "",
        label=parsed.label or _label_from_input(parsed.filename),
        project=parsed.project or "",
        project_path=project_path,
        project_sha256=project_sha256,
        input_path=_absolute_input_path(workdir, parsed.filename),
        method=parsed.functional or method,
        basis=parsed.basis or basis,
        execution_mode=(
            "test_fake"
            if "--fake" in argv or "--test" in argv
            else "hpc"
            if parsed.action == "sub"
            else "local"
        ),
        semantic_verdict=context.semantic_verdict,
        intent_verdict=context.intent_verdict,
    )
    stdout_path = run_dir / "stdout.log"
    stderr_path = run_dir / "stderr.log"
    run.stdout_path = str(stdout_path)
    run.stderr_path = str(stderr_path)
    outputs_before = _output_fingerprints(workdir, run.input_path)
    started_monotonic = time.monotonic()

    _emit(run, "validating", store, event_sink)
    run.status = CalculationStatus.STARTING.value
    run.stage = "Starting chemsmart command"
    run.started_at = datetime.now(UTC).isoformat()
    run.updated_at = run.started_at
    _emit(run, "starting", store, event_sink)

    timed_out = False
    with stdout_path.open("w", encoding="utf-8") as stdout_handle, stderr_path.open(
        "w", encoding="utf-8"
    ) as stderr_handle:
        process = subprocess.Popen(
            argv,
            cwd=str(workdir),
            text=True,
            stdout=stdout_handle,
            stderr=stderr_handle,
            start_new_session=True,
        )
        run.pid = process.pid
        run.status = CalculationStatus.RUNNING.value
        run.stage = "Process running"
        _register_process(run.run_id, process)
        _emit(run, "started", store, event_sink)
        last_progress = ""
        last_emit = 0.0
        try:
            while process.poll() is None:
                elapsed = time.monotonic() - started_monotonic
                run.elapsed_s = elapsed
                if elapsed > timeout_s:
                    timed_out = True
                    _terminate_process(process)
                    break
                output_path = _discover_output(
                    workdir, run.input_path, outputs_before
                )
                if output_path is not None:
                    run.output_path = str(output_path)
                progress = _progress_summary(
                    Path(run.output_path) if run.output_path else stdout_path,
                    run.program,
                    run.kind,
                )
                now = time.monotonic()
                if progress and progress != last_progress:
                    run.stage = progress
                    last_progress = progress
                    _emit(run, "progress", store, event_sink)
                    last_emit = now
                elif now - last_emit >= 2.0:
                    _emit(run, "heartbeat", store, event_sink)
                    last_emit = now
                time.sleep(0.2)
            process.wait()
        finally:
            _unregister_process(run.run_id)

    run.returncode = process.returncode
    run.elapsed_s = time.monotonic() - started_monotonic
    output_path = _discover_output(workdir, run.input_path, outputs_before)
    if output_path is not None:
        run.output_path = str(output_path)
    elif process.returncode == 0:
        reported_output = _reported_output_path(
            stdout_path,
            stderr_path,
            cwd=workdir,
            program=run.program,
        )
        if reported_output is not None:
            run.output_path = str(reported_output)
            run.reused_output = str(reported_output) in outputs_before
    chemistry = inspect_output(
        run.output_path,
        program=run.program,
        kind=run.kind,
    )
    _apply_chemistry_summary(run, chemistry)
    cancelled = _consume_cancelled(run.run_id)
    if timed_out:
        run.status = CalculationStatus.TIMEOUT.value
        run.stage = f"Timed out after {timeout_s}s"
        run.error = run.stage
    elif cancelled:
        run.status = CalculationStatus.CANCELLED.value
        run.stage = "Cancelled by user"
    elif process.returncode != 0:
        run.status = CalculationStatus.PROCESS_FAILED.value
        run.stage = f"Process exited with code {process.returncode}"
        run.error = _last_error(stderr_path, stdout_path)
    elif parsed.action == "sub":
        run.status = CalculationStatus.COMPLETED.value
        run.stage = "Submission command completed"
    elif run.execution_mode == "test_fake":
        run.status = CalculationStatus.COMPLETED.value
        run.stage = "Safe test command completed"
    elif not run.output_path:
        run.status = CalculationStatus.CHEMISTRY_FAILED.value
        run.stage = "No Gaussian/ORCA output file was produced"
        run.error = _last_error(stderr_path, stdout_path)
    elif run.output_path and run.normal_termination is not True:
        run.status = CalculationStatus.CHEMISTRY_FAILED.value
        run.stage = f"{run.program.upper() or 'Chemistry'} did not terminate normally"
        run.error = _last_error(Path(run.output_path), stderr_path)
    else:
        run.status = CalculationStatus.COMPLETED.value
        run.stage = (
            "Existing completed calculation reused"
            if run.reused_output
            else "Calculation completed"
        )
    run.ended_at = datetime.now(UTC).isoformat()
    run.updated_at = run.ended_at
    _emit(run, "completed" if run.status == "completed" else "failed", store, event_sink)
    return {
        "calculation": run.to_dict(),
        "stdout_tail": _tail_file(stdout_path),
        "stderr_tail": _tail_file(stderr_path),
    }


def cancel_calculation(run_id: str, pid: int | None = None) -> bool:
    """Cancel an active local calculation without touching scheduler jobs."""

    with _ACTIVE_LOCK:
        process = _ACTIVE_PROCESSES.get(run_id)
        _CANCELLED_RUNS.add(run_id)
    if process is not None and process.poll() is None:
        _terminate_process(process)
        return True
    if pid:
        try:
            os.killpg(pid, signal.SIGTERM)
        except (ProcessLookupError, PermissionError):
            return False
        return True
    return False


def inspect_calculation(
    run_id: str = "latest",
    output_path: str | None = None,
    session_root: str | None = None,
) -> dict[str, Any]:
    """Inspect a recent local calculation using deterministic output parsers."""

    run = _find_run(run_id, session_root=session_root)
    selected_output = output_path or (run.output_path if run else "")
    if not selected_output and run is None:
        selected_output = _latest_output_in_workspace(Path.cwd())
    if not selected_output:
        return {
            "ok": False,
            "status": "not_found",
            "error": "No calculation receipt or Gaussian/ORCA output was found in this workspace.",
        }
    program = run.program if run else ""
    kind = run.kind if run else ""
    summary = inspect_output(selected_output, program=program, kind=kind)
    if run is not None:
        payload = run.to_dict()
        payload.update({key: value for key, value in summary.items() if value is not None})
    else:
        payload = {
            "run_id": "workspace-output",
            "output_path": str(Path(selected_output).resolve()),
            "program": program,
            "kind": kind,
            **summary,
        }
    return {
        "ok": True,
        "status": payload.get("status") or "parsed",
        "calculation": payload,
    }


def inspect_output(
    output_path: str | Path,
    *,
    program: str = "",
    kind: str = "",
) -> dict[str, Any]:
    path = Path(output_path) if output_path else Path()
    if not output_path or not path.is_file():
        return {}
    text = _read_output_window(path)
    detected = program.lower() or _detect_program(text, path)
    result: dict[str, Any] = {
        "output_path": str(path.resolve()),
        "program": detected,
        "kind": kind,
    }
    if detected == "orca":
        energy_matches = re.findall(
            r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text
        )
        scf_matches = re.findall(
            r"SCF CONVERGED AFTER\s+(\d+)\s+CYCLES", text, re.IGNORECASE
        )
        optimization_cycles = re.findall(
            r"GEOMETRY OPTIMIZATION CYCLE\s+(\d+)", text
        )
        scan_points = re.findall(
            r"(?:RELAXED SURFACE SCAN (?:STEP|POINT)|SCAN POINT)\s+(\d+)",
            text,
            re.IGNORECASE,
        )
        scan_totals = re.findall(
            r"SCAN POINT\s+\d+\s+(?:OF|/)\s+(\d+)",
            text,
            re.IGNORECASE,
        )
        neb_images = re.findall(
            r"Number of images \(incl\. end points\)\s+\.{2,}\s+(\d+)",
            text,
        )
        qmmm_total = re.findall(
            r"Size of QMMM System\s+\.{2,}\s+(\d+)", text
        )
        qmmm_mm = re.findall(
            r"Point charges in QM calc\. from MM atoms\.{2,}\s+(\d+)",
            text,
        )
        result.update(
            energy=float(energy_matches[-1]) if energy_matches else None,
            scf_cycles=int(scf_matches[-1]) if scf_matches else None,
            optimization_cycles=(
                int(optimization_cycles[-1])
                if optimization_cycles
                else None
            ),
            optimization_converged=(
                "THE OPTIMIZATION HAS CONVERGED" in text
                if optimization_cycles
                else None
            ),
            scan_points_completed=(
                max(map(int, scan_points)) if scan_points else None
            ),
            scan_points_total=(
                int(scan_totals[-1]) if scan_totals else None
            ),
            neb_images=int(neb_images[-1]) if neb_images else None,
            qmmm_total_atoms=(
                int(qmmm_total[-1]) if qmmm_total else None
            ),
            qmmm_mm_atoms=int(qmmm_mm[-1]) if qmmm_mm else None,
            normal_termination="ORCA TERMINATED NORMALLY" in text,
        )
        if result.get("qmmm_total_atoms") is not None and result.get(
            "qmmm_mm_atoms"
        ) is not None:
            result["qmmm_qm_atoms"] = (
                int(result["qmmm_total_atoms"])
                - int(result["qmmm_mm_atoms"])
            )
        frequency_text = _last_section(text, "VIBRATIONAL FREQUENCIES")
        orca_frequencies = [
            float(value)
            for value in re.findall(
                r"^\s*\d+:\s+(-?\d+(?:\.\d+)?)\s+cm\*\*-1",
                frequency_text,
                re.MULTILINE,
            )
            if abs(float(value)) > 1e-6
        ]
        _add_frequency_summary(result, orca_frequencies)
    else:
        energy_matches = re.findall(
            r"SCF Done:\s+E\([^)]*\)\s*=\s*(-?\d+\.\d+)", text
        )
        scf_matches = re.findall(
            r"SCF Done:.*?after\s+(\d+)\s+cycles", text
        )
        cycle_matches = re.findall(r"Step number\s+(\d+)", text)
        if not cycle_matches:
            cycle_matches = re.findall(r"Cycle\s+(\d+)", text)
        scan_points = re.findall(
            r"Scan point\s+(\d+)(?:\s+out of\s+(\d+))?",
            text,
            re.IGNORECASE,
        )
        result.update(
            energy=float(energy_matches[-1]) if energy_matches else None,
            scf_cycles=int(scf_matches[-1]) if scf_matches else None,
            optimization_cycles=(
                int(cycle_matches[-1]) if cycle_matches else None
            ),
            optimization_converged=(
                "Stationary point found" in text
                if cycle_matches
                else None
            ),
            scan_points_completed=(
                max(int(point) for point, _ in scan_points)
                if scan_points
                else None
            ),
            scan_points_total=(
                max(int(total) for _, total in scan_points if total)
                if any(total for _, total in scan_points)
                else None
            ),
            normal_termination="Normal termination of Gaussian" in text,
        )
        frequency_text = _last_section(text, "Harmonic frequencies")
        frequencies = [
            float(value)
            for line in re.findall(
                r"Frequencies\s+--\s+([^\n]+)", frequency_text
            )
            for value in re.findall(r"-?\d+(?:\.\d+)?", line)
        ]
        _add_frequency_summary(result, frequencies)
    if energy_matches:
        energies = [float(value) for value in energy_matches]
        result["energy_min"] = min(energies)
        result["energy_max"] = max(energies)
    result["path_direction"] = _path_direction(text)
    result["chemistry_elapsed_s"] = _program_runtime_seconds(text, detected)
    _merge_cclib_summary(result, path)
    return result


def load_calculation_runs(session_root: str | Path) -> list[CalculationRun]:
    """Load persisted calculations across sessions without requiring state.json."""

    root = Path(session_root)
    runs: list[CalculationRun] = []
    if not root.exists():
        return runs
    for receipt in root.glob("**/calculations/*/receipt.json"):
        try:
            payload = json.loads(receipt.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        run = CalculationRun.from_dict(payload)
        if _refresh_recovered_run(run):
            temporary = receipt.with_suffix(".json.tmp")
            try:
                temporary.write_text(
                    json.dumps(run.to_dict(), indent=2, sort_keys=True),
                    encoding="utf-8",
                )
                temporary.replace(receipt)
            except OSError:
                temporary.unlink(missing_ok=True)
        runs.append(run)
    return sorted(runs, key=lambda item: item.updated_at or item.started_at)


def _emit(
    run: CalculationRun,
    kind: str,
    store: CalculationStore | None,
    sink: EventSink | None,
) -> None:
    run.updated_at = datetime.now(UTC).isoformat()
    event = CalculationEvent(kind=kind, run=run)
    if store is not None:
        store.record(event)
    if sink is not None:
        sink(event)


def _register_process(run_id: str, process: subprocess.Popen[str]) -> None:
    with _ACTIVE_LOCK:
        _ACTIVE_PROCESSES[run_id] = process


def _unregister_process(run_id: str) -> None:
    with _ACTIVE_LOCK:
        _ACTIVE_PROCESSES.pop(run_id, None)


def _consume_cancelled(run_id: str) -> bool:
    with _ACTIVE_LOCK:
        if run_id not in _CANCELLED_RUNS:
            return False
        _CANCELLED_RUNS.remove(run_id)
        return True


def _terminate_process(process: subprocess.Popen[str]) -> None:
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
        process.wait(timeout=3)
    except (ProcessLookupError, subprocess.TimeoutExpired):
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass


def _output_fingerprints(cwd: Path, input_path: str) -> dict[str, tuple[int, int]]:
    return {
        str(path.resolve()): (path.stat().st_mtime_ns, path.stat().st_size)
        for path in _output_candidates(cwd, input_path)
        if path.is_file()
    }


def _output_candidates(cwd: Path, input_path: str) -> Iterable[Path]:
    roots = {cwd}
    if input_path:
        roots.add(Path(input_path).resolve().parent)
    seen: set[Path] = set()
    for root in roots:
        if not root.is_dir():
            continue
        for suffix in ("*.out", "*.log"):
            for path in root.glob(suffix):
                resolved = path.resolve()
                if resolved not in seen:
                    seen.add(resolved)
                    yield resolved


def _discover_output(
    cwd: Path,
    input_path: str,
    before: dict[str, tuple[int, int]],
) -> Path | None:
    changed: list[Path] = []
    for path in _output_candidates(cwd, input_path):
        stat = path.stat()
        if before.get(str(path)) != (stat.st_mtime_ns, stat.st_size):
            changed.append(path)
    return max(changed, key=lambda path: path.stat().st_mtime_ns) if changed else None


def _reported_output_path(
    stdout_path: Path,
    stderr_path: Path,
    *,
    cwd: Path,
    program: str,
) -> Path | None:
    text = _read_tail(stdout_path, max_chars=120_000) + "\n" + _read_tail(
        stderr_path, max_chars=40_000
    )
    folders = re.findall(r"\bfolder=([^,>\n]+)", text)
    labels = re.findall(r"\blabel[=:]\s*([A-Za-z0-9_.-]+)", text)
    if not labels:
        return None
    roots = [Path(value.strip()).expanduser() for value in folders[-3:]]
    roots.append(cwd)
    suffixes = (".out",) if program == "orca" else (".log", ".out")
    for label in reversed(labels):
        for root in reversed(roots):
            for suffix in suffixes:
                candidate = (root / f"{label}{suffix}").resolve()
                if candidate.is_file():
                    return candidate
    return None


def _progress_summary(path: Path, program: str, kind: str) -> str:
    if not path.is_file():
        return "Process running"
    text = _read_tail(path, max_chars=120_000)
    if program == "orca":
        cycles = re.findall(r"^\s*(\d+)\s+-?\d+\.\d+", text, re.MULTILINE)
        if "ORCA TERMINATED NORMALLY" in text:
            return "ORCA terminated normally"
        if cycles:
            return f"SCF cycle {cycles[-1]}"
        if "VIBRATIONAL FREQUENCIES" in text:
            return "Frequency analysis"
        if "GEOMETRY OPTIMIZATION CYCLE" in text:
            matches = re.findall(r"GEOMETRY OPTIMIZATION CYCLE\s+(\d+)", text)
            return f"Optimization cycle {matches[-1]}" if matches else "Geometry optimization"
    else:
        if "Normal termination of Gaussian" in text:
            return "Gaussian terminated normally"
        cycles = re.findall(r"Cycle\s+(\d+)", text)
        if cycles:
            return f"Optimization cycle {cycles[-1]}"
        if "Frequencies --" in text:
            return "Frequency analysis"
        if "SCF Done:" in text:
            return "SCF converged"
    if "scan" in kind.lower():
        return "Coordinate scan running"
    return "Process running"


def _apply_chemistry_summary(run: CalculationRun, summary: dict[str, Any]) -> None:
    if not summary:
        return
    run.normal_termination = summary.get("normal_termination")
    run.energy = summary.get("energy")
    run.scf_cycles = summary.get("scf_cycles")
    run.optimization_cycles = summary.get("optimization_cycles")
    run.optimization_converged = summary.get("optimization_converged")
    run.frequencies = list(summary.get("frequencies") or [])
    run.imag_freqs = list(summary.get("imag_freqs") or [])
    run.frequency_count = summary.get("frequency_count")
    run.imaginary_frequency_count = summary.get(
        "imaginary_frequency_count"
    )
    run.scan_points_completed = summary.get("scan_points_completed")
    run.scan_points_total = summary.get("scan_points_total")
    run.energy_min = summary.get("energy_min")
    run.energy_max = summary.get("energy_max")
    run.neb_images = summary.get("neb_images")
    run.qmmm_qm_atoms = summary.get("qmmm_qm_atoms")
    run.qmmm_mm_atoms = summary.get("qmmm_mm_atoms")
    run.qmmm_total_atoms = summary.get("qmmm_total_atoms")
    run.path_direction = str(summary.get("path_direction") or "")
    run.chemistry_elapsed_s = summary.get("chemistry_elapsed_s")
    run.parser_backend = str(summary.get("parser_backend") or "native")
    run.parser_warnings = list(summary.get("parser_warnings") or [])
    run.atom_count = summary.get("atom_count")
    run.parsed_charge = summary.get("parsed_charge")
    run.parsed_multiplicity = summary.get("parsed_multiplicity")


def _project_metadata(
    cwd: Path,
    program: str,
    project: str,
) -> tuple[str, str, str, str]:
    if not program or not project:
        return "", "", "", ""
    path = cwd / ".chemsmart" / program / f"{project}.yaml"
    if not path.is_file():
        return str(path), "", "", ""
    try:
        import hashlib

        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        document = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError):
        return str(path), "", "", ""
    block = document.get("gas") if isinstance(document, dict) else None
    if not isinstance(block, dict):
        block = document if isinstance(document, dict) else {}
    return (
        str(path.resolve()),
        digest,
        str(block.get("functional") or block.get("method") or ""),
        str(block.get("basis") or ""),
    )


def _absolute_input_path(cwd: Path, value: str | None) -> str:
    if not value:
        return ""
    path = Path(value).expanduser()
    return str((cwd / path).resolve() if not path.is_absolute() else path.resolve())


def _label_from_input(value: str | None) -> str:
    return Path(value).stem if value else "calculation"


def _tail_file(path: Path, max_chars: int = 12_000) -> str:
    return _read_tail(path, max_chars=max_chars).strip()


def _read_output_window(path: Path) -> str:
    """Read bounded metadata from the start and results from the end."""

    try:
        size = path.stat().st_size
    except OSError:
        return ""
    tail = _read_tail(path, max_chars=2_000_000)
    if size <= 8_000_000:
        return tail
    try:
        with path.open("rb") as handle:
            head = handle.read(800_000)
    except OSError:
        return tail
    sections = _read_last_marker_windows(
        path,
        (b"VIBRATIONAL FREQUENCIES", b"Harmonic frequencies"),
    )
    return (
        head.decode("utf-8", errors="replace")
        + "\n"
        + "\n".join(sections)
        + "\n"
        + tail
    )


def _read_tail(path: Path, *, max_chars: int) -> str:
    try:
        with path.open("rb") as handle:
            handle.seek(0, os.SEEK_END)
            size = handle.tell()
            handle.seek(max(0, size - max_chars * 4))
            payload = handle.read()
    except OSError:
        return ""
    return payload.decode("utf-8", errors="replace")[-max_chars:]


def _read_last_marker_windows(
    path: Path,
    markers: tuple[bytes, ...],
    *,
    after_bytes: int = 2_000_000,
) -> list[str]:
    windows: list[str] = []
    try:
        with path.open("rb") as handle, mmap.mmap(
            handle.fileno(), 0, access=mmap.ACCESS_READ
        ) as mapped:
            for marker in markers:
                start = mapped.rfind(marker)
                if start < 0:
                    continue
                payload = mapped[start : min(len(mapped), start + after_bytes)]
                windows.append(payload.decode("utf-8", errors="replace"))
    except (OSError, ValueError):
        return []
    return windows


def _last_section(text: str, marker: str) -> str:
    start = text.rfind(marker)
    return text[start:] if start >= 0 else text


def _add_frequency_summary(
    result: dict[str, Any], frequencies: list[float]
) -> None:
    if not frequencies:
        return
    imaginary = [value for value in frequencies if value < 0]
    result["frequencies"] = frequencies[:12]
    result["imag_freqs"] = imaginary[:6]
    result["frequency_count"] = len(frequencies)
    result["imaginary_frequency_count"] = len(imaginary)


def _merge_cclib_summary(result: dict[str, Any], path: Path) -> None:
    cclib_result = parse_cclib_output(path)
    status = cclib_result.get("status")
    if status != "ok":
        result["parser_backend"] = "native"
        warning = cclib_result.get("warning")
        result["parser_warnings"] = [str(warning)] if warning else []
        return

    result["parser_backend"] = "native+cclib"
    result["parser_warnings"] = []
    for key in (
        "energy",
        "scf_cycles",
        "optimization_converged",
        "frequency_count",
        "imaginary_frequency_count",
    ):
        if result.get(key) is None and cclib_result.get(key) is not None:
            result[key] = cclib_result[key]
    if not result.get("frequencies") and cclib_result.get("frequencies"):
        frequencies = list(cclib_result["frequencies"])
        result["frequencies"] = frequencies[:12]
        result["imag_freqs"] = [value for value in frequencies if value < 0][
            :6
        ]
    result["atom_count"] = cclib_result.get("atom_count")
    result["parsed_charge"] = cclib_result.get("charge")
    result["parsed_multiplicity"] = cclib_result.get("multiplicity")


def _path_direction(text: str) -> str:
    lowered = text.lower()
    has_forward = bool(
        re.search(r"\b(?:irc|qrc)?\s*forward direction\b", lowered)
    )
    has_reverse = bool(
        re.search(r"\b(?:irc|qrc)?\s*(?:reverse|backward) direction\b", lowered)
    )
    if has_forward and has_reverse:
        return "both"
    if has_forward:
        return "forward"
    if has_reverse:
        return "reverse"
    return ""


def _program_runtime_seconds(text: str, program: str) -> float | None:
    if program == "orca":
        matches = re.findall(
            r"TOTAL RUN TIME:\s+(\d+) days\s+(\d+) hours\s+(\d+) minutes\s+"
            r"(\d+) seconds\s+(\d+) msec",
            text,
        )
        if not matches:
            return None
        days, hours, minutes, seconds, milliseconds = map(int, matches[-1])
        return (
            days * 86_400
            + hours * 3_600
            + minutes * 60
            + seconds
            + milliseconds / 1_000
        )
    matches = re.findall(
        r"Elapsed time:\s+(\d+) days\s+(\d+) hours\s+(\d+) minutes\s+"
        r"([\d.]+) seconds",
        text,
    )
    if not matches:
        return None
    days, hours, minutes, seconds = matches[-1]
    return (
        int(days) * 86_400
        + int(hours) * 3_600
        + int(minutes) * 60
        + float(seconds)
    )


def _last_error(*paths: Path) -> str:
    for path in paths:
        lines = [line.strip() for line in _read_tail(path, max_chars=20_000).splitlines() if line.strip()]
        if lines:
            return "\n".join(lines[-40:])
    return "Calculation failed without diagnostic output."


def _detect_program(text: str, path: Path) -> str:
    if "O   R   C   A" in text or "ORCA TERMINATED" in text:
        return "orca"
    if "Gaussian" in text or path.suffix.lower() == ".log":
        return "gaussian"
    return ""


def _find_run(run_id: str, *, session_root: str | None) -> CalculationRun | None:
    root = Path(session_root or Path.home() / ".chemsmart" / "agent" / "sessions")
    runs = [run for run in load_calculation_runs(root) if Path(run.cwd).resolve() == Path.cwd().resolve()]
    if not runs:
        return None
    if run_id in {"", "latest"}:
        return runs[-1]
    return next((run for run in reversed(runs) if run.run_id == run_id), None)


def _refresh_recovered_run(run: CalculationRun) -> bool:
    if run.status not in {
        CalculationStatus.VALIDATING.value,
        CalculationStatus.STARTING.value,
        CalculationStatus.RUNNING.value,
    }:
        return False
    changed = False
    if run.pid and _pid_is_alive(run.pid):
        candidate = _recovered_output_candidate(run)
        if candidate and candidate != run.output_path:
            run.output_path = candidate
            changed = True
        progress_path = Path(run.output_path) if run.output_path else Path(run.stdout_path)
        progress = _progress_summary(progress_path, run.program, run.kind)
        if progress and progress != run.stage:
            run.stage = progress
            changed = True
        elapsed = _elapsed_since(run.started_at)
        if elapsed is not None and elapsed > run.elapsed_s:
            run.elapsed_s = elapsed
            changed = True
        if changed:
            run.updated_at = datetime.now(UTC).isoformat()
        return changed

    candidate = _recovered_output_candidate(run)
    if candidate:
        run.output_path = candidate
    summary = inspect_output(
        run.output_path,
        program=run.program,
        kind=run.kind,
    )
    _apply_chemistry_summary(run, summary)
    if run.normal_termination is True:
        run.status = CalculationStatus.COMPLETED.value
        run.stage = "Calculation completed while TUI was detached"
    elif run.output_path:
        run.status = CalculationStatus.CHEMISTRY_FAILED.value
        run.stage = f"{run.program.upper() or 'Chemistry'} did not terminate normally"
        run.error = _last_error(Path(run.output_path), Path(run.stderr_path))
    else:
        run.status = CalculationStatus.PROCESS_FAILED.value
        run.stage = "Detached calculation process is no longer running"
        run.error = _last_error(Path(run.stderr_path), Path(run.stdout_path))
    run.ended_at = datetime.now(UTC).isoformat()
    run.updated_at = run.ended_at
    return True


def _pid_is_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _recovered_output_candidate(run: CalculationRun) -> str:
    if run.output_path and Path(run.output_path).is_file():
        return str(Path(run.output_path).resolve())
    roots = {Path(run.cwd)}
    if run.input_path:
        roots.add(Path(run.input_path).parent)
    label = run.label.lower()
    candidates = [
        path
        for root in roots
        if root.is_dir()
        for suffix in ("*.out", "*.log")
        for path in root.glob(suffix)
        if path.is_file()
    ]
    if not candidates:
        return ""
    matching = [path for path in candidates if label and label in path.stem.lower()]
    selected = max(
        matching or candidates,
        key=lambda path: path.stat().st_mtime_ns,
    )
    return str(selected.resolve())


def _elapsed_since(value: str) -> float | None:
    if not value:
        return None
    try:
        started = datetime.fromisoformat(value)
    except ValueError:
        return None
    if started.tzinfo is None:
        started = started.replace(tzinfo=UTC)
    return max(0.0, (datetime.now(UTC) - started).total_seconds())


def _latest_output_in_workspace(cwd: Path) -> str:
    outputs = [path for suffix in ("*.out", "*.log") for path in cwd.glob(suffix) if path.is_file()]
    if not outputs:
        return ""
    return str(max(outputs, key=lambda path: path.stat().st_mtime_ns).resolve())


__all__ = [
    "CalculationContext",
    "CalculationEvent",
    "CalculationRun",
    "CalculationStatus",
    "CalculationStore",
    "cancel_calculation",
    "execute_observed_process",
    "inspect_calculation",
    "inspect_output",
    "load_calculation_runs",
    "reset_calculation_context",
    "set_calculation_context",
]
