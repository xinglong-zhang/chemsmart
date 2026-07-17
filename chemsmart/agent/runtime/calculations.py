"""Observable lifecycle and chemistry summaries for CLI calculations."""

from __future__ import annotations

import json
import os
import re
import signal
import subprocess
import threading
import time
import uuid
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

import yaml

from chemsmart.agent.model_command_parser import (
    ParsedModelCommand,
    parse_model_command,
)
from chemsmart.agent.runtime.calculation_models import (
    UTC,
    CalculationContext,
    CalculationEvent,
    CalculationRun,
    CalculationStatus,
    CalculationStore,
    EventSink,
    apply_chemistry_summary as _apply_chemistry_summary,
    current_calculation_context,
    reset_calculation_context,
    set_calculation_context,
)
from chemsmart.agent.runtime.output_io import read_tail as _read_tail
from chemsmart.agent.runtime.result_parsing import inspect_output


_ACTIVE_PROCESSES: dict[str, subprocess.Popen[str]] = {}
_ACTIVE_LOCK = threading.Lock()
_CANCELLED_RUNS: set[str] = set()


@dataclass(slots=True)
class _ExecutionSetup:
    run: CalculationRun
    parsed: ParsedModelCommand
    store: CalculationStore | None
    workdir: Path
    stdout_path: Path
    stderr_path: Path
    outputs_before: dict[str, tuple[int, int]]


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

    setup = _prepare_execution(
        argv,
        command=command,
        context=context,
        cwd=cwd,
    )
    run = setup.run
    started_monotonic = time.monotonic()
    _emit(run, "validating", setup.store, event_sink)
    run.status = CalculationStatus.STARTING.value
    run.stage = "Starting chemsmart command"
    run.started_at = datetime.now(UTC).isoformat()
    run.updated_at = run.started_at
    _emit(run, "starting", setup.store, event_sink)
    returncode, timed_out = _run_process(
        argv,
        setup,
        timeout_s=timeout_s,
        started_monotonic=started_monotonic,
        event_sink=event_sink,
    )
    _finalize_execution(
        setup,
        returncode=returncode,
        timed_out=timed_out,
        timeout_s=timeout_s,
        started_monotonic=started_monotonic,
    )
    run.ended_at = datetime.now(UTC).isoformat()
    run.updated_at = run.ended_at
    _emit(
        run,
        "completed" if run.status == "completed" else "failed",
        setup.store,
        event_sink,
    )
    return {
        "calculation": run.to_dict(),
        "stdout_tail": _tail_file(setup.stdout_path),
        "stderr_tail": _tail_file(setup.stderr_path),
    }


def _prepare_execution(
    argv: list[str],
    *,
    command: str,
    context: CalculationContext | None,
    cwd: str | Path | None,
) -> _ExecutionSetup:
    workdir = Path(cwd or Path.cwd()).resolve()
    active_context = (
        context or current_calculation_context() or CalculationContext()
    )
    parsed = parse_model_command(command, cwd=str(workdir))
    run_id = f"calc-{uuid.uuid4().hex[:10]}"
    store = (
        CalculationStore(active_context.session_dir)
        if active_context.session_dir
        else None
    )
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
        session_id=active_context.session_id,
        turn_id=active_context.turn_id,
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
        execution_mode=_execution_mode(argv, parsed),
        semantic_verdict=active_context.semantic_verdict,
        intent_verdict=active_context.intent_verdict,
    )
    stdout_path = run_dir / "stdout.log"
    stderr_path = run_dir / "stderr.log"
    run.stdout_path = str(stdout_path)
    run.stderr_path = str(stderr_path)
    return _ExecutionSetup(
        run=run,
        parsed=parsed,
        store=store,
        workdir=workdir,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        outputs_before=_output_fingerprints(workdir, run.input_path),
    )


def _execution_mode(argv: list[str], parsed: ParsedModelCommand) -> str:
    if "--fake" in argv or "--test" in argv:
        return "test_fake"
    return "hpc" if parsed.action == "sub" else "local"


def _run_process(
    argv: list[str],
    setup: _ExecutionSetup,
    *,
    timeout_s: int,
    started_monotonic: float,
    event_sink: EventSink | None,
) -> tuple[int, bool]:
    with (
        setup.stdout_path.open("w", encoding="utf-8") as stdout_handle,
        setup.stderr_path.open("w", encoding="utf-8") as stderr_handle,
    ):
        process = subprocess.Popen(
            argv,
            cwd=str(setup.workdir),
            text=True,
            stdout=stdout_handle,
            stderr=stderr_handle,
            start_new_session=True,
        )
        _mark_process_started(setup, process, event_sink)
        try:
            timed_out = _monitor_process(
                setup,
                process,
                timeout_s=timeout_s,
                started_monotonic=started_monotonic,
                event_sink=event_sink,
            )
            process.wait()
        finally:
            _unregister_process(setup.run.run_id)
    return int(process.returncode or 0), timed_out


def _mark_process_started(
    setup: _ExecutionSetup,
    process: subprocess.Popen[str],
    event_sink: EventSink | None,
) -> None:
    run = setup.run
    run.pid = process.pid
    run.status = CalculationStatus.RUNNING.value
    run.stage = "Process running"
    _register_process(run.run_id, process)
    _emit(run, "started", setup.store, event_sink)


def _monitor_process(
    setup: _ExecutionSetup,
    process: subprocess.Popen[str],
    *,
    timeout_s: int,
    started_monotonic: float,
    event_sink: EventSink | None,
) -> bool:
    last_progress = ""
    last_emit = 0.0
    while process.poll() is None:
        elapsed = time.monotonic() - started_monotonic
        setup.run.elapsed_s = elapsed
        if elapsed > timeout_s:
            _terminate_process(process)
            return True
        _refresh_output_path(setup)
        progress = _progress_summary(
            Path(setup.run.output_path)
            if setup.run.output_path
            else setup.stdout_path,
            setup.run.program,
            setup.run.kind,
        )
        now = time.monotonic()
        if progress and progress != last_progress:
            setup.run.stage = progress
            last_progress = progress
            _emit(setup.run, "progress", setup.store, event_sink)
            last_emit = now
        elif now - last_emit >= 2.0:
            _emit(setup.run, "heartbeat", setup.store, event_sink)
            last_emit = now
        time.sleep(0.2)
    return False


def _refresh_output_path(setup: _ExecutionSetup) -> None:
    output_path = _discover_output(
        setup.workdir,
        setup.run.input_path,
        setup.outputs_before,
    )
    if output_path is not None:
        setup.run.output_path = str(output_path)


def _finalize_execution(
    setup: _ExecutionSetup,
    *,
    returncode: int,
    timed_out: bool,
    timeout_s: int,
    started_monotonic: float,
) -> None:
    run = setup.run
    run.returncode = returncode
    run.elapsed_s = time.monotonic() - started_monotonic
    _refresh_output_path(setup)
    if not run.output_path and returncode == 0:
        reported_output = _reported_output_path(
            setup.stdout_path,
            setup.stderr_path,
            cwd=setup.workdir,
            program=run.program,
        )
        if reported_output is not None:
            run.output_path = str(reported_output)
            run.reused_output = str(reported_output) in setup.outputs_before
    _apply_chemistry_summary(
        run,
        inspect_output(run.output_path, program=run.program, kind=run.kind),
    )
    _classify_completion(setup, timed_out=timed_out, timeout_s=timeout_s)


def _classify_completion(
    setup: _ExecutionSetup,
    *,
    timed_out: bool,
    timeout_s: int,
) -> None:
    run = setup.run
    if timed_out:
        run.status = CalculationStatus.TIMEOUT.value
        run.stage = f"Timed out after {timeout_s}s"
        run.error = run.stage
    elif _consume_cancelled(run.run_id):
        run.status = CalculationStatus.CANCELLED.value
        run.stage = "Cancelled by user"
    elif run.returncode != 0:
        run.status = CalculationStatus.PROCESS_FAILED.value
        run.stage = f"Process exited with code {run.returncode}"
        run.error = _last_error(setup.stderr_path, setup.stdout_path)
    elif setup.parsed.action == "sub":
        run.status = CalculationStatus.COMPLETED.value
        run.stage = "Submission command completed"
    elif run.execution_mode == "test_fake":
        run.status = CalculationStatus.COMPLETED.value
        run.stage = "Safe test command completed"
    elif not run.output_path:
        run.status = CalculationStatus.CHEMISTRY_FAILED.value
        run.stage = "No Gaussian/ORCA output file was produced"
        run.error = _last_error(setup.stderr_path, setup.stdout_path)
    elif run.normal_termination is not True:
        run.status = CalculationStatus.CHEMISTRY_FAILED.value
        run.stage = (
            f"{run.program.upper() or 'Chemistry'} did not terminate normally"
        )
        run.error = _last_error(Path(run.output_path), setup.stderr_path)
    else:
        run.status = CalculationStatus.COMPLETED.value
        run.stage = (
            "Existing completed calculation reused"
            if run.reused_output
            else "Calculation completed"
        )


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
        payload.update(
            {key: value for key, value in summary.items() if value is not None}
        )
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


def _output_fingerprints(
    cwd: Path, input_path: str
) -> dict[str, tuple[int, int]]:
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
    return (
        max(changed, key=lambda path: path.stat().st_mtime_ns)
        if changed
        else None
    )


def _reported_output_path(
    stdout_path: Path,
    stderr_path: Path,
    *,
    cwd: Path,
    program: str,
) -> Path | None:
    text = (
        _read_tail(stdout_path, max_chars=120_000)
        + "\n"
        + _read_tail(stderr_path, max_chars=40_000)
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
    progress = (
        _orca_progress(text) if program == "orca" else _gaussian_progress(text)
    )
    if progress:
        return progress
    if "scan" in kind.lower():
        return "Coordinate scan running"
    return "Process running"


def _orca_progress(text: str) -> str:
    cycles = re.findall(r"^\s*(\d+)\s+-?\d+\.\d+", text, re.MULTILINE)
    if "ORCA TERMINATED NORMALLY" in text:
        return "ORCA terminated normally"
    if cycles:
        return f"SCF cycle {cycles[-1]}"
    if "VIBRATIONAL FREQUENCIES" in text:
        return "Frequency analysis"
    if "GEOMETRY OPTIMIZATION CYCLE" in text:
        matches = re.findall(r"GEOMETRY OPTIMIZATION CYCLE\s+(\d+)", text)
        return (
            f"Optimization cycle {matches[-1]}"
            if matches
            else "Geometry optimization"
        )
    return ""


def _gaussian_progress(text: str) -> str:
    if "Normal termination of Gaussian" in text:
        return "Gaussian terminated normally"
    cycles = re.findall(r"Cycle\s+(\d+)", text)
    if cycles:
        return f"Optimization cycle {cycles[-1]}"
    if "Frequencies --" in text:
        return "Frequency analysis"
    return "SCF converged" if "SCF Done:" in text else ""


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
    return str(
        (cwd / path).resolve() if not path.is_absolute() else path.resolve()
    )


def _label_from_input(value: str | None) -> str:
    return Path(value).stem if value else "calculation"


def _tail_file(path: Path, max_chars: int = 12_000) -> str:
    return _read_tail(path, max_chars=max_chars).strip()


def _last_error(*paths: Path) -> str:
    for path in paths:
        lines = [
            line.strip()
            for line in _read_tail(path, max_chars=20_000).splitlines()
            if line.strip()
        ]
        if lines:
            return "\n".join(lines[-40:])
    return "Calculation failed without diagnostic output."


def _find_run(
    run_id: str, *, session_root: str | None
) -> CalculationRun | None:
    root = Path(
        session_root or Path.home() / ".chemsmart" / "agent" / "sessions"
    )
    runs = [
        run
        for run in load_calculation_runs(root)
        if Path(run.cwd).resolve() == Path.cwd().resolve()
    ]
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
        progress_path = (
            Path(run.output_path) if run.output_path else Path(run.stdout_path)
        )
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
        run.stage = (
            f"{run.program.upper() or 'Chemistry'} did not terminate normally"
        )
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
    matching = [
        path for path in candidates if label and label in path.stem.lower()
    ]
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
    outputs = [
        path
        for suffix in ("*.out", "*.log")
        for path in cwd.glob(suffix)
        if path.is_file()
    ]
    if not outputs:
        return ""
    return str(
        max(outputs, key=lambda path: path.stat().st_mtime_ns).resolve()
    )


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
