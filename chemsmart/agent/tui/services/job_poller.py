"""Background job polling and session-derived job helpers."""

from __future__ import annotations

import asyncio
import json
import os
import subprocess
import traceback
from concurrent.futures import Future
from dataclasses import dataclass, field
from datetime import datetime, timezone

UTC = timezone.utc
from pathlib import Path
from threading import Lock, Thread
from typing import Any

from textual import work
from textual.message import Message

from chemsmart.agent.core import (
    DecisionLog,
    SessionState,
    _resolve_refs,
    _restore_json_result,
)
from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.io.orca.output import ORCAOutput
from chemsmart.jobs.job import Job
from chemsmart.settings.server import Server
from chemsmart.settings.user import ChemsmartUserSettings
from chemsmart.utils.cluster import ClusterHelper
from chemsmart.utils.io import get_program_type_from_file

_JOB_SNAPSHOT_CACHE_LOCK = Lock()


@dataclass
class _JobSnapshotCacheState:
    snapshot: dict[str, dict[str, Any]] = field(default_factory=dict)
    last_error: Exception | None = None
    future: Future | None = None


_JOB_SNAPSHOT_CACHE: dict[str, _JobSnapshotCacheState] = {}


class JobStatusUpdated(Message):
    """Incremental job-status update emitted by the background poller."""

    def __init__(self, job_id: str, fields: dict[str, Any]) -> None:
        super().__init__()
        self.job_id = job_id
        self.fields = fields


class JobPollerError(Message):
    """Background polling failure surfaced to the transcript."""

    def __init__(self, summary: str, details: str) -> None:
        super().__init__()
        self.summary = summary
        self.details = details


class JobPollerMixin:
    session_root: Path

    @work(
        exclusive=True,
        exit_on_error=False,
        group="job-poller",
        name="job-poller",
    )
    async def run_job_poller(self, poll_interval: float = 5.0) -> None:
        previous: dict[str, dict[str, Any]] = {}
        failures = 0
        last_error_signature: str | None = None
        while self.is_mounted:
            try:
                snapshot = await asyncio.to_thread(
                    refresh_job_snapshot_cache,
                    self.session_root,
                )
                failures = 0
                last_error_signature = None
                stale_job_ids = sorted(set(previous) - set(snapshot))
                for stale_job_id in stale_job_ids:
                    self._emit_job_error(
                        "Removed stale job snapshot entry",
                        f"{stale_job_id} no longer exists in the current snapshot.",
                    )
                for job_id, current in snapshot.items():
                    changed = {
                        key: value
                        for key, value in current.items()
                        if previous.get(job_id, {}).get(key) != value
                    }
                    if changed:
                        self._emit_job_update(job_id, changed)
                previous = snapshot
                await asyncio.sleep(poll_interval)
            except Exception as exc:  # pragma: no cover - defensive backoff
                signature = f"{exc.__class__.__name__}:{exc}"
                if signature != last_error_signature:
                    self._emit_job_error(
                        "Job poller failed", traceback.format_exc()
                    )
                    last_error_signature = signature
                failures += 1
                await asyncio.sleep(min(poll_interval * (2**failures), 30.0))

    def _emit_job_error(self, summary: str, details: str) -> None:
        self.post_message(JobPollerError(summary, details))


class JobStateReader:
    """Compatibility loader for agent session state files."""

    @staticmethod
    def load(session_dir: Path) -> SessionState | None:
        for name in ("state.json", "session.json"):
            path = session_dir / name
            if path.exists():
                try:
                    return SessionState.load(path)
                except Exception as exc:
                    raise RuntimeError(
                        f"Malformed session state: {path}"
                    ) from exc
        return None


def available_server_names() -> list[str]:
    return sorted(ChemsmartUserSettings().all_available_servers)


def get_cached_job_snapshot(
    session_root: str | os.PathLike[str],
) -> dict[str, dict[str, Any]]:
    key = _job_snapshot_cache_key(session_root)
    with _JOB_SNAPSHOT_CACHE_LOCK:
        state = _JOB_SNAPSHOT_CACHE.get(key)
        if state is None:
            return {}
        return _copy_job_snapshot(state.snapshot)


def request_job_snapshot_refresh(
    session_root: str | os.PathLike[str],
) -> Future | None:
    root_path = Path(session_root)
    key = _job_snapshot_cache_key(root_path)
    with _JOB_SNAPSHOT_CACHE_LOCK:
        state = _JOB_SNAPSHOT_CACHE.setdefault(key, _JobSnapshotCacheState())
        future = state.future
        if future is not None and not future.done():
            return future
        future = Future()
        state.future = future
        future.add_done_callback(
            lambda done, cache_key=key: _clear_job_snapshot_future(
                cache_key, done
            )
        )
        Thread(
            target=_run_job_snapshot_refresh,
            args=(future, root_path),
            daemon=True,
            name="chemsmart-job-snapshot",
        ).start()
        return future


def refresh_job_snapshot_cache(
    session_root: str | os.PathLike[str],
) -> dict[str, dict[str, Any]]:
    root_path = Path(session_root)
    key = _job_snapshot_cache_key(root_path)
    try:
        snapshot = collect_job_snapshot(root_path)
    except Exception as exc:
        with _JOB_SNAPSHOT_CACHE_LOCK:
            state = _JOB_SNAPSHOT_CACHE.setdefault(
                key, _JobSnapshotCacheState()
            )
            state.last_error = exc
        raise
    with _JOB_SNAPSHOT_CACHE_LOCK:
        state = _JOB_SNAPSHOT_CACHE.setdefault(key, _JobSnapshotCacheState())
        state.snapshot = _copy_job_snapshot(snapshot)
        state.last_error = None
    return snapshot


def _job_snapshot_cache_key(session_root: str | os.PathLike[str]) -> str:
    return str(Path(session_root).resolve())


def _copy_job_snapshot(
    snapshot: dict[str, dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    return {
        job_id: dict(fields) if isinstance(fields, dict) else fields
        for job_id, fields in snapshot.items()
    }


def _clear_job_snapshot_future(cache_key: str, future: Future) -> None:
    with _JOB_SNAPSHOT_CACHE_LOCK:
        state = _JOB_SNAPSHOT_CACHE.get(cache_key)
        if state is not None and state.future is future:
            state.future = None


def _run_job_snapshot_refresh(
    future: Future,
    session_root: Path,
) -> None:
    if not future.set_running_or_notify_cancel():
        return
    try:
        result = refresh_job_snapshot_cache(session_root)
    except Exception as exc:
        future.set_exception(exc)
    else:
        future.set_result(result)


def collect_job_snapshot(session_root: Path) -> dict[str, dict[str, Any]]:
    session_root = Path(session_root)
    cluster_running_ids, cluster_running_names = _cluster_running_jobs()
    snapshot: dict[str, dict[str, Any]] = {}
    if not session_root.exists():
        return snapshot

    for session_dir in sorted(session_root.iterdir(), reverse=True):
        if not session_dir.is_dir():
            continue
        session_jobs = _session_jobs(
            session_dir,
            cluster_running_ids=cluster_running_ids,
            cluster_running_names=cluster_running_names,
        )
        snapshot.update(session_jobs)
    return snapshot


def queue_snapshot(
    jobs: dict[str, dict[str, Any]],
    *,
    server_name: str | None = None,
) -> list[dict[str, Any]]:
    rows = [
        job
        for job in jobs.values()
        if job.get("status") in {"queued", "running"}
        and (
            server_name is None
            or job.get("server_name") == server_name
            or job.get("host") == server_name
        )
    ]
    return sorted(rows, key=_job_sort_key)


def cancel_job(job: dict[str, Any]) -> dict[str, Any]:
    scheduler = str(job.get("scheduler") or "local").lower()
    job_id = str(job.get("job_id") or "").strip()
    if not job_id:
        raise ValueError("Job has no scheduler job id.")
    if scheduler == "local":
        raise ValueError("Local jobs cannot be cancelled via a scheduler.")
    if scheduler == "slurm":
        command = ["scancel", _scheduler_job_token(job_id)]
    elif scheduler == "pbs":
        command = ["qdel", _scheduler_job_token(job_id)]
    else:
        raise ValueError(f"Unsupported scheduler: {scheduler}")
    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        message = (result.stderr or result.stdout or "cancel failed").strip()
        raise RuntimeError(message)
    return {
        "command": command,
        "stdout": (result.stdout or "").strip(),
        "stderr": (result.stderr or "").strip(),
    }


def extract_run_result(
    value: str,
    jobs: dict[str, dict[str, Any]] | None = None,
) -> dict[str, Any]:
    output_path = _resolve_result_output_path(value, jobs or {})
    parser = _parser_for_output(output_path)
    molecule = parser.get_molecule(index="-1")
    frequencies = list(getattr(parser, "vibrational_frequencies", []) or [])
    imag_freqs = [float(freq) for freq in frequencies if float(freq) < 0]
    real_freqs = sorted(
        float(freq) for freq in frequencies if float(freq) >= 0
    )
    energies = list(getattr(parser, "energies", []) or [])
    energy = float(energies[-1]) if energies else None
    gibbs = getattr(parser, "gibbs_free_energy", None)
    if gibbs is not None:
        gibbs = float(gibbs)
    return {
        "source": str(value),
        "output_path": str(output_path),
        "input_path": _guess_input_path(output_path),
        "molecule": molecule,
        "energy": energy,
        "delta_g": gibbs,
        "frequencies": real_freqs[:3],
        "imag_freqs": imag_freqs,
        "normal_termination": bool(
            getattr(parser, "normal_termination", False)
        ),
    }


def format_jobs_table(rows: list[dict[str, Any]]) -> str:
    if not rows:
        return "No jobs found."
    headers = [
        "job_id",
        "name",
        "scheduler",
        "status",
        "started",
        "runtime",
        "host",
    ]
    widths = {
        header: max(
            [
                len(header),
                *[len(_display_job_value(row, header)) for row in rows],
            ]
        )
        for header in headers
    }
    line = " ".join(header.ljust(widths[header]) for header in headers)
    sep = " ".join("-" * widths[header] for header in headers)
    body = [
        " ".join(
            _display_job_value(row, header).ljust(widths[header])
            for header in headers
        )
        for row in rows
    ]
    return "\n".join([line, sep, *body])


def _session_jobs(
    session_dir: Path,
    *,
    cluster_running_ids: set[int],
    cluster_running_names: set[str],
) -> dict[str, dict[str, Any]]:
    state = JobStateReader.load(session_dir)
    if state is None or state.plan is None:
        return {}

    entries = DecisionLog(session_dir / "decision_log.jsonl").read_all()
    step_calls, step_results, step_errors = _step_events(entries)
    prior_results = _load_results(session_dir, state.current_step_index)
    jobs: dict[str, dict[str, Any]] = {}

    for step_index, step in enumerate(state.plan.steps, start=1):
        if step.tool not in {"run_local", "submit_hpc"}:
            continue
        has_call = step_index in step_calls
        result_entry = step_results.get(step_index)
        if (
            step.tool == "submit_hpc"
            and result_entry is not None
            and result_entry.get("payload", {}).get("from_preview")
        ):
            continue
        if not (
            has_call or result_entry is not None or step_index in step_errors
        ):
            continue
        try:
            resolved = _resolve_refs(
                step.args, prior_results[: step_index - 1]
            )
        except Exception:
            resolved = dict(step.args)
        job = resolved.get("job") if isinstance(resolved, dict) else None
        if not isinstance(job, Job):
            continue
        server_name = _resolved_server_name(resolved.get("server"))
        scheduler = _scheduler_name(step.tool, server_name)
        record = _job_record(
            session_dir=session_dir,
            state=state,
            step_index=step_index,
            job=job,
            server_name=server_name,
            scheduler=scheduler,
            call_entry=step_calls.get(step_index),
            result_entry=result_entry,
            error_entry=step_errors.get(step_index),
            cluster_running_ids=cluster_running_ids,
            cluster_running_names=cluster_running_names,
        )
        jobs[record["job_id"]] = record
    return jobs


def _job_record(
    *,
    session_dir: Path,
    state: SessionState,
    step_index: int,
    job: Job,
    server_name: str | None,
    scheduler: str,
    call_entry: dict[str, Any] | None,
    result_entry: dict[str, Any] | None,
    error_entry: dict[str, Any] | None,
    cluster_running_ids: set[int],
    cluster_running_names: set[str],
) -> dict[str, Any]:
    result_payload = (result_entry or {}).get("payload") or {}
    error_payload = (error_entry or {}).get("payload") or {}
    job_id = str(
        result_payload.get("payload", {}).get("job_id")
        or f"local:{state.session_id}:{step_index}"
    )
    output_path = _job_output_path(job)
    summary = _output_summary(output_path)
    started = (
        (call_entry or {}).get("payload", {}).get("ts_start")
        or (call_entry or {}).get("ts")
        or state.started_at
    )
    runtime_ms = (
        result_payload.get("step_wall_time_ms")
        or error_payload.get("step_wall_time_ms")
        or _runtime_ms_from_started(started)
    )
    status = _infer_status(
        scheduler=scheduler,
        job_id=job_id,
        job_name=str(job.label),
        result_payload=result_payload.get("payload") or {},
        error_payload=error_payload,
        summary=summary,
        cluster_running_ids=cluster_running_ids,
        cluster_running_names=cluster_running_names,
    )
    return {
        "job_id": job_id,
        "session_id": state.session_id,
        "name": str(job.label),
        "scheduler": scheduler,
        "status": status,
        "started": _format_started(started),
        "runtime": _format_runtime(runtime_ms),
        "host": _job_host(server_name, result_payload.get("payload") or {}),
        "server_name": server_name,
        "session_dir": str(session_dir),
        "cwd": state.cwd,
        "input_path": getattr(job, "inputfile", None),
        "output_path": output_path,
        "raw_started": started,
        "runtime_ms": runtime_ms,
    }


def _step_events(entries: list[dict[str, Any]]):
    step_calls: dict[int, dict[str, Any]] = {}
    step_results: dict[int, dict[str, Any]] = {}
    step_errors: dict[int, dict[str, Any]] = {}
    for entry in entries:
        payload = entry.get("payload") or {}
        step_index = int(payload.get("step_index") or -1)
        tool = payload.get("tool")
        if step_index < 1 or tool not in {"run_local", "submit_hpc"}:
            continue
        kind = entry.get("kind")
        if kind == "tool_call":
            step_calls[step_index] = entry
        elif kind == "tool_result":
            step_results[step_index] = entry
        elif kind == "tool_error":
            step_errors[step_index] = entry
    return step_calls, step_results, step_errors


def _load_results(session_dir: Path, count: int) -> list[Any]:
    results: list[Any] = []
    for step_index in range(1, count + 1):
        artifact = session_dir / f"step_{step_index:02d}.json"
        if not artifact.exists():
            break
        try:
            with artifact.open(encoding="utf-8") as handle:
                results.append(_restore_json_result(json.load(handle)))
        except Exception as exc:
            raise RuntimeError(
                f"Malformed result artifact: {artifact}"
            ) from exc
    return results


def _cluster_running_jobs() -> tuple[set[int], set[str]]:
    try:
        running_ids, running_names = (
            ClusterHelper().get_gaussian_running_jobs()
        )
    except Exception as exc:
        raise RuntimeError("Scheduler queue lookup failed") from exc
    return set(int(job_id) for job_id in running_ids), set(running_names)


def _resolved_server_name(value: Any) -> str | None:
    if isinstance(value, Server):
        return value.name
    if isinstance(value, str) and value.strip():
        return value.strip()
    return None


def _scheduler_name(tool_name: str, server_name: str | None) -> str:
    if tool_name == "run_local":
        return "local"
    if not server_name:
        return "local"
    try:
        return str(
            Server.from_servername(server_name).scheduler or "local"
        ).lower()
    except Exception:
        return "local"


def _job_output_path(job: Job) -> str | None:
    outputfile = getattr(job, "outputfile", None)
    if isinstance(outputfile, str) and outputfile:
        return os.path.abspath(outputfile)
    folder = getattr(job, "folder", None)
    label = getattr(job, "label", None)
    if not folder or not label:
        return None
    for suffix in (".log", ".out"):
        path = os.path.join(folder, f"{label}{suffix}")
        if os.path.exists(path):
            return os.path.abspath(path)
    return os.path.abspath(os.path.join(folder, f"{label}.log"))


def _output_summary(output_path: str | None) -> dict[str, Any]:
    if not output_path or not os.path.exists(output_path):
        return {}
    try:
        parser = _parser_for_output(Path(output_path))
        frequencies = list(
            getattr(parser, "vibrational_frequencies", []) or []
        )
        energies = list(getattr(parser, "energies", []) or [])
        return {
            "normal_termination": bool(
                getattr(parser, "normal_termination", False)
            ),
            "energy": float(energies[-1]) if energies else None,
            "imag_freqs": [
                float(freq) for freq in frequencies if float(freq) < 0
            ],
        }
    except Exception:
        return {}


def _infer_status(
    *,
    scheduler: str,
    job_id: str,
    job_name: str,
    result_payload: dict[str, Any],
    error_payload: dict[str, Any],
    summary: dict[str, Any],
    cluster_running_ids: set[int],
    cluster_running_names: set[str],
) -> str:
    if error_payload:
        return "failed"
    if summary.get("normal_termination"):
        return "done"
    if summary and summary.get("normal_termination") is False:
        return "failed"
    if scheduler == "local":
        return "done" if result_payload.get("ok", True) else "failed"
    if (
        _scheduler_job_number(job_id) in cluster_running_ids
        or job_name in cluster_running_names
    ):
        return "running"
    if result_payload.get("job_id"):
        return "queued"
    return "running"


def _job_host(server_name: str | None, result_payload: dict[str, Any]) -> str:
    if server_name:
        return server_name
    command = str(result_payload.get("command_executed") or "")
    if command.startswith("ssh "):
        parts = command.split()
        if len(parts) > 1:
            return parts[1]
    return "localhost"


def _runtime_ms_from_started(started: str | None) -> int | None:
    if not started:
        return None
    try:
        then = datetime.fromisoformat(started.replace("Z", "+00:00"))
    except ValueError:
        return None
    return int(
        (datetime.now(UTC) - then.astimezone(UTC)).total_seconds() * 1000
    )


def _format_started(value: str | None) -> str:
    if not value:
        return ""
    try:
        started = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return value
    return started.astimezone().strftime("%Y-%m-%d %H:%M")


def _format_runtime(runtime_ms: int | None) -> str:
    if runtime_ms is None:
        return ""
    total_seconds = max(0, int(runtime_ms // 1000))
    hours, remainder = divmod(total_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    if hours:
        return f"{hours}h {minutes:02d}m"
    if minutes:
        return f"{minutes}m {seconds:02d}s"
    return f"{seconds}s"


def _display_job_value(row: dict[str, Any], header: str) -> str:
    value = row.get(header) if isinstance(row, dict) else None
    return str(value) if value is not None and value != "" else "—"


def _job_sort_key(job: dict[str, Any]) -> tuple[int, str, str]:
    status_order = {
        "running": 0,
        "queued": 1,
        "failed": 2,
        "cancelled": 3,
        "done": 4,
    }
    return (
        status_order.get(str(job.get("status")), 9),
        str(job.get("raw_started") or ""),
        str(job.get("job_id") or ""),
    )


def _resolve_result_output_path(
    value: str,
    jobs: dict[str, dict[str, Any]],
) -> Path:
    value = value.strip()
    if value in jobs and jobs[value].get("output_path"):
        return Path(str(jobs[value]["output_path"]))
    candidate = Path(value).expanduser()
    if candidate.exists():
        if candidate.suffix.lower() in {".log", ".out"}:
            return candidate
        guessed = _guess_output_path(candidate)
        if guessed.exists():
            return guessed
    raise FileNotFoundError(f"Could not resolve a result file for {value!r}")


def _guess_output_path(path: Path) -> Path:
    if path.suffix.lower() in {".com", ".gjf"}:
        return path.with_suffix(".log")
    if path.suffix.lower() == ".inp":
        return path.with_suffix(".out")
    return path


def _guess_input_path(output_path: Path) -> str | None:
    stem = output_path.with_suffix("")
    for suffix in (".com", ".gjf", ".inp"):
        candidate = stem.with_suffix(suffix)
        if candidate.exists():
            return str(candidate)
    return None


def _parser_for_output(output_path: Path):
    suffix = output_path.suffix.lower()
    if suffix == ".log":
        return Gaussian16Output(str(output_path))
    if suffix == ".out":
        program = get_program_type_from_file(str(output_path))
        if program == "orca":
            return ORCAOutput(str(output_path))
        return Gaussian16Output(str(output_path))
    raise ValueError(f"Unsupported result file: {output_path}")


def _scheduler_job_number(job_id: str) -> int | None:
    token = _scheduler_job_token(job_id)
    return int(token) if token.isdigit() else None


def _scheduler_job_token(job_id: str) -> str:
    return str(job_id).split(".")[0]
