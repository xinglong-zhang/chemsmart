"""Bounded subprocess execution with public resource evidence.

The helper in this module owns no chemistry semantics.  It observes an
already-launched process, samples the resident memory of the process tree, and
terminates owned process groups when a time or memory boundary is crossed.
It deliberately returns an explicit ambiguous state when the complete tree
cannot be shown to have stopped.

Only typed ``Popen`` invocations are supported by callers; this module never
accepts or evaluates shell text.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
import os
import shutil
import signal
import subprocess
import time
from typing import Any


PROCESS_OBSERVATION_SCHEMA_VERSION = "chemsmart.process-observation.v1"


def _canonical_sha256(payload: dict[str, Any]) -> str:
    body = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    )
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class ProcessObservationV1:
    """Observable launch, resource, timeout, and termination facts."""

    schema_version: str
    state: str
    pid: int | None
    process_group_id: int | None
    process_group_owned: bool
    timeout_seconds: float
    memory_limit_mb: float | None
    returncode: int | None
    timed_out: bool
    memory_limit_exceeded: bool
    termination_requested: bool
    termination_signal: str
    termination_confirmed: bool | None
    termination_scope: str
    remaining_process_ids: tuple[int, ...]
    memory_observation_state: str
    last_rss_mb: float | None
    peak_rss_mb: float | None
    memory_sample_count: int
    wall_seconds: float
    findings: tuple[str, ...]
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != PROCESS_OBSERVATION_SCHEMA_VERSION:
            raise ValueError("unsupported process observation schema")
        if self.state not in {
            "exited",
            "launch_failed",
            "timed_out_terminated",
            "timed_out_ambiguous",
            "memory_limit_exceeded_terminated",
            "memory_limit_exceeded_ambiguous",
        }:
            raise ValueError("unsupported process observation state")
        if not math.isfinite(self.timeout_seconds) or self.timeout_seconds <= 0:
            raise ValueError("timeout_seconds must be finite and positive")
        if self.memory_limit_mb is not None and (
            not math.isfinite(self.memory_limit_mb)
            or self.memory_limit_mb <= 0
        ):
            raise ValueError("memory_limit_mb must be finite and positive")
        if self.memory_observation_state not in {
            "observed",
            "unavailable",
        }:
            raise ValueError("unsupported memory observation state")
        if self.memory_sample_count < 0:
            raise ValueError("memory_sample_count cannot be negative")
        if not math.isfinite(self.wall_seconds) or self.wall_seconds < 0:
            raise ValueError("wall_seconds must be finite and non-negative")
        if self.findings != tuple(sorted(set(self.findings))):
            raise ValueError("process findings must be sorted and unique")
        body = self.as_dict(include_digest=False)
        if self.receipt_sha256 != _canonical_sha256(body):
            raise ValueError("process observation digest mismatch")

    def as_dict(self, *, include_digest: bool = True) -> dict[str, Any]:
        payload = asdict(self)
        if not include_digest:
            payload.pop("receipt_sha256", None)
        return payload


@dataclass(frozen=True)
class ObservedProcessResult:
    """Transient controller result plus its durable public observation."""

    observation: ProcessObservationV1
    stdout: str | bytes | None
    stderr: str | bytes | None


@dataclass(frozen=True)
class _ProcessRow:
    pid: int
    parent_pid: int
    process_group_id: int
    rss_kib: int
    state: str


def _process_table() -> dict[int, _ProcessRow] | None:
    """Return one portable ``ps`` snapshot, or ``None`` if unavailable."""

    executable = shutil.which("ps")
    if not executable:
        return None
    try:
        completed = subprocess.run(
            [
                executable,
                "-axo",
                "pid=,ppid=,pgid=,rss=,state=",
            ],
            capture_output=True,
            text=True,
            timeout=2,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if completed.returncode != 0:
        return None
    rows: dict[int, _ProcessRow] = {}
    for raw in completed.stdout.splitlines():
        fields = raw.split(None, 4)
        if len(fields) != 5:
            continue
        try:
            pid, parent_pid, process_group_id, rss_kib = (
                int(fields[0]),
                int(fields[1]),
                int(fields[2]),
                int(fields[3]),
            )
        except ValueError:
            continue
        rows[pid] = _ProcessRow(
            pid=pid,
            parent_pid=parent_pid,
            process_group_id=process_group_id,
            rss_kib=max(0, rss_kib),
            state=fields[4].strip(),
        )
    return rows


def _descendant_ids(
    table: dict[int, _ProcessRow], root_pid: int
) -> set[int]:
    descendants = {root_pid}
    changed = True
    while changed:
        changed = False
        for row in table.values():
            if row.parent_pid in descendants and row.pid not in descendants:
                descendants.add(row.pid)
                changed = True
    return descendants


def _active_process_ids(
    table: dict[int, _ProcessRow] | None, process_ids: set[int]
) -> set[int] | None:
    if table is None:
        return None
    return {
        pid
        for pid in process_ids
        if pid in table and not table[pid].state.upper().startswith("Z")
    }


def _send_signal_to_tree(
    *,
    process: subprocess.Popen,
    known_process_ids: set[int],
    known_process_group_ids: set[int],
    signum: signal.Signals,
) -> str:
    """Signal owned groups first, then remaining known process IDs."""

    controller_group = os.getpgrp() if hasattr(os, "getpgrp") else None
    signalled_group = False
    if hasattr(os, "killpg"):
        for process_group_id in sorted(known_process_group_ids, reverse=True):
            if process_group_id <= 0 or process_group_id == controller_group:
                continue
            try:
                os.killpg(process_group_id, signum)
                signalled_group = True
            except ProcessLookupError:
                continue
            except PermissionError:
                continue
    for pid in sorted(known_process_ids, reverse=True):
        if pid <= 0 or pid == os.getpid():
            continue
        try:
            os.kill(pid, signum)
        except (ProcessLookupError, PermissionError):
            continue
    if process.poll() is None:
        try:
            if signum == signal.SIGTERM:
                process.terminate()
            else:
                process.kill()
        except (OSError, ProcessLookupError):
            pass
    return "process_group_tree" if signalled_group else "pid_tree"


def _terminate_process_tree(
    process: subprocess.Popen,
    *,
    known_process_ids: set[int],
    known_process_group_ids: set[int],
    grace_seconds: float,
) -> tuple[bool, str, tuple[int, ...], str]:
    scope = _send_signal_to_tree(
        process=process,
        known_process_ids=known_process_ids,
        known_process_group_ids=known_process_group_ids,
        signum=signal.SIGTERM,
    )
    deadline = time.monotonic() + max(0.0, grace_seconds)
    remaining: set[int] | None = set(known_process_ids)
    while time.monotonic() < deadline:
        table = _process_table()
        active = _active_process_ids(table, known_process_ids)
        if process.poll() is not None and active == set():
            return True, scope, (), "SIGTERM"
        if active is not None:
            remaining = active
        time.sleep(0.02)

    table = _process_table()
    if table is not None:
        known_process_ids.update(_descendant_ids(table, process.pid))
        known_process_group_ids.update(
            table[pid].process_group_id
            for pid in known_process_ids
            if pid in table
        )
    scope = _send_signal_to_tree(
        process=process,
        known_process_ids=known_process_ids,
        known_process_group_ids=known_process_group_ids,
        signum=signal.SIGKILL,
    )
    kill_deadline = time.monotonic() + max(0.2, grace_seconds)
    while time.monotonic() < kill_deadline:
        table = _process_table()
        active = _active_process_ids(table, known_process_ids)
        if process.poll() is not None and active == set():
            return True, scope, (), "SIGKILL"
        if active is not None:
            remaining = active
        time.sleep(0.02)
    active = _active_process_ids(_process_table(), known_process_ids)
    if active is not None:
        remaining = active
    confirmed = process.poll() is not None and remaining == set()
    return confirmed, scope, tuple(sorted(remaining or ())), "SIGKILL"


def _build_observation(**values: Any) -> ProcessObservationV1:
    body = {
        "schema_version": PROCESS_OBSERVATION_SCHEMA_VERSION,
        **values,
    }
    return ProcessObservationV1(
        **body,
        receipt_sha256=_canonical_sha256(body),
    )


def launch_failure_observation(
    *,
    timeout_seconds: float,
    memory_limit_mb: float | None,
    error_type: str,
    wall_seconds: float = 0.0,
) -> ProcessObservationV1:
    """Represent a controller failure before a child PID was returned."""

    finding = "process.launch_failed"
    if error_type:
        finding += "." + "".join(
            character.lower()
            if character.isalnum()
            else "-"
            for character in error_type
        ).strip("-")
    return _build_observation(
        state="launch_failed",
        pid=None,
        process_group_id=None,
        process_group_owned=False,
        timeout_seconds=float(timeout_seconds),
        memory_limit_mb=(
            None if memory_limit_mb is None else float(memory_limit_mb)
        ),
        returncode=None,
        timed_out=False,
        memory_limit_exceeded=False,
        termination_requested=False,
        termination_signal="",
        termination_confirmed=None,
        termination_scope="not_started",
        remaining_process_ids=(),
        memory_observation_state="unavailable",
        last_rss_mb=None,
        peak_rss_mb=None,
        memory_sample_count=0,
        wall_seconds=float(wall_seconds),
        findings=(finding,),
    )


def observe_process(
    process: subprocess.Popen,
    *,
    timeout_seconds: float,
    memory_limit_mb: float | None = None,
    sample_interval_seconds: float = 0.1,
    termination_grace_seconds: float = 1.0,
) -> ObservedProcessResult:
    """Wait for one process while enforcing time and observed-RSS bounds.

    Memory is sampled for the root and all descendants visible in a portable
    ``ps`` snapshot.  Crossing the configured memory boundary is a failed
    observation and terminates the process tree; the limit is never increased.
    If RSS cannot be observed, the receipt says so rather than fabricating a
    value.
    """

    timeout_seconds = float(timeout_seconds)
    if not math.isfinite(timeout_seconds) or timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be finite and positive")
    if memory_limit_mb is not None:
        memory_limit_mb = float(memory_limit_mb)
        if not math.isfinite(memory_limit_mb) or memory_limit_mb <= 0:
            raise ValueError("memory_limit_mb must be finite and positive")
    sample_interval_seconds = float(sample_interval_seconds)
    if (
        not math.isfinite(sample_interval_seconds)
        or sample_interval_seconds <= 0
    ):
        raise ValueError("sample interval must be finite and positive")

    started = time.monotonic()
    deadline = started + timeout_seconds
    known_process_ids = {int(process.pid)}
    known_process_group_ids: set[int] = set()
    try:
        process_group_id = os.getpgid(process.pid)
    except (AttributeError, OSError, ProcessLookupError):
        process_group_id = None
    if process_group_id is not None:
        known_process_group_ids.add(process_group_id)
    controller_group = os.getpgrp() if hasattr(os, "getpgrp") else None
    process_group_owned = bool(
        process_group_id is not None and process_group_id != controller_group
    )

    memory_sample_count = 0
    last_rss_mb: float | None = None
    peak_rss_mb: float | None = None
    limit_reason = ""
    stdout: str | bytes | None = None
    stderr: str | bytes | None = None

    while True:
        table = _process_table()
        if table is not None:
            descendants = _descendant_ids(table, process.pid)
            known_process_ids.update(descendants)
            known_process_group_ids.update(
                table[pid].process_group_id
                for pid in descendants
                if pid in table
            )
            active = _active_process_ids(table, known_process_ids) or set()
            rss_kib = sum(table[pid].rss_kib for pid in active if pid in table)
            if active:
                last_rss_mb = rss_kib / 1024.0
                peak_rss_mb = max(peak_rss_mb or 0.0, last_rss_mb)
                memory_sample_count += 1
                if (
                    memory_limit_mb is not None
                    and peak_rss_mb > memory_limit_mb
                ):
                    limit_reason = "memory_limit_exceeded"
                    break

        now = time.monotonic()
        remaining = deadline - now
        if remaining <= 0:
            limit_reason = "timed_out"
            break
        try:
            stdout, stderr = process.communicate(
                timeout=min(sample_interval_seconds, remaining)
            )
            break
        except subprocess.TimeoutExpired:
            continue

    termination_requested = bool(limit_reason)
    termination_confirmed: bool | None = None
    termination_scope = "not_required"
    termination_signal = ""
    remaining_process_ids: tuple[int, ...] = ()
    if termination_requested:
        (
            termination_confirmed,
            termination_scope,
            remaining_process_ids,
            termination_signal,
        ) = _terminate_process_tree(
            process,
            known_process_ids=known_process_ids,
            known_process_group_ids=known_process_group_ids,
            grace_seconds=termination_grace_seconds,
        )
        try:
            stdout, stderr = process.communicate(timeout=1)
        except subprocess.TimeoutExpired:
            termination_confirmed = False

    returncode = process.poll()
    wall_seconds = max(0.0, time.monotonic() - started)
    findings: list[str] = []
    if limit_reason == "timed_out":
        findings.append("process.timeout")
    elif limit_reason == "memory_limit_exceeded":
        findings.append("process.memory_limit_exceeded")
    if termination_requested and not termination_confirmed:
        findings.append("process.termination_ambiguous")
    if not limit_reason:
        state = "exited"
    else:
        state = limit_reason + (
            "_terminated" if termination_confirmed else "_ambiguous"
        )
    observation = _build_observation(
        state=state,
        pid=int(process.pid),
        process_group_id=process_group_id,
        process_group_owned=process_group_owned,
        timeout_seconds=timeout_seconds,
        memory_limit_mb=memory_limit_mb,
        returncode=returncode,
        timed_out=limit_reason == "timed_out",
        memory_limit_exceeded=limit_reason == "memory_limit_exceeded",
        termination_requested=termination_requested,
        termination_signal=termination_signal,
        termination_confirmed=termination_confirmed,
        termination_scope=termination_scope,
        remaining_process_ids=remaining_process_ids,
        memory_observation_state=(
            "observed" if memory_sample_count else "unavailable"
        ),
        last_rss_mb=last_rss_mb,
        peak_rss_mb=peak_rss_mb,
        memory_sample_count=memory_sample_count,
        wall_seconds=wall_seconds,
        findings=tuple(sorted(set(findings))),
    )
    return ObservedProcessResult(
        observation=observation,
        stdout=stdout,
        stderr=stderr,
    )


__all__ = [
    "ObservedProcessResult",
    "PROCESS_OBSERVATION_SCHEMA_VERSION",
    "ProcessObservationV1",
    "launch_failure_observation",
    "observe_process",
]
