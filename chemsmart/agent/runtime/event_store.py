"""Secure, locked, append-only Runtime V2 JSONL event store."""

from __future__ import annotations

from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import stat
from typing import Any, Iterable, Iterator, Mapping, TextIO

try:  # POSIX advisory locking.
    import fcntl as _fcntl
except ImportError:  # pragma: no cover - Windows import boundary
    _fcntl = None

try:  # Windows byte-range locking.
    import msvcrt as _msvcrt
except ImportError:  # pragma: no cover - POSIX import boundary
    _msvcrt = None

from chemsmart.agent._contracts import (
    ContractError,
    canonical_data,
    canonical_sha256,
)
from chemsmart.agent.permissions import (
    ApprovalResolutionV1,
    PermissionReceiptV1,
    PermissionRequestV1,
    _evaluate_permission,
)
from chemsmart.agent.runtime.events import (
    CAPABILITY_QUERIED,
    COMMAND_COMPILED,
    COMMAND_INSPECTED,
    ENGINE_BOUND,
    ENVIRONMENT_QUERIED,
    OPTIMIZED_GEOMETRY_HANDED_OFF,
    PERMISSION_RESOLVED,
    PROGRAM_BOUND,
    PROGRAM_EXECUTED,
    PROGRAM_PREFLIGHTED,
    PROJECT_VALIDATED,
    PROJECT_PROMOTED,
    RESULT_VERIFIED,
    SCIENTIFIC_DECISION_RECORDED,
    RUNTIME_TERMINATED,
    SAFE_PREVIEWED,
    SUBSTITUTION_ASSESSED,
    VALIDATOR_OBSERVED,
    RuntimeEvent,
)
from chemsmart.agent.runtime.reducer import RuntimeState, replay_events


class RuntimeEventStore:
    """One crash-stable stream with atomic one-shot approval consumption."""

    def __init__(self, path: str | Path, *, session_id: str) -> None:
        # abspath deliberately does not follow the final symlink; secure open
        # rejects links before and during open.
        self.path = Path(os.path.abspath(os.fspath(path)))
        self.session_id = str(session_id)
        if not self.session_id:
            raise ContractError("event store session_id must not be empty")
        self._prepare_private_parent()
        self._lock_path = self.path.with_name("." + self.path.name + ".lock")

    def read_events(self) -> tuple[RuntimeEvent, ...]:
        with self._locked_handle(exclusive=False) as handle:
            return self._read_locked(handle)

    def state(self) -> RuntimeState:
        return replay_events(self.read_events())

    def append(
        self,
        *,
        turn_id: str,
        kind: str,
        payload: Mapping[str, Any] | None = None,
        idempotency_key: str = "",
    ) -> RuntimeEvent:
        normalized_kind = str(getattr(kind, "value", kind))
        if normalized_kind == RUNTIME_TERMINATED and (
            payload or {}
        ).get("terminal_state") in {"complete", "planned"}:
            raise ContractError(
                "host-derived termination must pass RuntimeEventStore.terminate"
            )
        with self._locked_handle(exclusive=True) as handle:
            return self._append_locked(
                handle,
                turn_id=turn_id,
                kind=normalized_kind,
                payload=dict(payload or {}),
                idempotency_key=idempotency_key,
            )

    def consume_approval(
        self,
        *,
        turn_id: str,
        request: PermissionRequestV1,
        approval: ApprovalResolutionV1,
    ) -> tuple[PermissionReceiptV1, RuntimeEvent]:
        """Evaluate and consume an exact one-shot approval under one lock."""

        with self._locked_handle(exclusive=True) as handle:
            events = self._read_locked(handle)
            consumed = any(
                event.kind == PERMISSION_RESOLVED
                and event.payload.get("approval_id") == approval.approval_id
                and event.payload.get("decision") == "allow_once"
                for event in events
            )
            if consumed:
                raise ContractError("approval has already been consumed")
            receipt = _evaluate_permission(request, approval=approval)
            payload = {
                "receipt_sha256": receipt.receipt_sha256,
                "permission_request_sha256": request.request_sha256,
                "approval_id": approval.approval_id,
                "approval_resolution_sha256": approval.resolution_sha256,
                "decision": receipt.decision.value,
                "reason": receipt.reason,
            }
            event = self._append_locked(
                handle,
                turn_id=turn_id,
                kind=PERMISSION_RESOLVED,
                payload=payload,
                idempotency_key=(
                    "approval:" + approval.approval_id + ":" + request.request_sha256
                ),
                existing_events=events,
            )
            return receipt, event

    def terminate(
        self,
        *,
        turn_id: str,
        terminal_state: str,
        reason: str,
        required_receipt_sha256s: Iterable[str] = (),
    ) -> RuntimeEvent:
        """Create the sole terminal event; green status is host-derived."""

        with self._locked_handle(exclusive=True) as handle:
            events = self._read_locked(handle)
            state = replay_events(events)
            if state.terminal_state:
                raise ContractError("runtime is already terminal")
            payload: dict[str, Any] = {
                "terminal_state": str(terminal_state),
                "reason": str(reason),
            }
            if terminal_state == "complete":
                required = tuple(sorted(set(required_receipt_sha256s)))
                observed = _observed_receipts(state)
                if not required or not set(required).issubset(observed):
                    raise ContractError(
                        "completion requires receipts already observed in stream"
                    )
                green = tuple(
                    digest
                    for digest in required
                    if _receipt_is_green(events, digest)
                )
                if green != required:
                    raise ContractError("a required completion gate is red")
                if not state.preflight_receipts or (
                    state.preflight_receipts[-1] not in required
                ):
                    raise ContractError(
                        "completion requires the latest node preflight"
                    )
                gate = canonical_sha256(
                    {
                        "session_id": self.session_id,
                        "turn_id": turn_id,
                        "required_receipt_sha256s": required,
                        "green_receipt_sha256s": green,
                        "previous_hash": events[-1].event_hash if events else "",
                    }
                )
                payload.update(
                    {
                        "required_receipt_sha256s": required,
                        "green_receipt_sha256s": green,
                        "completion_gate_sha256": gate,
                    }
                )
            elif terminal_state == "planned":
                required = tuple(sorted(set(required_receipt_sha256s)))
                observed = _observed_receipts(state)
                if (
                    len(required) != 1
                    or required[0] not in observed
                    or not state.workflow_receipts
                    or state.workflow_receipts[-1] != required[0]
                ):
                    raise ContractError(
                        "planned termination requires the latest workflow draft"
                    )
                payload["plan_receipt_sha256"] = required[0]
            return self._append_locked(
                handle,
                turn_id=turn_id,
                kind=RUNTIME_TERMINATED,
                payload=payload,
                idempotency_key="runtime-terminal",
                existing_events=events,
            )

    def persist_public_transcript(
        self, *, turn_id: str, transcript: Iterable[Mapping[str, Any]]
    ) -> dict[str, Any]:
        """Atomically persist the exact sanitized visible transcript."""

        normalized = tuple(canonical_data(dict(item)) for item in transcript)
        transcript_sha256 = canonical_sha256(normalized)
        body = {
            "schema_version": "chemsmart.public-transcript.v1",
            "session_id": self.session_id,
            "turn_id": str(turn_id),
            "transcript_sha256": transcript_sha256,
            "transcript": normalized,
        }
        encoded = (
            json.dumps(
                body,
                sort_keys=True,
                indent=2,
                ensure_ascii=False,
            )
            + "\n"
        ).encode("utf-8")
        artifact_sha256 = hashlib.sha256(encoded).hexdigest()
        suffix = canonical_sha256({"turn_id": str(turn_id)})[:16]
        destination = self.path.with_name(
            f"public-transcript-{suffix}.json"
        )
        temporary = destination.with_name(
            "." + destination.name + f".tmp-{os.getpid()}"
        )
        with self._locked_handle(exclusive=True):
            if os.path.lexists(destination):
                existing = _secure_open_text(destination)
                try:
                    existing.seek(0)
                    observed = existing.read().encode("utf-8")
                finally:
                    existing.close()
                if observed != encoded:
                    raise ContractError("public transcript artifact conflicts")
            else:
                flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
                if hasattr(os, "O_NOFOLLOW"):
                    flags |= os.O_NOFOLLOW
                descriptor = os.open(temporary, flags, 0o600)
                try:
                    with os.fdopen(descriptor, "wb") as handle:
                        descriptor = -1
                        handle.write(encoded)
                        handle.flush()
                        os.fsync(handle.fileno())
                    os.replace(temporary, destination)
                    os.chmod(destination, 0o600)
                finally:
                    if descriptor >= 0:
                        os.close(descriptor)
                    if temporary.exists():
                        temporary.unlink()
        return {
            "schema_version": "chemsmart.public-transcript-artifact.v1",
            "artifact_id": "public_transcript." + suffix,
            "path": str(destination),
            "artifact_sha256": artifact_sha256,
            "transcript_sha256": transcript_sha256,
        }

    def _prepare_private_parent(self) -> None:
        self.path.parent.mkdir(parents=True, mode=0o700, exist_ok=True)
        if self.path.parent.is_symlink() or not self.path.parent.is_dir():
            raise ContractError("event-store parent must be a real directory")
        try:
            os.chmod(self.path.parent, 0o700)
        except OSError as exc:
            raise ContractError("event-store parent cannot be made private") from exc

    @contextmanager
    def _locked_handle(self, *, exclusive: bool) -> Iterator[TextIO]:
        lock_handle = _secure_open_text(self._lock_path)
        try:
            _acquire_lock(lock_handle, exclusive=exclusive)
            data_handle = _secure_open_text(self.path)
            try:
                yield data_handle
            finally:
                data_handle.close()
        finally:
            _release_lock(lock_handle)
            lock_handle.close()

    def _read_locked(self, handle: TextIO) -> tuple[RuntimeEvent, ...]:
        handle.seek(0)
        events = []
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                raw = json.loads(line)
                events.append(RuntimeEvent.from_dict(raw))
            except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
                raise ContractError(
                    f"invalid runtime event at JSONL line {line_number}"
                ) from exc
        state = replay_events(events)
        if state.session_id and state.session_id != self.session_id:
            raise ContractError("event store contains another session")
        return tuple(events)

    def _append_locked(
        self,
        handle: TextIO,
        *,
        turn_id: str,
        kind: str,
        payload: dict[str, Any],
        idempotency_key: str,
        existing_events: tuple[RuntimeEvent, ...] | None = None,
    ) -> RuntimeEvent:
        events = (
            existing_events
            if existing_events is not None
            else self._read_locked(handle)
        )
        if events and replay_events(events).terminal_state:
            raise ContractError("runtime terminal state is absorbing")
        if idempotency_key:
            identity = canonical_sha256({"kind": kind, "payload": payload})
            for existing in events:
                if existing.idempotency_key != idempotency_key:
                    continue
                existing_identity = canonical_sha256(
                    {"kind": existing.kind, "payload": existing.payload}
                )
                if existing_identity != identity:
                    raise ContractError(
                        "idempotency key conflicts with persisted action"
                    )
                return existing
        event = RuntimeEvent.create(
            sequence=len(events) + 1,
            session_id=self.session_id,
            turn_id=str(turn_id),
            kind=kind,
            payload=payload,
            previous_hash=events[-1].event_hash if events else "",
            idempotency_key=idempotency_key,
        )
        handle.seek(0, os.SEEK_END)
        handle.write(
            json.dumps(
                event.to_dict(),
                sort_keys=True,
                separators=(",", ":"),
                ensure_ascii=False,
            )
            + "\n"
        )
        handle.flush()
        os.fsync(handle.fileno())
        return event


def _secure_open_text(path: Path) -> TextIO:
    if os.path.lexists(path) and path.is_symlink():
        raise ContractError("runtime event paths must not be symbolic links")
    flags = os.O_RDWR | os.O_CREAT | os.O_APPEND
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    if hasattr(os, "O_BINARY"):
        flags |= os.O_BINARY
    try:
        descriptor = os.open(path, flags, 0o600)
    except OSError as exc:
        raise ContractError("runtime event file cannot be opened securely") from exc
    try:
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            raise ContractError("runtime event path must be a regular file")
        os.fchmod(descriptor, 0o600)
        return os.fdopen(descriptor, "a+", encoding="utf-8", newline="\n")
    except Exception:
        os.close(descriptor)
        raise


def _acquire_lock(handle: TextIO, *, exclusive: bool) -> None:
    if _fcntl is not None:
        operation = _fcntl.LOCK_EX if exclusive else _fcntl.LOCK_SH
        _fcntl.flock(handle.fileno(), operation)
        return
    if _msvcrt is None:  # pragma: no cover - unsupported platform
        raise ContractError("no portable file-lock implementation is available")
    handle.seek(0, os.SEEK_END)
    if handle.tell() == 0:
        handle.write("\0")
        handle.flush()
        os.fsync(handle.fileno())
    handle.seek(0)
    operation = _msvcrt.LK_LOCK if exclusive else _msvcrt.LK_RLCK
    _msvcrt.locking(handle.fileno(), operation, 1)


def _release_lock(handle: TextIO) -> None:
    if _fcntl is not None:
        _fcntl.flock(handle.fileno(), _fcntl.LOCK_UN)
        return
    if _msvcrt is not None:  # pragma: no branch - platform boundary
        handle.seek(0)
        _msvcrt.locking(handle.fileno(), _msvcrt.LK_UNLCK, 1)


def _receipt_is_green(
    events: tuple[RuntimeEvent, ...], digest: str, _seen: frozenset[str] = frozenset()
) -> bool:
    if digest in _seen:
        return False
    event = next(
        (
            item
            for item in reversed(events)
            if str(
                item.payload.get("receipt_sha256")
                or item.payload.get("binding_sha256")
                or ""
            )
            == digest
        ),
        None,
    )
    if event is None:
        return False
    value = event.payload
    if event.kind == CAPABILITY_QUERIED:
        return value.get("status") in {"supported", "preview_only"}
    if event.kind == ENVIRONMENT_QUERIED:
        return value.get("status") == "available"
    if event.kind in {PROGRAM_BOUND, ENGINE_BOUND}:
        return value.get("state") in {"resolved", "preview_only"}
    if event.kind == PROJECT_VALIDATED:
        return value.get("status") == "valid"
    if event.kind == PROJECT_PROMOTED:
        return value.get("status") == "validated"
    if event.kind == SCIENTIFIC_DECISION_RECORDED:
        return value.get("status") == "recorded"
    if event.kind == COMMAND_COMPILED:
        return value.get("status") == "compiled"
    if event.kind == COMMAND_INSPECTED:
        return value.get("status") == "valid"
    if event.kind == SAFE_PREVIEWED:
        return (
            value.get("status") == "previewed"
            and value.get("program_validation_status") == "valid"
            and _finding_count(value) == 0
        )
    if event.kind == VALIDATOR_OBSERVED:
        return value.get("status") == "valid" and _finding_count(value) == 0
    if event.kind == PROGRAM_PREFLIGHTED:
        preview_digest = str(value.get("safe_preview_receipt_sha256") or "")
        return (
            value.get("plan_state") == "previewed"
            and _finding_count(value) == 0
            and _receipt_is_green(
                events, preview_digest, _seen | frozenset({digest})
            )
        )
    if event.kind == RESULT_VERIFIED:
        return value.get("status") == "valid"
    if event.kind == PROGRAM_EXECUTED:
        return (
            value.get("execution_state") == "engine_complete"
            and value.get("engine_complete") is True
            and value.get("validated") is True
        )
    if event.kind == OPTIMIZED_GEOMETRY_HANDED_OFF:
        return value.get("status") == "validated_handoff"
    if event.kind == SUBSTITUTION_ASSESSED:
        return value.get("decision") in {"exact", "approved"}
    if event.kind == PERMISSION_RESOLVED:
        return value.get("decision") in {"auto_allow", "allow_once"}
    return False


def _finding_count(payload: Mapping[str, Any]) -> int:
    value = payload.get("critical_finding_count")
    return value if isinstance(value, int) and not isinstance(value, bool) else -1


def _observed_receipts(state: RuntimeState) -> set[str]:
    collections = (
        state.capability_receipts,
        state.environment_receipts,
        state.program_bindings,
        state.engine_bindings,
        state.project_receipts,
        state.project_promotion_receipts,
        state.scientific_decision_receipts,
        state.workflow_receipts,
        state.command_invocations,
        state.command_inspections,
        state.safe_preview_receipts,
        state.validator_receipts,
        state.preflight_receipts,
        state.substitution_receipts,
        state.result_verification_receipts,
        state.execution_receipts,
        state.optimized_geometry_handoff_receipts,
        state.permission_receipts,
        state.provider_turn_receipts,
        state.api_attempt_receipts,
    )
    return {digest for collection in collections for digest in collection if digest}


__all__ = ["RuntimeEventStore"]
