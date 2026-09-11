"""ORCA's own input check, run by the host on the materialised input.

Two REACH-1 cycles died at ORCA's input check under green previews
(2026-09-06): RIJK beside an analytical Hessian, and a bare RI with a
hybrid functional. ChemSmart's preview compiles the project; only ORCA
knows ORCA's rules, and it states them in the first half-second of a
run, before any parallel process is spawned. The owner ruled (R2,
2026-09-06) for a bounded probe launch at preflight: ORCA on the
host-materialised input, capped, stopped the moment its check passes,
a typed receipt on the review, and never an engine call -- engine
calls derive from execution receipts alone, and a probe mints none.
"""

from __future__ import annotations

import os
import re
import shutil
import signal
import subprocess
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from chemsmart.agent._contracts import (
    ContractError,
    canonical_sha256,
    file_sha256,
)
from chemsmart.io.native_failure import (
    ORCA_INPUT_CHECK_ABORT,
    summarize_orca_native_failure,
)

SCHEMA_VERSION = "chemsmart.input-check-probe-receipt.v1"
STATUSES = ("passed", "aborted", "not_run")

#: ORCA echoes the input under this banner once its checks have passed;
#: an aborted check never reaches it (read from real outputs on this
#: host, 2026-09-06).
INPUT_CHECK_PASSED = re.compile(r"^\s*INPUT FILE\s*$")


@dataclass(frozen=True)
class InputCheckProbeReceiptV1:
    schema_version: str
    node_id: str
    program: str
    input_sha256: str
    executable_sha256: str
    status: str
    reason: str
    engine_lines: tuple[str, ...]
    wall_seconds: float
    cap_seconds: float
    output_sha256: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != SCHEMA_VERSION:
            raise ContractError("unsupported input-check probe receipt schema")
        if self.status not in STATUSES:
            raise ContractError("invalid input-check probe status")
        if self.receipt_sha256 != canonical_sha256(_body(self)):
            raise ContractError("input-check probe receipt digest mismatch")

    def public_record(self) -> dict[str, Any]:
        return {**_body(self), "receipt_sha256": self.receipt_sha256}


def _body(receipt: Any) -> dict[str, Any]:
    return {
        "schema_version": receipt.schema_version,
        "node_id": receipt.node_id,
        "program": receipt.program,
        "input_sha256": receipt.input_sha256,
        "executable_sha256": receipt.executable_sha256,
        "status": receipt.status,
        "reason": receipt.reason,
        "engine_lines": tuple(receipt.engine_lines),
        "wall_seconds": round(float(receipt.wall_seconds), 3),
        "cap_seconds": float(receipt.cap_seconds),
        "output_sha256": receipt.output_sha256,
        "charged": False,
    }


def build_input_check_probe_receipt(
    *,
    node_id: str,
    program: str,
    input_sha256: str,
    executable_sha256: str,
    status: str,
    reason: str,
    engine_lines: tuple[str, ...] = (),
    wall_seconds: float = 0.0,
    cap_seconds: float,
    output_sha256: str = "",
) -> InputCheckProbeReceiptV1:
    fields = {
        "schema_version": SCHEMA_VERSION,
        "node_id": str(node_id),
        "program": str(program),
        "input_sha256": str(input_sha256),
        "executable_sha256": str(executable_sha256),
        "status": str(status),
        "reason": str(reason),
        "engine_lines": tuple(str(line) for line in engine_lines),
        "wall_seconds": round(float(wall_seconds), 3),
        "cap_seconds": float(cap_seconds),
        "output_sha256": str(output_sha256)
        or canonical_sha256({"output": ""}),
    }

    class _Draft:
        pass

    draft = _Draft()
    for key, value in fields.items():
        setattr(draft, key, value)
    return InputCheckProbeReceiptV1(
        **fields, receipt_sha256=canonical_sha256(_body(draft))
    )


def not_run_receipt(
    *,
    node_id: str,
    program: str,
    input_sha256: str,
    reason: str,
    cap_seconds: float,
) -> InputCheckProbeReceiptV1:
    """A probe the host decided not to launch, with the reason on it."""

    return build_input_check_probe_receipt(
        node_id=node_id,
        program=program,
        input_sha256=input_sha256,
        executable_sha256="",
        status="not_run",
        reason=reason,
        cap_seconds=cap_seconds,
    )


def input_check_passed(text: str) -> bool:
    return any(INPUT_CHECK_PASSED.match(line) for line in text.splitlines())


def probe_orca_input_check(
    *,
    node_id: str,
    input_path: Path,
    executable: Path,
    env: Mapping[str, str] | None = None,
    cap_seconds: float = 20.0,
    poll_seconds: float = 0.1,
) -> InputCheckProbeReceiptV1:
    """Run ORCA on one materialised input until its check concludes.

    The run happens in a throwaway directory in its own session, is
    stopped with SIGTERM then SIGKILL the moment the passed banner
    appears or the cap is reached, and leaves nothing behind but the
    receipt. The words are ORCA's: an abort is summarised by the same
    native-failure reader a real run's death would be.
    """

    input_path = Path(input_path)
    executable = Path(executable)
    if float(cap_seconds) <= 0.0:
        raise ContractError("an input-check probe needs a positive cap")
    input_sha256 = file_sha256(input_path)
    executable_sha256 = file_sha256(executable)
    work = Path(tempfile.mkdtemp(prefix="chemsmart-input-check-"))
    started = time.monotonic()
    status = "not_run"
    reason = ""
    engine_lines: tuple[str, ...] = ()
    text = ""
    try:
        staged = work / "probe.inp"
        shutil.copyfile(input_path, staged)
        out_path = work / "probe.out"
        with out_path.open("w", encoding="utf-8") as out:
            process = subprocess.Popen(
                [str(executable), str(staged)],
                stdout=out,
                stderr=subprocess.STDOUT,
                cwd=str(work),
                env=dict(env) if env else None,
                start_new_session=True,
            )
            try:
                while True:
                    exited = process.poll() is not None
                    text = out_path.read_text(
                        encoding="utf-8", errors="replace"
                    )
                    if input_check_passed(text):
                        status = "passed"
                        break
                    if exited:
                        break
                    if time.monotonic() - started > float(cap_seconds):
                        reason = (
                            "the check did not conclude within "
                            f"{float(cap_seconds):g} s"
                        )
                        break
                    time.sleep(float(poll_seconds))
            finally:
                _stop(process)
        text = out_path.read_text(encoding="utf-8", errors="replace")
        lines = text.splitlines()
        if status == "passed":
            reason = "the input passed ORCA's check; the run was stopped there"
        elif not reason:
            status = "aborted"
            aborted = any(
                ORCA_INPUT_CHECK_ABORT.search(line) for line in lines
            )
            reason = (
                "ORCA aborted the run at its input check"
                if aborted
                else "ORCA exited before its input check passed"
            )
            summary = summarize_orca_native_failure(lines)
            engine_lines = tuple(getattr(summary, "engine_lines", ()) or ())
            if not engine_lines:
                engine_lines = tuple(
                    line.rstrip() for line in lines[-8:] if line.strip()
                )
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return build_input_check_probe_receipt(
        node_id=node_id,
        program="orca",
        input_sha256=input_sha256,
        executable_sha256=executable_sha256,
        status=status,
        reason=reason,
        engine_lines=engine_lines,
        wall_seconds=time.monotonic() - started,
        cap_seconds=cap_seconds,
        output_sha256=canonical_sha256({"output": text}),
    )


def _stop(process: subprocess.Popen) -> None:
    if process.poll() is not None:
        return
    for signum, grace in ((signal.SIGTERM, 1.0), (signal.SIGKILL, 5.0)):
        try:
            os.killpg(process.pid, signum)
        except ProcessLookupError:
            return
        try:
            process.wait(timeout=grace)
            return
        except subprocess.TimeoutExpired:
            continue


def probe_observation_lines(
    receipt: InputCheckProbeReceiptV1,
) -> tuple[str, ...]:
    """The probe's word as review observations, never a verdict."""

    head = (
        f"input-check probe: {receipt.status} in "
        f"{receipt.wall_seconds:.1f} s of {receipt.cap_seconds:g} s, "
        f"not charged -- {receipt.reason}"
    )
    return (head, *(f"  orca: {line}" for line in receipt.engine_lines[:6]))
