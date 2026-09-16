"""ORCA's input check runs at preflight, bounded, and is never charged.

Two REACH-1 cycles died at ORCA's input check under green previews
(2026-09-06): RIJK beside an analytical Hessian, and a bare RI with a
hybrid. A preview is ChemSmart's compile; the check is ORCA's, stated
in the first half-second of a run. The owner ruled (R2) for a bounded
probe launch: the real program on the retained input, stopped at the
passed banner or the cap, its word an observation on the review, and
never an engine call because engine calls derive from execution
receipts alone.
"""

from __future__ import annotations

import json
import stat
from pathlib import Path
from types import SimpleNamespace

import pytest

from chemsmart.agent._contracts import file_sha256
from chemsmart.agent.input_check import (
    InputCheckProbeReceiptV1,
    probe_observation_lines,
    probe_orca_input_check,
)
from chemsmart.agent.runtime.events import EventKind
from chemsmart.agent.terminal_states import derive_run_outcome, read_run_events

from .test_a_guide_opens_when_something_asks import _host
from .test_a_terminal_state_is_derived_not_grepped import _reserve

_ABORT_TAIL = """\
WARNING: Analytical Hessian not available with RIJK approximation
       : If you want approximation for Coulomb AND Exchange please choose RIJCOSX!
  ===> : Skipping actual calculation

[file orca_main/main_input_check.cpp, line 8637]: Error (ORCA_MAIN): ... aborting the run
"""


def _fake_orca(tmp_path, body: str) -> Path:
    script = tmp_path / "orca"
    script.write_text("#!/bin/bash\n" + body)
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    return script


def _input(tmp_path) -> Path:
    path = tmp_path / "probe.inp"
    path.write_text(
        "! b3lyp def2-svp opt\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n"
    )
    return path


@pytest.mark.capability("rule:compile.the_probe_is_the_programs_check")
def test_a_passed_check_is_stopped_at_the_banner(tmp_path):
    orca = _fake_orca(
        tmp_path,
        "echo '                    WARNINGS'\n"
        "echo '                   INPUT FILE'\n"
        "sleep 30\n",
    )
    receipt = probe_orca_input_check(
        node_id="opt",
        input_path=_input(tmp_path),
        executable=orca,
        cap_seconds=10.0,
    )
    assert isinstance(receipt, InputCheckProbeReceiptV1)
    assert receipt.status == "passed"
    assert receipt.wall_seconds < 10.0
    assert receipt.engine_lines == ()
    assert receipt.public_record()["charged"] is False
    assert "not charged" in probe_observation_lines(receipt)[0]


@pytest.mark.capability("rule:compile.the_probe_is_the_programs_check")
def test_an_aborted_check_carries_orcas_own_lines(tmp_path):
    orca = _fake_orca(
        tmp_path,
        "echo '                    WARNINGS'\n"
        f"cat <<'EOF'\n{_ABORT_TAIL}EOF\nexit 1\n",
    )
    receipt = probe_orca_input_check(
        node_id="opt",
        input_path=_input(tmp_path),
        executable=orca,
        cap_seconds=10.0,
    )
    assert receipt.status == "aborted"
    assert "aborted the run at its input check" in receipt.reason
    assert any("RIJK" in line for line in receipt.engine_lines)
    lines = probe_observation_lines(receipt)
    assert lines[0].startswith("input-check probe: aborted")
    assert any("RIJK" in line for line in lines[1:])


@pytest.mark.capability("rule:compile.the_probe_is_the_programs_check")
def test_a_silent_program_hits_the_cap_and_is_not_run(tmp_path):
    orca = _fake_orca(tmp_path, "sleep 30\n")
    receipt = probe_orca_input_check(
        node_id="opt",
        input_path=_input(tmp_path),
        executable=orca,
        cap_seconds=0.5,
    )
    assert receipt.status == "not_run"
    assert "did not conclude within 0.5 s" in receipt.reason
    assert receipt.wall_seconds < 5.0


def _previewed(tmp_path, retention: Path):
    source = _input(tmp_path)
    digest = file_sha256(source)
    retention.mkdir()
    (retention / digest).write_bytes(source.read_bytes())
    return SimpleNamespace(
        status="previewed",
        artifacts=(
            SimpleNamespace(relative_path="job/probe.inp", sha256=digest),
        ),
    )


@pytest.mark.capability("rule:compile.the_probe_is_the_programs_check")
def test_the_host_probes_the_retained_input_and_records_its_word(
    tmp_path, monkeypatch
):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    monkeypatch.delenv("PBS_JOBID", raising=False)
    retention = tmp_path / "previews"
    preview = _previewed(tmp_path, retention)
    orca = _fake_orca(tmp_path, "echo '   INPUT FILE'\nsleep 30\n")
    host = _host(
        tmp_path,
        preview_retention_root=retention,
        input_check_executable=orca,
        input_check_cap_seconds=10.0,
    )
    receipt = host._probe_input_check(
        "t1", node_id="opt", program="orca", safe_preview=preview
    )
    assert receipt.status == "passed"
    assert host._input_check_by_node["opt"] is receipt
    events = [
        json.loads(line)
        for line in (tmp_path / "events.jsonl").read_text().splitlines()
    ]
    (probed,) = [
        e for e in events if e["kind"] == EventKind.INPUT_CHECK_PROBED.value
    ]
    assert probed["payload"]["status"] == "passed"
    assert probed["payload"]["charged"] is False
    assert (
        host._probe_input_check(
            "t2", node_id="sp", program="xtb", safe_preview=preview
        )
        is None
    )


@pytest.mark.capability("rule:compile.the_probe_is_the_programs_check")
def test_inside_an_allocation_the_probe_is_not_run(tmp_path, monkeypatch):
    monkeypatch.setenv("SLURM_JOB_ID", "12345")
    retention = tmp_path / "previews"
    preview = _previewed(tmp_path, retention)
    orca = _fake_orca(tmp_path, "echo '   INPUT FILE'\n")
    host = _host(
        tmp_path,
        preview_retention_root=retention,
        input_check_executable=orca,
    )
    receipt = host._probe_input_check(
        "t1", node_id="opt", program="orca", safe_preview=preview
    )
    assert receipt.status == "not_run"
    assert "scheduler allocation" in receipt.reason


@pytest.mark.capability("rule:compile.the_probe_is_the_programs_check")
def test_a_probe_event_is_never_an_engine_call(tmp_path):
    """What the probe costs is the difference it makes, not a zero.

    This asserted an absolute zero over a stream that also carried the
    launch reservation its scaffolding needed, so the number it read was
    the scaffolding's and the probe could have cost anything. It is now
    the same stream twice, with the probe event and without it.
    """

    from chemsmart.agent.runtime.event_store import RuntimeEventStore

    def _engine_calls(directory, *, probe: bool) -> int:
        directory.mkdir(parents=True, exist_ok=True)
        store = RuntimeEventStore(directory / "events.jsonl", session_id="s")
        _reserve(store, directory)
        if probe:
            store.append(
                turn_id="t1",
                kind=EventKind.INPUT_CHECK_PROBED.value,
                payload={
                    "receipt_sha256": "f" * 64,
                    "node_id": "opt",
                    "program": "orca",
                    "status": "aborted",
                    "reason": "ORCA aborted the run at its input check",
                    "input_sha256": "a" * 64,
                    "engine_lines": ("... aborting the run",),
                    "wall_seconds": 0.4,
                    "cap_seconds": 20.0,
                    "charged": False,
                },
                idempotency_key="input_check_probed:" + "f" * 64,
            )
        return derive_run_outcome(
            read_run_events(directory / "events.jsonl")
        ).engine_calls_consumed

    with_probe = _engine_calls(tmp_path / "probed", probe=True)
    without_probe = _engine_calls(tmp_path / "bare", probe=False)
    assert with_probe == without_probe
