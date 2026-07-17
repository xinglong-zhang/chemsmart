from __future__ import annotations

import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import pytest

from chemsmart.agent.runtime.calculations import (
    CalculationContext,
    CalculationEvent,
    CalculationRun,
    CalculationStatus,
    CalculationStore,
    cancel_calculation,
    execute_observed_process,
    inspect_calculation,
    inspect_output,
    load_calculation_runs,
)

COMMAND = "chemsmart run orca -p demo -f h2o.xyz -c 0 -m 1 sp"
REPO_ROOT = Path(__file__).parents[3]


def _workspace(tmp_path: Path) -> None:
    (tmp_path / "h2o.xyz").write_text(
        "3\nwater\nO 0 0 0\nH 0.75 0 0.5\nH -0.75 0 0.5\n",
        encoding="utf-8",
    )
    project = tmp_path / ".chemsmart" / "orca" / "demo.yaml"
    project.parent.mkdir(parents=True)
    project.write_text(
        "gas:\n  functional: pbe0\n  basis: def2-svp\n",
        encoding="utf-8",
    )


def _python_writer(source: str) -> list[str]:
    return [sys.executable, "-c", source]


def test_orca_process_lifecycle_persists_chemistry_receipt(
    tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    _workspace(tmp_path)
    session_dir = tmp_path / "session"
    events: list[CalculationEvent] = []
    source = """
from pathlib import Path
import time
path = Path('h2o.out')
path.write_text('SCF CONVERGED AFTER 8 CYCLES\\n')
time.sleep(0.25)
with path.open('a') as handle:
    handle.write('FINAL SINGLE POINT ENERGY     -76.269459371830\\n')
    handle.write('ORCA TERMINATED NORMALLY\\n')
"""

    result = execute_observed_process(
        _python_writer(source),
        command=COMMAND,
        timeout_s=5,
        context=CalculationContext(
            session_dir=session_dir,
            session_id="session-1",
            turn_id="turn-1",
            semantic_verdict="ok",
            intent_verdict="ok",
        ),
        event_sink=events.append,
    )

    run = result["calculation"]
    assert run["status"] == CalculationStatus.COMPLETED.value
    assert run["energy"] == pytest.approx(-76.269459371830)
    assert run["scf_cycles"] == 8
    assert run["normal_termination"] is True
    assert run["method"] == "pbe0"
    assert run["basis"] == "def2-svp"
    assert run["turn_id"] == "turn-1"
    assert [event.kind for event in events[:3]] == [
        "validating",
        "starting",
        "started",
    ]
    assert events[-1].kind == "completed"
    receipt = session_dir / "calculations" / run["run_id"] / "receipt.json"
    assert receipt.is_file()
    assert load_calculation_runs(tmp_path)[-1].run_id == run["run_id"]

    inspected = inspect_calculation(run["run_id"], session_root=str(tmp_path))
    assert inspected["ok"] is True
    assert inspected["calculation"]["energy"] == pytest.approx(
        -76.269459371830
    )


def test_exit_zero_without_output_is_chemistry_failure(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _workspace(tmp_path)

    result = execute_observed_process(
        _python_writer("print('process completed without chemistry output')"),
        command=COMMAND,
        timeout_s=5,
        cwd=tmp_path,
    )

    run = result["calculation"]
    assert run["returncode"] == 0
    assert run["status"] == CalculationStatus.CHEMISTRY_FAILED.value
    assert "No Gaussian/ORCA output" in run["stage"]


def test_fake_run_can_complete_without_chemistry_output(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _workspace(tmp_path)

    result = execute_observed_process(
        _python_writer("print('generated input')") + ["--fake"],
        command=COMMAND,
        timeout_s=5,
        cwd=tmp_path,
    )

    run = result["calculation"]
    assert run["status"] == CalculationStatus.COMPLETED.value
    assert run["execution_mode"] == "test_fake"


def test_existing_completed_output_is_reused_without_false_failure(
    tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    _workspace(tmp_path)
    output = tmp_path / "h2o_sp.out"
    output.write_text(
        "SCF CONVERGED AFTER 8 CYCLES\n"
        "FINAL SINGLE POINT ENERGY -76.269459371830\n"
        "ORCA TERMINATED NORMALLY\n",
        encoding="utf-8",
    )
    source = (
        "print('ORCASinglePointJob<folder="
        + str(tmp_path)
        + ", label=h2o_sp, jobrunner=ORCAJobRunner> is already complete, not running.')"
    )

    run = execute_observed_process(
        _python_writer(source),
        command=COMMAND,
        timeout_s=5,
        cwd=tmp_path,
    )["calculation"]

    assert run["status"] == CalculationStatus.COMPLETED.value
    assert run["reused_output"] is True
    assert run["output_path"] == str(output.resolve())
    assert run["normal_termination"] is True
    assert run["energy"] == pytest.approx(-76.269459371830)


def test_process_failure_and_timeout_are_distinct(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _workspace(tmp_path)

    failed = execute_observed_process(
        _python_writer(
            "import sys; print('ORCA error', file=sys.stderr); sys.exit(3)"
        ),
        command=COMMAND,
        timeout_s=5,
        cwd=tmp_path,
    )["calculation"]
    timed_out = execute_observed_process(
        _python_writer("import time; time.sleep(5)"),
        command=COMMAND,
        timeout_s=0,
        cwd=tmp_path,
    )["calculation"]

    assert failed["status"] == CalculationStatus.PROCESS_FAILED.value
    assert failed["returncode"] == 3
    assert "ORCA error" in failed["error"]
    assert timed_out["status"] == CalculationStatus.TIMEOUT.value


def test_started_event_can_cancel_process(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _workspace(tmp_path)

    def cancel_on_start(event: CalculationEvent) -> None:
        if event.kind == "started":
            assert cancel_calculation(event.run.run_id) is True

    run = execute_observed_process(
        _python_writer("import time; time.sleep(30)"),
        command=COMMAND,
        timeout_s=5,
        cwd=tmp_path,
        event_sink=cancel_on_start,
    )["calculation"]

    assert run["status"] == CalculationStatus.CANCELLED.value


def test_large_stdout_is_streamed_and_tool_tail_is_bounded(
    tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    _workspace(tmp_path)
    session_dir = tmp_path / "session"
    result = execute_observed_process(
        _python_writer("import sys; sys.stdout.write('x' * 10_000_000)")
        + ["--fake"],
        command=COMMAND,
        timeout_s=10,
        context=CalculationContext(session_dir=session_dir),
        cwd=tmp_path,
    )

    run = result["calculation"]
    assert Path(run["stdout_path"]).stat().st_size == 10_000_000
    assert len(result["stdout_tail"]) == 12_000
    assert run["status"] == CalculationStatus.COMPLETED.value


def test_inspect_orca_output_fixture(tmp_path):
    output = tmp_path / "water.out"
    output.write_text(
        "SCF CONVERGED AFTER 8 CYCLES\n"
        "FINAL SINGLE POINT ENERGY     -76.269459371830\n"
        "TOTAL RUN TIME: 0 days 0 hours 0 minutes 0 seconds 563 msec\n"
        "ORCA TERMINATED NORMALLY\n",
        encoding="utf-8",
    )

    summary = inspect_output(output, program="orca", kind="sp")

    assert summary["normal_termination"] is True
    assert summary["scf_cycles"] == 8
    assert summary["energy"] == pytest.approx(-76.269459371830)
    assert summary["chemistry_elapsed_s"] == pytest.approx(0.563)


def test_job_specific_orca_receipts_use_real_outputs():
    output_root = REPO_ROOT / "tests" / "data" / "ORCATests" / "outputs"

    optimization = inspect_output(
        output_root / "water_opt.out", program="orca", kind="opt"
    )
    neb = inspect_output(
        output_root / "neb_R-TS1-Si.out", program="orca", kind="neb"
    )
    qmmm = inspect_output(
        output_root / "methanol_ethane_qmmm.out",
        program="orca",
        kind="opt.qmmm",
    )

    assert optimization["optimization_cycles"] == 5
    assert optimization["optimization_converged"] is True
    assert optimization["imaginary_frequency_count"] == 0
    assert optimization["frequency_count"] == 3
    assert neb["neb_images"] == 10
    assert neb["imaginary_frequency_count"] == 1
    assert neb["imag_freqs"] == [-821.98]
    assert qmmm["qmmm_total_atoms"] == 36
    assert qmmm["qmmm_mm_atoms"] == 24
    assert qmmm["qmmm_qm_atoms"] == 12


def test_gaussian_ts_receipt_reports_imaginary_mode():
    output = (
        REPO_ROOT
        / "tests"
        / "data"
        / "GaussianTests"
        / "outputs"
        / "Pd_insertion_ts_r.log"
    )

    summary = inspect_output(output, program="gaussian", kind="ts")

    assert summary["normal_termination"] is True
    assert summary["imaginary_frequency_count"] == 1
    assert summary["imag_freqs"] == [-112.225]
    assert summary["optimization_converged"] is True


def test_restart_recovers_live_pid_and_completed_detached_run(
    tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    runtime_session = tmp_path / "sessions" / ".runtime" / "tui-one"
    store = CalculationStore(runtime_session)
    live_output = tmp_path / "live.out"
    live_output.write_text("  7  -76.2\n", encoding="utf-8")
    live = CalculationRun(
        run_id="calc-live",
        command=COMMAND,
        cwd=str(tmp_path),
        program="orca",
        kind="sp",
        label="live",
        status="running",
        stage="Process running",
        pid=os.getpid(),
        started_at=datetime.now(timezone.utc).isoformat(),
        output_path=str(live_output),
    )
    store.write_run(live)

    done_output = tmp_path / "done.out"
    done_output.write_text(
        "FINAL SINGLE POINT ENERGY -76.269459371830\n"
        "ORCA TERMINATED NORMALLY\n",
        encoding="utf-8",
    )
    done = CalculationRun(
        run_id="calc-detached",
        command=COMMAND,
        cwd=str(tmp_path),
        program="orca",
        kind="sp",
        label="done",
        status="running",
        stage="Process running",
        pid=999_999_999,
        started_at=datetime.now(timezone.utc).isoformat(),
        output_path=str(done_output),
    )
    store.write_run(done)

    recovered = {
        run.run_id: run for run in load_calculation_runs(tmp_path / "sessions")
    }

    assert recovered["calc-live"].status == "running"
    assert recovered["calc-live"].stage == "SCF cycle 7"
    assert recovered["calc-detached"].status == "completed"
    assert recovered["calc-detached"].normal_termination is True
    assert recovered["calc-detached"].energy == pytest.approx(-76.269459371830)
