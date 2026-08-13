"""Deterministic, chemistry-free checks for bounded process execution."""

import os
from pathlib import Path
import signal
import subprocess
import sys
import textwrap
import threading
import time

import pytest

from chemsmart.utils.process_observation import (
    launch_failure_observation,
    observe_process,
)


def _python_process(source):
    return subprocess.Popen(
        [sys.executable, "-c", source],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )


def _write_script(path: Path, source: str) -> Path:
    path.write_text(textwrap.dedent(source), encoding="utf-8")
    return path


def _wait_for_file(path: Path, timeout_seconds: float = 5.0) -> None:
    deadline = time.monotonic() + timeout_seconds
    while time.monotonic() < deadline:
        if path.is_file():
            return
        time.sleep(0.02)
    raise AssertionError(f"timed out waiting for {path}")


def _pid_exists(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    return True


def _wait_for_pid_exit(pid: int, timeout_seconds: float = 5.0) -> None:
    deadline = time.monotonic() + timeout_seconds
    while _pid_exists(pid) and time.monotonic() < deadline:
        time.sleep(0.02)
    assert not _pid_exists(pid), f"process {pid} remains present"


def test_process_observation_records_normal_exit_and_memory():
    process = _python_process(
        "import time; value = bytearray(2 * 1024 * 1024); "
        "print(len(value), flush=True); time.sleep(0.15)"
    )

    result = observe_process(
        process,
        timeout_seconds=2,
        memory_limit_mb=256,
        sample_interval_seconds=0.02,
    )

    assert result.stdout.strip() == "2097152"
    assert result.observation.state == "exited"
    assert result.observation.returncode == 0
    assert result.observation.process_group_owned is True
    assert result.observation.timed_out is False
    assert result.observation.memory_limit_exceeded is False
    assert result.observation.memory_observation_state == "observed"
    assert result.observation.peak_rss_mb > 0
    assert result.observation.termination_requested is False
    assert result.observation.findings == ()


def test_timeout_terminates_the_owned_process_group():
    process = _python_process(
        "import subprocess, sys, time; "
        "child = subprocess.Popen([sys.executable, '-c', "
        "'import time; time.sleep(30)']); "
        "print(child.pid, flush=True); time.sleep(30)"
    )

    result = observe_process(
        process,
        timeout_seconds=0.2,
        memory_limit_mb=256,
        sample_interval_seconds=0.02,
        termination_grace_seconds=0.2,
    )

    child_pid = int(result.stdout.strip())
    assert result.observation.state == "timed_out_terminated"
    assert result.observation.timed_out is True
    assert result.observation.termination_requested is True
    assert result.observation.termination_confirmed is True
    assert result.observation.termination_scope == "process_group_tree"
    assert result.observation.remaining_process_ids == ()
    assert "process.timeout" in result.observation.findings
    try:
        os.kill(child_pid, 0)
    except ProcessLookupError:
        pass
    else:  # pragma: no cover - assertion includes the observed child PID
        raise AssertionError("timed-out descendant remains alive")


def test_observed_memory_limit_is_not_silently_relaxed():
    process = _python_process(
        "import time; value = bytearray(8 * 1024 * 1024); time.sleep(30)"
    )

    result = observe_process(
        process,
        timeout_seconds=2,
        memory_limit_mb=1,
        sample_interval_seconds=0.02,
        termination_grace_seconds=0.2,
    )

    assert result.observation.state == "memory_limit_exceeded_terminated"
    assert result.observation.memory_limit_mb == 1
    assert result.observation.memory_limit_exceeded is True
    assert result.observation.peak_rss_mb > 1
    assert result.observation.termination_confirmed is True
    assert "process.memory_limit_exceeded" in result.observation.findings


def test_launch_failure_is_explicit_and_digest_bound():
    observation = launch_failure_observation(
        timeout_seconds=600,
        memory_limit_mb=4096,
        error_type="FileNotFoundError",
    )

    assert observation.state == "launch_failed"
    assert observation.pid is None
    assert observation.timeout_seconds == 600
    assert observation.memory_limit_mb == 4096
    assert observation.termination_confirmed is None
    assert observation.findings == ("process.launch_failed.filenotfounderror",)
    assert len(observation.receipt_sha256) == 64


def test_external_sigterm_reaps_nested_engine_session_before_controller_exit(
    tmp_path,
):
    """Model benchmark controller -> agent -> engine session -> child."""

    engine_ready = tmp_path / "engine-ready"
    agent_ready = tmp_path / "agent-ready"
    engine = _write_script(
        tmp_path / "engine.py",
        """
        import os
        from pathlib import Path
        import signal
        import subprocess
        import sys
        import time

        ready = Path(sys.argv[1])
        received = None

        def capture(signum, _frame):
            global received
            if received is None:
                received = signum

        signal.signal(signal.SIGTERM, capture)
        child = subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(60)"]
        )
        ready.write_text(
            f"{os.getpid()} {child.pid}\\n", encoding="utf-8"
        )
        while received is None:
            time.sleep(0.02)
        try:
            child.wait(timeout=2)
        except subprocess.TimeoutExpired:
            child.kill()
            child.wait()
        signal.signal(received, signal.SIG_DFL)
        os.kill(os.getpid(), received)
        """,
    )
    agent = _write_script(
        tmp_path / "agent.py",
        """
        from pathlib import Path
        import subprocess
        import sys
        import time

        from chemsmart.utils.process_observation import (
            ProcessSignalGuard,
            observe_process,
        )

        with ProcessSignalGuard() as signal_guard:
            process = subprocess.Popen(
                [sys.executable, sys.argv[1], sys.argv[2]],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                start_new_session=True,
            )
            Path(sys.argv[3]).write_text("launched\\n", encoding="utf-8")
            while signal_guard.external_signal is None:
                time.sleep(0.02)
            observe_process(
                process,
                timeout_seconds=60,
                memory_limit_mb=256,
                sample_interval_seconds=0.02,
                termination_grace_seconds=1,
                signal_guard=signal_guard,
            )
        """,
    )
    controller = _write_script(
        tmp_path / "controller.py",
        """
        import os
        from pathlib import Path
        import signal
        import subprocess
        import sys
        import time

        process = subprocess.Popen(
            [
                sys.executable,
                sys.argv[1],
                sys.argv[2],
                sys.argv[3],
                sys.argv[4],
            ],
            start_new_session=True,
        )
        Path(sys.argv[5]).write_text(f"{process.pid}\\n", encoding="utf-8")
        external_signal = None

        def forward(signum, _frame):
            global external_signal
            if external_signal is None:
                external_signal = signum
            try:
                os.killpg(process.pid, signum)
            except ProcessLookupError:
                pass

        prior = signal.getsignal(signal.SIGTERM)
        signal.signal(signal.SIGTERM, forward)
        try:
            while process.poll() is None:
                time.sleep(0.02)
            process.wait()
        finally:
            signal.signal(signal.SIGTERM, prior)
        if external_signal is not None:
            os.kill(os.getpid(), external_signal)
        """,
    )
    agent_pid_file = tmp_path / "agent-pid"
    env = os.environ.copy()
    repo_root = str(Path(__file__).parents[1])
    env["PYTHONPATH"] = os.pathsep.join(
        value for value in (repo_root, env.get("PYTHONPATH", "")) if value
    )
    process = subprocess.Popen(
        [
            sys.executable,
            str(controller),
            str(agent),
            str(engine),
            str(engine_ready),
            str(agent_ready),
            str(agent_pid_file),
        ],
        env=env,
    )
    engine_pid = None
    grandchild_pid = None
    agent_pid = None
    try:
        _wait_for_file(engine_ready)
        _wait_for_file(agent_ready)
        _wait_for_file(agent_pid_file)
        engine_pid, grandchild_pid = (
            int(value)
            for value in engine_ready.read_text(encoding="utf-8").split()
        )
        agent_pid = int(agent_pid_file.read_text(encoding="utf-8"))
        time.sleep(0.1)

        process.send_signal(signal.SIGTERM)

        assert process.wait(timeout=8) == -signal.SIGTERM
        _wait_for_pid_exit(engine_pid)
        _wait_for_pid_exit(grandchild_pid)
        _wait_for_pid_exit(agent_pid)
    finally:
        if process.poll() is None:
            process.kill()
            process.wait()
        for process_group in (engine_pid, agent_pid):
            if process_group is None:
                continue
            try:
                os.killpg(process_group, signal.SIGKILL)
            except ProcessLookupError:
                pass


def test_external_signal_restores_and_invokes_custom_prior_handler():
    seen = []
    previous = signal.getsignal(signal.SIGTERM)

    def custom_handler(signum, _frame):
        seen.append(signum)

    signal.signal(signal.SIGTERM, custom_handler)
    process = _python_process("import time; time.sleep(60)")
    sender = threading.Thread(
        target=lambda: (time.sleep(0.1), os.kill(os.getpid(), signal.SIGTERM))
    )
    try:
        sender.start()
        with pytest.raises(
            InterruptedError,
            match="process observation interrupted by external signal",
        ):
            observe_process(
                process,
                timeout_seconds=60,
                memory_limit_mb=256,
                sample_interval_seconds=0.02,
                termination_grace_seconds=0.2,
            )
        sender.join(timeout=2)

        assert seen == [signal.SIGTERM]
        assert signal.getsignal(signal.SIGTERM) is custom_handler
        assert process.poll() is not None
    finally:
        signal.signal(signal.SIGTERM, previous)
        if process.poll() is None:
            process.kill()
            process.wait()
