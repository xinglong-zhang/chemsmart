"""N elements write one run directory at once, and none of them dies of it.

Measured live on CUHK (goal `butane-wave-1`, array 2135192, 2026-09-16):
a wave of four started four processes in the same second, all four found
`execution-server.yaml` absent, all four called `os.open(..., O_EXCL)`,
one won and three died with `FileExistsError`. Three approved
calculations never ran, and the wave could never satisfy its own barrier.

This is §2 row 13 of the plan -- recorded as a falsified serial
assumption and never repaired, which is worse than not having noticed.

The file is a per-*run* fact: it carries the allocation the scheduler
granted this array, which is one `#SBATCH` shape for every element. So
identical bytes from a racing writer are the expected case and must
succeed; **different** bytes remain a contract error, because then two
elements disagree about what the host granted.
"""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.live_session import _write_private_exact


def test_concurrent_writers_of_identical_bytes_all_succeed(tmp_path):
    path = tmp_path / "execution-server.yaml"
    payload = b"SERVER:\n  NUM_CORES: 16\n"
    barrier = __import__("threading").Barrier(8)

    def write():
        barrier.wait()
        _write_private_exact(path, payload)
        return "ok"

    with ThreadPoolExecutor(max_workers=8) as pool:
        results = [f.result() for f in [pool.submit(write) for _ in range(8)]]

    assert results == ["ok"] * 8, (
        "an element died writing a file another element had just written "
        "with the same bytes; three approved calculations were lost to "
        "this on CUHK"
    )
    assert path.read_bytes() == payload
    assert oct(path.stat().st_mode)[-3:] == "600"


def test_a_writer_with_different_bytes_is_still_refused(tmp_path):
    """Not a weakening: two elements disagreeing is still an error."""

    path = tmp_path / "execution-server.yaml"
    _write_private_exact(path, b"SERVER:\n  NUM_CORES: 16\n")
    with pytest.raises(ContractError):
        _write_private_exact(path, b"SERVER:\n  NUM_CORES: 64\n")


def test_the_loser_of_the_race_sees_the_winners_bytes(tmp_path):
    """The race is resolved by reading, not by assuming."""

    path = tmp_path / "execution-server.yaml"
    payload = b"SERVER:\n  NUM_CORES: 16\n"

    real_open = os.open
    fired = {"once": False}

    def racing_open(target, flags, mode=0o777, **kw):
        # The window: the file appears between `exists()` and `open`.
        if (
            not fired["once"]
            and str(target) == str(path)
            and flags & os.O_EXCL
        ):
            fired["once"] = True
            Path(target).write_bytes(payload)
        return real_open(target, flags, mode, **kw)

    os.open = racing_open
    try:
        _write_private_exact(path, payload)
    finally:
        os.open = real_open

    assert fired["once"], "the race window was never entered"
    assert path.read_bytes() == payload
