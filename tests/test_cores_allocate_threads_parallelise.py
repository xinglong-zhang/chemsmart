"""Cores are allocated; threads are what a process runs. Not one number.

``NUM_CORES`` answered three different questions at once: how many CPUs the
scheduler allocates, how many MPI processes ORCA launches, and how many
OpenMP threads Gaussian, PySCF and xTB run. ``NUM_THREADS`` sat in every
server YAML the wizard has ever written, was defined as a property on
``Server``, and was read by nothing.

The distinction matters per program and is not a matter of taste:

- **ORCA** is MPI. ``%pal nprocs`` is a count of *processes*, one per
  allocated core, and ``%maxcore`` is memory *per rank* computed by
  dividing the total by that same count. Handing ORCA a smaller
  threads-per-process number would launch fewer ranks than the allocation
  paid for **and** give each rank too much memory, because the denominator
  would shrink while the total stayed. ORCA's runner sets no OpenMP
  variable at all. So ORCA reads the allocation, and must not be forced
  into the threads abstraction.
- **Gaussian, PySCF, xTB** are shared-memory. ``%nprocshared`` and the
  OMP/MKL/OpenBLAS variables are threads inside one process, which is
  exactly what a threads-per-process control means.

Owner ruling, 2026-09-16: keep NUM_CORES as the allocation authority, wire
NUM_THREADS as threads-per-process, and do not force ORCA into it.
"""

from __future__ import annotations

from chemsmart.settings.server import Server


def _server(**overrides):
    base = dict(
        SCHEDULER="SLURM",
        SUBMIT_COMMAND="sbatch",
        NUM_CORES=64,
        MEM_GB=160,
        NUM_HOURS=24,
        NUM_GPUS=0,
        QUEUE_NAME="chpc",
    )
    base.update(overrides)
    return Server("canned", **base)


def test_threads_default_to_the_allocation_so_nothing_changes_today():
    """Every profile ever written must behave exactly as it did.

    The old default was a bare 16, unrelated to the machine: a profile
    with 64 cores and no NUM_THREADS would have run 16 threads the moment
    anything started reading the key.
    """

    server = _server()
    assert "NUM_THREADS" not in server.kwargs
    assert server.num_threads == 64, (
        "threads must fall back to the allocation, not to an unrelated "
        f"constant; got {server.num_threads}"
    )


def test_an_operator_can_run_fewer_threads_than_cores():
    server = _server(NUM_THREADS=16)
    assert (server.num_cores, server.num_threads) == (64, 16)


def test_the_runner_carries_both_numbers():
    from chemsmart.jobs.runner import JobRunner

    runner = JobRunner(server=_server(NUM_THREADS=16), scratch=False)
    assert runner.num_cores == 64
    assert runner.num_threads == 16


def test_orca_launches_one_rank_per_allocated_core(tmp_path):
    """ORCA reads the allocation, never the threads-per-process number."""

    from chemsmart.jobs.runner import JobRunner

    runner = JobRunner(server=_server(NUM_THREADS=16), scratch=False)
    assert runner.num_cores == 64

    import io

    from chemsmart.jobs.orca.writer import ORCAInputWriter

    writer = ORCAInputWriter.__new__(ORCAInputWriter)
    writer.jobrunner = runner
    buffer = io.StringIO()
    writer._write_processors(buffer)
    assert "%pal nprocs 64 end" in buffer.getvalue(), buffer.getvalue()

    memory = io.StringIO()
    writer._write_memory(memory)
    # Memory per RANK: the total divided by the rank count, not by the
    # thread count. A smaller denominator would over-commit the node.
    expected = int((160 * 1000) / 64 * 0.75)
    assert f"%maxcore {expected}" in memory.getvalue(), memory.getvalue()


def test_gaussian_shares_memory_across_threads(tmp_path):
    from chemsmart.jobs.runner import JobRunner

    runner = JobRunner(server=_server(NUM_THREADS=16), scratch=False)
    import io
    from types import SimpleNamespace

    from chemsmart.jobs.gaussian.writer import GaussianInputWriter

    writer = GaussianInputWriter.__new__(GaussianInputWriter)
    writer.jobrunner = runner
    buffer = io.StringIO()
    writer.settings = SimpleNamespace(chk=False)
    writer.job = SimpleNamespace(label="j")
    writer._write_gaussian_header(buffer)
    text = buffer.getvalue()
    assert "%nprocshared=16" in text, text


def test_the_threaded_programs_bind_openmp_to_threads_not_cores():
    """PySCF and xTB bind OpenMP, which is threads inside one process."""

    import ast
    from pathlib import Path

    import chemsmart.jobs.pyscf.runner as pyscf_runner
    import chemsmart.jobs.xtb.runner as xtb_runner

    for module in (pyscf_runner, xtb_runner):
        source = Path(module.__file__).read_text(encoding="utf-8")
        tree = ast.parse(source)
        omp_lines = [
            node.lineno
            for node in ast.walk(tree)
            if isinstance(node, ast.Constant)
            and isinstance(node.value, str)
            and node.value.endswith("OMP_NUM_THREADS")
        ]
        assert omp_lines, f"{module.__name__} binds no OpenMP variable"
        assert "self.num_threads" in source, (
            f"{module.__name__} binds OpenMP from something other than "
            "the threads-per-process number"
        )
        assert (
            "str(self.num_cores)" not in source
        ), f"{module.__name__} still binds OpenMP to the allocation"


def test_orca_binds_no_openmp_variable_at_all():
    """Which is why ORCA must not be given a threads number.

    ORCA parallelises with MPI ranks and sets no OpenMP variable, so a
    threads-per-process value would have no consumer there -- and using
    it for %pal nprocs would under-fill the allocation while %maxcore
    over-committed the node.
    """

    from pathlib import Path

    import chemsmart.jobs.orca.runner as orca_runner

    source = Path(orca_runner.__file__).read_text(encoding="utf-8")
    assert "OMP_NUM_THREADS" not in source, (
        "ORCA binds an OpenMP variable, so the threads-per-process "
        "number would have a consumer here after all"
    )
    # It mentions num_threads only in a debug line; what matters is that
    # nothing ORCA writes is derived from it.
    from chemsmart.jobs.orca import writer as orca_writer

    writer_source = Path(orca_writer.__file__).read_text(encoding="utf-8")
    assert "num_threads" not in writer_source, (
        "an ORCA input field is derived from threads-per-process; %pal "
        "nprocs and %maxcore are both allocation quantities"
    )
