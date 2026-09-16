"""What the scheduler granted is what the engine is told it may use.

The execution profile written into the run directory carries the granted
allocation. The argv beside it carried `--num-cores`, `--num-gpus` and
`--mem-gb` taken from the approval bundle, and `JobRunner` prefers an
explicit CLI value over the profile -- so the profile was correct and
ignored.

The consequence is not a slow job. Reproduced by an adversarial review of
this round: a bundle asking 64 cores / 300 GB inside a 32-core / 154 GB
allocation produced `%pal nprocs 64 end` and a `%maxcore` sized for 300 GB.
That is an OOM kill or 64 ranks thrashing 32 cores, and it arrives as a
native program error -- a repairable terminal state -- so the woken
session is offered a *scientific* repair menu for a host arithmetic
decision it never made.

One number, one place: the profile. The argv stops carrying a second
opinion.
"""

from __future__ import annotations


from chemsmart.agent.execution import (
    build_execution_resource_spec,
    build_real_execution_argv,
)


def _resources(cores, memory_gb, gpus=0):
    return build_execution_resource_spec(
        execution_target="run",
        cores=cores,
        memory_gb=memory_gb,
        gpu_count=gpus,
        scratch_policy="server",
        node_timeout_seconds=900,
    )


def _argv(resources, server="/run/execution-server.yaml"):
    return build_real_execution_argv(
        compiled_argv=(
            "run",
            "-p",
            "proj.yaml",
            "orca",
            "opt",
        ),
        command_path=("run", "orca", "opt"),
        resources=resources,
        server=server,
    )


def test_the_argv_states_no_resource_the_profile_already_owns():
    argv = _argv(_resources(64, 300))
    for flag in ("--num-cores", "--mem-gb", "--num-gpus"):
        assert flag not in argv, (
            f"{flag} on the argv overrides the execution profile in "
            "JobRunner, so the granted allocation is written down and "
            f"then ignored: {argv}"
        )
    assert "--server" in argv, argv


def test_the_engine_reads_the_granted_allocation_through_the_runner(
    tmp_path,
):
    """Drive the real consumer: profile -> Server -> JobRunner -> writer.

    Asserting the argv alone would only prove a spelling. This builds the
    profile the executor writes, resolves it exactly as `chemsmart run`
    does, and reads what ORCA is actually told.
    """

    import io

    from chemsmart.agent.live_session import _write_execution_server_profile
    from chemsmart.jobs.orca.writer import ORCAInputWriter
    from chemsmart.jobs.runner import JobRunner
    from chemsmart.settings.server import Server

    run_directory = tmp_path / "run"
    run_directory.mkdir()
    # What the scheduler granted.
    granted = _resources(32, 154)
    profile = _write_execution_server_profile(run_directory, granted)

    # What the approval asked for, which is what the argv used to carry.
    argv = _argv(_resources(64, 300), server=str(profile))
    assert "--num-cores" not in argv

    server = Server.from_yaml(str(profile))
    runner = JobRunner(server=server, scratch=False)
    assert runner.num_cores == 32
    assert float(runner.mem_gb) == 154.0

    writer = ORCAInputWriter.__new__(ORCAInputWriter)
    writer.jobrunner = runner
    processors = io.StringIO()
    writer._write_processors(processors)
    assert (
        "%pal nprocs 32 end" in processors.getvalue()
    ), "ORCA was told to launch more MPI ranks than the allocation has"
    memory = io.StringIO()
    writer._write_memory(memory)
    expected = int((154 * 1000) / 32 * 0.75)
    assert f"%maxcore {expected}" in memory.getvalue(), (
        "per-rank memory was computed from an allocation the job does "
        "not have"
    )


def test_the_reviewed_argv_and_the_launched_argv_still_agree():
    """The launch is compared against the review before it runs.

    Both sides now omit the resource flags, so removing them keeps that
    comparison meaningful rather than silently diverging.
    """

    reviewed = _argv(_resources(64, 300))
    launched = _argv(_resources(32, 154))
    assert reviewed == launched, (
        "the review and the launch disagree about the command, which is "
        "the comparison that guards an approved execution"
    )


def test_a_scratch_policy_is_still_carried(tmp_path):
    """Only the resource numbers move to the profile; scratch is a
    behaviour the profile does not express."""

    assert "--scratch" in _argv(_resources(4, 8))
    none_scratch = build_execution_resource_spec(
        execution_target="run",
        cores=4,
        memory_gb=8,
        gpu_count=0,
        scratch_policy="none",
        node_timeout_seconds=900,
    )
    assert "--no-scratch" in _argv(none_scratch)
