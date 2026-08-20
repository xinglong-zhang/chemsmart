"""The verification table reports what it observed, and never raises.

Adapted from aiida's `verdi computer test`: each check is independent,
failures are counted rather than fatal, and a noisy login shell -- the
classic corrupter of parsed scheduler replies -- is named explicitly.
"""

from __future__ import annotations

from types import SimpleNamespace

from chemsmart.settings.wizard import ServerChoicesV1, run_verification


def _choices(tmp_path, scheduler="SLURM"):
    return ServerChoicesV1(
        scheduler=scheduler,
        queue_name="compute",
        num_hours=24,
        mem_gb=52,
        num_cores=6,
        num_gpus=0,
        num_threads=6,
        submit_command="/opt/slurm/current/bin/sbatch" if scheduler else None,
        scratch_dir=str(tmp_path / "scratch"),
    )


def _runner(responses):
    def run(command, **_kwargs):
        key = " ".join(command)
        returncode, stdout, stderr = responses.get(key, (0, "", ""))
        return SimpleNamespace(
            returncode=returncode, stdout=stdout, stderr=stderr
        )

    return run


def test_a_clean_host_passes_every_applicable_check(tmp_path):
    checks = run_verification(
        _choices(tmp_path),
        runner=_runner({"sinfo --version": (0, "slurm 26.05.2\n", "")}),
        which=lambda _name: "/usr/bin/sbatch",
    )
    by_name = {check.name: check for check in checks}

    assert by_name["clean login shell"].status == "ok"
    assert by_name["scheduler responds"].status == "ok"
    assert by_name["identity"].status == "ok"
    assert by_name["scratch round trip"].status == "ok"
    assert by_name["submit command present"].status == "ok"


def test_a_noisy_login_shell_is_a_named_failure(tmp_path):
    checks = run_verification(
        _choices(tmp_path),
        runner=_runner(
            {
                "bash -lc echo -n": (0, "module loaded: gaussian\n", ""),
                "sinfo --version": (0, "ok", ""),
            }
        ),
        which=lambda _name: "/usr/bin/sbatch",
    )
    shell = next(c for c in checks if c.name == "clean login shell")

    assert shell.status == "fail"
    assert "corrupt parsed scheduler replies" in shell.detail
    assert "module loaded" in shell.detail


def test_a_dead_scheduler_fails_that_check_only(tmp_path):
    checks = run_verification(
        _choices(tmp_path),
        runner=_runner({"sinfo --version": (1, "", "cannot contact")}),
        which=lambda _name: "/usr/bin/sbatch",
    )
    by_name = {check.name: check for check in checks}

    assert by_name["scheduler responds"].status == "fail"
    assert by_name["scratch round trip"].status == "ok"


def test_local_choices_skip_scheduler_and_submit_checks(tmp_path):
    choices = ServerChoicesV1(
        scheduler=None,
        queue_name=None,
        num_hours=1,
        mem_gb=4,
        num_cores=4,
        num_gpus=0,
        num_threads=4,
        submit_command=None,
        scratch_dir=None,
    )
    checks = run_verification(
        choices, runner=_runner({}), which=lambda _n: None
    )
    by_name = {check.name: check for check in checks}

    assert by_name["scheduler responds"].status == "skipped"
    assert by_name["submit command present"].status == "skipped"
    assert by_name["scratch round trip"].status == "skipped"
