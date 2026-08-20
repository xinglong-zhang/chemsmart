"""Detection believes evidence, not the mere existence of a binary.

The legacy detector treats an installed squeue as proof of SLURM even when
the controller is down; the wizard's detector requires the query command to
answer, records its evidence, and tells the two qsub-owning schedulers
apart before naming either.
"""

from __future__ import annotations

from types import SimpleNamespace

from chemsmart.settings.probe import detect_scheduler
from chemsmart.settings.probe.localhost import scratch_candidates


def _runner(responses):
    def run(command, **_kwargs):
        key = " ".join(command)
        returncode, stdout = responses.get(key, (127, ""))
        return SimpleNamespace(returncode=returncode, stdout=stdout, stderr="")

    return run


def test_a_responding_sinfo_is_slurm_with_named_evidence():
    detection = detect_scheduler(
        env={},
        which=lambda name: (
            f"/usr/bin/{name}" if name in {"sbatch", "sinfo"} else None
        ),
        runner=_runner({"sinfo --version": (0, "slurm 26.05.2\n")}),
    )
    assert detection.scheduler == "SLURM"
    assert any("slurm 26.05.2" in item for item in detection.evidence)


def test_an_installed_but_dead_slurm_is_not_claimed():
    detection = detect_scheduler(
        env={},
        which=lambda name: (
            f"/usr/bin/{name}" if name in {"sbatch", "sinfo"} else None
        ),
        runner=_runner({"sinfo --version": (1, "")}),
    )
    assert detection.scheduler is None
    assert any("controller may be down" in item for item in detection.evidence)


def test_qsub_with_sge_root_is_sge_not_pbs():
    detection = detect_scheduler(
        env={"SGE_ROOT": "/opt/sge"},
        which=lambda name: "/usr/bin/qsub" if name == "qsub" else None,
        runner=_runner({}),
    )
    assert detection.scheduler == "SGE"


def test_qsub_with_pbsnodes_is_pbs():
    detection = detect_scheduler(
        env={},
        which=lambda name: (
            f"/usr/bin/{name}" if name in {"qsub", "pbsnodes"} else None
        ),
        runner=_runner({"qstat --version": (0, "pbs_version = 2024.1\n")}),
    )
    assert detection.scheduler == "PBS"
    assert any("pbs_version" in item for item in detection.evidence)


def test_a_bare_host_is_local_with_the_reason_stated():
    detection = detect_scheduler(
        env={}, which=lambda _n: None, runner=_runner({})
    )
    assert detection.scheduler is None
    assert any("local execution only" in item for item in detection.evidence)


def test_scratch_candidates_prefer_the_user_scratch(tmp_path):
    (tmp_path / "scratch" / "alice").mkdir(parents=True)
    (tmp_path / "tmp").mkdir()
    env = {"SCRATCH": str(tmp_path / "scratch" / "alice")}

    import chemsmart.settings.probe.localhost as localhost

    found = scratch_candidates(env, user="alice")
    assert found[0] == str(tmp_path / "scratch" / "alice")
    assert localhost is not None
