from __future__ import annotations

from datetime import datetime, timedelta, timezone

UTC = timezone.utc

from chemsmart.agent.wizard.cache import CacheEntry, write_cache
from chemsmart.agent.wizard.parsers import QueueFacts
from chemsmart.agent.wizard.project import ProjectFinding
from chemsmart.agent.wizard.refresh import opportunistic_refresh, refresh_cache
from chemsmart.agent.wizard.scratch import ScratchFinding
from chemsmart.agent.wizard.software import (
    ModuleSystem,
    ProgramFinding,
    SoftwareSurvey,
)
from chemsmart.agent.wizard.survey import ScheduleSurvey


class StubRunner:
    pass


def _entry(*, status: str = "fresh", hours_ago: int = 0) -> CacheEntry:
    probed_at = (datetime.now(UTC) - timedelta(hours=hours_ago)).isoformat()
    return CacheEntry(
        server_name="perlmutter",
        host="login.cluster",
        mode="ssh",
        scheduler="SLURM",
        probed_at=probed_at.replace("+00:00", "Z"),
        source_commands={"scheduler": "sinfo --json"},
        partitions=[{"name": "debug", "default": True}],
        node_summary={"selected_queue": "debug", "cpu": 16},
        program_candidates={"gaussian": {"source": "path"}},
        status=status,
        last_error=None,
    )


def _write_server_yaml(tmp_path) -> None:
    path = tmp_path / ".chemsmart" / "server" / "perlmutter.yaml"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        """SERVER:
  SCHEDULER: SLURM
  HOST: login.cluster
""",
        encoding="utf-8",
    )


def _schedule() -> ScheduleSurvey:
    return ScheduleSurvey(
        scheduler="SLURM",
        submit_command="sbatch",
        queues=[
            QueueFacts(
                name="debug",
                default=True,
                max_walltime_hours=12,
                default_walltime_hours=4,
                default_mem_gb=64,
                default_cores=16,
                gpus_per_node=0,
                enabled=True,
                started=True,
            )
        ],
        chosen_queue="debug",
        evidence={"sinfo --json": "parsed"},
    )


def _software() -> SoftwareSurvey:
    return SoftwareSurvey(
        module_system=ModuleSystem(kind="lmod", version="8.7"),
        programs={
            "gaussian": ProgramFinding(
                program="gaussian",
                exefolder="/apps/gaussian",
                source="path",
                module_candidates=[],
                on_path=True,
            ),
            "orca": ProgramFinding(
                program="orca",
                exefolder=None,
                source="none",
                module_candidates=[],
                on_path=False,
            ),
            "nciplot": ProgramFinding(
                program="nciplot",
                exefolder=None,
                source="module",
                module_candidates=["nciplot/1.0"],
                on_path=False,
            ),
        },
        conda_base=None,
        conda_env=None,
    )


def test_refresh_cache_returns_existing_when_fresh(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))
    existing = _entry()
    write_cache(existing)

    monkeypatch.setattr(
        "chemsmart.agent.wizard.refresh.run_schedule_survey",
        lambda runner, topology: (_ for _ in ()).throw(
            AssertionError("unexpected probe")
        ),
    )

    refreshed = refresh_cache(StubRunner(), "perlmutter", force=False)

    assert refreshed == existing


def test_refresh_cache_stale_entry_reprobes_successfully(
    monkeypatch, tmp_path
):
    monkeypatch.setenv("HOME", str(tmp_path))
    _write_server_yaml(tmp_path)
    write_cache(_entry(hours_ago=30, status="stale"))
    monkeypatch.setattr(
        "chemsmart.agent.wizard.refresh.run_schedule_survey",
        lambda runner, topology: _schedule(),
    )
    monkeypatch.setattr(
        "chemsmart.agent.wizard.refresh.run_software_survey",
        lambda runner, topology: _software(),
    )
    monkeypatch.setattr(
        "chemsmart.agent.wizard.refresh.discover_scratch",
        lambda runner, topology: ScratchFinding(
            path="/scratch/user",
            source="env:SCRATCH",
            writable=True,
            candidates=[("env:SCRATCH", "/scratch/user")],
        ),
    )
    monkeypatch.setattr(
        "chemsmart.agent.wizard.refresh.discover_project",
        lambda runner, topology, scheduler: ProjectFinding(
            project="chem-123",
            source="groups",
            candidates=["chem-123"],
        ),
    )

    refreshed = refresh_cache(StubRunner(), "perlmutter", force=False)

    assert refreshed.status == "fresh"
    assert refreshed.scheduler == "SLURM"
    assert refreshed.node_summary["selected_queue"] == "debug"
    assert (
        refreshed.program_candidates["gaussian"]["exefolder"]
        == "/apps/gaussian"
    )


def test_opportunistic_refresh_preserves_last_good_cache_on_failure(
    monkeypatch,
    tmp_path,
):
    monkeypatch.setenv("HOME", str(tmp_path))
    write_cache(_entry(hours_ago=30))
    monkeypatch.setattr(
        "chemsmart.agent.wizard.refresh.refresh_cache",
        lambda runner, name, force=False: (_ for _ in ()).throw(
            RuntimeError("boom")
        ),
    )

    refreshed = opportunistic_refresh(StubRunner(), "perlmutter", ttl_hours=24)

    assert refreshed is not None
    assert refreshed.status == "stale"
    assert refreshed.last_error == "boom"


def test_opportunistic_refresh_returns_error_entry_without_prior_cache(
    monkeypatch,
    tmp_path,
):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setattr(
        "chemsmart.agent.wizard.refresh.refresh_cache",
        lambda runner, name, force=False: (_ for _ in ()).throw(
            RuntimeError("boom")
        ),
    )

    refreshed = opportunistic_refresh(StubRunner(), "perlmutter", ttl_hours=24)

    assert refreshed is not None
    assert refreshed.status == "error"
    assert refreshed.last_error == "boom"
