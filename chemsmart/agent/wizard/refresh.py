"""On-demand wizard node-refresh cache orchestration."""

from __future__ import annotations

import logging
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path

import yaml

from chemsmart.agent.wizard.cache import (
    CacheEntry,
    is_stale,
    load_cache,
    mark_status,
    write_cache,
)
from chemsmart.agent.wizard.paths import server_yaml_path, validate_server_name
from chemsmart.agent.wizard.project import discover_project
from chemsmart.agent.wizard.scratch import discover_scratch
from chemsmart.agent.wizard.software import run_software_survey
from chemsmart.agent.wizard.survey import run_schedule_survey
from chemsmart.agent.wizard.topology import Topology

logger = logging.getLogger(__name__)

_LOCAL_HOSTS = {"", "localhost", "local"}


def refresh_cache(runner, name: str, force: bool = False) -> CacheEntry:
    """Refresh a node sidecar cache entry for a server configuration."""

    validate_server_name(name)
    existing = load_cache(name)
    if existing is not None and not force and not is_stale(existing):
        return existing

    context = _load_server_context(name)
    topology = _topology_from_context(context)
    schedule_survey = run_schedule_survey(runner, topology)
    software_survey = run_software_survey(runner, topology)
    scratch_finding = discover_scratch(runner, topology)
    project_finding = discover_project(
        runner,
        topology,
        schedule_survey.scheduler,
    )

    entry = CacheEntry(
        server_name=name,
        host=context["host"],
        mode=context["mode"],
        scheduler=schedule_survey.scheduler,
        probed_at=_utc_now_isoformat(),
        source_commands=_build_source_commands(schedule_survey.scheduler),
        partitions=[asdict(queue) for queue in schedule_survey.queues],
        node_summary=_build_node_summary(
            schedule_survey,
            scratch_path=scratch_finding.path,
            scratch_writable=scratch_finding.writable,
            project=project_finding.project,
        ),
        program_candidates=_build_program_candidates(software_survey),
        status="fresh",
        last_error=None,
    )
    write_cache(entry)
    return entry


def opportunistic_refresh(
    runner,
    name: str,
    ttl_hours: int = 24,
) -> CacheEntry | None:
    """Refresh stale cache data without raising on failures."""

    validate_server_name(name)
    existing = load_cache(name)
    if existing is not None and not is_stale(existing, ttl_hours=ttl_hours):
        return existing

    try:
        return refresh_cache(runner, name, force=True)
    except Exception as exc:  # pragma: no cover - defensive wrapper
        logger.warning(
            "Wizard cache refresh failed for %s: %s",
            name,
            exc,
        )
        return _recover_cache_entry(name, exc)


def _recover_cache_entry(name: str, exc: Exception) -> CacheEntry:
    existing = load_cache(name)
    if existing is not None:
        entry = mark_status(existing, "stale", last_error=str(exc))
    else:
        entry = _build_error_entry(name, exc)
    write_cache(entry)
    return entry


def _build_error_entry(name: str, exc: Exception) -> CacheEntry:
    context = _load_server_context_or_default(name)
    return CacheEntry(
        server_name=name,
        host=context["host"],
        mode=context["mode"],
        scheduler=context["scheduler"],
        probed_at=_utc_now_isoformat(),
        source_commands={},
        partitions=[],
        node_summary={},
        program_candidates={},
        status="error",
        last_error=str(exc),
    )


def _load_server_context(name: str) -> dict[str, str | None]:
    path = _resolve_server_yaml_path(name)
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError(f"Missing YAML mapping in {path}")
    server_block = raw.get("SERVER")
    if not isinstance(server_block, dict):
        raise ValueError(f"Missing SERVER mapping in {path}")

    raw_host = server_block.get("HOST")
    host = raw_host.strip() if isinstance(raw_host, str) else None
    if host is None or host in _LOCAL_HOSTS:
        mode = "local"
        host = None
    else:
        mode = "ssh"
    return {
        "host": host,
        "mode": mode,
        "scheduler": _clean_optional_str(server_block.get("SCHEDULER")),
    }


def _load_server_context_or_default(name: str) -> dict[str, str | None]:
    try:
        return _load_server_context(name)
    except Exception:
        return {"host": None, "mode": "unknown", "scheduler": None}


def _resolve_server_yaml_path(server_name: str) -> Path:
    validate_server_name(server_name)
    return server_yaml_path(server_name)


def _topology_from_context(context: dict[str, str | None]) -> Topology:
    mode = str(context.get("mode") or "unknown")
    host = context.get("host")
    if mode == "ssh" and host:
        return Topology(mode="B", host=host, evidence=[f"yaml:HOST={host}"])
    return Topology(
        mode="A", host="localhost", evidence=["yaml:HOST=localhost"]
    )


def _build_source_commands(scheduler: str) -> dict[str, str]:
    scheduler_name = scheduler.upper()
    commands = {
        "topology": "printenv; sinfo; qstat; bqueues; qconf",
        "software": "type module; module --version; module -t avail; command -v; which; readlink -f; printenv CONDA_PREFIX; conda info --base",
        "scratch": "printenv SCRATCH/WORK/TMPDIR/HOME; test -d <home>/scratch -a -w <home>/scratch; test -w <candidate>",
        "project": "printenv SBATCH_ACCOUNT/SLURM_ACCOUNT/PBS_ACCOUNT/USER; sacctmgr show user ...; groups",
    }
    scheduler_commands = {
        "SLURM": "sinfo --json; scontrol show partition --oneliner",
        "PBS": "qstat -Q -f -F json; qstat -Q -f",
        "LSF": "bqueues -o 'queue_name status runlimit memlimit prolimit' -json; bqueues -l",
        "SGE": "qconf -sql; qconf -sq <queue>",
    }
    if scheduler_name in scheduler_commands:
        commands["scheduler"] = scheduler_commands[scheduler_name]
    return commands


def _build_node_summary(
    schedule_survey,
    *,
    scratch_path: str | None,
    scratch_writable: bool,
    project: str | None,
) -> dict:
    selected_queue = next(
        (
            queue
            for queue in schedule_survey.queues
            if queue.name == schedule_survey.chosen_queue
        ),
        None,
    )
    return {
        "selected_queue": schedule_survey.chosen_queue,
        "queue_count": len(schedule_survey.queues),
        "enabled_queue_count": sum(
            1 for queue in schedule_survey.queues if queue.enabled
        ),
        "started_queue_count": sum(
            1 for queue in schedule_survey.queues if queue.started
        ),
        "gpu_queue_count": sum(
            1
            for queue in schedule_survey.queues
            if (queue.gpus_per_node or 0) > 0
        ),
        "cpu": (
            None if selected_queue is None else selected_queue.default_cores
        ),
        "mem_gb": (
            None if selected_queue is None else selected_queue.default_mem_gb
        ),
        "gpu": (
            None
            if selected_queue is None
            else selected_queue.gpus_per_node or 0
        ),
        "state": {
            "enabled": (
                None if selected_queue is None else selected_queue.enabled
            ),
            "started": (
                None if selected_queue is None else selected_queue.started
            ),
        },
        "scratch_dir": scratch_path,
        "scratch_writable": scratch_writable,
        "project": project,
    }


def _build_program_candidates(software_survey) -> dict:
    return {
        "module_system": {
            "kind": software_survey.module_system.kind,
            "version": software_survey.module_system.version,
        },
        "conda": {
            "base": software_survey.conda_base,
            "env": software_survey.conda_env,
        },
        **{
            name: {
                "exefolder": finding.exefolder,
                "source": finding.source,
                "module_candidates": list(finding.module_candidates),
                "on_path": finding.on_path,
            }
            for name, finding in software_survey.programs.items()
        },
    }


def _clean_optional_str(value: object) -> str | None:
    if not isinstance(value, str):
        return None
    cleaned = value.strip()
    return cleaned or None


def _utc_now_isoformat() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")
