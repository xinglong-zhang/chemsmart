"""Wizard tool wrappers registered with the agent tool registry."""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
from typing import Any

from chemsmart.agent.wizard.cache import cache_path
from chemsmart.agent.wizard.orchestrator import run_wizard
from chemsmart.agent.wizard.probe import ProbeRunner
from chemsmart.agent.wizard.refresh import _recover_cache_entry, refresh_cache
from chemsmart.agent.wizard.verify import verify_server_yaml
from chemsmart.agent.wizard.write import write_server_yaml


def wizard_probe(
    server_name: str,
    ssh_host_hint: str | None = None,
) -> dict[str, Any]:
    """Probe a local or remote HPC target and render a candidate server YAML."""

    outcome = run_wizard(
        ProbeRunner,
        server_name=server_name,
        ssh_host_hint=ssh_host_hint,
        write=False,
    )
    return {
        "ok": outcome.validation.ok,
        "server_name": server_name,
        "ssh_host_hint": ssh_host_hint,
        "risk": "remote_probe" if ssh_host_hint else "local_probe",
        "yaml_text": outcome.plan.text,
        "plan": _json_safe(outcome.plan),
        "validation": _json_safe(outcome.validation),
        "surveys": _json_safe(outcome.surveys),
        "written_path": outcome.written_path,
    }


def wizard_write(
    server_name: str,
    yaml_text: str,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Write a wizard-rendered server YAML into ``~/.chemsmart/server``."""

    path = write_server_yaml(
        name=server_name,
        yaml_text=yaml_text,
        overwrite=overwrite,
    )
    return {
        "ok": True,
        "server_name": server_name,
        "written_path": path,
        "overwrite": overwrite,
    }


def wizard_verify(server_name: str) -> dict[str, Any]:
    """Verify wizard/server transport wiring for an existing YAML."""

    return _json_safe(verify_server_yaml(server_name))


def wizard_refresh(server_name: str, force: bool = False) -> dict[str, Any]:
    """Refresh or reuse the wizard node sidecar cache for a server."""

    try:
        entry = refresh_cache(ProbeRunner, server_name, force=force)
    except Exception as exc:  # pragma: no cover - bridged to cache recovery
        entry = _recover_cache_entry(server_name, exc)
    result = _json_safe(entry)
    result["ok"] = entry.status == "fresh"
    result["cache_path"] = str(cache_path(server_name))
    result["force"] = force
    return result


def _json_safe(value: Any) -> Any:
    if is_dataclass(value):
        return {key: _json_safe(item) for key, item in asdict(value).items()}
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_safe(item) for item in value]
    if isinstance(value, tuple):
        return [_json_safe(item) for item in value]
    return value
