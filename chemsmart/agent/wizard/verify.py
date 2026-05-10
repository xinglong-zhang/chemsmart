"""Read-only verification of wizard/server transport wiring."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

from chemsmart.agent.transport import build_submit_invocation
from chemsmart.agent.wizard.paths import server_yaml_path, validate_server_name
from chemsmart.settings.server import Server

_LOCAL_HOSTS = {"", "localhost", "local"}
_VERIFY_SCRIPT_PATH = "/tmp/chemsmart-wizard-verify.sh"
_VERIFY_WORKING_DIR = "/tmp"


@dataclass(frozen=True)
class VerifyResult:
    server_name: str
    host: str | None
    mode: str
    would_submit_via: str
    transport_invocation: list[str] | None
    warnings: list[str]
    errors: list[str]


def verify_server_yaml(server_name: str) -> VerifyResult:
    """Verify wizard/server transport wiring without submitting anything."""

    warnings: list[str] = []
    errors: list[str] = []
    try:
        path = _resolve_server_yaml_path(server_name)
    except ValueError as exc:
        errors.append(str(exc))
        return VerifyResult(
            server_name=server_name,
            host=None,
            mode="unknown",
            would_submit_via="unknown",
            transport_invocation=None,
            warnings=warnings,
            errors=errors,
        )

    try:
        raw_contents = yaml.safe_load(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        errors.append(f"Server YAML not found: {path}")
        return VerifyResult(
            server_name=server_name,
            host=None,
            mode="unknown",
            would_submit_via="unknown",
            transport_invocation=None,
            warnings=warnings,
            errors=errors,
        )
    except yaml.YAMLError as exc:
        errors.append(f"Invalid YAML in {path}: {exc}")
        return VerifyResult(
            server_name=server_name,
            host=None,
            mode="unknown",
            would_submit_via="unknown",
            transport_invocation=None,
            warnings=warnings,
            errors=errors,
        )

    server_block = (
        raw_contents.get("SERVER") if isinstance(raw_contents, dict) else None
    )
    if not isinstance(server_block, dict):
        errors.append(f"Missing SERVER mapping in {path}")
        return VerifyResult(
            server_name=server_name,
            host=None,
            mode="unknown",
            would_submit_via="unknown",
            transport_invocation=None,
            warnings=warnings,
            errors=errors,
        )

    if "HOST" not in server_block:
        warnings.append(
            "SERVER.HOST is missing; transport may fall back to the YAML "
            "filename stem."
        )

    try:
        server = Server.from_yaml(str(path))
    except Exception as exc:  # pragma: no cover - defensive normalization
        errors.append(f"Failed to load Server from YAML: {exc}")
        return VerifyResult(
            server_name=server_name,
            host=None,
            mode="unknown",
            would_submit_via="unknown",
            transport_invocation=None,
            warnings=warnings,
            errors=errors,
        )

    host = _normalize_host(server.kwargs.get("HOST"))
    mode = "local" if host in {None, "", "localhost", "local"} else "ssh"

    try:
        transport_invocation, _ = build_submit_invocation(
            script_path=_VERIFY_SCRIPT_PATH,
            working_dir=_VERIFY_WORKING_DIR,
            server=server,
        )
    except Exception as exc:  # pragma: no cover - defensive normalization
        errors.append(f"Failed to build transport invocation: {exc}")
        transport_invocation = None

    would_submit_via = _submit_mode_from_invocation(transport_invocation)
    if transport_invocation is not None and mode != would_submit_via:
        warnings.append(
            "HOST implies local mode, but transport would submit via ssh. "
            "This usually means a legacy YAML is relying on filename-stem "
            "fallback."
        )

    return VerifyResult(
        server_name=server_name,
        host=host,
        mode=mode,
        would_submit_via=would_submit_via,
        transport_invocation=transport_invocation,
        warnings=warnings,
        errors=errors,
    )


def _resolve_server_yaml_path(server_name: str) -> Path:
    validate_server_name(server_name)
    return server_yaml_path(server_name)


def _normalize_host(host: object) -> str | None:
    if not isinstance(host, str):
        return None
    value = host.strip()
    if not value:
        return None
    return value


def _submit_mode_from_invocation(invocation: list[str] | None) -> str:
    if not invocation:
        return "unknown"
    return "ssh" if invocation[0] == "ssh" else "local"
