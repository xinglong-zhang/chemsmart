"""HPC read-only tool wrappers."""

from __future__ import annotations

from typing import Any

from chemsmart.agent.transport import ExecTransport, SshExecTransport
from chemsmart.agent.wizard.probe import ALL_PROBE_SPECS, ProbeError

_OUTPUT_LIMIT = 4096
_TRANSPORT_FACTORY = SshExecTransport


def ssh_probe(
    server: str,
    probe_name: str,
    timeout_s: int = 15,
) -> dict[str, Any]:
    """Run a predefined read-only probe on a remote HPC server."""

    if probe_name not in ALL_PROBE_SPECS:
        return {
            "error": "unknown_probe",
            "available": sorted(ALL_PROBE_SPECS),
        }

    spec = ALL_PROBE_SPECS[probe_name]
    transport = _make_transport()
    try:
        result = transport.exec(spec.command, server, timeout_s)
    except TimeoutError:
        return {
            "error": "timeout",
            "server": server,
            "probe": probe_name,
        }
    except ProbeError as exc:
        return {
            "error": "probe_requires_slots",
            "message": str(exc),
            "server": server,
            "probe": probe_name,
        }

    parsed = None
    if spec.parser is not None and result.stdout.strip():
        try:
            parsed = spec.parser(result.stdout)
        except Exception:
            parsed = None
    return {
        "server": server,
        "probe": probe_name,
        "returncode": result.returncode,
        "stdout_truncated": _truncate(result.stdout),
        "stderr_truncated": _truncate(result.stderr),
        "parsed": parsed,
        "duration_s": result.duration_s,
    }


def _make_transport() -> ExecTransport:
    return _TRANSPORT_FACTORY()


def _truncate(text: str) -> str:
    return text[:_OUTPUT_LIMIT]
