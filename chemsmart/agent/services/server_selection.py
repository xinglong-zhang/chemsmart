"""Shared server resolution for validation and submission services."""

from __future__ import annotations

from chemsmart.settings.server import Server
from chemsmart.settings.user import ChemsmartUserSettings


def coerce_server(server) -> Server:
    """Resolve a configured server object from a name or existing object."""

    if server is None:
        server = default_submit_server_name()
    if isinstance(server, Server):
        return server
    if isinstance(server, str):
        return Server.from_servername(server)
    raise TypeError("server must be a chemsmart.settings.server.Server or str")


def default_submit_server_name() -> str:
    """Return the sole configured server or require an explicit selection."""

    settings = ChemsmartUserSettings()
    available = sorted(settings.all_available_servers)
    if len(available) == 1:
        return available[0]
    if not available:
        raise ValueError(
            "submit_hpc requires server when no configured servers are available"
        )
    raise ValueError(
        "submit_hpc requires server when multiple configured servers are "
        f"available: {', '.join(available)}"
    )


__all__ = ["coerce_server", "default_submit_server_name"]
