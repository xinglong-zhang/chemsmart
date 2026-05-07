"""
chemsmart.agent — AI-scientist agent layer (Wave 1 scaffold).

AgentSession and run_agent are placeholders; full implementation
arrives in subsequent waves (see bin/plan.md).
"""

from chemsmart.agent.tools import (
    build_gaussian_settings,
    build_job,
    build_molecule,
    build_orca_settings,
    dry_run_input,
    recommend_method,
)


class AgentSession:
    def __init__(self, *args, **kwargs):
        raise NotImplementedError


def run_agent(*args, **kwargs):
    raise NotImplementedError


__all__ = [
    "AgentSession",
    "run_agent",
    "build_molecule",
    "build_gaussian_settings",
    "build_orca_settings",
    "build_job",
    "dry_run_input",
    "recommend_method",
]
