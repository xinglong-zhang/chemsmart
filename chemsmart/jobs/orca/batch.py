"""
ORCA batch job implementation.

This module provides the ORCA-specific thin batch wrapper around the shared
``BatchJob`` orchestration layer.
"""

from typing import Any

from chemsmart.jobs.batch import BatchJob


class OrcaBatchJob(BatchJob):
    """
    Generic ORCA batch job class for running a collection of ORCA jobs.

    Shared batch orchestration is inherited from ``BatchJob``; this subclass
    only applies ORCA runner node-pinning behavior when needed.
    """

    PROGRAM = "orca"

    def _configure_runner_for_node(
        self,
        runner: Any,
        node: str,
        job: Any,
    ) -> Any:
        """Apply ORCA runner node pinning through shared batch helper."""
        return self._wrap_runner_command_for_node(runner=runner, node=node)


ORCABatchJob = OrcaBatchJob
