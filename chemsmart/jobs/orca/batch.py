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

    Engine-specific input writing and command execution remain delegated to the
    child ORCA jobs and their associated job runners.
    """

    PROGRAM = "orca"

    def _configure_runner_for_node(
        self,
        runner: Any,
        node: str,
        job: Any,
    ) -> Any:
        return runner


ORCABatchJob = OrcaBatchJob
