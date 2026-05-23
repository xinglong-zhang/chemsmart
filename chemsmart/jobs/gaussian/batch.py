"""
Gaussian batch job implementation.

This module provides the Gaussian-specific ``GaussianBatchJob`` thin wrapper
around the shared batch orchestration provided by ``BatchJob``.
"""

import logging
from typing import Any

from chemsmart.jobs.batch import BatchJob

logger = logging.getLogger(__name__)


class GaussianBatchJob(BatchJob):
    """
    Generic Gaussian job class for running a batch of calculations.

    Shared batch orchestration is inherited from ``BatchJob``; this subclass
    only applies Gaussian-specific node-pinning behavior when needed.
    """

    PROGRAM = "gaussian"

    def _configure_runner_for_node(
        self,
        runner: Any,
        node: str,
        job: Any,
    ) -> Any:
        """Apply Gaussian runner node pinning through shared batch helper."""
        return self._wrap_runner_command_for_node(runner=runner, node=node)
