"""Gaussian ``BatchJob`` for running a collection of Gaussian jobs."""

import logging
from typing import Any

from chemsmart.jobs.batch import BatchJob

logger = logging.getLogger(__name__)


class GaussianBatchJob(BatchJob):
    """Batch controller for a collection of Gaussian child jobs."""

    PROGRAM = "gaussian"

    def _configure_runner_for_node(
        self,
        runner: Any,
        node: str,
        job: Any,
    ) -> Any:
        """Pin child commands to a SLURM node via srun --nodelist."""
        return self._wrap_runner_command_for_node(runner=runner, node=node)
