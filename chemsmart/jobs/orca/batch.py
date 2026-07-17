"""ORCA ``BatchJob`` for running a collection of ORCA jobs."""

from typing import Any

from chemsmart.jobs.batch import BatchJob


class OrcaBatchJob(BatchJob):
    """Batch controller for a collection of ORCA child jobs."""

    PROGRAM = "orca"

    def _configure_runner_for_node(
        self,
        runner: Any,
        node: str,
        job: Any,
    ) -> Any:
        """Pin child commands to a SLURM node via srun --nodelist."""
        return self._wrap_runner_command_for_node(runner=runner, node=node)


ORCABatchJob = OrcaBatchJob
