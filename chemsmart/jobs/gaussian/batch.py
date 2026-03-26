"""
Gaussian batch job implementation.

This module provides the Gaussian-specific ``GaussianBatchJob`` thin wrapper
around the shared batch orchestration provided by ``BatchJob``.
"""

import logging
import os
import types

from chemsmart.jobs.batch import BatchJob

logger = logging.getLogger(__name__)


class GaussianBatchJob(BatchJob):
    """
    Generic Gaussian job class for running a batch of calculations.

    Shared batch orchestration is inherited from ``BatchJob``; this subclass
    only applies Gaussian-specific node-pinning behavior when needed.
    """

    PROGRAM = "gaussian"

    def _configure_runner_for_node(self, runner, node, job):
        if not os.environ.get("SLURM_JOB_NODELIST"):
            return runner

        original_get_command = runner._get_command

        def patched_get_command_slurm(self_runner, job_obj):
            command = original_get_command(job_obj)
            prefix = f"srun --nodelist={node} --exclusive -N1 -n1 "
            return prefix + command

        runner._get_command = types.MethodType(
            patched_get_command_slurm,
            runner,
        )
        return runner
