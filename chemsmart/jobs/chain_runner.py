"""Job runners for standalone chain workflows."""

import logging

from chemsmart.jobs.runner import JobRunner

logger = logging.getLogger(__name__)


class ChainJobRunner(JobRunner):
    """Runner that identifies ``ChainJob`` for ``JobRunner.from_job``.

    Chain jobs orchestrate child phases in ``ChainMixin._run``. This runner
    does not launch an external program.
    """

    JOBTYPES = ["chain"]
    PROGRAM = "chain"
    FAKE = False
    SCRATCH = None

    def __init__(
        self, server, scratch=None, fake=False, scratch_dir=None, **kwargs
    ):
        if scratch is None:
            scratch = self.SCRATCH
        super().__init__(
            server=server,
            scratch=scratch,
            scratch_dir=scratch_dir,
            fake=fake,
            **kwargs,
        )
        logger.debug("ChainJobRunner initialized")
        logger.debug(f"Jobrunner server: {self.server}")
        logger.debug(f"Jobrunner scratch: {self.scratch}")
        logger.debug(f"Jobrunner fake mode: {self.fake}")

    @property
    def executable(self):
        return None

    def _get_command(self, job):
        return None

    def _create_process(self, job, command, env):
        return None

    def run(self, job, **kwargs):
        """Chain jobs orchestrate child phases; no external process is run."""
        logger.debug(
            "ChainJobRunner does not execute %s; phases run via ChainMixin.",
            job,
        )


class FakeChainJobRunner(ChainJobRunner):
    """Fake chain runner for CLI ``--fake`` and tests."""

    FAKE = True

    def __init__(self, server, scratch=None, fake=True, **kwargs):
        super().__init__(server=server, scratch=scratch, fake=fake, **kwargs)
