"""ORCA TDDFT / TDA excited-state job implementation."""

import logging

from chemsmart.jobs.orca.job import ORCAJob

logger = logging.getLogger(__name__)


class ORCATDDFTJob(ORCAJob):
    """ORCA TDDFT / TDA excited-state calculation job.

    Uses :class:`ORCATDDFTJobSettings` and delegates all input-file
    generation to :class:`ORCAInputWriter`, which emits the ``%tddft`` and
    optional ``%rel`` blocks based on the settings.  A bare TD job is
    vertical; excited-state ``Opt``/``Freq`` are opt-in via the settings
    flags (``opt_excited`` / ``freq`` / ``numfreq``) or, in the CLI, via
    ``-r Opt`` / ``-r Freq`` on the ``orca`` group.
    """

    TYPE = "orcatd"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
