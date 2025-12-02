import logging
import os
from typing import Type

from chemsmart.analysis.iterate import IterateAnalyzer, SkeletonPreprocessor
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.job import Job
from chemsmart.jobs.runner import JobRunner
from chemsmart.jobs.iterate.settings import IterateJobSettings

logger = logging.getLogger(__name__)


class IterateJob(Job):


    PROGRAM = "Iterate"
    TYPE = "iterate"

    def __init__(
            self, 
            settings,
            label, 
            jobrunner, 
            **kwargs,
            ):
        logger.debug("Initializing IterateJob")
        super().__init__(
            settings=settings, jobrunner=jobrunner, **kwargs
        )
        