"""
ORCA job module initialization.

This module provides ORCA job functionality including various job types,
settings, and runners for ORCA quantum chemistry calculations.
"""

from .batch import ORCABatchJob, OrcaBatchJob
from .irc import ORCAIRCJob
from .job import ORCAGeneralJob, ORCAInpJob, ORCAJob
from .modred import ORCAModredJob
from .opt import ORCAOptJob
from .pka import ORCApKaJob
from .qmmm import ORCAQMMMJob
from .qrc import ORCAQRCJob
from .runner import ORCAJobRunner
from .scan import ORCAScanJob
from .singlepoint import ORCASinglePointJob
from .ts import ORCATSJob

# Get all available ORCA job subclasses
jobs = ORCAJob.subclasses()


__all__ = [
    "ORCAOptJob",
    "ORCAIRCJob",
    "OrcaBatchJob",
    "ORCABatchJob",
    "ORCAJob",
    "ORCAInpJob",
    "ORCAGeneralJob",
    "ORCAModredJob",
    "ORCAJobRunner",
    "ORCAQRCJob",
    "ORCAScanJob",
    "ORCASinglePointJob",
    "ORCATSJob",
    "ORCAQMMMJob",
    "ORCApKaJob",
    "jobs",
]
