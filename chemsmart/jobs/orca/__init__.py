"""
ORCA job module initialization.

This module provides ORCA job functionality including various job types,
settings, and runners for ORCA quantum chemistry calculations.
"""

from .irc import ORCAIRCJob
from .job import ORCAGeneralJob, ORCAInpJob, ORCAJob
from .modred import ORCAModredJob
from .opt import ORCAOptJob
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
    "ORCAJob",
    "ORCAInpJob",
    "ORCAGeneralJob",
    "ORCAModredJob",
    "ORCAJobRunner",
    "ORCAQRCJob",
    "ORCAScanJob",
    "ORCASinglePointJob",
    "ORCATSJob",
    "jobs",
]
