"""
ORCA job module initialization.

This module provides ORCA job functionality including various job types,
settings, and runners for ORCA quantum chemistry calculations.
"""

from .fukui import ORCAFukuiJob
from .irc import ORCAIRCJob
from .job import ORCAGeneralJob, ORCAInpJob, ORCAJob
from .modred import ORCAModredJob
from .opt import ORCAOptJob
from .qmmm import ORCAQMMMJob
from .qrc import ORCAQRCJob
from .redox import ORCARedoxJob
from .runner import ORCAJobRunner
from .scan import ORCAScanJob
from .singlepoint import ORCASinglePointJob
from .ts import ORCATSJob

# Get all available ORCA job subclasses
jobs = ORCAJob.subclasses()


__all__ = [
    "ORCAFukuiJob",
    "ORCAOptJob",
    "ORCAIRCJob",
    "ORCAJob",
    "ORCAInpJob",
    "ORCAGeneralJob",
    "ORCAModredJob",
    "ORCAJobRunner",
    "ORCAQRCJob",
    "ORCARedoxJob",
    "ORCAScanJob",
    "ORCASinglePointJob",
    "ORCATSJob",
    "ORCAQMMMJob",
    "jobs",
]
