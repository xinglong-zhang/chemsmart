from .irc import ORCAIRCJob
from .job import ORCAGeneralJob, ORCAInpJob, ORCAJob
from .modred import ORCAModredJob
from .opt import ORCAOptJob
from .qmmm import ORCAQMMMJob
from .runner import ORCAJobRunner
from .scan import ORCAScanJob
from .singlepoint import ORCASinglePointJob
from .ts import ORCATSJob

jobs = ORCAJob.subclasses()


__all__ = [
    "ORCAOptJob",
    "ORCAIRCJob",
    "ORCAJob",
    "ORCAInpJob",
    "ORCAGeneralJob",
    "ORCAModredJob",
    "ORCAJobRunner",
    "ORCAScanJob",
    "ORCASinglePointJob",
    "ORCATSJob",
    "ORCAQMMMJob",
    "jobs",
]
