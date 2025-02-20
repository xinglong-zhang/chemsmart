from .job import XTBJob
from .opt import XTBOptJob

jobs = XTBJob.subclasses()

__all__ = [
    "XTBOptJob",
    "jobs",
]
