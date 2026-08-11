"""
PySCF job module initialization.

Importing this package registers the PySCF job and jobrunner subclasses in
the ``RegistryMixin`` registries. ``JobRunner.from_job`` can only select a
runner whose module has been imported, so these imports are load-bearing
rather than convenience re-exports.
"""

from .hess import PySCFHessJob
from .job import PySCFGeneralJob, PySCFJob
from .opt import PySCFOptJob
from .runner import (
    FakePySCFJobRunner,
    PySCFJobRunner,
    PySCFPreflightError,
    PySCFResultValidationError,
)
from .settings import PySCFJobSettings
from .singlepoint import PySCFSinglePointJob
from .td import PySCFTDJob

# Get all available PySCF job subclasses
jobs = PySCFJob.subclasses()

__all__ = [
    "FakePySCFJobRunner",
    "PySCFGeneralJob",
    "PySCFHessJob",
    "PySCFJob",
    "PySCFJobRunner",
    "PySCFPreflightError",
    "PySCFResultValidationError",
    "PySCFJobSettings",
    "PySCFOptJob",
    "PySCFSinglePointJob",
    "PySCFTDJob",
    "jobs",
]
