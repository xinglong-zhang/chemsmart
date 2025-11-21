"""
Gaussian computational chemistry job classes and utilities.

This module provides a comprehensive collection of Gaussian job types
for various computational chemistry calculations including optimization,
frequency analysis, single point calculations, IRC, NCI analysis,
and more. All job classes inherit from the base GaussianJob class
and provide specialized functionality for different calculation types.

The module also includes job runners and utilities for managing
Gaussian calculations in computational workflows.
"""

from .crest import GaussianCrestJob
from .custom import GaussianCustomJob
from .dias import GaussianDIASJob
from .irc import GaussianIRCJob
from .job import GaussianComJob, GaussianGeneralJob, GaussianJob
from .link import GaussianLinkJob
from .modred import GaussianModredJob
from .nci import GaussianNCIJob
from .opt import GaussianOptJob
from .qrc import GaussianQRCJob
from .resp import GaussianRESPJob
from .runner import GaussianJobRunner
from .scan import GaussianScanJob
from .singlepoint import GaussianSinglePointJob
from .tddft import GaussianTDDFTJob
from .traj import GaussianTrajJob
from .ts import GaussianTSJob
from .uvvis import GaussianUVVISJob
from .wbi import GaussianWBIJob

# Dynamically collect all registered Gaussian job subclasses
jobs = GaussianJob.subclasses()

__all__ = [
    "GaussianCrestJob",
    "GaussianCustomJob",
    "GaussianDIASJob",
    "GaussianIRCJob",
    "GaussianComJob",
    "GaussianGeneralJob",
    "GaussianJob",
    "GaussianLinkJob",
    "GaussianModredJob",
    "GaussianNCIJob",
    "GaussianOptJob",
    "GaussianQRCJob",
    "GaussianRESPJob",
    "GaussianJobRunner",
    "GaussianTrajJob",
    "GaussianScanJob",
    "GaussianSinglePointJob",
    "GaussianTDDFTJob",
    "GaussianTSJob",
    "GaussianUVVISJob",
    "GaussianWBIJob",
    "jobs",
]

# signals to the linter these imports are intentional
# imports as explicitly used

# If I comment all these out, I get the following error during run:
#   File "/Users/xinglongzhang/bin/chemsmart/chemsmart/jobs/runner.py", line 192, in from_job
#     raise ValueError(
# ValueError: Could not find any runners for job:
# GaussianOptJob<folder=<run/folder>, label=final_prd_opt_scan_gas_opt_opt>.
# Runners in registry: [].
#  Fake: True
#
