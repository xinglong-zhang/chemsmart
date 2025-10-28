"""
NCIPLOT job module initialization.

This module provides NCIPLOT job functionality including job classes and
runners for non-covalent interaction analysis.
"""

from .job import NCIPLOTJob
from .runner import NCIPLOTJobRunner

# Get all available NCIPLOT job subclasses
jobs = NCIPLOTJob.subclasses()

__all__ = [
    "NCIPLOTJob",
    "NCIPLOTJobRunner",
    "jobs",
]
