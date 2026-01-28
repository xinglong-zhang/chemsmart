"""
Grouper Job Management Module.

This module provides comprehensive grouper job classes and runners for
molecular structure clustering tasks including RMSD-based grouping,
fingerprint similarity, torsion fingerprint deviation, and more.

This follows the same pattern as the mol module.
"""

from .connectivity import ConnectivityGrouperJob
from .formula import FormulaGrouperJob
from .hrmsd import HRMSDGrouperJob
from .irmsd import IRMSDGrouperJob
from .isomorphism import IsomorphismGrouperJob

# Base job class
from .job import GrouperJob
from .pymolrmsd import PymolRMSDGrouperJob

# Strategy-specific job classes (entry points)
from .rmsd import RMSDGrouperJob

# Grouper algorithm classes (from runner.py)
from .runner import (  # Job Runner; Helper functions; Config; Base classes; RMSD-based groupers; Fingerprint groupers; Graph-based groupers; Factory (deprecated, kept for compatibility)
    BasicRMSDGrouper,
    ConnectivityGrouper,
    ConnectivityGrouperSharedMemory,
    FormulaGrouper,
    GrouperJobRunner,
    HungarianRMSDGrouper,
    IRMSDGrouper,
    MoleculeGrouper,
    PymolRMSDGrouper,
    RDKitIsomorphismGrouper,
    RMSDGrouper,
    RMSDGrouperSharedMemory,
    SpyRMSDGrouper,
    StructureGrouperConfig,
    StructureGrouperFactory,
    TanimotoSimilarityGrouper,
    TorsionFingerprintGrouper,
    to_graph_wrapper,
)
from .spyrmsd import SpyRMSDGrouperJob
from .tanimoto import TanimotoGrouperJob
from .tfd import TFDGrouperJob

__all__ = [
    # Job Runner
    "GrouperJobRunner",
    # Base job class
    "GrouperJob",
    # Strategy-specific job classes
    "RMSDGrouperJob",
    "HRMSDGrouperJob",
    "SpyRMSDGrouperJob",
    "IRMSDGrouperJob",
    "PymolRMSDGrouperJob",
    "TanimotoGrouperJob",
    "TFDGrouperJob",
    "ConnectivityGrouperJob",
    "FormulaGrouperJob",
    "IsomorphismGrouperJob",
    # Helper functions
    "to_graph_wrapper",
    # Config
    "StructureGrouperConfig",
    # Base classes
    "MoleculeGrouper",
    "RMSDGrouper",
    # RMSD-based groupers
    "BasicRMSDGrouper",
    "HungarianRMSDGrouper",
    "SpyRMSDGrouper",
    "IRMSDGrouper",
    "PymolRMSDGrouper",
    "RMSDGrouperSharedMemory",
    # Fingerprint groupers
    "TanimotoSimilarityGrouper",
    "TorsionFingerprintGrouper",
    # Graph-based groupers
    "RDKitIsomorphismGrouper",
    "FormulaGrouper",
    "ConnectivityGrouper",
    "ConnectivityGrouperSharedMemory",
    # Factory
    "StructureGrouperFactory",
]
