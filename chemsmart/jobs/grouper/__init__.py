"""
Grouper Job Management Module.

This module provides grouper job classes and runners for molecular structure
clustering tasks including RMSD-based grouping, fingerprint similarity,
torsion fingerprint deviation, and more.

Main classes:
- GrouperJob: Base job class for all grouping strategies
- GrouperJobRunner: Job runner for executing grouper jobs

For the factory interface to create groupers directly, use:
    from chemsmart.utils.grouper import StructureGrouperFactory
"""

# Base job class
from .job import GrouperJob

# Grouper algorithm classes (for advanced users who need direct access)
# Job runner
from .runner import (  # Base classes; RMSD-based groupers; Fingerprint groupers; Graph-based groupers; Helper
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
    TanimotoSimilarityGrouper,
    TorsionFingerprintGrouper,
    to_graph_wrapper,
)

__all__ = [
    # Main interface
    "GrouperJob",
    "GrouperJobRunner",
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
    # Helper
    "StructureGrouperConfig",
    "to_graph_wrapper",
]
