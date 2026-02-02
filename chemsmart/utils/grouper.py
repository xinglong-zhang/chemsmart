"""
Molecular structure grouping factory.

Provides a factory interface for creating grouper instances used by
crest and traj jobs. The actual grouper implementations are in
chemsmart.jobs.grouper.

Available strategies:
- rmsd: Basic RMSD grouping
- hrmsd: Hungarian RMSD grouping
- spyrmsd: spyrmsd-based RMSD grouping
- irmsd: Invariant RMSD grouping (symmetry-aware)
- pymolrmsd: PyMOL-based RMSD alignment
- tanimoto: Fingerprint-based Tanimoto similarity
- torsion: Torsion Fingerprint Deviation (TFD)
- isomorphism: RDKit graph isomorphism
- formula: Chemical formula grouping
- connectivity: Molecular connectivity grouping
- energy: Energy difference grouping
"""

import logging

from chemsmart.jobs.grouper.connectivity import ConnectivityGrouper
from chemsmart.jobs.grouper.energy import EnergyGrouper
from chemsmart.jobs.grouper.formula import FormulaGrouper
from chemsmart.jobs.grouper.isomorphism import RDKitIsomorphismGrouper
from chemsmart.jobs.grouper.rmsd import (
    BasicRMSDGrouper,
    HungarianRMSDGrouper,
    IRMSDGrouper,
    PymolRMSDGrouper,
    SpyRMSDGrouper,
)
from chemsmart.jobs.grouper.tanimoto import TanimotoSimilarityGrouper
from chemsmart.jobs.grouper.tfd import TorsionFingerprintGrouper

logger = logging.getLogger(__name__)

# Strategy name to grouper class mapping
GROUPER_CLASSES = {
    "rmsd": BasicRMSDGrouper,
    "hrmsd": HungarianRMSDGrouper,
    "spyrmsd": SpyRMSDGrouper,
    "irmsd": IRMSDGrouper,
    "pymolrmsd": PymolRMSDGrouper,
    "tanimoto": TanimotoSimilarityGrouper,
    "torsion": TorsionFingerprintGrouper,
    "isomorphism": RDKitIsomorphismGrouper,
    "formula": FormulaGrouper,
    "connectivity": ConnectivityGrouper,
    "energy": EnergyGrouper,
}

# Strategies that support threshold parameter (grouping threshold)
THRESHOLD_SUPPORTED = {
    "rmsd",
    "hrmsd",
    "spyrmsd",
    "irmsd",
    "pymolrmsd",
    "tanimoto",
    "torsion",
    "energy",
}

# Strategies that support ignore_hydrogens parameter
IGNORE_HYDROGENS_SUPPORTED = {
    "rmsd",
    "hrmsd",
    "spyrmsd",
    "irmsd",
    "pymolrmsd",
    "connectivity",
    "torsion",
}


class StructureGrouperFactory:
    """
    Factory for creating molecular grouper instances.

    Provides a unified entry point to construct groupers by strategy name.
    Used by GaussianCrestJob and GaussianTrajJob to group conformers.

    Supported strategies:
    - "rmsd": BasicRMSDGrouper
    - "hrmsd": HungarianRMSDGrouper
    - "spyrmsd": SpyRMSDGrouper
    - "irmsd": IRMSDGrouper (symmetry-aware)
    - "pymolrmsd": PymolRMSDGrouper
    - "tanimoto": TanimotoSimilarityGrouper
    - "torsion": TorsionFingerprintGrouper
    - "isomorphism": RDKitIsomorphismGrouper
    - "formula": FormulaGrouper
    - "connectivity": ConnectivityGrouper
    - "energy": EnergyGrouper (energy difference grouping)
    """

    @staticmethod
    def create(
        structures,
        strategy="rmsd",
        num_procs=1,
        threshold=None,
        num_groups=None,
        ignore_hydrogens=False,
        label=None,
        **kwargs,
    ):
        """
        Create a molecular grouper instance by strategy name.

        Args:
            structures: Iterable of Molecule objects to group.
            strategy (str): Grouping strategy name. Defaults to "rmsd".
            num_procs (int): Number of workers for parallel computation.
            threshold (float): Threshold for grouping (strategy dependent).
            num_groups (int): Number of groups to create (alternative to threshold).
            ignore_hydrogens (bool): Whether to ignore hydrogen atoms.
            label (str): Label/name for output files.
            **kwargs: Extra options forwarded to the grouper constructor.

        Returns:
            MoleculeGrouper: An instance of the selected grouper subclass.

        Raises:
            ValueError: If strategy is not a supported name.
        """
        strategy = strategy.lower() if strategy else "rmsd"

        if strategy not in GROUPER_CLASSES:
            raise ValueError(
                f"Unknown grouping strategy: '{strategy}'. "
                f"Supported strategies: {list(GROUPER_CLASSES.keys())}"
            )

        grouper_cls = GROUPER_CLASSES[strategy]
        logger.info(f"Using {strategy} grouping strategy.")

        # Build kwargs for the grouper
        grouper_kwargs = {
            "molecules": structures,
            "num_procs": num_procs,
            "label": label,
        }

        # Add threshold/num_groups for strategies that support it
        if strategy in THRESHOLD_SUPPORTED:
            grouper_kwargs["threshold"] = threshold
            if strategy in {
                "rmsd",
                "hrmsd",
                "spyrmsd",
                "irmsd",
                "pymolrmsd",
                "tanimoto",
                "torsion",
                "energy",
            }:
                grouper_kwargs["num_groups"] = num_groups

        # Add ignore_hydrogens for strategies that support it
        if strategy in IGNORE_HYDROGENS_SUPPORTED:
            grouper_kwargs["ignore_hydrogens"] = ignore_hydrogens

        # Add any additional kwargs
        grouper_kwargs.update(kwargs)

        return grouper_cls(**grouper_kwargs)
