import numpy as np

from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.grouper import (
    ConnectivityGrouper,
    FormulaGrouper,
    GeometryGrouper,
    HybridMoleculeGrouper,
    PymatgenMoleculeGrouper,
    RCMAdjacencyGrouper,
    RCMSimilarityGrouper,
    RDKitFingerprintGrouper,
    RDKitIsomorphismGrouper,
    RMSDGrouper,
    StructureGrouperFactory,
)

# molecules for testing
# methanol
methanol = Molecule.from_pubchem(identifier="CO")
# f = open("methanol.xyz", "w")
# methanol.write_coordinates(f)


# rotated methanol
ase_atoms = methanol.to_ase()
ase_atoms.rotate(90, [0, 0, 1])
methanol_rot1 = Molecule.from_ase_atoms(ase_atoms)
# g = open("methanol_rot1.xyz", "w")
# methanol_rot1.write_coordinates(g)

# ethanol
ethanol = Molecule.from_pubchem(identifier="CCO")


class TestGrouper:
    def test_geometry_grouper(self):
        # subclass RMSDGrouper
        assert np.any(
            methanol.positions != methanol_rot1.positions
        ), "Rotated molecule should have different positions."
        molecules = list([methanol, methanol_rot1])
        grouper = GeometryGrouper(molecules)
        groups, group_indices = grouper.group()
        assert (
            len(groups) == 1
        ), "Molecules should form one group based on geometry."
        assert (
            len(group_indices) == 1
        ), "Molecules should form one group based on geometry."
        rmsd = grouper._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 0.0), "RMSD should be close to zero."
        unique_structures = grouper.unique()
        assert (
            len(unique_structures) == 1
        ), "Molecules should form one group based on geometry."

    def test_formula_grouper(self):
        molecules = list([methanol, ethanol])
        grouper = FormulaGrouper(molecules)
        groups = grouper.group()
        unique_structures = grouper.unique()
        assert (
            len(groups) == 2
        ), "Molecules should form two groups based on formula."
        assert (
            len(unique_structures) == 2
        ), "Molecules should form two groups based on formula."

    def test_connectivity_grouper(self):
        grouper = ConnectivityGrouper()
        assert grouper is not None

    def test_rcm_adjacency_grouper(self):
        grouper = RCMAdjacencyGrouper()
        assert grouper is not None

    def test_rcm_similarity_grouper(self):
        grouper = RCMSimilarityGrouper()
        assert grouper is not None

    def test_hybrid_molecule_grouper(self):
        grouper = HybridMoleculeGrouper()
        assert grouper is not None

    def test_pymatgen_molecule_grouper(self):
        grouper = PymatgenMoleculeGrouper()
        assert grouper is not None

    def test_rmsd_molecule_grouper(self):
        grouper = RMSDGrouper()
        assert grouper is not None

    def test_structure_grouper_factory(self):
        factory = StructureGrouperFactory()
        assert factory is not None

    def test_rdkit_isomorphism_grouper(self):
        grouper = RDKitIsomorphismGrouper()
        assert grouper is not None

    def test_rdkit_fingerprint_grouper(self):
        grouper = RDKitFingerprintGrouper()
        assert grouper is not None
