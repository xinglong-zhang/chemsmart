import numpy as np

from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.grouper import (
    ConnectivityGrouper,
    FormulaGrouper,
    GeometryGrouper,
    HybridMoleculeGrouper,
    PymatgenMoleculeGrouper,
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

# ethanol
ethanol = Molecule.from_pubchem(identifier="CCO")

methanol_molecules = [methanol, methanol_rot1]
methanol_and_ethanol = [methanol, ethanol]


class TestGrouper:
    def test_geometry_grouper(self):
        # subclass RMSDGrouper
        assert np.any(
            methanol.positions != methanol_rot1.positions
        ), "Rotated molecule should have different positions."
        grouper = GeometryGrouper(methanol_molecules)
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

        grouper2 = GeometryGrouper(methanol_and_ethanol)
        groups, group_indices = grouper2.group()
        assert (
            len(groups) == 2
        ), "Molecules should form two groups based on geometry."
        assert (
            len(group_indices) == 2
        ), "Molecules should form two groups based on geometry."
        rmsd = grouper2._calculate_rmsd((0, 1))
        assert (
            rmsd is np.inf
        ), "RMSD is set to be infinity for different molecules."
        unique_structures = grouper2.unique()
        assert (
            len(unique_structures) == 2
        ), "Molecules should form two groups based on geometry."

    def test_formula_grouper(self):
        grouper = FormulaGrouper(methanol_molecules)
        groups, group_indices = grouper.group()
        unique_structures = grouper.unique()
        assert (
            len(groups) == 1
        ), "Molecules should form one group based on formula."

        assert (
            len(unique_structures) == 1
        ), "Molecules should form one group based on formula."

        grouper2 = FormulaGrouper(methanol_and_ethanol)
        groups, group_indices = grouper2.group()
        unique_structures = grouper2.unique()
        assert (
            len(groups) == 2
        ), "Molecules should form two groups based on formula."
        assert (
            len(unique_structures) == 2
        ), "Molecules should form two groups based on formula."

    def test_connectivity_grouper(self):
        grouper = ConnectivityGrouper(methanol_molecules)
        groups, group_indices = grouper.group()
        assert (
            len(groups) == 1
        ), "Molecules should form one group based on connectivity."
        assert (
            len(group_indices) == 1
        ), "Molecules should form one group based on connectivity."
        unique_structures = grouper.unique()
        assert (
            len(unique_structures) == 1
        ), "Molecules should form one group based on connectivity."

        grouper2 = ConnectivityGrouper(methanol_and_ethanol)
        groups, group_indices = grouper2.group()
        assert (
            len(groups) == 2
        ), "Molecules should form two groups based on connectivity."
        assert (
            len(group_indices) == 2
        ), "Molecules should form two groups based on connectivity."
        unique_structures = grouper2.unique()
        assert (
            len(unique_structures) == 2
        ), "Molecules should form two groups based on connectivity."

    def test_rcm_similarity_grouper(self):
        grouper = RCMSimilarityGrouper(methanol_molecules)
        groups, group_indices = grouper.group()
        assert (
            len(groups) == 1
        ), "Molecules should form one group based on RCM similarity."
        assert (
            len(group_indices) == 1
        ), "Molecules should form one group based on RCM similarity."
        unique_structures = grouper.unique()
        assert (
            len(unique_structures) == 1
        ), "Molecules should form one group based on RCM similarity."
        grouper2 = RCMSimilarityGrouper(methanol_and_ethanol)
        groups, group_indices = grouper2.group()
        assert (
            len(groups) == 2
        ), "Molecules should form two groups based on RCM similarity."
        assert (
            len(group_indices) == 2
        ), "Molecules should form two groups based on RCM similarity."
        unique_structures = grouper2.unique()
        assert (
            len(unique_structures) == 2
        ), "Molecules should form two groups based on RCM similarity."

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
