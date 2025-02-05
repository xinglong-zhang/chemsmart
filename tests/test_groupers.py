import numpy as np
import pytest

from chemsmart.utils.grouper import (
    ConnectivityGrouper,
    FormulaGrouper,
    RCMSimilarityGrouper,
    RDKitFingerprintGrouper,
    RDKitIsomorphismGrouper,
    RMSDGrouper,
    StructureGrouperFactory,
)


class TestGrouper:
    NUM_PROCS = 10

    def test_rmsd_grouper(self, methanol_molecules, methanol_and_ethanol):
        methanol = methanol_molecules[0]
        methanol_rot1 = methanol_molecules[1]
        assert np.any(
            methanol.positions != methanol_rot1.positions
        ), "Rotated molecule should have different positions."
        grouper = RMSDGrouper(methanol_molecules)
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

        grouper2 = RMSDGrouper(methanol_and_ethanol)
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

    def test_rmsd_grouper_for_large_number_of_mols(
        self, conformers_from_rdkit
    ):

        grouper = RMSDGrouper(conformers_from_rdkit, num_procs=self.NUM_PROCS)
        groups, group_indices = grouper.group()
        assert len(groups) == 296
        assert len(group_indices) == 296
        rmsd = grouper._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 1.72375)
        unique_structures = grouper.unique()
        assert len(unique_structures) == 296

        grouper2 = RMSDGrouper(
            conformers_from_rdkit, rmsd_threshold=1.0, num_procs=self.NUM_PROCS
        )
        # increased threshold, so should have less distinct groups
        groups, group_indices = grouper2.group()
        assert len(groups) == 160
        assert len(group_indices) == 160
        rmsd = grouper2._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 1.72375)
        unique_structures = grouper2.unique()
        assert len(unique_structures) == 160

        grouper3 = RMSDGrouper(
            conformers_from_rdkit, rmsd_threshold=1.5, num_procs=self.NUM_PROCS
        )
        # even greater threshold, so should have even less distinct groups
        groups, group_indices = grouper3.group()
        assert len(groups) == 10
        assert len(group_indices) == 10
        rmsd = grouper2._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 1.72375)
        unique_structures = grouper3.unique()
        assert len(unique_structures) == 10

        grouper4 = RMSDGrouper(
            conformers_from_rdkit, rmsd_threshold=2.0, num_procs=self.NUM_PROCS
        )
        # even greater threshold, so should have even less distinct groups
        groups, group_indices = grouper4.group()
        assert len(groups) == 3
        assert len(group_indices) == 3
        rmsd = grouper2._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 1.72375)
        unique_structures = grouper4.unique()
        assert len(unique_structures) == 3

    def test_formula_grouper(
        self, methanol_molecules, methanol_and_ethanol, conformers_from_rdkit
    ):
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

        grouper3 = FormulaGrouper(conformers_from_rdkit)
        # based on Formula, should all be the same even for 300 conformers
        groups, group_indices = grouper3.group()
        unique_structures = grouper3.unique()
        assert len(groups) == 1
        assert len(unique_structures) == 1

    @pytest.mark.slow
    def test_connectivity_grouper(
        self, methanol_molecules, methanol_and_ethanol
    ):
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

    @pytest.mark.slow
    def test_connectivity_grouper_for_large_number_of_mols(
        self, conformers_from_rdkit
    ):
        grouper3 = ConnectivityGrouper(
            conformers_from_rdkit, num_procs=self.NUM_PROCS
        )
        # based on connectivity, should all be the same even for 300 conformers
        groups, group_indices = grouper3.group()
        unique_structures = grouper3.unique()
        assert len(groups) == 1
        assert len(unique_structures) == 1

    def test_rcm_similarity_grouper(
        self, methanol_molecules, methanol_and_ethanol
    ):
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

    def test_rcm_similarity_grouper_for_large_number_of_mols(
        self, conformers_from_rdkit
    ):

        grouper = RCMSimilarityGrouper(
            conformers_from_rdkit, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 296
        assert len(group_indices) == 296
        rmsd = grouper._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 1.72375)
        unique_structures = grouper.unique()
        assert len(unique_structures) == 296

        grouper2 = RCMSimilarityGrouper(
            conformers_from_rdkit,
            similarity_threshold=1.0,
            num_procs=self.NUM_PROCS,
        )
        # increased threshold, so should have less distinct groups
        groups, group_indices = grouper2.group()
        assert len(groups) == 160
        assert len(group_indices) == 160
        rmsd = grouper2._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 1.72375)
        unique_structures = grouper2.unique()
        assert len(unique_structures) == 160

        grouper3 = RCMSimilarityGrouper(
            conformers_from_rdkit, rmsd_threshold=1.5, num_procs=10
        )
        # even greater threshold, so should have even less distinct groups
        groups, group_indices = grouper3.group()
        assert len(groups) == 10
        assert len(group_indices) == 10
        rmsd = grouper2._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 1.72375)
        unique_structures = grouper3.unique()
        assert len(unique_structures) == 10

        grouper4 = RCMSimilarityGrouper(
            conformers_from_rdkit, rmsd_threshold=2.0, num_procs=10
        )
        # even greater threshold, so should have even less distinct groups
        groups, group_indices = grouper4.group()
        assert len(groups) == 3
        assert len(group_indices) == 3
        rmsd = grouper2._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 1.72375)
        unique_structures = grouper4.unique()
        assert len(unique_structures) == 3

    def test_structure_grouper_factory(self):
        factory = StructureGrouperFactory()
        assert factory is not None

    def test_rdkit_isomorphism_grouper(self):
        grouper = RDKitIsomorphismGrouper()
        assert grouper is not None

    def test_rdkit_fingerprint_grouper(self):
        grouper = RDKitFingerprintGrouper()
        assert grouper is not None
