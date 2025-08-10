import os

import numpy as np
import pytest

from chemsmart.io.xyz.file import XYZFile
from chemsmart.utils.grouper import (
    ConnectivityGrouper,
    FormulaGrouper,
    RDKitIsomorphismGrouper,
    RMSDGrouper,
    StructureGrouperFactory,
    TanimotoSimilarityGrouper,
)
from chemsmart.utils.utils import kabsch_align


class TestGrouper:
    NUM_PROCS = 4

    def test_rmsd_grouper(self, methanol_molecules, methanol_and_ethanol):
        methanol = methanol_molecules[0]
        methanol_rot1 = methanol_molecules[1]
        assert np.any(
            methanol.positions != methanol_rot1.positions
        ), "Rotated molecule should have different positions."
        grouper = RMSDGrouper(methanol_molecules)
        groups, group_indices = grouper.group()
        rmsd = grouper._calculate_rmsd((0, 1))
        assert np.isclose(
            rmsd, 0.0, rtol=1e-3
        ), "RMSD should be close to zero."
        rmsd = grouper._calculate_rmsd((1, 2))
        assert np.isclose(
            rmsd, 0.0, rtol=1e-3
        ), "RMSD should be close to zero."
        _, _, _, _, rmsd_kabsch = kabsch_align(
            methanol.positions, methanol_rot1.positions
        )
        assert np.isclose(rmsd_kabsch, 0.0, rtol=1e-3)
        assert (
            len(groups) == 1
        ), "Molecules should form one group based on geometry."
        assert (
            len(group_indices) == 1
        ), "Molecules should form one group based on geometry."
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

    def test_rmsd_grouper_for_crest_conformers(
        self, multiple_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)

        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = RMSDGrouper(
            molecules, threshold=0.2, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 17
        assert len(group_indices) == 17
        unique_structures = grouper.unique()
        assert len(unique_structures) == 17

        # rmsd calculation from grouper
        rmsd = grouper._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 0.301, rtol=1e-3)

        # rmsd calculation from Kabsh alignment
        from chemsmart.utils.utils import kabsch_align

        _, _, _, _, rmsd = kabsch_align(
            molecules[0].positions, molecules[1].positions
        )
        assert np.isclose(rmsd, 0.409, rtol=1e-3)

        grouper2 = RMSDGrouper(
            molecules, threshold=0.5, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper2.group()
        assert len(groups) == 12
        assert len(group_indices) == 12
        unique_structures = grouper2.unique()
        assert len(unique_structures) == 12

        grouper3 = RMSDGrouper(
            molecules, threshold=1.0, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper3.group()
        assert len(groups) == 7
        assert len(group_indices) == 7
        unique_structures = grouper3.unique()
        assert len(unique_structures) == 7

        grouper4 = RMSDGrouper(
            molecules, threshold=2.0, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper4.group()
        assert len(groups) == 1
        assert len(group_indices) == 1
        unique_structures = grouper4.unique()
        assert len(unique_structures) == 1

    def test_rmsd_grouper_for_crest_conformers_ignore_Hs(
        self, multiple_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)

        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = RMSDGrouper(
            molecules,
            threshold=0.2,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=True,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 17
        assert len(group_indices) == 17
        unique_structures = grouper.unique()
        assert len(unique_structures) == 17

        # rmsd calculation from grouper
        rmsd = grouper._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 0.301, rtol=1e-3)  # removed H atoms

        # rmsd calculation from Kabsh alignment
        from chemsmart.utils.utils import kabsch_align

        _, _, _, _, rmsd = kabsch_align(
            molecules[0].positions, molecules[1].positions
        )
        assert np.isclose(rmsd, 0.409, rtol=1e-3)  # did not remove H atoms

        grouper2 = RMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=True,
        )
        groups, group_indices = grouper2.group()
        assert len(groups) == 12
        assert len(group_indices) == 12
        unique_structures = grouper2.unique()
        assert len(unique_structures) == 12

        grouper3 = RMSDGrouper(
            molecules,
            threshold=1.0,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=True,
        )
        groups, group_indices = grouper3.group()
        assert len(groups) == 7
        assert len(group_indices) == 7
        unique_structures = grouper3.unique()
        assert len(unique_structures) == 7

        grouper4 = RMSDGrouper(
            molecules,
            threshold=1.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=True,
        )
        groups, group_indices = grouper4.group()
        assert len(groups) == 1
        assert len(group_indices) == 1
        unique_structures = grouper4.unique()
        assert len(unique_structures) == 1

    def test_rmsd_grouper_for_Ticatalysis_conformers(
        self, conformers_test_directory
    ):
        xyz_files = sorted(
            [
                os.path.join(conformers_test_directory, f)
                for f in os.listdir(conformers_test_directory)
                if f.endswith(".xyz")
            ]
        )
        molecules = []
        for path in xyz_files:
            mol = XYZFile(filename=path).get_molecules(return_list=False)
            molecules.append(mol)

        assert len(molecules) == 10

        grouper1 = RMSDGrouper(
            molecules,
            threshold=0.2,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper1.group()

        assert len(groups) == 4
        assert len(group_indices) == 4
        unique_structures = grouper1.unique()
        assert len(unique_structures) == 4

        grouper2 = RMSDGrouper(
            molecules,
            threshold=1,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper2.group()

        assert len(groups) == 1
        assert len(group_indices) == 1
        unique_structures = grouper2.unique()
        assert len(unique_structures) == 1

        grouper3 = RMSDGrouper(
            molecules,
            threshold=0.1,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=True,
        )
        groups, group_indices = grouper3.group()

        assert len(groups) == 4
        assert len(group_indices) == 4
        unique_structures = grouper3.unique()
        assert len(unique_structures) == 4

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

    def test_connectivity_grouper_for_crest_conformers(
        self, multiple_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)

        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = ConnectivityGrouper(
            molecules, num_procs=self.NUM_PROCS, threshold=0.2
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 3
        assert len(group_indices) == 3
        unique_structures = grouper.unique()
        assert len(unique_structures) == 3

        grouper2 = ConnectivityGrouper(
            molecules, num_procs=self.NUM_PROCS, threshold=0.5
        )
        groups, group_indices = grouper2.group()
        assert len(groups) == 1
        assert len(group_indices) == 1
        unique_structures = grouper2.unique()
        assert len(unique_structures) == 1

    def test_tanimoto_similarity_grouper(
        self, methanol_molecules, methanol_and_ethanol
    ):
        grouper = TanimotoSimilarityGrouper(methanol_molecules)
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
        grouper2 = TanimotoSimilarityGrouper(methanol_and_ethanol)
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

    def test_tanimoto_grouper_for_crest_conformers(
        self, multiple_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)

        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = TanimotoSimilarityGrouper(
            molecules, threshold=0.95, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 1
        assert len(group_indices) == 1
        unique_structures = grouper.unique()
        assert len(unique_structures) == 1

        grouper2 = TanimotoSimilarityGrouper(
            molecules, threshold=0.8, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper2.group()
        assert len(groups) == 1
        assert len(group_indices) == 1
        unique_structures = grouper2.unique()
        assert len(unique_structures) == 1

        grouper3 = TanimotoSimilarityGrouper(
            molecules, threshold=0.5, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper3.group()
        assert len(groups) == 1
        assert len(group_indices) == 1
        unique_structures = grouper3.unique()
        assert len(unique_structures) == 1

    def test_rdkit_isomorphism_grouper(
        self, methanol_molecules, methanol_and_ethanol
    ):
        grouper = RDKitIsomorphismGrouper(methanol_molecules)
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
        grouper2 = RDKitIsomorphismGrouper(methanol_and_ethanol)
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

    def test_structure_grouper_factory(self, methanol_molecules):
        factory = StructureGrouperFactory()
        rmsd_grouper = factory.create(methanol_molecules, strategy="rmsd")
        assert isinstance(rmsd_grouper, RMSDGrouper)
        formula_grouper = factory.create(
            methanol_molecules, strategy="formula"
        )
        assert isinstance(formula_grouper, FormulaGrouper)
        connectivity_grouper = factory.create(
            methanol_molecules, strategy="connectivity"
        )
        assert isinstance(connectivity_grouper, ConnectivityGrouper)
        tanimoto_grouper = factory.create(
            methanol_molecules, strategy="tanimoto"
        )
        assert isinstance(tanimoto_grouper, TanimotoSimilarityGrouper)
        rdkit_isomorphism_grouper = factory.create(
            methanol_molecules, strategy="isomorphism"
        )
        assert isinstance(rdkit_isomorphism_grouper, RDKitIsomorphismGrouper)
