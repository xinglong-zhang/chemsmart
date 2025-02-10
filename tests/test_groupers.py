import numpy as np
import pytest

from chemsmart.io.molecules.structure import XYZFile
from chemsmart.utils.grouper import (
    ConnectivityGrouper,
    ConnectivityGrouperSharedMemory,
    FormulaGrouper,
    RDKitIsomorphismGrouper,
    RMSDGrouper,
    RMSDGrouperSharedMemory,
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

        molecules = xyz_file.get_molecule(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = RMSDGrouper(
            molecules, num_procs=self.NUM_PROCS, rmsd_threshold=0.2
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 18
        assert len(group_indices) == 18
        unique_structures = grouper.unique()
        assert len(unique_structures) == 18

        # rmsd calculation from grouper
        rmsd = grouper._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 0.409, rtol=1e-3)

        # rmsd calculation from Kabsh alignment
        from chemsmart.utils.utils import kabsch_align

        _, _, _, _, rmsd = kabsch_align(
            molecules[0].positions, molecules[1].positions
        )
        assert np.isclose(rmsd, 0.409, rtol=1e-3)

        grouper2 = RMSDGrouper(
            molecules, num_procs=self.NUM_PROCS, rmsd_threshold=0.5
        )
        groups, group_indices = grouper2.group()
        assert len(groups) == 12
        assert len(group_indices) == 12
        unique_structures = grouper2.unique()
        assert len(unique_structures) == 12

        grouper3 = RMSDGrouper(
            molecules, num_procs=self.NUM_PROCS, rmsd_threshold=1.0
        )
        groups, group_indices = grouper3.group()
        assert len(groups) == 10
        assert len(group_indices) == 10
        unique_structures = grouper3.unique()
        assert len(unique_structures) == 10

        grouper4 = RMSDGrouper(
            molecules, num_procs=self.NUM_PROCS, rmsd_threshold=1.5
        )
        groups, group_indices = grouper4.group()
        assert len(groups) == 6
        assert len(group_indices) == 6
        unique_structures = grouper4.unique()
        assert len(unique_structures) == 6

        grouper5 = RMSDGrouper(
            molecules, num_procs=self.NUM_PROCS, rmsd_threshold=2.0
        )
        groups, group_indices = grouper5.group()
        assert len(groups) == 4
        assert len(group_indices) == 4
        unique_structures = grouper5.unique()
        assert len(unique_structures) == 4

        grouper6 = RMSDGrouper(
            molecules, num_procs=self.NUM_PROCS, rmsd_threshold=2.5
        )
        groups, group_indices = grouper6.group()
        assert len(groups) == 3
        assert len(group_indices) == 3
        unique_structures = grouper6.unique()
        assert len(unique_structures) == 3

    @pytest.mark.slow
    def test_rmsd_grouper_for_large_number_of_mols(
        self, conformers_from_rdkit
    ):
        # conformers_from_rdkit = conformers_from_rdkit[:12]

        grouper = RMSDGrouper(conformers_from_rdkit, num_procs=self.NUM_PROCS)
        groups, group_indices = grouper.group()
        # groupers seems sensitive to numpy version
        assert len(groups) == 297
        assert len(group_indices) == 297
        mol1 = conformers_from_rdkit[0]

        mol2 = conformers_from_rdkit[1]

        ## from pymol alignment
        #  ExecutiveRMS: 4 atoms rejected during cycle 1 (RMSD=1.67).
        #  ExecutiveRMS: 2 atoms rejected during cycle 2 (RMSD=1.32).
        #  ExecutiveRMS: 4 atoms rejected during cycle 3 (RMSD=1.20).
        #  ExecutiveRMS: 2 atoms rejected during cycle 4 (RMSD=1.00).
        #  ExecutiveRMS: 1 atoms rejected during cycle 5 (RMSD=0.90).
        #  Executive: RMSD =    0.860 (40 to 40 atoms)

        # rmsd calculation from Kabsh alignment
        from chemsmart.utils.utils import kabsch_align

        _, _, _, _, rmsd = kabsch_align(mol1.positions, mol2.positions)
        assert np.isclose(rmsd, 1.670, rtol=1e-3)

        # rmsd calculation from grouper
        rmsd = grouper._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 1.670, rtol=1e-3)
        unique_structures = grouper.unique()
        assert len(unique_structures) == 297

        grouper2 = RMSDGrouper(
            conformers_from_rdkit, rmsd_threshold=1.0, num_procs=self.NUM_PROCS
        )
        # increased threshold, so should have less distinct groups
        groups, group_indices = grouper2.group()
        assert len(groups) == 163
        assert len(group_indices) == 163
        rmsd = grouper2._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 1.670, rtol=1e-3)
        unique_structures = grouper2.unique()
        assert len(unique_structures) == 163

        grouper3 = RMSDGrouper(
            conformers_from_rdkit, rmsd_threshold=1.5, num_procs=self.NUM_PROCS
        )
        # even greater threshold, so should have even less distinct groups
        groups, group_indices = grouper3.group()
        assert len(groups) == 6
        assert len(group_indices) == 6
        rmsd = grouper2._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 1.670, rtol=1e-3)
        unique_structures = grouper3.unique()
        assert len(unique_structures) == 6

        grouper4 = RMSDGrouper(
            conformers_from_rdkit, rmsd_threshold=2.0, num_procs=self.NUM_PROCS
        )
        # even greater threshold, so should have even less distinct groups
        groups, group_indices = grouper4.group()
        assert len(groups) == 4
        assert len(group_indices) == 4
        rmsd = grouper2._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 1.670, rtol=1e-3)
        unique_structures = grouper4.unique()
        assert len(unique_structures) == 4
        # took 18 sec 848 ms

    @pytest.mark.slow
    def test_rmsd_grouper_shared_memory_for_large_number_of_mols(
        self, conformers_from_rdkit
    ):

        grouper = RMSDGrouperSharedMemory(
            conformers_from_rdkit, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 297
        assert len(group_indices) == 297
        unique_structures = grouper.unique()
        assert len(unique_structures) == 297

        grouper2 = RMSDGrouperSharedMemory(
            conformers_from_rdkit, rmsd_threshold=1.0, num_procs=self.NUM_PROCS
        )
        # increased threshold, so should have less distinct groups
        groups, group_indices = grouper2.group()
        assert len(groups) == 163
        assert len(group_indices) == 163
        unique_structures = grouper2.unique()
        assert len(unique_structures) == 163

        grouper3 = RMSDGrouperSharedMemory(
            conformers_from_rdkit, rmsd_threshold=1.5, num_procs=self.NUM_PROCS
        )
        # even greater threshold, so should have even less distinct groups
        groups, group_indices = grouper3.group()
        assert len(groups) == 6
        assert len(group_indices) == 6
        unique_structures = grouper3.unique()
        assert len(unique_structures) == 6

        grouper4 = RMSDGrouperSharedMemory(
            conformers_from_rdkit, rmsd_threshold=2.0, num_procs=self.NUM_PROCS
        )
        # even greater threshold, so should have even less distinct groups
        groups, group_indices = grouper4.group()
        assert len(groups) == 4
        assert len(group_indices) == 4
        unique_structures = grouper4.unique()
        assert len(unique_structures) == 4
        # took 19 sec 247 ms

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

        molecules = xyz_file.get_molecule(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = ConnectivityGrouper(
            molecules, num_procs=self.NUM_PROCS, bond_cutoff_buffer=0.2
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 3
        assert len(group_indices) == 3
        unique_structures = grouper.unique()
        assert len(unique_structures) == 3

        grouper2 = ConnectivityGrouper(
            molecules, num_procs=self.NUM_PROCS, bond_cutoff_buffer=0.5
        )
        groups, group_indices = grouper2.group()
        assert len(groups) == 1
        assert len(group_indices) == 1
        unique_structures = grouper2.unique()
        assert len(unique_structures) == 1

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
        assert len(groups) == 209
        assert len(unique_structures) == 209

    @pytest.mark.slow
    def test_connectivity_grouper_shared_for_large_number_of_mols(
        self, conformers_from_rdkit
    ):
        grouper3 = ConnectivityGrouperSharedMemory(
            conformers_from_rdkit, num_procs=self.NUM_PROCS
        )
        # based on connectivity, should all be the same even for 300 conformers
        groups, group_indices = grouper3.group()
        unique_structures = grouper3.unique()
        assert len(groups) == 209
        assert len(unique_structures) == 209
        # uses 1 min 50 sec; most calls to time.sleep()

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

    @pytest.mark.slow
    def test_tanimoto_similarity_grouper_for_large_number_of_mols(
        self, conformers_from_rdkit
    ):

        grouper = TanimotoSimilarityGrouper(
            conformers_from_rdkit, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 2
        assert len(group_indices) == 2
        unique_structures = grouper.unique()
        assert len(unique_structures) == 2

        grouper2 = TanimotoSimilarityGrouper(
            conformers_from_rdkit,
            similarity_threshold=0.95,  # very strict threshold
            num_procs=self.NUM_PROCS,
        )
        # increased threshold, so should have less distinct groups
        groups, group_indices = grouper2.group()
        assert len(groups) == 92
        assert len(group_indices) == 92
        unique_structures = grouper2.unique()
        assert len(unique_structures) == 92

        grouper3 = TanimotoSimilarityGrouper(
            conformers_from_rdkit, similarity_threshold=0.8, num_procs=10
        )
        # even greater threshold, so should have even less distinct groups
        groups, group_indices = grouper3.group()
        assert len(groups) == 1
        assert len(group_indices) == 1
        unique_structures = grouper3.unique()
        assert len(unique_structures) == 1

        grouper4 = TanimotoSimilarityGrouper(
            conformers_from_rdkit, similarity_threshold=0.5, num_procs=10
        )
        # even greater threshold, so should have even less distinct groups
        groups, group_indices = grouper4.group()
        assert len(groups) == 1
        assert len(group_indices) == 1
        unique_structures = grouper4.unique()
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

    @pytest.mark.slow
    def test_rdkit_isomorphism_grouper_for_large_number_of_mols(
        self, conformers_from_rdkit
    ):
        conformers_from_rdkit = conformers_from_rdkit[:5]
        # seems very slow for 300 conformers, even for 20 confs

        grouper = RDKitIsomorphismGrouper(
            conformers_from_rdkit, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 1
        assert len(group_indices) == 1
        unique_structures = grouper.unique()
        assert len(unique_structures) == 1

        grouper2 = RDKitIsomorphismGrouper(
            conformers_from_rdkit,
            use_stereochemistry=False,
            num_procs=self.NUM_PROCS,
        )
        # increased threshold, so should have less distinct groups
        groups, group_indices = grouper2.group()
        assert len(groups) == 1
        assert len(group_indices) == 1
        unique_structures = grouper2.unique()
        assert len(unique_structures) == 1

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
