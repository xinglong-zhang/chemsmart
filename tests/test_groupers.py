import os
import shutil
import tempfile

import numpy as np
import pytest

from chemsmart.io.xyz.xyzfile import XYZFile
from chemsmart.jobs.grouper.runner import (
    BasicRMSDGrouper,
    ConnectivityGrouper,
    FormulaGrouper,
    HungarianRMSDGrouper,
    IRMSDGrouper,
    PymolRMSDGrouper,
    RDKitIsomorphismGrouper,
    RMSDGrouper,
    SpyRMSDGrouper,
    TanimotoSimilarityGrouper,
    TorsionFingerprintGrouper,
)
from chemsmart.utils.grouper import StructureGrouperFactory
from chemsmart.utils.utils import kabsch_align


@pytest.fixture
def temp_working_dir():
    """
    Pytest fixture to create a temporary directory and change to it for testing.
    This prevents group_result folders from being created in the project directory.
    """
    original_dir = os.getcwd()
    temp_dir = tempfile.mkdtemp()

    try:
        os.chdir(temp_dir)
        yield temp_dir
    finally:
        os.chdir(original_dir)
        # Clean up the temporary directory
        shutil.rmtree(temp_dir, ignore_errors=True)


class Test_BasicRMSD_grouper_and_basic_functionality:
    NUM_PROCS = 4

    def test_rmsd_grouper(
        self, methanol_molecules, methanol_and_ethanol, temp_working_dir
    ):
        methanol = methanol_molecules[0]
        methanol_rot1 = methanol_molecules[1]
        assert np.any(
            methanol.positions != methanol_rot1.positions
        ), "Rotated molecule should have different positions."
        grouper = BasicRMSDGrouper(methanol_molecules)
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

        grouper2 = BasicRMSDGrouper(methanol_and_ethanol)
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
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)

        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = BasicRMSDGrouper(
            molecules, threshold=0.2, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 18
        assert len(group_indices) == 18
        unique_structures = grouper.unique()
        assert len(unique_structures) == 18

        # rmsd calculation from grouper
        rmsd = grouper._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 0.409, rtol=1e-3)

        _, _, _, _, rmsd = kabsch_align(
            molecules[0].positions, molecules[1].positions
        )
        assert np.isclose(rmsd, 0.409, rtol=1e-3)

        grouper2 = BasicRMSDGrouper(
            molecules, threshold=0.5, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper2.group()
        assert len(groups) == 12
        assert len(group_indices) == 12
        unique_structures = grouper2.unique()
        assert len(unique_structures) == 12

        grouper3 = BasicRMSDGrouper(
            molecules, threshold=1.0, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper3.group()
        assert len(groups) == 11
        assert len(group_indices) == 11
        unique_structures = grouper3.unique()
        assert len(unique_structures) == 11

        grouper4 = BasicRMSDGrouper(
            molecules, threshold=1.5, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper4.group()
        assert len(groups) == 9
        assert len(group_indices) == 9
        unique_structures = grouper4.unique()
        assert len(unique_structures) == 9

        grouper5 = BasicRMSDGrouper(
            molecules, threshold=2.0, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper5.group()
        assert len(groups) == 8
        assert len(group_indices) == 8
        unique_structures = grouper5.unique()
        assert len(unique_structures) == 8

        grouper6 = BasicRMSDGrouper(
            molecules, threshold=2.5, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper6.group()
        assert len(groups) == 4
        assert len(group_indices) == 4
        unique_structures = grouper6.unique()
        assert len(unique_structures) == 4

    def test_num_groups_parameter(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)

        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = BasicRMSDGrouper(
            molecules, threshold=None, num_groups=17, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 17
        assert len(group_indices) == 17
        unique_structures = grouper.unique()
        assert len(unique_structures) == 17

    def test_pick_the_lowestenergy_conformers(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)

        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = BasicRMSDGrouper(
            molecules, threshold=None, num_groups=3, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 3
        assert len(group_indices) == 3
        unique_structures = grouper.unique()
        assert len(unique_structures) == 3

        energies = [mol.energy for mol in unique_structures]
        assert -126.2575508 in energies
        assert -126.25153216 in energies
        assert -126.24909661 in energies

    def test_rmsd_grouper_for_crest_conformers_ignore_Hs(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)

        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = BasicRMSDGrouper(
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

        _, _, _, _, rmsd = kabsch_align(
            molecules[0].positions, molecules[1].positions
        )
        assert np.isclose(rmsd, 0.409, rtol=1e-3)  # did not remove H atoms

        grouper2 = BasicRMSDGrouper(
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

        grouper3 = BasicRMSDGrouper(
            molecules,
            threshold=1.0,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=True,
        )
        groups, group_indices = grouper3.group()
        assert len(groups) == 10
        assert len(group_indices) == 10
        unique_structures = grouper3.unique()
        assert len(unique_structures) == 10

        grouper4 = BasicRMSDGrouper(
            molecules,
            threshold=1.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=True,
        )
        groups, group_indices = grouper4.group()
        assert len(groups) == 7
        assert len(group_indices) == 7
        unique_structures = grouper4.unique()
        assert len(unique_structures) == 7

        grouper5 = BasicRMSDGrouper(
            molecules,
            threshold=2.0,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=True,
        )
        groups, group_indices = grouper5.group()
        assert len(groups) == 5
        assert len(group_indices) == 5
        unique_structures = grouper5.unique()
        assert len(unique_structures) == 5

        grouper6 = BasicRMSDGrouper(
            molecules,
            threshold=2.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=True,
        )
        groups, group_indices = grouper6.group()
        assert len(groups) == 5
        assert len(group_indices) == 5
        unique_structures = grouper6.unique()
        assert len(unique_structures) == 5

    def test_rmsd_grouper_for_rotated_molecules(
        self, two_rotated_molecules_xyz_file, temp_working_dir
    ):
        xyz_file = XYZFile(filename=two_rotated_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 2
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 2
        assert len(group_indices) == 2
        unique_structures = grouper.unique()
        assert len(unique_structures) == 2

        # rmsd calculation from grouper
        rmsd = grouper._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 0.611, rtol=1e-3)


class Test_HungarianRMSD_grouper:
    NUM_PROCS = 4

    def test_hrmsd_grouper_for_rotated_molecules(
        self, two_rotated_molecules_xyz_file, temp_working_dir
    ):
        xyz_file = XYZFile(filename=two_rotated_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 2
        grouper = HungarianRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 1
        assert len(group_indices) == 1
        unique_structures = grouper.unique()
        assert len(unique_structures) == 1

        # rmsd calculation from grouper
        rmsd = grouper._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 0.2294, rtol=1e-3)

    def test_hrmsd_grouper_for_crest_molecules(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = HungarianRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 12
        assert len(group_indices) == 12
        unique_structures = grouper.unique()
        assert len(unique_structures) == 12

        # rmsd calculation from grouper
        rmsd = grouper._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 0.4091, rtol=1e-3)
        rmsd = grouper._calculate_rmsd((0, 2))
        assert np.isclose(rmsd, 0.5899, rtol=1e-3)
        rmsd = grouper._calculate_rmsd((0, 3))
        assert np.isclose(rmsd, 1.8891, rtol=1e-3)


class Test_SpyRMSD_grouper:
    NUM_PROCS = 4

    def test_spyrmsd_grouper_for_rotated_molecules(
        self, two_rotated_molecules_xyz_file, temp_working_dir
    ):
        xyz_file = XYZFile(filename=two_rotated_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 2
        grouper = SpyRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 1
        assert len(group_indices) == 1
        unique_structures = grouper.unique()
        assert len(unique_structures) == 1

        # rmsd calculation from grouper
        rmsd = grouper._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 0.2125, rtol=1e-3)

    def test_spyrmsd_grouper_for_crest_molecules(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = SpyRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 12
        assert len(group_indices) == 12
        unique_structures = grouper.unique()
        assert len(unique_structures) == 12

        # rmsd calculation from grouper
        rmsd = grouper._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 0.4091, rtol=1e-3)
        rmsd = grouper._calculate_rmsd((0, 2))
        assert np.isclose(rmsd, 1.3925, rtol=1e-3)
        rmsd = grouper._calculate_rmsd((0, 3))
        assert np.isclose(rmsd, 2.1789, rtol=1e-3)
        rmsd = grouper._calculate_rmsd((0, 4))
        assert np.isclose(rmsd, 1.8202, rtol=1e-3)
        rmsd = grouper._calculate_rmsd((0, 5))
        assert np.isclose(rmsd, 2.0029, rtol=1e-3)


class Test_IRMSD_grouper:
    NUM_PROCS = 4

    def test_irmsd_grouper_for_rotated_molecules(
        self, two_rotated_molecules_xyz_file, temp_working_dir
    ):
        xyz_file = XYZFile(filename=two_rotated_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 2
        grouper = IRMSDGrouper(
            molecules,
            threshold=0.125,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 2
        assert len(group_indices) == 2
        unique_structures = grouper.unique()
        assert len(unique_structures) == 2

        # rmsd calculation from grouper
        rmsd = grouper._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 0.21249, rtol=1e-3)

    def test_irmsd_grouper_for_crest_molecules(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = IRMSDGrouper(
            molecules,
            threshold=0.125,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 18
        assert len(group_indices) == 18
        unique_structures = grouper.unique()
        assert len(unique_structures) == 18

        # rmsd calculation from grouper
        rmsd = grouper._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 0.4091, rtol=1e-3)
        rmsd = grouper._calculate_rmsd((0, 4))
        assert np.isclose(rmsd, 1.8202, rtol=1e-3)
        rmsd = grouper._calculate_rmsd((1, 2))
        assert np.isclose(rmsd, 1.4626, rtol=1e-3)
        rmsd = grouper._calculate_rmsd((11, 17))
        assert np.isclose(rmsd, 0.4497, rtol=1e-3)
        rmsd = grouper._calculate_rmsd((16, 17))
        assert np.isclose(rmsd, 3.1502, rtol=1e-3)


class Test_PymolRMSD_grouper:
    NUM_PROCS = 4

    @classmethod
    def setup_class(cls):
        """Initialize PyMOL once for all tests in this class."""
        try:
            import pymol
            from pymol import cmd

            pymol.finish_launching(["pymol", "-qc"])
            cmd.reinitialize()
        except Exception:
            pass

    def teardown_method(self, method):
        """Clean up PyMOL objects after each test method to prevent slowdown."""
        try:
            from pymol import cmd

            # Delete all objects instead of reinitialize to keep session alive
            cmd.delete("all")
        except Exception:
            pass

    def test_pymolrmsd_grouper_for_rotated_molecules(
        self, two_rotated_molecules_xyz_file, temp_working_dir
    ):
        xyz_file = XYZFile(filename=two_rotated_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 2
        grouper = PymolRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 1
        assert len(group_indices) == 1
        unique_structures = grouper.unique()
        assert len(unique_structures) == 1

        # rmsd calculation from grouper
        rmsd = grouper._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 0.0000, rtol=1e-3)

        # Explicitly cleanup to prevent __del__ from calling quit()
        if hasattr(grouper, "_temp_dir") and grouper._temp_dir:
            import shutil

            shutil.rmtree(grouper._temp_dir, ignore_errors=True)
            grouper._temp_dir = None
        grouper.cmd = None  # Prevent __del__ from calling quit()

    def test_pymolrmsd_grouper_for_crest_molecules(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = PymolRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 7
        assert len(group_indices) == 7
        unique_structures = grouper.unique()
        assert len(unique_structures) == 7

        # rmsd calculation from grouper
        rmsd = grouper._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 0.074175, rtol=1e-3)
        rmsd = grouper._calculate_rmsd((0, 2))
        assert np.isclose(rmsd, 0.023745, rtol=1e-3)
        rmsd = grouper._calculate_rmsd((3, 4))
        assert np.isclose(rmsd, 1.725241, rtol=1e-3)
        rmsd = grouper._calculate_rmsd((7, 8))
        assert np.isclose(rmsd, 2.309451, rtol=1e-3)
        rmsd = grouper._calculate_rmsd((14, 16))
        assert np.isclose(rmsd, 1.837319, rtol=1e-3)

        # Explicitly cleanup to prevent __del__ from calling quit()
        if hasattr(grouper, "_temp_dir") and grouper._temp_dir:
            import shutil

            shutil.rmtree(grouper._temp_dir, ignore_errors=True)
            grouper._temp_dir = None
        grouper.cmd = None  # Prevent __del__ from calling quit()

    @classmethod
    def teardown_class(cls):
        """Clean up PyMOL after all tests in this class are done."""
        try:
            from pymol import cmd

            cmd.reinitialize()
        except Exception:
            pass


class Test_Tanimoto_similarity_grouper:
    NUM_PROCS = 4

    def test_tanimoto_similarity_grouper(
        self, methanol_molecules, methanol_and_ethanol, temp_working_dir
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
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)

        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = TanimotoSimilarityGrouper(
            molecules,
            threshold=0.98,
            fingerprint_type="usrcat",
            num_procs=self.NUM_PROCS,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 4
        assert len(group_indices) == 4
        unique_structures = grouper.unique()
        assert len(unique_structures) == 4

        grouper = TanimotoSimilarityGrouper(
            molecules,
            threshold=0.999,
            fingerprint_type="usrcat",
            num_procs=self.NUM_PROCS,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 14
        assert len(group_indices) == 14
        unique_structures = grouper.unique()
        assert len(unique_structures) == 14


class Test_TorsionFingerprint_grouper:
    NUM_PROCS = 4

    def test_torsionfingerprint_grouper_for_rotated_molecules(
        self, two_rotated_molecules_xyz_file, temp_working_dir
    ):
        xyz_file = XYZFile(filename=two_rotated_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 2
        grouper = TorsionFingerprintGrouper(
            molecules,
            threshold=0.1,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 1
        assert len(group_indices) == 1
        unique_structures = grouper.unique()
        assert len(unique_structures) == 1

    def test_torsionfingerprint_grouper_for_crest_molecules(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = TorsionFingerprintGrouper(
            molecules,
            threshold=0.05,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 10
        assert len(group_indices) == 10
        unique_structures = grouper.unique()
        assert len(unique_structures) == 10

    def test_use_weights_parameter(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = TorsionFingerprintGrouper(
            molecules,
            threshold=0.05,
            num_procs=self.NUM_PROCS,
            use_weights=False,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 15
        assert len(group_indices) == 15
        unique_structures = grouper.unique()
        assert len(unique_structures) == 15


class Test_other_groupers:
    NUM_PROCS = 4

    def test_formula_grouper(
        self,
        methanol_molecules,
        methanol_and_ethanol,
        conformers_from_rdkit,
        temp_working_dir,
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
        self, methanol_molecules, methanol_and_ethanol, temp_working_dir
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
        self, multiple_molecules_xyz_file, temp_working_dir
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

    def test_rdkit_isomorphism_grouper(
        self, methanol_molecules, methanol_and_ethanol, temp_working_dir
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


class Testfactory:

    def test_structure_grouper_factory(
        self, methanol_molecules, temp_working_dir
    ):
        factory = StructureGrouperFactory()
        rmsd_grouper = factory.create(methanol_molecules, strategy="rmsd")
        assert isinstance(rmsd_grouper, RMSDGrouper)
        hrmsd_grouper = factory.create(methanol_molecules, strategy="hrmsd")
        assert isinstance(hrmsd_grouper, RMSDGrouper)
        spyrmsd_grouper = factory.create(
            methanol_molecules, strategy="spyrmsd"
        )
        assert isinstance(spyrmsd_grouper, RMSDGrouper)
        irmsd_grouper = factory.create(methanol_molecules, strategy="irmsd")
        assert isinstance(irmsd_grouper, RMSDGrouper)
        pymolrmsd_grouper = factory.create(
            methanol_molecules, strategy="pymolrmsd"
        )
        assert isinstance(pymolrmsd_grouper, RMSDGrouper)

        # Explicitly cleanup pymolrmsd_grouper to prevent __del__ from calling quit()
        if (
            hasattr(pymolrmsd_grouper, "_temp_dir")
            and pymolrmsd_grouper._temp_dir
        ):
            import shutil

            shutil.rmtree(pymolrmsd_grouper._temp_dir, ignore_errors=True)
            pymolrmsd_grouper._temp_dir = None
        pymolrmsd_grouper.cmd = None  # Prevent __del__ from calling quit()

        torsion_grouper = factory.create(
            methanol_molecules, strategy="torsion"
        )
        assert isinstance(torsion_grouper, TorsionFingerprintGrouper)
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

    @classmethod
    def teardown_class(cls):
        """Clean up PyMOL after all tests in this class are done."""
        try:
            from pymol import cmd

            cmd.reinitialize()
        except Exception:
            pass
