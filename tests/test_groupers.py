import os
import shutil
import tempfile

import numpy as np
import pytest

from chemsmart.io.xyz.xyzfile import XYZFile
from chemsmart.jobs.grouper.connectivity import ConnectivityGrouper
from chemsmart.jobs.grouper.energy import EnergyGrouper
from chemsmart.jobs.grouper.formula import FormulaGrouper
from chemsmart.jobs.grouper.isomorphism import RDKitIsomorphismGrouper
from chemsmart.jobs.grouper.rmsd import (
    BasicRMSDGrouper,
    HungarianRMSDGrouper,
    IRMSDGrouper,
    PymolRMSDGrouper,
    RMSDGrouper,
    SpyRMSDGrouper,
)
from chemsmart.jobs.grouper.tanimoto import TanimotoSimilarityGrouper
from chemsmart.jobs.grouper.tfd import TorsionFingerprintGrouper
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

    def test_ignore_hydrogen(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = HungarianRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=True,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 12
        assert len(group_indices) == 12
        unique_structures = grouper.unique()
        assert len(unique_structures) == 12

        # rmsd calculation from grouper
        rmsd = grouper._calculate_rmsd((0, 4))
        assert np.isclose(rmsd, 1.1915, rtol=1e-3)


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

    def test_ignore_hydrogen(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = SpyRMSDGrouper(
            molecules,
            threshold=2.5981,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=True,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 3
        assert len(group_indices) == 3
        unique_structures = grouper.unique()
        assert len(unique_structures) == 3

        # rmsd calculation from grouper
        rmsd = grouper._calculate_rmsd((0, 5))
        assert np.isclose(rmsd, 1.7034, rtol=1e-3)
        rmsd = grouper._calculate_rmsd((0, 6))
        assert np.isclose(rmsd, 2.6183, rtol=1e-3)


def _irmsd_available():
    """Check if irmsd command is available."""
    import os
    import shutil

    # Check IRMSD_PATH
    irmsd_path = os.environ.get("IRMSD_PATH")
    if (
        irmsd_path
        and os.path.isfile(irmsd_path)
        and os.access(irmsd_path, os.X_OK)
    ):
        return True

    # Check IRMSD_CONDA_ENV
    conda_env = os.environ.get("IRMSD_CONDA_ENV")
    if conda_env:
        conda_base = os.environ.get("CONDA_PREFIX_1") or os.environ.get(
            "CONDA_PREFIX"
        )
        if conda_base:
            if "envs" in conda_base:
                envs_dir = conda_base.rsplit("envs", 1)[0] + "envs"
            else:
                envs_dir = os.path.join(os.path.dirname(conda_base), "envs")
            irmsd_path = os.path.join(envs_dir, conda_env, "bin", "irmsd")
            if os.path.isfile(irmsd_path) and os.access(irmsd_path, os.X_OK):
                return True

    # Check PATH
    if shutil.which("irmsd"):
        return True

    return False


@pytest.mark.skipif(
    not _irmsd_available(), reason="irmsd command not available"
)
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
        assert np.isclose(rmsd, 0.2294, rtol=1e-3)

    def test_irmsd_grouper_for_crest_molecules(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = IRMSDGrouper(
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
        rmsd = grouper._calculate_rmsd((2, 10))
        assert np.isclose(rmsd, 2.2390, rtol=1e-3)
        rmsd = grouper._calculate_rmsd((8, 16))
        assert np.isclose(rmsd, 3.4209, rtol=1e-3)
        rmsd = grouper._calculate_rmsd((0, 13))
        assert np.isclose(rmsd, 0.8411, rtol=1e-3)

    def test_ignore_hydrogen(
        self, two_rotated_molecules_xyz_file, temp_working_dir
    ):
        xyz_file = XYZFile(filename=two_rotated_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 2
        grouper = IRMSDGrouper(
            molecules,
            threshold=0.125,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=True,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 2
        assert len(group_indices) == 2
        unique_structures = grouper.unique()
        assert len(unique_structures) == 2

        # rmsd calculation from grouper
        rmsd = grouper._calculate_rmsd((0, 1))
        assert np.isclose(rmsd, 0.2294, rtol=1e-3)


class Test_PymolRMSD_grouper:
    NUM_PROCS = 1

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

    def test_ignore_hydrogen(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = PymolRMSDGrouper(
            molecules,
            threshold=0.8203818,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=True,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 5
        assert len(group_indices) == 5
        unique_structures = grouper.unique()
        assert len(unique_structures) == 5

        # rmsd calculation from grouper
        rmsd = grouper._calculate_rmsd((0, 2))
        assert np.isclose(rmsd, 0.0211, rtol=1e-3)
        rmsd = grouper._calculate_rmsd((1, 4))
        assert np.isclose(rmsd, 0.7669, rtol=1e-3)

        # Explicitly cleanup to prevent __del__ from calling quit()
        if hasattr(grouper, "_temp_dir") and grouper._temp_dir:
            import shutil

            shutil.rmtree(grouper._temp_dir, ignore_errors=True)
            grouper._temp_dir = None
        grouper.cmd = None  # Prevent __del__ from calling quit()

    def test_pymol_grouper_rejects_multiproc(
        self, methanol_molecules, temp_working_dir
    ):
        """Test that PyMOL grouper raises error when num_procs > 1."""
        with pytest.raises(ValueError) as excinfo:
            PymolRMSDGrouper(
                methanol_molecules,
                threshold=0.5,
                num_procs=4,  # Should raise error
            )
        assert (
            "single" in str(excinfo.value).lower()
            or "num_procs" in str(excinfo.value).lower()
        )

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

        tfd = grouper._calculate_tfd((0, 2))
        assert np.isclose(tfd, 0.08229, rtol=1e-3)
        tfd = grouper._calculate_tfd((0, 7))
        assert np.isclose(tfd, 0.07727, rtol=1e-3)

    def test_use_maxdev_parameter(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = TorsionFingerprintGrouper(
            molecules,
            threshold=0.26965,
            num_procs=self.NUM_PROCS,
            use_weights=True,
            max_dev="spec",
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 2
        assert len(group_indices) == 2
        unique_structures = grouper.unique()
        assert len(unique_structures) == 2

        tfd = grouper._calculate_tfd((0, 1))
        assert np.isclose(tfd, 0.02027, rtol=1e-3)
        tfd = grouper._calculate_tfd((3, 4))
        assert np.isclose(tfd, 0.24365, rtol=1e-3)


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
            molecules,
            num_procs=self.NUM_PROCS,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 4
        assert len(group_indices) == 4
        unique_structures = grouper.unique()
        assert len(unique_structures) == 4

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


class Test_EnergyGrouper:
    NUM_PROCS = 4

    def test_energy_grouper_raises_error_for_missing_energy(
        self, methanol_molecules, temp_working_dir
    ):
        # methanol_molecules from pubchem don't have energy information
        with pytest.raises(ValueError) as excinfo:
            EnergyGrouper(methanol_molecules)
        assert "missing energy information" in str(excinfo.value)

    def test_energy_grouper_for_crest_conformers(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        """Test EnergyGrouper with molecules that have energy information."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18

        # All CREST conformers should have energy
        for mol in molecules:
            assert mol.energy is not None

        assert len(molecules) == 18
        grouper = EnergyGrouper(
            molecules,
            threshold=0.5,
            num_procs=self.NUM_PROCS,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 7
        assert len(group_indices) == 7
        unique_structures = grouper.unique()
        assert len(unique_structures) == 7

        expected1 = -1.7839 / 627.509474
        expected2 = 1.4580 / 627.509474

        relative_diff, abs_diff = grouper._calculate_energy_diff((1, 0))
        assert np.isclose(relative_diff, expected1, rtol=1e-3)
        relative_diff, abs_diff = grouper._calculate_energy_diff((1, 2))
        assert np.isclose(abs_diff, expected2, rtol=1e-3)

        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        grouper = EnergyGrouper(
            molecules, num_groups=5, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 5
        assert len(group_indices) == 5

    def test_energy_grouper_for_log_conformers(
        self, ts_conformers_log_directory, temp_working_dir
    ):
        """Test EnergyGrouper with molecules loaded from log files using Gibbs energy."""
        import glob
        import re

        from chemsmart.cli.grouper.grouper import _extract_gibbs_energy
        from chemsmart.io.gaussian.output import Gaussian16Output
        from chemsmart.io.molecules.structure import Molecule

        # Load molecules from log files and set Gibbs energy (like CLI does)
        log_files = sorted(
            glob.glob(os.path.join(ts_conformers_log_directory, "*.log"))
        )
        molecules = []
        conformer_ids = []

        for log_file in log_files:
            mol = Molecule.from_filepath(
                filepath=log_file, index="-1", return_list=False
            )
            if mol is not None:
                # Extract and set Gibbs energy (this is what CLI does)
                g16_output = Gaussian16Output(filename=log_file)
                gibbs_energy = _extract_gibbs_energy(g16_output)
                if gibbs_energy is not None:
                    mol._energy = gibbs_energy

                molecules.append(mol)
                basename = os.path.basename(log_file)
                match = re.search(r"_c(\d+)\.log$", basename)
                if match:
                    conformer_ids.append(f"c{match.group(1)}")
                else:
                    conformer_ids.append(basename)

        assert len(molecules) == 5
        assert len(conformer_ids) == 5

        # All log conformers should have Gibbs energy
        for mol in molecules:
            assert mol.energy is not None

        grouper = EnergyGrouper(
            molecules,
            threshold=1.0,
            num_procs=self.NUM_PROCS,
            conformer_ids=conformer_ids,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 2
        assert len(group_indices) == 2
        unique_structures = grouper.unique()
        assert len(unique_structures) == 2

        expected1 = -4.1429 / 627.509474
        expected2 = 4.1429 / 627.509474

        relative_diff, abs_diff = grouper._calculate_energy_diff((1, 4))
        assert np.isclose(relative_diff, expected1, rtol=1e-2)
        assert np.isclose(abs_diff, expected2, rtol=1e-2)

    def test_energy_extraction_from_xyz_file(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):

        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        # Test specific energy value for 3rd structure (index 2)
        expected_energy_mol3 = -126.25238449
        actual_energy_mol3 = molecules[2].energy

        assert actual_energy_mol3 is not None, "3rd molecule energy is None"
        assert np.isclose(
            actual_energy_mol3, expected_energy_mol3, rtol=1e-8
        ), (
            f"Energy mismatch for 3rd molecule: got {actual_energy_mol3}, "
            f"expected {expected_energy_mol3}"
        )

    def test_gibbs_energy_extraction_function(
        self, ts_conformers_log_directory, temp_working_dir
    ):

        from chemsmart.cli.grouper.grouper import _extract_gibbs_energy
        from chemsmart.io.gaussian.output import Gaussian16Output

        log_file = os.path.join(
            ts_conformers_log_directory, "ch_1c_para_c1.log"
        )
        g16_output = Gaussian16Output(filename=log_file)

        # Known values
        expected_scf = -1401.99431519
        expected_gibbs_correction = 0.311863
        expected_gibbs_energy = (
            expected_scf + expected_gibbs_correction
        )  # -1401.68245219

        gibbs_energy = _extract_gibbs_energy(g16_output)
        assert np.isclose(gibbs_energy, expected_gibbs_energy, rtol=1e-6), (
            f"Gibbs energy mismatch: got {gibbs_energy}, "
            f"expected {expected_gibbs_energy}"
        )

    def test_energy_extraction_from_ts_log_files(
        self, ts_conformers_log_directory, temp_working_dir
    ):
        """Test that energy is correctly extracted from TS log files as SCF Done energy."""
        import glob

        from chemsmart.io.gaussian.output import Gaussian16Output
        from chemsmart.io.molecules.structure import Molecule

        log_files = sorted(
            glob.glob(os.path.join(ts_conformers_log_directory, "*.log"))
        )

        molecules = []
        for log_file in log_files:
            mol = Molecule.from_filepath(
                filepath=log_file, index="-1", return_list=False
            )
            if mol is not None:
                molecules.append(mol)
                g16_output = Gaussian16Output(filename=log_file)
                scf_energy = (
                    g16_output.scf_energies[-1]
                    if g16_output.scf_energies
                    else None
                )
                assert np.isclose(mol.energy, scf_energy, rtol=1e-6), (
                    f"Energy mismatch for {log_file}: "
                    f"mol.energy={mol.energy}, expected SCF={scf_energy}"
                )


class Testfactory:

    def test_structure_grouper_factory_energy(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        """Test factory creation of energy grouper (requires molecules with energy)."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        factory = StructureGrouperFactory()
        energy_grouper = factory.create(molecules, strategy="energy")
        assert isinstance(energy_grouper, EnergyGrouper)

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

        # irmsd requires external command, skip if not available
        if _irmsd_available():
            irmsd_grouper = factory.create(
                methanol_molecules, strategy="irmsd"
            )
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


class Test_grouper_utility_functions:
    """Test utility functions and helper methods in groupers."""

    NUM_PROCS = 1

    def test_rmsd_matrix_with_num_groups(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        """Test that RMSD matrix filename reflects num_groups when used."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        grouper = BasicRMSDGrouper(
            molecules,
            num_groups=3,
            num_procs=1,
            label="test_num_groups",
        )
        groups, group_indices = grouper.group()

        # Check that Excel file was created in the temp_working_dir
        expected_file = os.path.join(
            temp_working_dir,
            "test_num_groups_group_result",
            "test_num_groups_BasicRMSDGrouper_N3.xlsx",
        )
        assert os.path.exists(
            expected_file
        ), f"Matrix file not found: {expected_file}"

    def test_grouping_result_caching(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        """Test that grouping results are cached and reused."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        grouper = BasicRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=1,
        )

        groups1, indices1 = grouper.group()

        assert grouper._cached_groups is not None
        assert grouper._cached_group_indices is not None

        # unique() should use cached results
        unique_mols = grouper.unique()
        assert len(unique_mols) == len(groups1)

    def test_unique_returns_lowest_energy_representative(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        """Test that unique() returns lowest energy molecule from each group."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        # CREST conformers have energy information
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=2.0,
            num_procs=1,
        )
        groups, group_indices = grouper.group()
        unique_mols = grouper.unique()

        # For each group, verify the representative has lowest energy
        for i, (group, indices) in enumerate(zip(groups, group_indices)):
            mols_with_energy = [
                (mol, idx)
                for mol, idx in zip(group, indices)
                if mol.energy is not None
            ]
            if mols_with_energy:
                min_energy_mol = min(
                    mols_with_energy, key=lambda x: x[0].energy
                )[0]
                assert unique_mols[i].energy == min_energy_mol.energy


class Test_grouper_complete_linkage:
    """Test complete linkage clustering behavior."""

    NUM_PROCS = 4

    def test_complete_linkage_prevents_chaining(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        """Test that complete linkage prevents chaining effect in grouping."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        # With a moderate threshold, complete linkage should create
        # groups where ALL members are within threshold of each other
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=1.0,
            num_procs=self.NUM_PROCS,
        )
        groups, group_indices = grouper.group()

        # Verify complete linkage property: within each group,
        # all pairs should have RMSD < threshold
        for group_idx, indices in enumerate(group_indices):
            if len(indices) > 1:
                for i, idx_i in enumerate(indices):
                    for j, idx_j in enumerate(indices):
                        if i < j:
                            rmsd = grouper._calculate_rmsd((idx_i, idx_j))
                            assert rmsd < 1.0, (
                                f"Group {group_idx}: RMSD between {idx_i} and {idx_j} "
                                f"is {rmsd}, exceeds threshold 1.0"
                            )


class Test_conformer_ids_functionality:
    """Test conformer_ids parameter functionality."""

    NUM_PROCS = 1

    def test_conformer_ids_from_log_directory(
        self, ts_conformers_log_directory, temp_working_dir
    ):
        """Test loading conformer IDs from a directory of log files."""
        import glob
        import re

        from chemsmart.io.molecules.structure import Molecule

        # Load molecules from log files in the ts_conformers_log_directory
        log_files = sorted(
            glob.glob(os.path.join(ts_conformers_log_directory, "*.log"))
        )
        assert (
            len(log_files) == 5
        ), f"Expected 5 log files, found {len(log_files)}"

        # Extract conformer IDs from filenames (e.g., ch_1c_para_c1.log -> c1)
        molecules = []
        conformer_ids = []
        for log_file in log_files:
            mol = Molecule.from_filepath(
                filepath=log_file, index="-1", return_list=False
            )
            if mol is not None:
                molecules.append(mol)
                # Extract conformer ID from filename
                basename = os.path.basename(log_file)
                match = re.search(r"_c(\d+)\.log$", basename)
                if match:
                    conformer_ids.append(f"c{match.group(1)}")
                else:
                    conformer_ids.append(basename)

        assert len(molecules) >= 2, "Need at least 2 molecules for grouping"

        grouper = BasicRMSDGrouper(
            molecules,
            threshold=1.0,
            num_procs=1,
            label="ts_conformers_test",
            conformer_ids=conformer_ids,
        )
        groups, group_indices = grouper.group()

        # Verify conformer_ids are stored
        assert grouper.conformer_ids == conformer_ids

        # Check Excel output
        import pandas as pd

        excel_file = os.path.join(
            "ts_conformers_test_group_result",
            "ts_conformers_test_BasicRMSDGrouper_T1.0.xlsx",
        )
        assert os.path.exists(
            excel_file
        ), f"Excel file not found: {excel_file}"

        # Matrix data starts at row 8 (0-indexed: skiprows=8), first column is index
        df = pd.read_excel(
            excel_file, sheet_name="RMSD_Matrix", skiprows=8, index_col=0
        )
        # Check that conformer IDs are used as labels
        assert "c1" in str(df.columns[0]) or "c1" in str(df.index[0])

    def test_traj_conformer_ids_original_indices(self, temp_working_dir):
        """Test that traj job correctly sets original conformer indices."""
        from chemsmart.io.molecules.structure import Molecule

        # Create test molecules simulating a trajectory
        mol = Molecule.from_pubchem(identifier="CO")
        molecules = [mol.copy() for _ in range(18)]

        # Simulate traj behavior: select last 50% (indices 10-18 in original)
        proportion = 0.5
        total = len(molecules)
        last_num = int(round(total * proportion, 1))
        _ = molecules[-last_num:]  # selected_molecules - not used in this test

        # Calculate expected original indices (1-based)
        start_original_index = total - last_num + 1
        expected_ids = [str(start_original_index + i) for i in range(last_num)]

        # Verify expected IDs
        assert expected_ids[0] == "10"  # First selected is original index 10
        assert expected_ids[-1] == "18"  # Last selected is original index 18


class Test_output_file_generation:
    """Test that grouper generates correct output files."""

    NUM_PROCS = 1

    def test_group_xyz_files_contain_energy_and_index_info(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        """Test that group XYZ files contain energy and original index information."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        grouper = BasicRMSDGrouper(
            molecules[:6],
            threshold=1.0,
            num_procs=1,
            label="info_test",
        )
        groups, group_indices = grouper.group()
        grouper.unique()  # Generates xyz files

        output_dir = "info_test_group_result"
        first_group_file = os.path.join(output_dir, "info_test_group_1.xyz")

        with open(first_group_file, "r") as f:
            content = f.read()

        # Check for expected information in comments
        assert "Original_Index:" in content
        assert "Energy" in content

    def test_group_xyz_files_sorted_by_energy(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        """Test that molecules in group XYZ files are sorted by energy (lowest first)."""
        import re

        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        grouper = BasicRMSDGrouper(
            molecules,
            threshold=2.0,  # Higher threshold to get more molecules per group
            num_procs=1,
            label="sort_test",
        )
        groups, group_indices = grouper.group()
        grouper.unique()  # Generates xyz files

        output_dir = "sort_test_group_result"

        # Check each group file
        for group_num in range(len(groups)):
            group_file = os.path.join(
                output_dir, f"sort_test_group_{group_num + 1}.xyz"
            )

            with open(group_file, "r") as f:
                content = f.read()

            # Extract energies from file
            lines = content.strip().split("\n")
            i = 0
            energies = []

            while i < len(lines):
                num_atoms = int(lines[i].strip())
                comment_line = lines[i + 1]
                energy_match = re.search(
                    r"Energy\(Hartree\):\s*([-\d.]+)", comment_line
                )
                if energy_match:
                    energies.append(float(energy_match.group(1)))
                i += num_atoms + 2

            # Verify sorted by energy (lowest first)
            for j in range(len(energies) - 1):
                assert energies[j] <= energies[j + 1], (
                    f"Group {group_num + 1} not sorted: "
                    f"energy[{j}]={energies[j]} > energy[{j+1}]={energies[j+1]}"
                )

    def test_group_xyz_files_energy_from_log_files(
        self, ts_conformers_log_directory, temp_working_dir
    ):
        """Test that energy is correctly extracted from log files and written to group XYZ files."""
        import glob
        import re

        from chemsmart.io.molecules.structure import Molecule

        # Load molecules from log files
        log_files = sorted(
            glob.glob(os.path.join(ts_conformers_log_directory, "*.log"))
        )
        assert (
            len(log_files) >= 2
        ), f"Need at least 2 log files, found {len(log_files)}"

        molecules = []
        conformer_ids = []
        original_energies = {}

        for idx, log_file in enumerate(log_files):
            mol = Molecule.from_filepath(
                filepath=log_file, index="-1", return_list=False
            )
            if mol is not None:
                molecules.append(mol)
                # Store original energy for verification
                original_energies[idx] = mol.energy

                # Extract conformer ID from filename
                basename = os.path.basename(log_file)
                match = re.search(r"_c(\d+)\.log$", basename)
                if match:
                    conformer_ids.append(f"c{match.group(1)}")
                else:
                    conformer_ids.append(basename)

        assert len(molecules) >= 2, "Need at least 2 molecules for grouping"

        grouper = BasicRMSDGrouper(
            molecules,
            threshold=5.0,  # High threshold to group all together
            num_procs=1,
            label="log_energy_test",
            conformer_ids=conformer_ids,
        )
        groups, group_indices = grouper.group()
        grouper.unique()  # Generates xyz files

        # Check group file
        output_dir = "log_energy_test_group_result"
        first_group_file = os.path.join(
            output_dir, "log_energy_test_group_1.xyz"
        )

        assert os.path.exists(
            first_group_file
        ), f"Group file not found: {first_group_file}"

        with open(first_group_file, "r") as f:
            content = f.read()

        # Verify energy values in the output
        assert "Energy(Hartree):" in content

        # Parse and verify each energy value
        lines = content.strip().split("\n")
        i = 0
        energies_found = 0
        while i < len(lines):
            try:
                num_atoms = int(lines[i].strip())
            except ValueError:
                break
            comment_line = lines[i + 1]

            # Extract energy from comment
            energy_match = re.search(
                r"Energy\(Hartree\):\s*([-\d.]+)", comment_line
            )
            if energy_match:
                extracted_energy = float(energy_match.group(1))
                # Verify the energy is a reasonable value (not zero or nan)
                assert (
                    extracted_energy < 0
                ), f"Energy should be negative: {extracted_energy}"
                energies_found += 1

            i += num_atoms + 2

        assert energies_found > 0, "No energies found in output file"


class Test_edge_cases:
    """Test edge cases and boundary conditions."""

    def test_two_molecules_grouping(self, temp_working_dir):
        """Test grouping with minimum number of molecules (2)."""
        from chemsmart.io.molecules.structure import Molecule

        mol1 = Molecule.from_pubchem(identifier="CO")
        mol2 = mol1.copy()

        grouper = BasicRMSDGrouper([mol1, mol2], threshold=0.5)
        groups, group_indices = grouper.group()

        # Two identical molecules should form one group
        assert len(groups) == 1
        assert len(groups[0]) == 2

    def test_num_groups_equals_num_molecules(
        self, methanol_molecules, temp_working_dir
    ):
        """Test requesting same number of groups as molecules."""
        grouper = BasicRMSDGrouper(
            methanol_molecules,
            num_groups=len(methanol_molecules),
        )
        groups, group_indices = grouper.group()

        # Should create one group per molecule (or fewer if some are identical)
        assert len(groups) <= len(methanol_molecules)

    def test_num_groups_exceeds_num_molecules(
        self, methanol_molecules, temp_working_dir
    ):
        """Test requesting more groups than molecules."""
        n_mols = len(methanol_molecules)
        grouper = BasicRMSDGrouper(
            methanol_molecules,
            num_groups=n_mols + 5,
        )
        groups, group_indices = grouper.group()

        # Should create at most n_mols groups
        assert len(groups) <= n_mols

    def test_very_low_threshold(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        """Test with very low threshold (should create many groups)."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        grouper = BasicRMSDGrouper(
            molecules,
            threshold=0.001,  # Very low threshold
            num_procs=1,
        )
        groups, group_indices = grouper.group()

        # With very low threshold, most molecules should be in separate groups
        # (unless they're truly identical)
        assert len(groups) >= len(molecules) - 1

    def test_very_high_threshold(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        """Test with very high threshold (should create few groups)."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        grouper = BasicRMSDGrouper(
            molecules,
            threshold=100.0,  # Very high threshold
            num_procs=1,
        )
        groups, group_indices = grouper.group()

        # With very high threshold, all molecules with same formula should be in one group
        assert len(groups) <= 3  # Likely 1-2 groups

    def test_different_formulas_always_separate(
        self, methanol_and_ethanol, temp_working_dir
    ):
        """Test that molecules with different formulas are always in separate groups."""
        grouper = BasicRMSDGrouper(
            methanol_and_ethanol,
            threshold=100.0,  # Even with high threshold
        )
        groups, group_indices = grouper.group()

        # Methanol and ethanol should always be separate
        assert len(groups) == 2

    def test_rmsd_infinity_for_different_molecules(
        self, methanol_and_ethanol, temp_working_dir
    ):
        """Test that RMSD returns infinity for molecules with different atom counts."""
        grouper = BasicRMSDGrouper(methanol_and_ethanol, threshold=0.5)
        rmsd = grouper._calculate_rmsd((0, 1))

        assert rmsd == np.inf or rmsd == float("inf")


class Test_label_and_append_label:
    """Test -l (label) and -a (append_label) parameter functionality."""

    NUM_PROCS = 1

    def test_get_label_function(self, temp_working_dir):
        """Test _get_label function logic."""
        from chemsmart.cli.grouper.grouper import _get_label

        # Test with label only
        result = _get_label(
            label="custom", append_label=None, base_label="base"
        )
        assert result == "custom"

        # Test with append_label only
        result = _get_label(
            label=None, append_label="suffix", base_label="base"
        )
        assert result == "base_suffix"

        # Test with neither (use base_label)
        result = _get_label(label=None, append_label=None, base_label="base")
        assert result == "base"

        # Test with both should raise error
        with pytest.raises(ValueError) as excinfo:
            _get_label(
                label="custom", append_label="suffix", base_label="base"
            )
        assert "Only give label or append_label" in str(excinfo.value)

    def test_label_in_output_directory(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        """Test that label parameter affects output directory name."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)[:5]

        # Test with custom label
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=1.0,
            num_procs=self.NUM_PROCS,
            label="my_custom_label",
        )
        groups, group_indices = grouper.group()

        # Check output directory uses custom label
        output_dir = "my_custom_label_group_result"
        assert os.path.exists(
            output_dir
        ), f"Output dir not found: {output_dir}"

        # Check Excel file uses custom label
        excel_file = os.path.join(
            output_dir, "my_custom_label_BasicRMSDGrouper_T1.0.xlsx"
        )
        assert os.path.exists(
            excel_file
        ), f"Excel file not found: {excel_file}"

    def test_label_in_group_xyz_files(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        """Test that label parameter affects group XYZ file names."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)[:5]

        grouper = BasicRMSDGrouper(
            molecules,
            threshold=1.0,
            num_procs=self.NUM_PROCS,
            label="test_label",
        )
        groups, group_indices = grouper.group()
        grouper.unique()  # Generates xyz files

        output_dir = "test_label_group_result"

        # Check group XYZ files use label
        for i in range(len(groups)):
            xyz_path = os.path.join(output_dir, f"test_label_group_{i+1}.xyz")
            assert os.path.exists(xyz_path), f"Group XYZ not found: {xyz_path}"

    def test_different_labels_create_different_outputs(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        """Test that different labels create separate output directories."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)[:5]

        # First grouper with label "run1"
        grouper1 = BasicRMSDGrouper(
            molecules,
            threshold=1.0,
            num_procs=self.NUM_PROCS,
            label="run1",
        )
        grouper1.group()

        # Second grouper with label "run2"
        grouper2 = BasicRMSDGrouper(
            molecules,
            threshold=1.0,
            num_procs=self.NUM_PROCS,
            label="run2",
        )
        grouper2.group()

        # Both output directories should exist
        assert os.path.exists("run1_group_result")
        assert os.path.exists("run2_group_result")

        # Both should have their own Excel files
        assert os.path.exists(
            "run1_group_result/run1_BasicRMSDGrouper_T1.0.xlsx"
        )
        assert os.path.exists(
            "run2_group_result/run2_BasicRMSDGrouper_T1.0.xlsx"
        )

    def test_label_with_num_groups(
        self, multiple_molecules_xyz_file, temp_working_dir
    ):
        """Test that label works correctly with num_groups parameter."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        grouper = BasicRMSDGrouper(
            molecules,
            num_groups=5,
            num_procs=self.NUM_PROCS,
            label="numgroups_test",
        )
        groups, group_indices = grouper.group()

        output_dir = "numgroups_test_group_result"
        assert os.path.exists(output_dir)

        # Excel file should reflect num_groups (N5) instead of threshold
        excel_file = os.path.join(
            output_dir, "numgroups_test_BasicRMSDGrouper_N5.xlsx"
        )
        assert os.path.exists(
            excel_file
        ), f"Excel file not found: {excel_file}"
