"""
Direct unit tests for :class:`TanimotoSimilarityGrouper` in
``chemsmart.jobs.grouper.tanimoto`` that supplement the integration-style
tests in ``test_groupers.py`` (which only exercise the default "rdkit"
fingerprint under threshold-based grouping).

Covers constructor validation, ``_get_fingerprint`` for every supported
fingerprint type, the ``ignore_hydrogens`` path, the ``num_groups``-based
auto-threshold grouping strategy (including the "more groups requested
than molecules" and merge-to-target branches), and ``__repr__``.
"""

import os

import pytest

from chemsmart.io.xyz.xyzfile import XYZFile
from chemsmart.jobs.grouper.tanimoto import TanimotoSimilarityGrouper


@pytest.fixture()
def crest_molecules(xyz_directory):
    xyz_file = XYZFile(
        filename=os.path.join(xyz_directory, "crest_conformers.xyz")
    )
    return xyz_file.get_molecules(index=":", return_list=True)


@pytest.mark.usefixtures("temporary_working_dir")
class TestTanimotoGrouperConstruction:
    def test_threshold_and_num_groups_mutually_exclusive(
        self, methanol_molecules
    ):
        with pytest.raises(ValueError, match="Cannot specify both"):
            TanimotoSimilarityGrouper(
                methanol_molecules, threshold=0.9, num_groups=2
            )

    def test_default_threshold_is_point_nine(self, methanol_molecules):
        grouper = TanimotoSimilarityGrouper(methanol_molecules)
        assert grouper.threshold == 0.9
        assert grouper.num_groups is None

    def test_num_groups_leaves_threshold_none(self, methanol_molecules):
        grouper = TanimotoSimilarityGrouper(methanol_molecules, num_groups=1)
        assert grouper.threshold is None
        assert grouper.num_groups == 1

    def test_fingerprint_type_is_lowercased(self, methanol_molecules):
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules, fingerprint_type="MACCS"
        )
        assert grouper.fingerprint_type == "maccs"


@pytest.mark.usefixtures("temporary_working_dir")
class TestGetFingerprint:
    @pytest.mark.parametrize(
        "fingerprint_type",
        [
            "rdkit",
            "rdk",
            "morgan",
            "maccs",
            "atompair",
            "torsion",
            "usr",
            "usrcat",
        ],
    )
    def test_supported_fingerprint_types_produce_a_fingerprint(
        self, methanol_molecules, fingerprint_type
    ):
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules, fingerprint_type=fingerprint_type
        )
        fp = grouper._get_fingerprint(grouper.rdkit_molecules[0])
        assert fp is not None

    def test_unknown_fingerprint_type_falls_back_to_rdkit_default(
        self, methanol_molecules
    ):
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules, fingerprint_type="not-a-real-type"
        )
        fp = grouper._get_fingerprint(grouper.rdkit_molecules[0])
        assert fp is not None

    def test_fingerprint_generation_failure_returns_none(
        self, methanol_molecules
    ):
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules, fingerprint_type="morgan"
        )
        # Passing a non-Mol object makes the RDKit call raise internally;
        # the method should swallow it and return None rather than
        # propagating the exception.
        assert grouper._get_fingerprint(object()) is None


@pytest.mark.usefixtures("temporary_working_dir")
class TestIgnoreHydrogens:
    def test_ignore_hydrogens_still_groups_successfully(
        self, methanol_molecules
    ):
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules, ignore_hydrogens=True
        )
        assert len(grouper.rdkit_molecules) == len(methanol_molecules)
        groups, index_groups = grouper.group()
        assert len(groups) == 1


@pytest.mark.usefixtures("temporary_working_dir")
class TestNumGroupsGrouping:
    def test_num_groups_creates_requested_number_of_groups(
        self, crest_molecules
    ):
        grouper = TanimotoSimilarityGrouper(
            crest_molecules,
            num_groups=3,
            fingerprint_type="usrcat",
            num_procs=2,
        )
        groups, index_groups = grouper.group()
        assert len(groups) == 3
        assert len(index_groups) == 3
        assert grouper._auto_threshold is not None

    def test_num_groups_at_least_as_large_as_molecule_count_gives_singletons(
        self, methanol_molecules
    ):
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules, num_groups=10, fingerprint_type="usrcat"
        )
        groups, index_groups = grouper.group()
        assert len(groups) == len(methanol_molecules)
        assert all(len(g) == 1 for g in groups)

    def test_num_groups_one_merges_everything(self, crest_molecules):
        grouper = TanimotoSimilarityGrouper(
            crest_molecules,
            num_groups=1,
            fingerprint_type="usrcat",
            num_procs=2,
        )
        groups, index_groups = grouper.group()
        assert len(groups) == 1
        assert sum(len(g) for g in groups) == len(crest_molecules)


@pytest.mark.usefixtures("temporary_working_dir")
class TestInvalidMoleculeHandling:
    def test_molecule_with_no_rdkit_conversion_is_skipped(
        self, methanol_molecules, mocker
    ):
        from chemsmart.io.molecules.structure import Molecule

        # ``Molecule.to_rdkit`` is only used to build the RDKit molecule
        # for the 3D "usr"/"usrcat" descriptors; other fingerprint types
        # perceive bonds directly from the XYZ geometry instead.
        mocker.patch.object(Molecule, "to_rdkit", return_value=None)
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules, fingerprint_type="usr"
        )
        assert grouper.rdkit_molecules == []

    def test_ignore_hydrogens_falls_back_on_removehs_failure(
        self, methanol_molecules, mocker
    ):
        mocker.patch(
            "chemsmart.jobs.grouper.tanimoto.Chem.RemoveHs",
            side_effect=Exception("boom"),
        )
        # The with-H fallback on a failed RemoveHs only exists on the
        # "usr"/"usrcat" preparation path.
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules, ignore_hydrogens=True, fingerprint_type="usr"
        )
        # Falls back to the original (with-H) rdkit mol for every molecule.
        assert len(grouper.rdkit_molecules) == len(methanol_molecules)

    def test_group_returns_empty_when_no_valid_fingerprints(
        self, methanol_molecules, mocker
    ):
        grouper = TanimotoSimilarityGrouper(methanol_molecules)
        mocker.patch.object(grouper, "_get_fingerprint", return_value=None)
        groups, index_groups = grouper.group()
        assert groups == []
        assert index_groups == []


@pytest.mark.usefixtures("temporary_working_dir")
@pytest.mark.usefixtures("temporary_working_dir")
class TestRecordResultsDirect:
    def test_record_results_without_grouping_time_or_index_groups(
        self, methanol_molecules, tmp_path
    ):
        import numpy as np

        grouper = TanimotoSimilarityGrouper(
            methanol_molecules, threshold=0.9, label=str(tmp_path / "run")
        )
        n = len(methanol_molecules)
        matrix = np.eye(n, dtype=np.float32)
        # Should not raise even with grouping_time=None and
        # index_groups=None (both optional-argument False branches).
        grouper._record_results(
            tanimoto_matrix=matrix,
            valid_indices=list(range(n)),
            grouping_time=None,
            groups=None,
            index_groups=None,
        )


@pytest.mark.usefixtures("temporary_working_dir")
class TestRepr:
    def test_repr_with_threshold(self, methanol_molecules):
        grouper = TanimotoSimilarityGrouper(methanol_molecules, threshold=0.8)
        rep = repr(grouper)
        assert "threshold=0.8" in rep
        assert "TanimotoSimilarityGrouper" in rep

    def test_repr_with_num_groups(self, methanol_molecules):
        grouper = TanimotoSimilarityGrouper(methanol_molecules, num_groups=2)
        rep = repr(grouper)
        assert "num_groups=2" in rep
