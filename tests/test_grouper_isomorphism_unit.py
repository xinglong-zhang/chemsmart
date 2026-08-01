"""
Direct unit tests for chemsmart.jobs.grouper.isomorphism.RDKitIsomorphismGrouper,
covering branches the real end-to-end grouping test in test_groupers.py
does not reach: RDKit conversion failures, hash-computation failures,
the invalid-molecule own-group fallback in group(), and the
grouping_time-not-given branch in _record_results.
"""

from unittest.mock import MagicMock, patch

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.grouper.isomorphism import RDKitIsomorphismGrouper


def _make_water():
    return Molecule(
        symbols=["O", "H"],
        positions=[[0.0, 0.0, 0.0], [0.96, 0.0, 0.0]],
        charge=0,
        multiplicity=1,
    )


class TestMolToRdkit:
    def test_returns_none_when_rdkit_parse_fails(self):
        grouper = RDKitIsomorphismGrouper([_make_water()])
        with patch(
            "chemsmart.jobs.grouper.isomorphism.Chem.MolFromXYZBlock",
            return_value=None,
        ):
            assert grouper._mol_to_rdkit(_make_water()) is None

    def test_unexpected_exception_is_caught_and_returns_none(self):
        grouper = RDKitIsomorphismGrouper([_make_water()])
        with patch(
            "chemsmart.jobs.grouper.isomorphism.Chem.MolFromXYZBlock",
            side_effect=RuntimeError("boom"),
        ):
            assert grouper._mol_to_rdkit(_make_water()) is None

    def test_ignore_hydrogens_falls_back_when_removehs_fails(self):
        """If Chem.RemoveHs raises, the original (H-containing) molecule
        is re-parsed and returned instead of failing outright."""
        grouper = RDKitIsomorphismGrouper(
            [_make_water()], ignore_hydrogens=True
        )
        with patch(
            "chemsmart.jobs.grouper.isomorphism.Chem.RemoveHs",
            side_effect=RuntimeError("boom"),
        ):
            result = grouper._mol_to_rdkit(_make_water())
        assert result is not None

    def test_ignore_hydrogens_succeeds_and_resanitizes(self):
        """The normal (no-exception) path: hydrogens are actually
        removed and the molecule is re-sanitized without kekulization."""
        methanol = Molecule(
            symbols=["C", "H", "H", "H", "O", "H"],
            positions=[
                [0.0, 0.0, 0.0],
                [1.09, 0.0, 0.0],
                [-0.36, 1.03, 0.0],
                [-0.36, -0.51, 0.89],
                [-0.4, -0.7, -1.0],
                [-1.35, -0.7, -1.0],
            ],
            charge=0,
            multiplicity=1,
        )
        grouper = RDKitIsomorphismGrouper([methanol], ignore_hydrogens=True)
        result = grouper._mol_to_rdkit(methanol)
        assert result is not None

    def test_ignore_hydrogens_fallback_reconversion_also_fails(self):
        """Covers the "still None after the fallback re-conversion"
        arm: RemoveHs fails, and the retry MolFromXYZBlock call also
        fails to produce a molecule."""
        grouper = RDKitIsomorphismGrouper(
            [_make_water()], ignore_hydrogens=True
        )
        call_count = {"n": 0}
        from rdkit import Chem as real_chem

        real_from_xyz_block = real_chem.MolFromXYZBlock

        def fake_from_xyz_block(xyz_string):
            call_count["n"] += 1
            if call_count["n"] == 1:
                return real_from_xyz_block(xyz_string)
            return None

        with (
            patch(
                "chemsmart.jobs.grouper.isomorphism.Chem.MolFromXYZBlock",
                side_effect=fake_from_xyz_block,
            ),
            patch(
                "chemsmart.jobs.grouper.isomorphism.Chem.RemoveHs",
                side_effect=RuntimeError("boom"),
            ),
        ):
            result = grouper._mol_to_rdkit(_make_water())

        assert result is None
        assert call_count["n"] == 2


class TestGetMolHash:
    def test_none_input_returns_none(self):
        grouper = RDKitIsomorphismGrouper([_make_water()])
        assert grouper._get_mol_hash(None) is None

    def test_hash_failure_returns_none(self):
        grouper = RDKitIsomorphismGrouper([_make_water()])
        with patch(
            "chemsmart.jobs.grouper.isomorphism.rdMolHash.MolHash",
            side_effect=RuntimeError("boom"),
        ):
            assert grouper._get_mol_hash(MagicMock()) is None


class TestGroupInvalidMoleculeFallback:
    def test_molecule_with_no_hash_gets_its_own_group(
        self, tmp_path, monkeypatch
    ):
        """group() calls self.record(...) unconditionally, which writes
        output files relative to cwd when output_dir is unset -- chdir
        into a throwaway directory to avoid polluting the repo."""
        monkeypatch.chdir(tmp_path)
        mol = _make_water()
        grouper = RDKitIsomorphismGrouper([mol], num_procs=1)
        with patch.object(
            RDKitIsomorphismGrouper, "_mol_to_rdkit", return_value=None
        ):
            groups, index_groups = grouper.group()
        assert groups == [[mol]]
        assert index_groups == [[0]]


class TestRecordResultsWithoutGroupingTime:
    def test_grouping_time_none_skips_time_header(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        mol = _make_water()
        grouper = RDKitIsomorphismGrouper([mol], num_procs=1, label="test")

        # Should not raise despite grouping_time being omitted.
        grouper._record_results(
            hashes=["abc"],
            groups=[[mol]],
            index_groups=[[0]],
            grouping_time=None,
        )


class TestRepr:
    def test_repr_includes_num_procs(self):
        grouper = RDKitIsomorphismGrouper([_make_water()], num_procs=4)
        assert repr(grouper) == "RDKitIsomorphismGrouper(num_procs=4)"
