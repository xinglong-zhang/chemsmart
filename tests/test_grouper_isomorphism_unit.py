"""
Direct unit tests for chemsmart.jobs.grouper.isomorphism.RDKitIsomorphismGrouper,
covering branches the real end-to-end grouping test in test_groupers.py
does not reach: the invalid-molecule own-group fallback in group(), and
the grouping_time-not-given branch in _record_results.

RDKit conversion and hash-computation failures are now handled inside
``Molecule.get_rdkit_hash`` (chemsmart.io.molecules.structure) rather
than in this grouper, and are covered there in test_structures.py.
"""

from unittest.mock import patch

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.grouper.isomorphism import RDKitIsomorphismGrouper


def _make_water():
    return Molecule(
        symbols=["O", "H"],
        positions=[[0.0, 0.0, 0.0], [0.96, 0.0, 0.0]],
        charge=0,
        multiplicity=1,
    )


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
        with patch.object(Molecule, "get_rdkit_hash", return_value=None):
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
        # RDKitIsomorphismGrouper does not support multiprocessing, so
        # any requested num_procs > 1 is clamped down to 1.
        grouper = RDKitIsomorphismGrouper([_make_water()], num_procs=4)
        assert repr(grouper) == "RDKitIsomorphismGrouper(num_procs=1)"
