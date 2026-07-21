"""
Direct unit tests for :class:`SubstituentPreprocessor` in
``chemsmart.jobs.iterate.iterate``.
"""

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.iterate.iterate import SubstituentPreprocessor


class TestSubstituentPreprocessorEarlyReturn:
    def test_unsaturated_link_atom_returns_unchanged(self):
        ch2 = Molecule(
            symbols=["C", "H", "H"],
            positions=[
                [0.0, 0.0, 0.0],
                [1.09, 0.0, 0.0],
                [-0.36, 1.03, 0.0],
            ],
        )
        pre = SubstituentPreprocessor(molecule=ch2, link_index=1)
        result = pre.run()
        assert result is ch2
        assert pre._final_indices == [0, 1, 2]

    def test_fallback_indices_when_unsaturated(self):
        ch2 = Molecule(
            symbols=["C", "H", "H"],
            positions=[
                [0.0, 0.0, 0.0],
                [1.09, 0.0, 0.0],
                [-0.36, 1.03, 0.0],
            ],
        )
        pre = SubstituentPreprocessor(molecule=ch2, link_index=1)
        assert pre._final_indices is None
        assert pre.get_new_link_index() == 1


class TestSubstituentPreprocessorAutoDetect:
    def test_saturated_link_atom_removes_smallest_group(
        self, methanol_molecule
    ):
        pre = SubstituentPreprocessor(molecule=methanol_molecule, link_index=1)
        result = pre.run()
        assert len(result) == 5
        assert pre._final_indices is not None
