"""
Direct unit tests for the pure geometry/array helper methods of
:class:`IterateAnalyzer` in ``chemsmart.jobs.iterate.iterate``.

The full ``run()``/``_optimize_lagrange`` numerical optimization
pipeline is not exercised here (it requires realistic multi-atom
skeleton/substituent geometries and is inherently slow); this focuses
on the well-defined, independently-testable array/geometry
transformations and input-validation branches.
"""

import numpy as np
import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.iterate.iterate import IterateAnalyzer


class TestMoleculeToArray:
    def test_converts_symbols_and_positions(self, methanol_molecule):
        arr = IterateAnalyzer._molecule_to_array(methanol_molecule)
        assert arr.shape == (6, 4)
        # Carbon (Z=6), Oxygen (Z=8)
        assert arr[0, 0] == 6
        assert arr[1, 0] == 8
        np.testing.assert_allclose(arr[:, 1:4], methanol_molecule.positions)


class TestUpdateMoleculePositions:
    def test_returns_new_molecule_with_updated_positions(
        self, methanol_molecule
    ):
        new_positions = methanol_molecule.positions + 10.0
        updated = IterateAnalyzer._update_molecule_positions(
            methanol_molecule, new_positions
        )
        assert updated is not methanol_molecule
        np.testing.assert_allclose(updated.positions, new_positions)
        assert updated.symbols == methanol_molecule.symbols
        assert updated.charge == methanol_molecule.charge
        assert updated.multiplicity == methanol_molecule.multiplicity


class TestCombineMolecules:
    def test_concatenates_symbols_and_positions(
        self, methanol_molecule, ethanol_molecule
    ):
        combined = IterateAnalyzer._combine_molecules(
            methanol_molecule, ethanol_molecule
        )
        assert len(combined) == len(methanol_molecule) + len(ethanol_molecule)
        assert combined.chemical_symbols == (
            list(methanol_molecule.chemical_symbols)
            + list(ethanol_molecule.chemical_symbols)
        )
        np.testing.assert_allclose(
            combined.positions[: len(methanol_molecule)],
            methanol_molecule.positions,
        )

    def test_inherits_charge_and_multiplicity_from_first_molecule(
        self, methanol_molecule, ethanol_molecule
    ):
        combined = IterateAnalyzer._combine_molecules(
            methanol_molecule, ethanol_molecule
        )
        assert combined.charge == methanol_molecule.charge
        assert combined.multiplicity == methanol_molecule.multiplicity

    def test_merges_frozen_atoms_with_offset(self, methanol_molecule):
        mol1 = methanol_molecule
        mol2 = Molecule(
            symbols=["H", "H"],
            positions=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            frozen_atoms=[1, 0],
        )
        combined = IterateAnalyzer._combine_molecules(mol1, mol2)
        assert combined.frozen_atoms == [0] * len(mol1) + [1, 0]


class TestCalcRelativeCoords:
    def test_base_atom_becomes_origin(self):
        arr = np.array(
            [
                [6.0, 1.0, 1.0, 1.0],
                [8.0, 2.0, 3.0, 4.0],
                [1.0, 0.0, 0.0, 0.0],
            ]
        )
        result = IterateAnalyzer._calc_relative_coords(arr, base_index=0)
        np.testing.assert_allclose(result[0, 1:4], [0.0, 0.0, 0.0])
        np.testing.assert_allclose(result[1, 1:4], [1.0, 2.0, 3.0])
        np.testing.assert_allclose(result[2, 1:4], [-1.0, -1.0, -1.0])
        # atomic numbers unchanged
        np.testing.assert_allclose(result[:, 0], arr[:, 0])

    def test_raises_for_negative_base_index(self):
        arr = np.zeros((3, 4))
        with pytest.raises(IndexError, match="out of range"):
            IterateAnalyzer._calc_relative_coords(arr, base_index=-1)

    def test_raises_for_out_of_range_base_index(self):
        arr = np.zeros((3, 4))
        with pytest.raises(IndexError, match="out of range"):
            IterateAnalyzer._calc_relative_coords(arr, base_index=3)


class TestFibonacciSphere:
    def test_returns_correct_shape(self):
        points = IterateAnalyzer._fibonacci_sphere(50)
        assert points.shape == (50, 3)

    def test_points_are_unit_vectors(self):
        points = IterateAnalyzer._fibonacci_sphere(96)
        norms = np.linalg.norm(points, axis=1)
        np.testing.assert_allclose(norms, 1.0, atol=1e-10)

    def test_z_coordinates_span_full_range(self):
        points = IterateAnalyzer._fibonacci_sphere(100)
        z_coords = points[:, 2]
        assert z_coords.max() < 1.0
        assert z_coords.min() > -1.0
        assert z_coords.max() - z_coords.min() > 1.9  # spans most of [-1, 1]


class TestFindOptimalPositionValidation:
    def test_raises_for_out_of_range_skeleton_link_index(self):
        skeleton = np.array([[6.0, 0.0, 0.0, 0.0], [8.0, 1.0, 0.0, 0.0]])
        sub = np.array([[1.0, 0.0, 0.0, 0.0]])
        with pytest.raises(IndexError, match="skeleton_link_index"):
            IterateAnalyzer._find_optimal_position(
                skeleton, sub, skeleton_link_index=5, sub_link_index=0
            )

    def test_raises_for_out_of_range_sub_link_index(self):
        skeleton = np.array([[6.0, 0.0, 0.0, 0.0], [8.0, 1.0, 0.0, 0.0]])
        sub = np.array([[1.0, 0.0, 0.0, 0.0]])
        with pytest.raises(IndexError, match="sub_link_index"):
            IterateAnalyzer._find_optimal_position(
                skeleton, sub, skeleton_link_index=0, sub_link_index=5
            )

    def test_raises_for_unknown_method(self):
        skeleton = np.array([[6.0, 0.0, 0.0, 0.0], [8.0, 1.0, 0.0, 0.0]])
        sub = np.array([[1.0, 0.0, 0.0, 0.0]])
        with pytest.raises(ValueError, match="Unknown method"):
            IterateAnalyzer._find_optimal_position(
                skeleton,
                sub,
                skeleton_link_index=0,
                sub_link_index=0,
                method="not_a_real_method",
            )


class TestIterateAnalyzerConstruction:
    def test_converts_link_indices_to_zero_based(
        self, methanol_molecule, ethanol_molecule
    ):
        analyzer = IterateAnalyzer(
            skeleton=methanol_molecule,
            substituent=ethanol_molecule,
            skeleton_link_index=1,
            substituent_link_index=2,
        )
        assert analyzer.skeleton_link_index == 0
        assert analyzer.substituent_link_index == 1

    def test_default_parameters(self, methanol_molecule, ethanol_molecule):
        analyzer = IterateAnalyzer(
            skeleton=methanol_molecule,
            substituent=ethanol_molecule,
            skeleton_link_index=1,
            substituent_link_index=1,
        )
        assert analyzer.buffer == 0.3
        assert analyzer.method == "lagrange_multipliers"
        assert analyzer.sphere_direction_samples_num == 96
        assert analyzer.axial_rotations_sample_num == 6
