"""Tests for :class:`IterateAnalyzer` in ``chemsmart.jobs.iterate.iterate``:
geometry helpers and the Lagrange-multiplier optimisation used to place
substituents on skeletons.
"""

from unittest.mock import patch

import numpy as np
import pytest

from chemsmart.jobs.iterate.iterate import IterateAnalyzer
from tests.helpers import make_scripted_minimize


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

    def test_merges_frozen_atoms_with_offset(
        self, methanol_molecule, frozen_hydrogen_molecule
    ):
        combined = IterateAnalyzer._combine_molecules(
            methanol_molecule, frozen_hydrogen_molecule
        )
        assert combined.frozen_atoms == [0] * len(methanol_molecule) + [1, 0]


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


class TestIterateAnalyzerRun:
    def test_run_returns_combined_molecule_on_success(
        self,
        methanol_fluorine_iterate_analyzer,
        methanol_molecule,
        fluorine_atom_molecule,
    ):
        result = methanol_fluorine_iterate_analyzer.run()
        assert result is not None
        assert len(result) == len(methanol_molecule) + len(
            fluorine_atom_molecule
        )

    def test_run_returns_none_when_optimization_fails(
        self, methanol_fluorine_iterate_analyzer
    ):
        with patch.object(
            IterateAnalyzer, "_find_optimal_position", return_value=None
        ):
            assert methanol_fluorine_iterate_analyzer.run() is None


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


class TestOptimizeLagrangeHappyPath:
    def test_single_atom_substituent_returns_valid_position(
        self, iterate_skeleton_array
    ):
        # A 1-atom substituent has zero relative offset from its own
        # link atom, so sub_centroid_norm < 1e-3 (fallback principal
        # axis branch), and the sub link atom starts away from the
        # skeleton link (initial_link_dist >= 1e-6, normal branch).
        sub = np.array([[1.0, 5.0, 5.0, 5.0]])
        result = IterateAnalyzer._find_optimal_position(
            iterate_skeleton_array,
            sub,
            skeleton_link_index=0,
            sub_link_index=0,
            sphere_direction_samples_num=4,
            axial_rotations_sample_num=2,
        )
        assert result is not None
        assert result.shape == sub.shape
        assert result[0, 0] == 1.0  # atomic number preserved

    def test_multi_atom_substituent_with_nonzero_centroid(
        self, iterate_skeleton_array
    ):
        # A 2-atom substituent gives a non-trivial relative offset,
        # exercising the sub_centroid_norm >= 1e-3 branch.
        sub = np.array([[1.0, 5.0, 5.0, 5.0], [1.0, 5.5, 5.5, 5.5]])
        result = IterateAnalyzer._find_optimal_position(
            iterate_skeleton_array,
            sub,
            skeleton_link_index=0,
            sub_link_index=0,
            sphere_direction_samples_num=4,
            axial_rotations_sample_num=2,
        )
        assert result is not None
        assert result.shape == sub.shape

    def test_sub_link_starting_at_skeleton_link_position(
        self, iterate_skeleton_array
    ):
        # sub_coord's link atom placed exactly at skeleton_link_coords
        # forces initial_link_dist < 1e-6, taking the "default
        # direction" branch (skeleton has multiple atoms not centered
        # on the link, so dir_norm >= 1e-6 here).
        sub = np.array([[1.0, 0.0, 0.0, 0.0], [1.0, 0.5, 0.5, 0.5]])
        result = IterateAnalyzer._find_optimal_position(
            iterate_skeleton_array,
            sub,
            skeleton_link_index=0,
            sub_link_index=0,
            sphere_direction_samples_num=4,
            axial_rotations_sample_num=2,
        )
        assert result is not None

    def test_single_atom_skeleton_forces_dir_norm_fallback(
        self, iterate_single_atom_skeleton_array
    ):
        # A single-atom skeleton means skeleton_com == skeleton_link,
        # so `direction` has near-zero norm, forcing the [1, 0, 0]
        # fallback direction branch (in addition to initial_link_dist
        # < 1e-6, since the sub link also starts at the same point).
        sub = np.array([[1.0, 0.0, 0.0, 0.0]])
        result = IterateAnalyzer._find_optimal_position(
            iterate_single_atom_skeleton_array,
            sub,
            skeleton_link_index=0,
            sub_link_index=0,
            sphere_direction_samples_num=4,
            axial_rotations_sample_num=2,
        )
        assert result is not None


class TestOptimizeLagrangePreOptimizationFallback:
    def test_pre_optimization_failure_falls_back_to_projection(
        self,
        caplog,
        iterate_skeleton_array,
        failed_minimize_result,
    ):
        # First minimize() call (pre-optimization) fails; the rest use
        # the real solver.
        sub = np.array([[1.0, 5.0, 5.0, 5.0]])
        with patch(
            "chemsmart.jobs.iterate.iterate.minimize",
            side_effect=make_scripted_minimize(failed_minimize_result),
        ):
            with caplog.at_level("WARNING"):
                result = IterateAnalyzer._find_optimal_position(
                    iterate_skeleton_array,
                    sub,
                    skeleton_link_index=0,
                    sub_link_index=0,
                    sphere_direction_samples_num=4,
                    axial_rotations_sample_num=2,
                )
        assert "Initial position optimization failed" in caplog.text
        assert result is not None

    def test_pre_optimization_exception_falls_back_to_projection(
        self, caplog, iterate_skeleton_array
    ):
        # First minimize() call (pre-optimization) raises; the rest use
        # the real solver.
        sub = np.array([[1.0, 5.0, 5.0, 5.0]])
        with patch(
            "chemsmart.jobs.iterate.iterate.minimize",
            side_effect=make_scripted_minimize(RuntimeError("boom")),
        ):
            with caplog.at_level("WARNING"):
                result = IterateAnalyzer._find_optimal_position(
                    iterate_skeleton_array,
                    sub,
                    skeleton_link_index=0,
                    sub_link_index=0,
                    sphere_direction_samples_num=4,
                    axial_rotations_sample_num=2,
                )
        assert "Initial position optimization error" in caplog.text
        assert result is not None


class TestOptimizeLagrangeMainLoopFailures:
    def test_all_attempts_fail_returns_none(
        self,
        caplog,
        iterate_skeleton_array,
        failed_minimize_result,
    ):
        # Let the pre-optimization (first call) succeed via the real
        # solver so we reach the main loop; fail everything after that.
        sub = np.array([[1.0, 5.0, 5.0, 5.0]])
        with patch(
            "chemsmart.jobs.iterate.iterate.minimize",
            side_effect=make_scripted_minimize(
                None, default=failed_minimize_result
            ),
        ):
            with caplog.at_level("ERROR"):
                result = IterateAnalyzer._find_optimal_position(
                    iterate_skeleton_array,
                    sub,
                    skeleton_link_index=0,
                    sub_link_index=0,
                    sphere_direction_samples_num=4,
                    axial_rotations_sample_num=2,
                )
        assert result is None
        assert "All optimization attempts failed" in caplog.text

    def test_exception_in_main_loop_is_logged_and_skipped(
        self, caplog, iterate_skeleton_array
    ):
        # First call: real pre-optimization. Second call (first main-loop
        # attempt): raise. Remaining calls: real solver, so the run still
        # succeeds overall.
        sub = np.array([[1.0, 5.0, 5.0, 5.0]])
        with patch(
            "chemsmart.jobs.iterate.iterate.minimize",
            side_effect=make_scripted_minimize(
                None, RuntimeError("solver crashed")
            ),
        ):
            with caplog.at_level("WARNING"):
                result = IterateAnalyzer._find_optimal_position(
                    iterate_skeleton_array,
                    sub,
                    skeleton_link_index=0,
                    sub_link_index=0,
                    sphere_direction_samples_num=4,
                    axial_rotations_sample_num=2,
                )
        assert "Optimization attempt" in caplog.text
        assert "failed" in caplog.text
        assert result is not None
