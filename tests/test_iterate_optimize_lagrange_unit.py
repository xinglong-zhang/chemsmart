"""
Direct unit tests for IterateAnalyzer._optimize_lagrange (via the public
_find_optimal_position dispatcher), the core SLSQP-based substituent
placement optimization in chemsmart.jobs.iterate.iterate.

Previously entirely untested because the full pipeline was assumed to
be "inherently slow" -- but with small sphere_direction_samples_num /
axial_rotations_sample_num values (as used throughout this file), a
full optimization run completes in well under 100ms, making direct
testing (including error/fallback branches, via mocking
scipy.optimize.minimize) practical.
"""

from unittest.mock import patch

import numpy as np
from scipy.optimize import minimize as real_minimize

from chemsmart.jobs.iterate.iterate import IterateAnalyzer

# A simple non-linear skeleton (like a bent triatomic) with the link
# atom (index 0) off-center, so skeleton_com != skeleton_link_coords.
SKELETON = np.array(
    [
        [6.0, 0.0, 0.0, 0.0],
        [6.0, 1.5, 0.0, 0.0],
        [1.0, -0.5, 0.9, 0.0],
        [1.0, -0.5, -0.9, 0.0],
    ]
)

# Single-atom skeleton: skeleton_com coincides exactly with the (only)
# link atom, used to force the dir_norm < 1e-6 fallback branch.
SINGLE_ATOM_SKELETON = np.array([[6.0, 0.0, 0.0, 0.0]])


class TestOptimizeLagrangeHappyPath:
    def test_single_atom_substituent_returns_valid_position(self):
        # A 1-atom substituent has zero relative offset from its own
        # link atom, so sub_centroid_norm < 1e-3 (fallback principal
        # axis branch), and the sub link atom starts away from the
        # skeleton link (initial_link_dist >= 1e-6, normal branch).
        sub = np.array([[1.0, 5.0, 5.0, 5.0]])
        result = IterateAnalyzer._find_optimal_position(
            SKELETON,
            sub,
            skeleton_link_index=0,
            sub_link_index=0,
            sphere_direction_samples_num=4,
            axial_rotations_sample_num=2,
        )
        assert result is not None
        assert result.shape == sub.shape
        assert result[0, 0] == 1.0  # atomic number preserved

    def test_multi_atom_substituent_with_nonzero_centroid(self):
        # A 2-atom substituent gives a non-trivial relative offset,
        # exercising the sub_centroid_norm >= 1e-3 branch.
        sub = np.array([[1.0, 5.0, 5.0, 5.0], [1.0, 5.5, 5.5, 5.5]])
        result = IterateAnalyzer._find_optimal_position(
            SKELETON,
            sub,
            skeleton_link_index=0,
            sub_link_index=0,
            sphere_direction_samples_num=4,
            axial_rotations_sample_num=2,
        )
        assert result is not None
        assert result.shape == sub.shape

    def test_sub_link_starting_at_skeleton_link_position(self):
        # sub_coord's link atom placed exactly at skeleton_link_coords
        # forces initial_link_dist < 1e-6, taking the "default
        # direction" branch (skeleton has multiple atoms not centered
        # on the link, so dir_norm >= 1e-6 here).
        sub = np.array([[1.0, 0.0, 0.0, 0.0], [1.0, 0.5, 0.5, 0.5]])
        result = IterateAnalyzer._find_optimal_position(
            SKELETON,
            sub,
            skeleton_link_index=0,
            sub_link_index=0,
            sphere_direction_samples_num=4,
            axial_rotations_sample_num=2,
        )
        assert result is not None

    def test_single_atom_skeleton_forces_dir_norm_fallback(self):
        # A single-atom skeleton means skeleton_com == skeleton_link,
        # so `direction` has near-zero norm, forcing the [1, 0, 0]
        # fallback direction branch (in addition to initial_link_dist
        # < 1e-6, since the sub link also starts at the same point).
        sub = np.array([[1.0, 0.0, 0.0, 0.0]])
        result = IterateAnalyzer._find_optimal_position(
            SINGLE_ATOM_SKELETON,
            sub,
            skeleton_link_index=0,
            sub_link_index=0,
            sphere_direction_samples_num=4,
            axial_rotations_sample_num=2,
        )
        assert result is not None


class TestOptimizeLagrangePreOptimizationFallback:
    def _make_side_effect(self, first_call_result):
        """Return a minimize() side_effect: the first call returns
        first_call_result (a canned/failing outcome); every subsequent
        call delegates to the real scipy minimize."""
        call_count = {"n": 0}

        def _side_effect(*args, **kwargs):
            call_count["n"] += 1
            if call_count["n"] == 1:
                return first_call_result
            return real_minimize(*args, **kwargs)

        return _side_effect

    def test_pre_optimization_failure_falls_back_to_projection(self, caplog):
        class _FakeResult:
            success = False
            message = "did not converge"

        sub = np.array([[1.0, 5.0, 5.0, 5.0]])
        with patch(
            "chemsmart.jobs.iterate.iterate.minimize",
            side_effect=self._make_side_effect(_FakeResult()),
        ):
            with caplog.at_level("WARNING"):
                result = IterateAnalyzer._find_optimal_position(
                    SKELETON,
                    sub,
                    skeleton_link_index=0,
                    sub_link_index=0,
                    sphere_direction_samples_num=4,
                    axial_rotations_sample_num=2,
                )
        assert "Initial position optimization failed" in caplog.text
        assert result is not None

    def test_pre_optimization_exception_falls_back_to_projection(self, caplog):
        def _combined_side_effect(*args, **kwargs):
            _combined_side_effect.calls += 1
            if _combined_side_effect.calls == 1:
                raise RuntimeError("boom")
            return real_minimize(*args, **kwargs)

        _combined_side_effect.calls = 0

        sub = np.array([[1.0, 5.0, 5.0, 5.0]])
        with patch(
            "chemsmart.jobs.iterate.iterate.minimize",
            side_effect=_combined_side_effect,
        ):
            with caplog.at_level("WARNING"):
                result = IterateAnalyzer._find_optimal_position(
                    SKELETON,
                    sub,
                    skeleton_link_index=0,
                    sub_link_index=0,
                    sphere_direction_samples_num=4,
                    axial_rotations_sample_num=2,
                )
        assert "Initial position optimization error" in caplog.text
        assert result is not None


class TestOptimizeLagrangeMainLoopFailures:
    def test_all_attempts_fail_returns_none(self, caplog):
        class _FakeResult:
            success = False
            message = "did not converge"
            fun = float("inf")
            x = np.zeros(6)

        def _side_effect(*args, **kwargs):
            # Let the pre-optimization (first call) succeed via the
            # real solver so we reach the main loop; fail everything
            # after that.
            if _side_effect.calls == 0:
                _side_effect.calls += 1
                return real_minimize(*args, **kwargs)
            return _FakeResult()

        _side_effect.calls = 0

        sub = np.array([[1.0, 5.0, 5.0, 5.0]])
        with patch(
            "chemsmart.jobs.iterate.iterate.minimize",
            side_effect=_side_effect,
        ):
            with caplog.at_level("ERROR"):
                result = IterateAnalyzer._find_optimal_position(
                    SKELETON,
                    sub,
                    skeleton_link_index=0,
                    sub_link_index=0,
                    sphere_direction_samples_num=4,
                    axial_rotations_sample_num=2,
                )
        assert result is None
        assert "All optimization attempts failed" in caplog.text

    def test_exception_in_main_loop_is_logged_and_skipped(self, caplog):
        def _side_effect(*args, **kwargs):
            _side_effect.calls += 1
            # First call: real pre-optimization.
            if _side_effect.calls == 1:
                return real_minimize(*args, **kwargs)
            # Second call (first main-loop attempt): raise.
            if _side_effect.calls == 2:
                raise RuntimeError("solver crashed")
            # Remaining calls: real solver, so the run still succeeds
            # overall.
            return real_minimize(*args, **kwargs)

        _side_effect.calls = 0

        sub = np.array([[1.0, 5.0, 5.0, 5.0]])
        with patch(
            "chemsmart.jobs.iterate.iterate.minimize",
            side_effect=_side_effect,
        ):
            with caplog.at_level("WARNING"):
                result = IterateAnalyzer._find_optimal_position(
                    SKELETON,
                    sub,
                    skeleton_link_index=0,
                    sub_link_index=0,
                    sphere_direction_samples_num=4,
                    axial_rotations_sample_num=2,
                )
        assert "Optimization attempt" in caplog.text
        assert "failed" in caplog.text
        assert result is not None
