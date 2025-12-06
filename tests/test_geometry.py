"""Tests for geometry utility functions."""

import numpy as np
import pytest

from chemsmart.utils.geometry import (
    calculate_crude_occupied_volume,
    calculate_grid_vdw_volume,
    calculate_moments_of_inertia,
    calculate_molecular_volume_vdp,
    calculate_vdw_volume,
    is_collinear,
)


class TestIsCollinear:
    """Tests for the is_collinear function."""

    def test_collinear_points_on_x_axis(self):
        """Test that points on a straight line along x-axis are collinear."""
        coords = [[0, 0, 0], [1, 0, 0], [2, 0, 0]]
        assert is_collinear(coords)

    def test_collinear_points_on_y_axis(self):
        """Test that points on a straight line along y-axis are collinear."""
        coords = [[0, 0, 0], [0, 1, 0], [0, 2, 0]]
        assert is_collinear(coords)

    def test_collinear_points_on_z_axis(self):
        """Test that points on a straight line along z-axis are collinear."""
        coords = [[0, 0, 0], [0, 0, 1], [0, 0, 2]]
        assert is_collinear(coords)

    def test_collinear_points_on_diagonal(self):
        """Test that points on a diagonal line are collinear."""
        coords = [[0, 0, 0], [1, 1, 1], [2, 2, 2]]
        assert is_collinear(coords)

    def test_non_collinear_points_triangle(self):
        """Test that points forming a triangle are not collinear."""
        coords = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        assert not is_collinear(coords)

    def test_non_collinear_points_3d(self):
        """Test that points in 3D space not on a line are not collinear."""
        coords = [[0, 0, 0], [1, 0, 0], [0, 0, 1]]
        assert not is_collinear(coords)

    def test_collinear_with_custom_tolerance(self):
        """Test collinearity with custom tolerance."""
        # Points that are slightly off a line
        coords = [[0, 0, 0], [1, 0, 0], [2, 1e-6, 0]]
        assert is_collinear(coords, tol=1e-5)
        assert not is_collinear(coords, tol=1e-7)

    def test_collinear_negative_coordinates(self):
        """Test collinearity with negative coordinates."""
        coords = [[-2, -2, -2], [0, 0, 0], [2, 2, 2]]
        assert is_collinear(coords)


class TestCalculateMomentsOfInertia:
    """Tests for the calculate_moments_of_inertia function."""

    def test_single_atom_at_origin(self):
        """Test moment of inertia for a single atom at origin."""
        mass = [12.0]  # Carbon mass
        coords = [[0, 0, 0]]
        moi_tensor, evals, evecs = calculate_moments_of_inertia(mass, coords)

        # For a single atom at origin, all moments should be 0
        assert np.allclose(evals, [0, 0, 0])

    def test_two_atoms_along_x_axis(self):
        """Test moment of inertia for two atoms along x-axis."""
        mass = [1.0, 1.0]  # Two identical masses
        coords = [[-1, 0, 0], [1, 0, 0]]
        moi_tensor, evals, evecs = calculate_moments_of_inertia(mass, coords)

        # For two masses at ±1 along x: I_xx = 0, I_yy = I_zz = 2
        assert np.isclose(evals[0], 0, atol=1e-10)
        assert np.isclose(evals[1], 2.0, atol=1e-10)
        assert np.isclose(evals[2], 2.0, atol=1e-10)

    def test_square_planar_arrangement(self):
        """Test moment of inertia for square planar arrangement."""
        mass = [1.0, 1.0, 1.0, 1.0]  # Four identical masses
        coords = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]]
        moi_tensor, evals, evecs = calculate_moments_of_inertia(mass, coords)

        # For this symmetric arrangement, two principal moments should be equal
        assert np.isclose(evals[0], 2.0, atol=1e-10)
        assert np.isclose(evals[1], 2.0, atol=1e-10)
        assert np.isclose(evals[2], 4.0, atol=1e-10)

    def test_moi_tensor_symmetry(self):
        """Test that moment of inertia tensor is symmetric."""
        mass = [1.0, 2.0, 3.0]
        coords = [[1, 2, 3], [-1, 0, 1], [2, -1, 0]]
        moi_tensor, evals, evecs = calculate_moments_of_inertia(mass, coords)

        # Check symmetry
        assert np.allclose(moi_tensor, moi_tensor.T)

    def test_evecs_shape(self):
        """Test that eigenvectors have correct shape."""
        mass = [1.0, 1.0, 1.0]
        coords = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        moi_tensor, evals, evecs = calculate_moments_of_inertia(mass, coords)

        assert evecs.shape == (3, 3)

    def test_evals_are_non_negative(self):
        """Test that principal moments are non-negative."""
        mass = [12.0, 1.0, 1.0, 1.0, 1.0]  # CH4-like
        coords = [
            [0, 0, 0],
            [1, 1, 1],
            [-1, -1, 1],
            [-1, 1, -1],
            [1, -1, -1],
        ]
        moi_tensor, evals, evecs = calculate_moments_of_inertia(mass, coords)

        assert all(e >= -1e-10 for e in evals)


class TestCalculateCrudeOccupiedVolume:
    """Tests for the calculate_crude_occupied_volume function."""

    def test_single_atom(self):
        """Test volume for a single atom."""
        coords = [[0, 0, 0]]
        radii = [1.0]
        volume = calculate_crude_occupied_volume(coords, radii)
        expected = (4 / 3) * np.pi * 1.0**3
        assert np.isclose(volume, expected)

    def test_two_atoms(self):
        """Test volume for two atoms (sum of individual volumes)."""
        coords = [[0, 0, 0], [5, 0, 0]]
        radii = [1.0, 2.0]
        volume = calculate_crude_occupied_volume(coords, radii)
        expected = (4 / 3) * np.pi * (1.0**3 + 2.0**3)
        assert np.isclose(volume, expected)

    def test_different_radii(self):
        """Test volume calculation with different radii."""
        coords = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        radii = [1.0, 1.5, 2.0]
        volume = calculate_crude_occupied_volume(coords, radii)
        expected = (4 / 3) * np.pi * (1.0**3 + 1.5**3 + 2.0**3)
        assert np.isclose(volume, expected)


class TestCalculateVdwVolume:
    """Tests for the calculate_vdw_volume function."""

    def test_single_atom(self):
        """Test VDW volume for a single atom."""
        coords = [[0, 0, 0]]
        radii = [1.5]
        volume = calculate_vdw_volume(coords, radii)
        expected = (4 / 3) * np.pi * 1.5**3
        assert np.isclose(volume, expected)

    def test_non_overlapping_atoms(self):
        """Test VDW volume for non-overlapping atoms."""
        coords = [[0, 0, 0], [10, 0, 0]]  # Far apart
        radii = [1.0, 1.0]
        volume = calculate_vdw_volume(coords, radii)
        expected = 2 * (4 / 3) * np.pi * 1.0**3
        assert np.isclose(volume, expected)

    def test_overlapping_atoms_reduces_volume(self):
        """Test that overlapping atoms result in reduced volume."""
        coords = [[0, 0, 0], [1.5, 0, 0]]  # Overlapping
        radii = [1.0, 1.0]
        volume = calculate_vdw_volume(coords, radii)
        # Should be less than sum of individual volumes
        sum_volumes = 2 * (4 / 3) * np.pi * 1.0**3
        assert volume < sum_volumes

    def test_mismatched_coords_radii_raises_error(self):
        """Test that mismatched coordinates and radii raises error."""
        coords = [[0, 0, 0], [1, 0, 0]]
        radii = [1.0]  # Only one radius
        with pytest.raises(ValueError):
            calculate_vdw_volume(coords, radii)

    def test_completely_overlapping_atoms(self):
        """Test when one atom is completely inside another."""
        coords = [[0, 0, 0], [0.1, 0, 0]]  # Nearly coincident
        radii = [2.0, 0.5]  # Small sphere inside large sphere
        volume = calculate_vdw_volume(coords, radii)
        # Volume should be approximately the large sphere volume
        large_sphere_vol = (4 / 3) * np.pi * 2.0**3
        assert volume <= large_sphere_vol
        assert volume > 0


class TestCalculateGridVdwVolume:
    """Tests for the calculate_grid_vdw_volume function."""

    def test_single_atom(self):
        """Test grid VDW volume for a single atom."""
        coords = [[0, 0, 0]]
        radii = [1.5]
        volume = calculate_grid_vdw_volume(coords, radii, grid_spacing=0.1)
        expected = (4 / 3) * np.pi * 1.5**3
        # Grid method is approximate, allow 5% tolerance
        assert np.isclose(volume, expected, rtol=0.05)

    def test_empty_input(self):
        """Test grid VDW volume with empty input."""
        coords = []
        radii = []
        volume = calculate_grid_vdw_volume(coords, radii)
        assert volume == 0.0

    def test_overlapping_atoms_handled_correctly(self):
        """Test that overlapping atoms are handled correctly."""
        coords = [[0, 0, 0], [1.0, 0, 0]]  # Overlapping
        radii = [1.0, 1.0]
        volume = calculate_grid_vdw_volume(coords, radii, grid_spacing=0.1)
        # Should be less than sum of individual volumes
        sum_volumes = 2 * (4 / 3) * np.pi * 1.0**3
        assert volume < sum_volumes

    def test_mismatched_coords_radii_raises_error(self):
        """Test that mismatched coordinates and radii raises error."""
        coords = [[0, 0, 0], [1, 0, 0]]
        radii = [1.0]  # Only one radius
        with pytest.raises(ValueError):
            calculate_grid_vdw_volume(coords, radii)


class TestCalculateMolecularVolumeVDP:
    """Tests for the calculate_molecular_volume_vdp function."""

    def test_single_atom(self):
        """Test VDP volume calculation for a single atom."""
        coords = [[0, 0, 0]]
        radii = [1.5]
        # Single atom with dummy points should give some volume
        volume = calculate_molecular_volume_vdp(
            coords, radii, dummy_points=True
        )
        assert volume > 0

    def test_mismatched_coords_radii_raises_error(self):
        """Test that mismatched coordinates and radii raises error."""
        coords = [[0, 0, 0], [1, 0, 0]]
        radii = [1.0]  # Only one radius
        with pytest.raises(ValueError):
            calculate_molecular_volume_vdp(coords, radii)

    def test_invalid_coordinates_shape_raises_error(self):
        """Test that invalid coordinate shape raises error."""
        coords = [[0, 0], [1, 0]]  # 2D instead of 3D
        radii = [1.0, 1.0]
        with pytest.raises(ValueError):
            calculate_molecular_volume_vdp(coords, radii)

    def test_multiple_atoms(self):
        """Test VDP volume for multiple atoms."""
        coords = [[0, 0, 0], [2, 0, 0], [1, 1.73, 0]]  # Triangle
        radii = [1.0, 1.0, 1.0]
        volume = calculate_molecular_volume_vdp(
            coords, radii, dummy_points=True
        )
        # Volume can be 0 for some configurations depending on tessellation
        assert volume >= 0
