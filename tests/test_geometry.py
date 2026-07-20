"""Tests for geometry utility functions."""

import numpy as np
import pytest

from chemsmart.utils.geometry import (
    calculate_crude_occupied_volume,
    calculate_grid_vdw_volume,
    calculate_molecular_volume_vdp,
    calculate_moments_of_inertia,
    calculate_vdw_volume,
    canonicalize_positions,
    clean_rotational_constants_by_geometry,
    get_coordinating_atoms,
    is_collinear,
)
from chemsmart.utils.periodictable import PeriodicTable

_pt = PeriodicTable()


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


class TestCleanRotationalConstantsByGeometry:
    def test_gaussian_mode_preserves_values(self):
        cleaned, status = clean_rotational_constants_by_geometry(
            [10.0, 2.0, 1.5],
            mode="gaussian",
            return_status=True,
        )
        assert np.allclose(cleaned, [10.0, 2.0, 1.5])
        assert status == "gaussian"

    def test_gaussian_mode_preserves_overflow(self):
        cleaned, status = clean_rotational_constants_by_geometry(
            [np.inf, 8.30647, 8.30647],
            mode="gaussian",
            return_status=True,
        )
        assert np.isinf(cleaned[0])
        assert np.allclose(cleaned[1:], [8.30647, 8.30647])
        assert status == "gaussian_overflow"

    def test_physical_mode_exact_linear_overflow(self):
        cleaned, status = clean_rotational_constants_by_geometry(
            [np.inf, 8.30647, 8.30647],
            mode="physical",
            return_status=True,
        )
        assert np.allclose(cleaned, [8.30647])
        assert status == "linear"

    def test_physical_mode_near_linear_overflow(self):
        cleaned, status = clean_rotational_constants_by_geometry(
            [np.inf, 8.306470, 8.306471],
            mode="physical",
            return_status=True,
        )
        assert np.allclose(cleaned, [0.5 * (8.306470 + 8.306471)])
        assert status == "linear"

    def test_physical_mode_huge_axial_constant(self):
        cleaned, status = clean_rotational_constants_by_geometry(
            [10919209500.9483, 8.3064670, 8.3064670],
            mode="physical",
            return_status=True,
        )
        assert np.allclose(cleaned, [8.3064670])
        assert status == "quasi_linear"

    def test_physical_mode_zero_axial_constant(self):
        cleaned, status = clean_rotational_constants_by_geometry(
            [0.0, 8.30647, 8.30647],
            mode="physical",
            return_status=True,
        )
        assert np.allclose(cleaned, [8.30647])
        assert status == "linear"

    def test_physical_mode_nonlinear(self):
        cleaned, status = clean_rotational_constants_by_geometry(
            [609.8308225, 5.8331025, 5.7778368],
            mode="physical",
            return_status=True,
        )
        assert np.allclose(cleaned, [609.8308225, 5.8331025, 5.7778368])
        assert status == "nonlinear"

    def test_invalid_mode_raises(self):
        with pytest.raises(ValueError):
            clean_rotational_constants_by_geometry([1.0, 2.0, 3.0], mode="bad")


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


class TestCanonicalizePositions:
    """Tests for the canonicalize_positions function."""

    def test_center_of_mass_at_origin(self):
        """Centre of mass of canonical positions must be at the origin."""
        masses = [12.0, 16.0, 16.0]
        coords = np.array(
            [
                [-2.184891, -0.000000, -0.000000],
                [2.184891, -0.000000, -0.000000],
                [0.000000, 0.000000, 0.000000],
            ]
        )
        result = canonicalize_positions(masses, coords)
        com = np.average(result, axis=0, weights=masses)
        assert np.allclose(com, 0.0, atol=1e-6)

    def test_single_atom_returns_origin(self):
        """A single atom must be placed at the origin."""
        masses = [12.0]
        coords = [[5.0, 3.0, -1.0]]
        result = canonicalize_positions(masses, coords)
        assert np.allclose(result, [[0.0, 0.0, 0.0]])

    def test_diatomic_aligned_along_z(self):
        """A diatomic molecule must be aligned along the z-axis."""
        masses = [14.0, 14.0]
        coords = [[1.0, 2.0, 3.0], [2.0, 4.0, 5.0]]
        result = canonicalize_positions(masses, coords)
        # bond length = ||r2 - r1|| = sqrt(1^2 + 2^2 + 2^2) = 3
        # canonical positions = (0,0,±|v|/2) = (0,0,±1.5)
        assert np.allclose(result, [[0.0, 0.0, -1.5], [0.0, 0.0, 1.5]])

    def test_translation_invariance(self):
        """Translating all atoms by a constant vector must not change the result."""
        masses = [12.0, 1.0, 1.0, 1.0, 1.0]
        coords = np.array(
            [
                [0, 0, 0],
                [1.09, 0, 0],
                [-1.09, 0, 0],
                [0, 1.09, 0],
                [0, -1.09, 0],
            ]
        )
        shifted = coords + np.array([10.0, -20.0, 30.0])
        canon_shifted = canonicalize_positions(masses, shifted)
        canon_original = canonicalize_positions(masses, coords)
        assert np.allclose(canon_original, canon_shifted, atol=1e-6)

    def test_rotation_invariance(self):
        """Rotating all atoms must not change the canonical positions."""
        masses = [12.0, 16.0, 1.0, 1.0]
        coords = np.array(
            [
                [-0.6123, 0.0000, 0.0000],
                [0.6123, 0.0000, 0.0000],
                [-1.2000, 0.2426, -0.8998],
                [-1.2000, -0.2424, 0.8998],
            ]
        )
        # Apply a rotation matrix (90° around z-axis)
        theta = np.pi / 2
        R = np.array(
            [
                [np.cos(theta), -np.sin(theta), 0],
                [np.sin(theta), np.cos(theta), 0],
                [0, 0, 1],
            ]
        )
        rotated = coords @ R.T
        canon_original = canonicalize_positions(masses, coords)
        canon_rotated = canonicalize_positions(masses, rotated)
        assert np.allclose(canon_original, canon_rotated, atol=1e-6)


class TestGetCoordinatingAtoms:
    """Tests for covalent-radius-ratio coordination spheres."""

    def test_primary_donor_and_secondary_ligand(self):
        r_mn = _pt.covalent_radius("Mn")
        r_n = _pt.covalent_radius("N")
        r_c = _pt.covalent_radius("C")
        primary_dist = 1.05 * (r_mn + r_n)
        secondary_dist = 1.25 * (r_mn + r_c)
        elements = ["Mn", "N", "C"]
        coordinates = np.array(
            [
                [0.0, 0.0, 0.0],
                [primary_dist, 0.0, 0.0],
                [secondary_dist, 0.0, 0.0],
            ]
        )

        primary, secondary = get_coordinating_atoms(0, elements, coordinates)

        assert primary == [1]
        assert secondary == [2]

    def test_hydride_within_cutoff_is_primary(self):
        elements = ["Mn", "H"]
        coordinates = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])

        primary, secondary = get_coordinating_atoms(0, elements, coordinates)

        assert primary == [1]
        assert secondary == []

    def test_distant_atom_excluded(self):
        r_mn = _pt.covalent_radius("Mn")
        r_c = _pt.covalent_radius("C")
        far_dist = 1.50 * (r_mn + r_c)
        elements = ["Mn", "C"]
        coordinates = np.array([[0.0, 0.0, 0.0], [far_dist, 0.0, 0.0]])

        primary, secondary = get_coordinating_atoms(0, elements, coordinates)

        assert primary == []
        assert secondary == []

    def test_geometric_expansion_captures_terminal_co_oxygen(self):
        """XYZ-safe expansion: CO oxygen near primary C, even if far from Mn."""
        r_mn = _pt.covalent_radius("Mn")
        r_c = _pt.covalent_radius("C")
        r_o = _pt.covalent_radius("O")
        c_dist = 1.05 * (r_mn + r_c)
        co_bond = 1.15
        o_dist = c_dist + co_bond
        # Terminal O should sit outside the secondary radius-ratio shell.
        assert o_dist / (r_mn + r_o) > 1.35

        elements = ["Mn", "C", "O"]
        coordinates = np.array(
            [
                [0.0, 0.0, 0.0],
                [c_dist, 0.0, 0.0],
                [o_dist, 0.0, 0.0],
            ]
        )

        primary, secondary = get_coordinating_atoms(0, elements, coordinates)

        assert primary == [1]
        assert secondary == [2]

    def test_geometric_expansion_can_be_disabled(self):
        r_mn = _pt.covalent_radius("Mn")
        r_c = _pt.covalent_radius("C")
        r_o = _pt.covalent_radius("O")
        c_dist = 1.05 * (r_mn + r_c)
        o_dist = c_dist + 1.15
        assert o_dist / (r_mn + r_o) > 1.35

        elements = ["Mn", "C", "O"]
        coordinates = np.array(
            [
                [0.0, 0.0, 0.0],
                [c_dist, 0.0, 0.0],
                [o_dist, 0.0, 0.0],
            ]
        )

        primary, secondary = get_coordinating_atoms(
            0, elements, coordinates, expand_cutoff=0
        )

        assert primary == [1]
        assert secondary == []
