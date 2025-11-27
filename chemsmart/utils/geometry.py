"""
Geometric calculations and molecular structure analysis utilities.

This module provides functions for geometric analysis of molecular structures,
including collinearity testing and moment of inertia calculations. These
utilities are essential for computational chemistry applications involving
molecular geometry optimization and thermochemical property calculations.

Key functionality includes:
- Collinearity detection for three-point systems
- Moment of inertia tensor calculations for molecular systems
- Principal axis determination for rotational analysis
- Molecular volume calculations using various methods
"""

import math

import numpy as np


def is_collinear(coords, tol=1e-5):
    """
    Check if three points are collinear using cross product method.

    Determines whether three points lie on the same straight line by
    calculating the cross product of vectors formed by the points.
    Used for identifying linear molecular configurations.

    Args:
        coords (array-like): List or array of three coordinate points,
            each containing [x, y, z] coordinates.
        tol (float, optional): Tolerance for collinearity test.
            Defaults to 1e-5.

    Returns:
        bool: True if points are collinear within tolerance, False otherwise.
    """
    a = coords[0]
    b = coords[1]
    c = coords[2]
    vec1 = np.array(b) - np.array(a)
    vec2 = np.array(c) - np.array(a)
    cross_product = np.cross(vec1, vec2)
    return (
        np.linalg.norm(cross_product) < tol
    )  # If cross product is close to zero, points are collinear


# def calculate_moments_of_inertia(mass, coords):
#     """Calculate the principal moments of inertia from mass and coordinates.
#     Mass parameter is a list of atomic masses corresponding to each coordinate in coords."""
#     # Convert inputs to numpy arrays
#     mass = np.array(mass)
#     coords = np.array(coords)
#
#     # Step 1: Compute the center of mass (COM)
#     center_of_mass = np.average(coords, axis=0, weights=mass)
#     shifted_coords = coords - center_of_mass  # Shifted coordinates
#
#     # Step 2: Compute the moment of inertia tensor
#     moi_tensor = np.zeros((3, 3))
#
#     for k in range(len(mass)):  # Loop over atoms
#         x, y, z = shifted_coords[k]  # Relative to COM
#         moi_tensor[0, 0] += mass[k] * (y ** 2 + z ** 2)  # Ixx
#         moi_tensor[1, 1] += mass[k] * (x ** 2 + z ** 2)  # Iyy
#         moi_tensor[2, 2] += mass[k] * (x ** 2 + y ** 2)  # Izz
#
#         moi_tensor[0, 1] -= mass[k] * x * y  # Ixy
#         moi_tensor[0, 2] -= mass[k] * x * z  # Ixz
#         moi_tensor[1, 2] -= mass[k] * y * z  # Iyz
#
#     # Since the inertia tensor is symmetric
#     moi_tensor[1, 0] = moi_tensor[0, 1]
#     moi_tensor[2, 0] = moi_tensor[0, 2]
#     moi_tensor[2, 1] = moi_tensor[1, 2]
#
#     # Step 3: Compute principal moments of inertia (eigenvalues)
#     moments_of_inertia_principal_axes = np.linalg.eigvalsh(moi_tensor)  # Sorted eigenvalues
#
#     return moments_of_inertia_principal_axes, moi_tensor

# def calculate_moments_of_inertia(mass, coords):
#     """Calculate the moment of inertia tensor and principal moments of inertia.
#
#     Parameters:
#     - mass (list or np.array): Atomic masses corresponding to each coordinate.
#     - coords (list or np.array): Nx3 array of atomic coordinates.
#
#     Returns:
#     - moments_of_inertia_principal_axes (np.array): Sorted eigenvalues of moi_tensor.
#     - moi_tensor (np.array): 3x3 moment of inertia tensor.
#     """
#     # Convert inputs to NumPy arrays
#     mass = np.array(mass).reshape(-1, 1)  # Ensure column vector for broadcasting
#     coords = np.array(coords)
#
#     # Compute center of mass (COM)
#     center_of_mass = np.average(coords, axis=0, weights=mass.flatten())
#     shifted_coords = coords - center_of_mass  # Shift coordinates to COM frame
#
#     # Compute squared distances (x², y², z²)
#     r2 = np.sum(shifted_coords ** 2, axis=1, keepdims=True)  # Shape (N,1)
#
#     # Compute moment of inertia tensor using vectorized operations
#     moi_tensor = np.diag(np.sum(mass * (r2 - shifted_coords ** 2), axis=0)) - \
#                  (mass * shifted_coords).T @ shifted_coords
#
#     print(f'moi_tensor: {moi_tensor}')
#
#     evals, evecs = np.linalg.eigh(moi_tensor)  # Eigenvalues and eigenvectors
#     # Instead of returning evecs.T as in ASE, we return evecs directly, as in Gaussian
#     return moi_tensor, evals, evecs


def calculate_moments_of_inertia(mass, coords):
    """
    Calculate the moment of inertia tensor and principal moments of inertia.

    Computes the moment of inertia tensor for a molecular system and
    determines the principal moments of inertia through eigenvalue
    decomposition. Essential for rotational analysis and thermochemical
    property calculations.

    Parameters
    ----------
    mass : array-like
        Atomic masses corresponding to each coordinate. Must have same
        length as coords array.
    coords : array-like
        Nx3 array of atomic coordinates in Cartesian space.

    Returns
    -------
    moi_tensor : np.ndarray
        3x3 moment of inertia tensor about the center of mass.
    evals : np.ndarray
        Principal moments of inertia (eigenvalues) sorted in ascending order.
    evecs : np.ndarray
        Principal axes (eigenvectors) as row vectors, transposed from
        NumPy's column format for compatibility with ASE conventions.
    """
    # Convert inputs to NumPy arrays
    mass = np.array(mass)
    coords = np.array(coords)

    # Step 1: Compute the center of mass (COM)
    center_of_mass = np.average(coords, axis=0, weights=mass)
    shifted_coords = coords - center_of_mass  # Shift coordinates to COM frame

    # Step 2: Compute the moment of inertia tensor
    moi_tensor = np.zeros((3, 3))

    for k in range(len(mass)):  # Loop over atoms
        x, y, z = shifted_coords[k]  # Relative to COM
        moi_tensor[0, 0] += mass[k] * (y**2 + z**2)  # Ixx
        moi_tensor[1, 1] += mass[k] * (x**2 + z**2)  # Iyy
        moi_tensor[2, 2] += mass[k] * (x**2 + y**2)  # Izz

        moi_tensor[0, 1] -= mass[k] * x * y  # Ixy
        moi_tensor[0, 2] -= mass[k] * x * z  # Ixz
        moi_tensor[1, 2] -= mass[k] * y * z  # Iyz

    # Since the inertia tensor is symmetric
    moi_tensor[1, 0] = moi_tensor[0, 1]
    moi_tensor[2, 0] = moi_tensor[0, 2]
    moi_tensor[2, 1] = moi_tensor[1, 2]

    # Step 3: Compute principal moments of inertia (eigenvalues)
    evals, evecs = np.linalg.eigh(moi_tensor)
    # Eigenvalues and eigenvectors
    # NumPy’s np.linalg.eigh returns the eigenvectors as columns in a matrix
    # so transpose to row vectors as in ASE
    return moi_tensor, evals, evecs.transpose()


def calculate_voronoi_dirichlet_occupied_volume(coords, radii, dispersion):
    """
    Estimate the occupied volume of a molecule using Voronoi-Dirichlet method,
    scaled by atomic radii to ensure a physically reasonable result.

    Parameters:
    - coords (list or np.array): Nx3 array of atomic coordinates.
    - radii (list or np.array): Atomic radii corresponding to each coordinate.
    - dispersion (float): Size of blocks for the Voronoi grid (used by pyvoro).

    Returns:
    - occupied_volume (float): Estimated physically occupied volume.

    Raises:
    - ImportError: If pyvoro is not installed. Install with: pip install pyvoro
        Note: pyvoro requires Python < 3.12
    """
    try:
        import pyvoro
    except ImportError:
        raise ImportError(
            "pyvoro is required for Voronoi-Dirichlet volume calculation. "
            "Install with: pip install chemsmart[voronoi]. "
            "Note: pyvoro requires Python < 3.12."
        )

    coords = np.array(coords)
    radii = np.array(radii)

    if len(coords) != len(radii):
        raise ValueError("Number of coordinates must match number of radii.")
    if coords.shape[1] != 3:
        raise ValueError("Coordinates must be 3D (Nx3 array).")

    padding = np.max(radii)
    box_min = np.min(coords, axis=0) - padding
    box_max = np.max(coords, axis=0) + padding
    limits = [[box_min[i], box_max[i]] for i in range(3)]

    try:
        cells = pyvoro.compute_voronoi(
            points=coords,
            limits=limits,
            dispersion=dispersion,
            radii=radii,
            periodic=[False, False, False],
        )
    except Exception as e:
        raise RuntimeError(f"Voronoi-Dirichlet tessellation failed: {e}")

    occupied_volume = 0.0
    for i, cell in enumerate(cells):
        cell_volume = cell["volume"]
        atomic_volume = (4 / 3) * np.pi * radii[i] ** 3

        # We cannot occupy more than the atomic volume
        occupied_volume += min(atomic_volume, cell_volume)

    return occupied_volume


def calculate_molecular_volume_vdp(coordinates, vdw_radii, dummy_points=True):
    """
    Calculate the molecular volume using Voronoi-Dirichlet polyhedra (VDP),
    approximating volume within van der Waals spheres by tetrahedral clipping.

    Parameters:
    coordinates : array-like, shape (n_atoms, 3)
    vdw_radii : array-like, shape (n_atoms,)
    dummy_points : bool

    Returns:
    float : Molecular volume in Å³
    """
    import numpy as np
    from scipy.spatial import Delaunay, Voronoi

    coordinates = np.array(coordinates)
    vdw_radii = np.array(vdw_radii)

    if coordinates.shape[0] != vdw_radii.shape[0]:
        raise ValueError("Mismatch between coordinates and radii.")
    if coordinates.shape[1] != 3:
        raise ValueError("Coordinates must be 3D.")

    original_coordinates = np.array(coordinates)
    original_radii = np.array(vdw_radii)
    n_atoms = original_coordinates.shape[0]
    coordinates = original_coordinates.copy()
    vdw_radii = original_radii.copy()

    # Add dummy points if needed
    if n_atoms < 4 and dummy_points:
        center = np.mean(coordinates, axis=0)
        dummy_distance = 10.0
        dummy_offsets = dummy_distance * np.array(
            [[1, 1, 1], [-1, -1, 1], [-1, 1, -1], [1, -1, -1]]
        )
        dummy_coords = center + dummy_offsets
        coordinates = np.vstack([coordinates, dummy_coords])
        vdw_radii = np.append(vdw_radii, [1e10] * 4)

    try:
        vor = Voronoi(coordinates)
    except Exception as e:
        raise RuntimeError(f"Voronoi tessellation failed: {e}")

    molecular_volume = 0.0
    n_original_atoms = n_atoms

    for i in range(n_original_atoms):
        region_idx = vor.point_region[i]
        region = vor.regions[region_idx]

        if -1 in region or not region:
            continue

        vertices = vor.vertices[region]
        if len(vertices) < 4:
            continue  # Cannot form a volume

        try:
            tri = Delaunay(vertices)
        except ValueError:
            continue  # Skip malformed regions

        atom_center = coordinates[i]
        vdw_radius = vdw_radii[i]
        vdw_radius2 = vdw_radius**2

        cell_volume = 0.0
        for simplex in tri.simplices:
            tetra = vertices[simplex]
            centroid = np.mean(tetra, axis=0)
            if np.sum((centroid - atom_center) ** 2) <= vdw_radius2:
                v1, v2, v3, v4 = tetra
                vol = abs(np.dot(v1 - v4, np.cross(v2 - v4, v3 - v4))) / 6.0
                cell_volume += vol

        molecular_volume += cell_volume

    return molecular_volume


def calculate_crude_occupied_volume(coords, radii):
    """
    Calculate the occupied volume of a molecule as the sum of atomic volumes.

    Parameters:
    - coords (list or np.array): Nx3 array of atomic coordinates.
    - radii (list or np.array): Atomic radii corresponding to each coordinate.

    Returns:
    - occupied_volume (float): Total occupied volume of the molecule.
    Ignores overlaps between atoms and gives an upper bound to true occupied volumes.
    """
    coords = np.array(coords)
    radii = np.array(radii)

    # Volume of a sphere: (4/3) * pi * r^3
    volumes = (4 / 3) * np.pi * np.power(radii, 3)
    occupied_volume = np.sum(volumes)

    return occupied_volume


def calculate_vdw_volume(coords, radii):
    """
    Calculate VDW volume from atomic coordinates and VDW radii.

    Parameters:
    - coords: List of [x, y, z] coordinates in Ångstroms for each atom.
    - radii: List of VDW radii in Ångstroms for each atom.

    Returns:
    - Volume in cubic Ångstroms (Å³).
    """
    if len(coords) != len(radii):
        raise ValueError("Number of coordinates must match number of radii")

    # Calculate individual sphere volumes
    volume = 0.0
    for radius in radii:
        volume += (4 / 3) * math.pi * (radius**3)

    # Overlap correction
    overlap_volume = 0.0
    for i in range(len(coords)):
        for j in range(i + 1, len(coords)):
            # Calculate distance between atoms i and j
            x1, y1, z1 = coords[i]
            x2, y2, z2 = coords[j]
            distance = math.sqrt(
                (x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2
            )
            r_i, r_j = radii[i], radii[j]
            sum_radii = r_i + r_j

            # Check for overlap
            if distance < sum_radii:
                # Approximate overlap volume using spherical cap formula
                overlap = (
                    math.pi
                    * (sum_radii - distance) ** 2
                    * (distance + 2 * sum_radii)
                ) / (12 * distance)
                overlap_volume += overlap

    return volume - overlap_volume
