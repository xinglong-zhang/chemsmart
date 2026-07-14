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

import logging
import math

import numpy as np

from chemsmart.utils.periodictable import PeriodicTable

logger = logging.getLogger(__name__)
_pt = PeriodicTable()


def get_coordinating_atoms(
    metal_index,
    elements,
    coordinates,
    tau_primary=1.15,
    tau_secondary=1.35,
    expand_cutoff=1.6,
):
    """Categorize atoms coordinating to a metal by covalent-radius ratio.

    Direct ligands are assigned via metal–atom covalent-radius ratios.
    For XYZ inputs without connectivity, a purely geometric expansion then
    pulls in covalently bound partner atoms of those direct ligands (for
    example the O of CO, or H of H2O) that lie within ``expand_cutoff`` of
    any primary-sphere atom. No topology / residue metadata is used.

    Parameters
    ----------
    metal_index : int
        0-based index of the metal center.
    elements : sequence of str
        Element symbols for all atoms.
    coordinates : array-like
        Cartesian coordinates, shape ``(n_atoms, 3)``.
    tau_primary, tau_secondary : float
        Radius-ratio thresholds for the primary and secondary shells.
    expand_cutoff : float
        Distance (Å) used to expand from direct ligands to their bound
        partners. Set to ``0`` or ``None`` to disable expansion.

    Returns
    -------
    tuple[list[int], list[int]]
        ``(primary_sphere, secondary_sphere)`` as 0-based atom indices.
        Expanded partner atoms are placed in ``secondary_sphere`` when they
        are not already primary.
    """
    primary_sphere = []
    secondary_sphere = []

    coordinates = np.asarray(coordinates, dtype=float)
    xyz_m = coordinates[metal_index]
    element_m = elements[metal_index]
    r_m = _pt.covalent_radius(element_m)

    distances = np.linalg.norm(coordinates - xyz_m, axis=1)

    for idx, (element, dist) in enumerate(zip(elements, distances)):
        if idx == metal_index:
            continue

        # Hydride / hydrogen exception: tight catalytic hydrides.
        if element == "H" and dist <= 1.8:
            primary_sphere.append(idx)
            continue

        r_a = _pt.covalent_radius(element)
        ratio = dist / (r_m + r_a)

        if ratio <= tau_primary:
            primary_sphere.append(idx)
        elif ratio <= tau_secondary:
            secondary_sphere.append(idx)

    # Geometric partner expansion (XYZ-safe; no byres / neighbor topology).
    # Captures the second atom of diatomic ligands and similar small molecules.
    if primary_sphere and expand_cutoff:
        primary_set = set(primary_sphere)
        secondary_set = set(secondary_sphere)
        primary_coords = coordinates[primary_sphere]
        for idx in range(len(elements)):
            if (
                idx == metal_index
                or idx in primary_set
                or idx in secondary_set
            ):
                continue
            partner_dists = np.linalg.norm(
                primary_coords - coordinates[idx], axis=1
            )
            if np.min(partner_dists) <= expand_cutoff:
                secondary_sphere.append(idx)
                secondary_set.add(idx)

    return primary_sphere, secondary_sphere


def is_collinear(coords, tol=1e-2):
    """
    Check if three points are collinear using cross product method.

    Determines whether three points lie on the same straight line by
    calculating the cross product of vectors formed by the points.
    Used for identifying linear molecular configurations.

    Args:
        coords (array-like): List or array of three coordinate points,
            each containing [x, y, z] coordinates.
        tol (float, optional): Tolerance for collinearity test.
            Defaults to 1e-2.

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
# Mass parameter is a list of atomic masses
# corresponding to each coordinate in coords."""
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
# moments_of_inertia_principal_axes =
# np.linalg.eigvalsh(moi_tensor) # Sorted eigenvalues
#
#     return moments_of_inertia_principal_axes, moi_tensor

# def calculate_moments_of_inertia(mass, coords):
# """Calculate the moment of inertia
# tensor and principal moments of inertia.
#
#     Parameters:
# - mass (list or np.array): Atomic
# masses corresponding to each coordinate.
#     - coords (list or np.array): Nx3 array of atomic coordinates.
#
#     Returns:
# - moments_of_inertia_principal_axes
# (np.array): Sorted eigenvalues of moi_tensor.
#     - moi_tensor (np.array): 3x3 moment of inertia tensor.
#     """
#     # Convert inputs to NumPy arrays
# mass = np.array(mass).reshape(-1, 1) #
# Ensure column vector for broadcasting
#     coords = np.array(coords)
#
#     # Compute center of mass (COM)
#     center_of_mass = np.average(coords, axis=0, weights=mass.flatten())
# shifted_coords = coords - center_of_mass
# # Shift coordinates to COM frame
#
#     # Compute squared distances (x², y², z²)
#     r2 = np.sum(shifted_coords ** 2, axis=1, keepdims=True)  # Shape (N,1)
#
#     # Compute moment of inertia tensor using vectorized operations
# moi_tensor = np.diag(np.sum(mass * (r2
# - shifted_coords ** 2), axis=0)) - \
#                  (mass * shifted_coords).T @ shifted_coords
#
#     print(f'moi_tensor: {moi_tensor}')
#
#     evals, evecs = np.linalg.eigh(moi_tensor)  # Eigenvalues and eigenvectors
# # Instead of returning evecs.T as in ASE,
# we return evecs directly, as in Gaussian
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


def canonicalize_positions(masses, coords, moment_tol=1e-6):
    """Compute deterministic canonical atomic positions.

    The canonical frame is translation- and rotation-invariant:

    Algorithm
    ---------
    1. Translate the coordinates so that the centre of mass is at the origin.
    2. Rotate coordinates into the principal-axes frame obtained from the
       moment-of-inertia tensor (eigenvalues sorted in ascending order).
    3. Apply a deterministic sign convention to each principal axis so that the
       result is unique.
    4. Enforce a right-handed coordinate system.

    Notes
    -----
    - This is a practical heuristic for structure normalization, not a
      mathematically perfect canonicalization.
    - Near-degenerate or highly symmetric cases may remain numerically
      sensitive because the principal axes themselves may not be uniquely
      defined.
    - In the centre-of-mass frame, the first mass-weighted coordinate moment
      along each axis is theoretically zero, so in practice the third and
      fifth moments are the main sign-disambiguation criteria.

    Parameters
    ----------
    masses : array-like
        Atomic masses, shape (N,).
    coords : array-like
        Cartesian coordinates, shape (N, 3).
    moment_tol : float, optional
        Tolerance for treating a mass-weighted moment as zero (default 1e-6).

    Returns
    -------
    canonical : np.ndarray
        Canonical coordinates, shape (N, 3).
    """
    masses = np.asarray(masses, dtype=float)
    coords = np.asarray(coords, dtype=float)
    n_atoms = len(masses)

    # --- Step 1: translate to centre of mass ----------------------------
    com = np.average(coords, axis=0, weights=masses)
    shifted = coords - com

    # --- Edge case: single atom -----------------------------------------
    if n_atoms == 1:
        return shifted

    # --- Edge case: diatomic → align along z-axis -----------------------
    if n_atoms == 2:
        vec = shifted[1] - shifted[0]
        norm = np.linalg.norm(vec)
        if norm < 1e-14:
            return shifted
        z_hat = vec / norm
        # Build a right-handed frame with z_hat as the third axis.
        trial = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(z_hat, trial)) > 0.9:
            trial = np.array([0.0, 1.0, 0.0])
        x_hat = np.cross(z_hat, trial)
        x_hat /= np.linalg.norm(x_hat)
        y_hat = np.cross(z_hat, x_hat)
        y_hat /= np.linalg.norm(y_hat)
        R = np.vstack(
            [x_hat, y_hat, z_hat]
        )  # rotation matrix (rows = new axes)
        rotated = shifted @ R.T

        # Deterministic z-axis sign convention
        if masses[0] > masses[1]:
            rotated[:, 2] *= -1
        elif masses[0] == masses[1]:
            if rotated[0, 2] > 0:
                rotated[:, 2] *= -1
        return rotated

    # --- Step 2: compute principal-axes frame ---------------------------
    moi_tensor, evals, evecs_rows = calculate_moments_of_inertia(
        masses, coords
    )

    # evecs_rows has shape (3, 3), each row is a principal axis
    rotated = shifted @ evecs_rows.T

    # --- Step 3: deterministic sign convention --------------------------
    axis_scores = []
    axis_signs = np.ones(3, dtype=float)
    for i in range(3):
        # In COM frame this is theoretically zero, but we keep it as the
        # first numerical check.
        first_moment = np.dot(masses, rotated[:, i])

        if abs(first_moment) > moment_tol:
            if first_moment < 0:
                rotated[:, i] *= -1
                axis_signs[i] = -1
            score = abs(first_moment)
        else:
            third_moment = np.dot(masses, rotated[:, i] ** 3)
            if abs(third_moment) > moment_tol:
                if third_moment < 0:
                    rotated[:, i] *= -1
                    axis_signs[i] = -1
                score = abs(third_moment)
            else:
                fifth_moment = np.dot(masses, rotated[:, i] ** 5)
                if fifth_moment < -moment_tol:
                    rotated[:, i] *= -1
                    axis_signs[i] = -1
                score = abs(fifth_moment)

        axis_scores.append(score)

    # --- Step 4: ensure right-handed coordinate system ------------------
    # effective determinant after all sign flips
    effective_det = np.linalg.det(evecs_rows) * np.prod(axis_signs)
    if effective_det < 0:
        idx_flip = int(np.argmin(axis_scores))
        rotated[:, idx_flip] *= -1

    return rotated


def calculate_voronoi_dirichlet_occupied_volume(
    coords, radii, dispersion=None
):
    """
    Estimate the occupied volume of a molecule using Voronoi-Dirichlet method,
    scaled by atomic radii to ensure a physically reasonable result.

    Uses scipy for Voronoi tessellation. Mirror images of all atoms are added
    across the faces, edges, and corners of the bounding box to ensure all
    Voronoi cells are bounded. For each atom, the occupied contribution is the
    minimum of its atomic sphere volume and its Voronoi cell volume.

    Parameters:
    - coords (list or np.array): Nx3 array of atomic coordinates.
    - radii (list or np.array): Atomic radii corresponding to each coordinate.
    - dispersion (float, optional): Padding extent used to construct the
      bounding box for Voronoi tessellation. When provided, this sets how far
      the mirrored bounding box extends beyond the molecule; when omitted or
      None, the maximum atomic radius is used as the padding. Larger values
      produce bigger cells for edge atoms, which can reduce their contribution.

    Returns:
    - occupied_volume (float): Estimated physically occupied volume.
      Atoms whose Voronoi cells remain unbounded after mirroring are excluded
      from the sum; in practice this should not occur for typical molecules.
    """
    from scipy.spatial import ConvexHull, Voronoi

    coords = np.array(coords)
    radii = np.array(radii)

    if len(coords) != len(radii):
        raise ValueError("Number of coordinates must match number of radii.")
    if coords.shape[1] != 3:
        raise ValueError("Coordinates must be 3D (Nx3 array).")

    # Use dispersion as bounding box padding when provided; fall back to the
    # maximum atomic radius for a tight, physically meaningful box
    padding = dispersion if dispersion is not None else np.max(radii)
    box_min = np.min(coords, axis=0) - padding
    box_max = np.max(coords, axis=0) + padding

    # Add mirror images of all atoms across each face, edge, and corner of the
    # bounding box to ensure all Voronoi cells are bounded within the box
    all_coords = [coords]
    for di in [-1, 0, 1]:
        for dj in [-1, 0, 1]:
            for dk in [-1, 0, 1]:
                if di == 0 and dj == 0 and dk == 0:
                    continue
                mirror = coords.copy()
                for axis, (offset, box_lo, box_hi) in enumerate(
                    zip([di, dj, dk], box_min, box_max)
                ):
                    if offset == -1:
                        mirror[:, axis] = 2 * box_lo - mirror[:, axis]
                    elif offset == 1:
                        mirror[:, axis] = 2 * box_hi - mirror[:, axis]
                all_coords.append(mirror)

    all_points = np.vstack(all_coords)

    try:
        vor = Voronoi(all_points)
    except Exception as e:
        raise RuntimeError(f"Voronoi-Dirichlet tessellation failed: {e}")

    n = len(coords)
    occupied_volume = 0.0
    for i in range(n):
        region_idx = vor.point_region[i]
        region = vor.regions[region_idx]

        # Skip atoms whose cell is still open (contains vertex at infinity)
        if -1 in region or not region:
            continue

        vertices = vor.vertices[region]
        try:
            hull = ConvexHull(vertices)
            cell_volume = hull.volume
        except Exception:
            continue

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
    Ignores overlaps between atoms and gives
    an upper bound to true occupied volumes.
    """
    coords = np.array(coords)
    radii = np.array(radii)

    # Volume of a sphere: (4/3) * pi * r^3
    volumes = (4 / 3) * np.pi * np.power(radii, 3)
    occupied_volume = np.sum(volumes)

    return occupied_volume


def calculate_vdw_volume(coords, radii):
    """
    Calculate VDW volume from atomic coordinates and VDW radii using pairwise
    overlap correction.

    Uses the sum of spherical volumes minus pairwise overlap corrections.
    The overlap volume between two spheres with radii r_i and r_j at distance d
    is computed using the exact lens-shaped intersection formula derived from
    spherical cap volumes:

        V_overlap = π(r_i + r_j - d)² × [d² +
        2d(r_i + r_j) - 3(r_i - r_j)²] / (12d)

    Parameters
    ----------
    coords : list or array-like
        List of [x, y, z] coordinates in Ångstroms for each atom.
    radii : list or array-like
        List of VDW radii in Ångstroms for each atom.

    Returns
    -------
    float
        Volume in cubic Ångstroms (Å³).

    Limitations
    -----------
    This method only corrects for pairwise (two-body) overlaps and does not
    account for higher-order overlaps (three or more spheres overlapping at
    the same point). For densely packed molecules, this leads to overcounting
    of the overlap correction, resulting in underestimated volumes.

    For example, if spheres A, B, and C all overlap at a common region:
    - The pairwise method subtracts the A-B, B-C, and A-C overlaps
    - But the triple-overlap region (where all three meet) gets subtracted
      three times instead of once
    - This causes the calculated volume to be smaller than the true volume

    For more accurate volume calculations on complex molecules, consider using
    ``calculate_grid_vdw_volume()`` which uses numerical grid integration
    and correctly handles all overlap orders.

    See Also
    --------
    calculate_grid_vdw_volume : Grid-based volume calculation (more accurate)
    calculate_crude_occupied_volume : Sum of spheres without overlap correction
    """
    if len(coords) != len(radii):
        raise ValueError("Number of coordinates must match number of radii")

    # Calculate individual sphere volumes
    volume = 0.0
    for radius in radii:
        volume += (4 / 3) * math.pi * (radius**3)

    # Overlap correction using the exact lens-shaped intersection formula
    overlap_volume = 0.0
    for i in range(len(coords)):
        for j in range(i + 1, len(coords)):
            # Calculate distance between atoms i and j
            x1, y1, z1 = coords[i]
            x2, y2, z2 = coords[j]
            d = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
            r_i, r_j = radii[i], radii[j]
            sum_radii = r_i + r_j
            diff_radii = r_i - r_j

            # Check for overlap
            if d < sum_radii and d > 0:
                # Check if one sphere is completely inside the other
                if d <= abs(diff_radii):
                    # Smaller sphere is completely inside the larger
                    overlap = (4 / 3) * math.pi * min(r_i, r_j) ** 3
                else:
                    # Partial overlap: use exact
                    # lens-shaped intersection formula
                    # V = π(r_i + r_j - d)² × [d² + 2d(r_i
                    # + r_j) - 3(r_i - r_j)²] / (12d)
                    overlap = (
                        math.pi
                        * (sum_radii - d) ** 2
                        * (d**2 + 2 * d * sum_radii - 3 * diff_radii**2)
                    ) / (12 * d)
                overlap_volume += overlap

    return volume - overlap_volume


def calculate_grid_vdw_volume(coords, radii, grid_spacing=0.2):
    """
    Calculate VDW volume using grid-based numerical integration.

    This method places the molecule in a 3D grid and counts grid points
    that fall inside any atomic VDW sphere. This approach correctly handles
    all orders of atomic overlaps (pairwise, triple, etc.) and provides
    more accurate volume estimates for complex molecules.

    This implementation is similar in concept
    to RDKit's DoubleCubicLatticeVolume
    algorithm, which also uses grid-based integration for molecular volume
    calculation.

    Parameters
    ----------
    coords : list or array-like
        List of [x, y, z] coordinates in Ångstroms for each atom.
    radii : list or array-like
        List of VDW radii in Ångstroms for each atom.
    grid_spacing : float, optional
        Spacing between grid points in Ångstroms. Smaller values give more
        accurate results but increase computation time. Default is 0.2 Å,
        which provides a good balance between accuracy and speed.

    Returns
    -------
    float
        Volume in cubic Ångstroms (Å³).

    Notes
    -----
    The algorithm works as follows:
    1. Create a bounding box around the molecule with padding equal to the
       maximum VDW radius
    2. Generate a regular 3D grid of points within the bounding box
    3. For each grid point, check if it falls inside any atomic VDW sphere
    4. Count the number of points inside and multiply by the volume per
       grid point (grid_spacing³)

    The accuracy depends on grid_spacing:
    - 0.2 Å: Good balance of accuracy and speed (default)
    - 0.1 Å: Higher accuracy, ~8x slower
    - 0.5 Å: Faster but less accurate

    See Also
    --------
    calculate_vdw_volume : Analytical pairwise overlap method (faster but
        less accurate for complex molecules)

    Examples
    --------
    >>> coords = [[0, 0, 0], [1.5, 0, 0]]
    >>> radii = [1.7, 1.2]
    >>> vol = calculate_grid_vdw_volume(coords, radii)
    """
    coords = np.array(coords)
    radii = np.array(radii)

    if len(coords) != len(radii):
        raise ValueError("Number of coordinates must match number of radii")

    if len(coords) == 0:
        return 0.0

    # Determine bounding box with padding
    max_radius = np.max(radii)
    min_coords = np.min(coords, axis=0) - max_radius
    max_coords = np.max(coords, axis=0) + max_radius

    # Generate grid points using meshgrid for vectorized operations
    x_range = np.arange(
        min_coords[0], max_coords[0] + grid_spacing, grid_spacing
    )
    y_range = np.arange(
        min_coords[1], max_coords[1] + grid_spacing, grid_spacing
    )
    z_range = np.arange(
        min_coords[2], max_coords[2] + grid_spacing, grid_spacing
    )

    # Create 3D meshgrid and reshape to (N, 3) array of all grid points
    xx, yy, zz = np.meshgrid(x_range, y_range, z_range, indexing="ij")
    grid_points = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])

    # For each atom, compute squared distances from all grid points
    # and check if any grid point is inside the atom's VDW sphere
    radii_sq = radii**2
    inside_any = np.zeros(len(grid_points), dtype=bool)

    for i in range(len(coords)):
        # Squared distances from atom i to all grid points
        diff = grid_points - coords[i]
        dist_sq = np.sum(diff**2, axis=1)
        inside_any |= dist_sq <= radii_sq[i]

    # Count points inside and compute volume
    points_inside = np.sum(inside_any)
    volume_per_point = grid_spacing**3
    return points_inside * volume_per_point


def clean_rotational_constants_by_geometry(
    rotational_constants,
    mode="physical",
    linear_rel_tol=1e-4,
    linear_abs_tol=1e-6,
    zero_abs_tol=1e-12,
    quasi_linear_ratio=1e4,
    return_status=False,
):
    """Clean Gaussian-style rotational constants without changing units.

    Gaussian prints rotational constants as ``A, B, C`` on lines such as::

        Rotational constants (GHZ):   A   B   C

    For ordinary nonlinear molecules, three finite constants are printed.
    For linear and quasi-linear molecules, Gaussian may instead print an
    overflow token such as ``*************`` for the axial constant, a huge
    finite axial constant, or ``0.0`` for the axial constant while still
    printing two perpendicular constants. Physically, only one perpendicular
    rotational constant is meaningful for a linear rotor.

    This helper supports two modes:

    ``mode="gaussian"``
        Preserve Gaussian's printed values exactly, except that overflow
        tokens should already have been converted to ``np.inf`` by the parser.
        Here ``np.inf`` is only a sentinel for an overflowed/unknown printed
        value; it is not a true known physical infinity.

    ``mode="physical"``
        Collapse effectively linear or quasi-linear ``[A, B, C]`` triples to a
        single perpendicular constant ``[B_perp]`` while keeping nonlinear
        triples unchanged.

    Units are preserved. For Gaussian output this usually means GHz, but this
    utility performs no unit conversion and can be used with any consistent
    units.

    Parameters
    ----------
    rotational_constants : array-like
        Rotational constants in Gaussian order ``A, B, C`` or an already
        cleaned single-value linear-rotor representation.
    mode : {"gaussian", "physical"}, optional
        Cleanup mode. ``"gaussian"`` preserves the parsed values; ``"physical"``
        applies linear/quasi-linear cleanup rules.
    linear_rel_tol : float, optional
        Relative tolerance for deciding whether ``B`` and ``C`` are close.
    linear_abs_tol : float, optional
        Absolute tolerance for deciding whether ``B`` and ``C`` are close, in
        the same units as ``rotational_constants``.
    zero_abs_tol : float, optional
        Absolute tolerance for deciding whether the axial constant ``A`` is
        effectively zero.
    quasi_linear_ratio : float, optional
        If ``A / B_perp`` exceeds this value and ``B``/``C`` are collapsible,
        the rotor is treated as quasi-linear in physical mode.
    return_status : bool, optional
        If ``True``, also return a status string describing how the values were
        interpreted.

    Returns
    -------
    np.ndarray
        Cleaned rotational constants in the same units as the input.
    tuple[np.ndarray, str]

        Always returned as ``(cleaned_constants, status)``.
        The status is one of
        ``"gaussian"``, ``"gaussian_overflow"``, ``"linear"``,
        ``"quasi_linear"``, ``"nonlinear"``, or ``"unknown"``.

    """

    vals = np.asarray(rotational_constants, dtype=float)

    if mode == "gaussian":
        status = "gaussian_overflow" if np.isinf(vals).any() else "gaussian"
        cleaned = vals.copy()
    elif mode != "physical":
        raise ValueError(
            f"Unsupported rotational-constant cleanup mode: {mode!r}."
        )
    elif vals.size == 1:
        cleaned = vals.copy()
        status = "linear"
    elif vals.size != 3:
        cleaned = vals.copy()
        status = "unknown"
    else:
        A, B, C = vals
        if not (np.isfinite(B) and np.isfinite(C)):
            cleaned = vals.copy()
            status = "unknown"
        elif np.isinf(A):
            if B == C:
                cleaned = np.array([C], dtype=float)
                status = "linear"
            elif np.isclose(
                B,
                C,
                rtol=linear_rel_tol,
                atol=linear_abs_tol,
            ):
                cleaned = np.array([0.5 * (B + C)])
                status = "linear"
            else:
                cleaned = vals.copy()
                status = "nonlinear"
        else:
            if B == C:
                B_perp = C
                bc_collapsible = True
            elif np.isclose(
                B,
                C,
                rtol=linear_rel_tol,
                atol=linear_abs_tol,
            ):
                B_perp = 0.5 * (B + C)
                bc_collapsible = True
            else:
                cleaned = vals.copy()
                status = "nonlinear"
                bc_collapsible = False

            if bc_collapsible:
                axial_zero = np.isfinite(A) and abs(A) <= zero_abs_tol
                axial_huge = (
                    np.isfinite(A)
                    and A > 0.0
                    and B_perp > 0.0
                    and A / B_perp > quasi_linear_ratio
                )
                if axial_zero:
                    cleaned = np.array([B_perp], dtype=float)
                    status = "linear"
                elif axial_huge:
                    cleaned = np.array([B_perp], dtype=float)
                    status = "quasi_linear"
                else:
                    cleaned = vals.copy()
                    status = "nonlinear"

    logger.debug(f"Cleaned rotational constants: {cleaned}, status: {status}")
    if return_status:
        return cleaned, status
    return cleaned
