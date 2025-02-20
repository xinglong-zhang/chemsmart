import numpy as np


def is_collinear(coords, tol=1e-5):
    """Check if three points are collinear using cross product method."""
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
#     # instead of returning evecs.T as in ase, we return evecs directly, as in Gaussian
#     return moi_tensor, evals, evecs


def calculate_moments_of_inertia(mass, coords):
    """Calculate the moment of inertia tensor and principal moments of inertia.

    Parameters:
    - mass (list or np.array): Atomic masses corresponding to each coordinate.
    - coords (list or np.array): Nx3 array of atomic coordinates.

    Returns:
    - moments_of_inertia_principal_axes (np.array): Sorted eigenvalues of moi_tensor.
    - moi_tensor (np.array): 3x3 moment of inertia tensor.
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
    evals, evecs = np.linalg.eigh(moi_tensor)  # Eigenvalues and eigenvectors
    # instead of returning evecs.T as in ase, we return evecs directly, as in Gaussian
    return moi_tensor, evals, evecs
