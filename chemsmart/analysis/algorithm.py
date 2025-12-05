import numpy as np
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation
from ase.data import atomic_numbers, covalent_radii

# Default buffer for bond cutoff calculations in Angstroms
DEFAULT_BUFFER = 0.3

# Element symbol to atomic number mapping
ELEMENT_TO_ATOMIC_NUMBER = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
    'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
    'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26,
    'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34,
    'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42,
    'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
    'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58,
    'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66,
    'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74,
    'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82,
    'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90,
    'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98,
    'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103
}

# Atomic number to element symbol mapping (reverse of ELEMENT_TO_ATOMIC_NUMBER)
ATOMIC_NUMBER_TO_ELEMENT = {v: k for k, v in ELEMENT_TO_ATOMIC_NUMBER.items()}

def calc_pairwise_distances(mol_a: np.ndarray, mol_b: np.ndarray) -> np.ndarray:
    """
    Calculate pairwise distances between all atoms in two molecules.
    
    Parameters
    ----------
    mol_a : np.ndarray
        First molecule, shape (n, 4), each row is [element_index, x, y, z]
    mol_b : np.ndarray
        Second molecule, shape (m, 4), each row is [element_index, x, y, z]
    
    Returns
    -------
    np.ndarray
        Distance matrix of shape (n, m), where element [i, j] is the distance
        between atom i in mol_a and atom j in mol_b
    """
    coords_a = mol_a[:, 1:4]
    coords_b = mol_b[:, 1:4]
    
    # Use broadcasting to compute all pairwise distances
    # coords_a: (n, 3) -> (n, 1, 3)
    # coords_b: (m, 3) -> (1, m, 3)
    diff = coords_a[:, np.newaxis, :] - coords_b[np.newaxis, :, :]
    distances = np.linalg.norm(diff, axis=2)
    
    return distances


def calc_relative_coords(mol: np.ndarray, base_index: int) -> np.ndarray:
    """
    Calculate relative coordinates of all atoms with respect to a base atom.
    
    Parameters
    ----------
    mol : np.ndarray
        Molecule array, shape (n, 4), each row is [atomic_number, x, y, z]
    base_index : int
        Index of the base atom (0-based)
    
    Returns
    -------
    np.ndarray
        Array of shape (n, 4), same format as input, but with coordinates
        relative to the base atom. The base atom will have coordinates (0, 0, 0).
    
    Raises
    ------
    IndexError
        If base_index is out of range
    """
    n = mol.shape[0]
    if base_index < 0 or base_index >= n:
        raise IndexError(f"base_index {base_index} out of range [0, {n-1}]")
    
    # Copy the input array to avoid modifying the original
    result = mol.copy()
    
    # Get the base atom coordinates
    base_coords = mol[base_index, 1:4]
    
    # Subtract base coordinates from all atoms
    result[:, 1:4] = mol[:, 1:4] - base_coords
    
    return result


def find_optimal_position_for_single_atom(
    atom_coord: np.ndarray,
    mol: np.ndarray,
    link_index: int,
    buffer: float = DEFAULT_BUFFER
) -> np.ndarray:
    """
    Find the optimal position for an atom such that the sum of distances
    to all atoms in mol is maximized, subject to constraints:
    1. The distance to the link atom equals the max bond distance (equality constraint)
    2. The distance to all atoms >= min bond distance (inequality constraint)
    
    Uses Lagrange multiplier method with numerical optimization.
    
    Parameters
    ----------
    atom_coord : np.ndarray
        Atom info, shape (4,), [element_index, x, y, z]. Used to get element type
        and initial position.
    mol : np.ndarray
        Molecule array, shape (n, 4), each row is [atomic_number, x, y, z]
    link_index : int
        Index of the link atom in mol (0-based)
    buffer : float
        Buffer for max bond distance calculation (default: 0.3 Å)
    
    Returns
    -------
    np.ndarray
        Optimal position, shape (4,), [element_index, x, y, z]
    
    Raises
    ------
    IndexError
        If link_index is out of range
    """
    n_atoms = mol.shape[0]
    if link_index < 0 or link_index >= n_atoms:
        raise IndexError(f"link_index {link_index} out of range [0, {n_atoms-1}]")
    
    # Get element index of the atom
    atom_element = int(atom_coord[0])
    mol_elements = mol[:, 0].astype(int)
    
    # Calculate bond cutoff for link atom (equality constraint)
    # bond_dist = covalent_radius_A + covalent_radius_B + buffer (max bonding distance)
    bond_dist = covalent_radii[atom_element] + covalent_radii[mol_elements[link_index]] + buffer
    
    # Calculate minimum distance array (inequality constraint)
    # min_dist = covalent_radius_A + covalent_radius_B (no buffer, to prevent overlap)
    min_dist_array = np.array([
        covalent_radii[atom_element] + covalent_radii[z]
        for z in mol_elements
    ])
    
    # Get link atom coordinates
    link_coords = mol[link_index, 1:4]
    
    # Get all atom coordinates in mol
    mol_coords = mol[:, 1:4]
    
    # =========================================================================
    # Lagrange Multiplier Formulation:
    # 
    # Objective: max sum_i ||x - mol_i||  (equivalent to min -sum_i ||x - mol_i||)
    # 
    # Constraints:
    #   Equality:   ||x - link||^2 = bond_dist^2
    #   Inequality: ||x - mol_i||^2 >= min_dist_i^2  for all i
    #
    # We solve this using scipy.optimize.minimize with SLSQP.
    # =========================================================================
    
    def objective(x):
        """Negative sum of distances (to maximize distance sum)"""
        diff = x - mol_coords  # (n, 3)
        distances = np.linalg.norm(diff, axis=1)  # (n,)
        return -np.sum(distances)
    
    def objective_gradient(x):
        """Gradient of negative objective function"""
        diff = x - mol_coords  # (n, 3)
        distances = np.linalg.norm(diff, axis=1, keepdims=True)  # (n, 1)
        # Avoid division by zero
        distances = np.maximum(distances, 1e-10)
        grad = -np.sum(diff / distances, axis=0)  # (3,)
        return grad
    
    def eq_constraint(x):
        """Equality constraint: ||x - link||^2 - bond_dist^2 = 0"""
        diff = x - link_coords
        return np.dot(diff, diff) - bond_dist**2
    
    def eq_constraint_gradient(x):
        """Gradient of equality constraint"""
        return 2 * (x - link_coords)
    
    def ineq_constraint(x):
        """
        Inequality constraints: ||x - mol_i||^2 - min_dist_i^2 >= 0 for all i
        SLSQP requires constraints in form f(x) >= 0
        """
        diff = x - mol_coords  # (n, 3)
        dist_sq = np.sum(diff**2, axis=1)  # (n,)
        return dist_sq - min_dist_array**2
    
    def ineq_constraint_gradient(x):
        """Gradient of inequality constraints, shape (n, 3)"""
        return 2 * (x - mol_coords)
    
    # Initial guess: use the original atom position, projected onto the constraint sphere
    x_init = atom_coord[1:4].copy()
    
    # If initial point is too close to link, perturb it
    init_diff = x_init - link_coords
    init_dist = np.linalg.norm(init_diff)
    if init_dist < 1e-6:
        # Default direction
        x_init = link_coords + np.array([bond_dist, 0, 0])
    else:
        # Project onto the constraint sphere
        x_init = link_coords + (init_diff / init_dist) * bond_dist
    
    # Define constraints for scipy
    constraints = [
        {
            'type': 'eq',
            'fun': eq_constraint,
            'jac': eq_constraint_gradient
        },
        {
            'type': 'ineq',
            'fun': ineq_constraint,
            'jac': ineq_constraint_gradient
        }
    ]
    
    # Solve using SLSQP (Sequential Least Squares Programming)
    result = minimize(
        objective,
        x_init,
        method='SLSQP',
        jac=objective_gradient,
        constraints=constraints,
        options={'ftol': 1e-9, 'maxiter': 1000}
    )
    
    # Extract optimal position
    optimal_x = result.x
    
    # Build result
    optimal_coord = np.array([atom_coord[0], optimal_x[0], optimal_x[1], optimal_x[2]])
    
    return optimal_coord


def find_optimal_position(
    skeleton_coord: np.ndarray,
    sub_coord: np.ndarray,
    skeleton_link_index: int,
    sub_link_index: int,
    buffer: float = DEFAULT_BUFFER
) -> np.ndarray:
    """
    Find the optimal position for a substituent molecule (sub) to be attached to a 
    skeleton molecule, such that the total distance sum is maximized while avoiding
    atomic overlaps.
    
    The optimization is performed by finding the optimal position for the link atom
    in sub (sub_link), and then placing all other atoms in sub according to their
    relative positions to sub_link.
    
    Parameters
    ----------
    skeleton_coord : np.ndarray
        Skeleton molecule, shape (n, 4), each row is [atomic_number, x, y, z]
    sub_coord : np.ndarray
        Substituent molecule, shape (m, 4), each row is [atomic_number, x, y, z]
    skeleton_link_index : int
        Index of the link atom in skeleton (0-based)
    sub_link_index : int
        Index of the link atom in sub (0-based)
    buffer : float
        Buffer for max bond distance calculation (default: 0.3 Å)
    
    Returns
    -------
    np.ndarray
        Optimal position of the substituent molecule, shape (m, 4),
        each row is [atomic_number, x, y, z]
    
    Constraints
    -----------
    1. All atoms in sub maintain their relative positions to sub_link (rigid body)
    2. Distance(sub_link, skeleton_link) = bond_dist (equality constraint)
    3. Distance(sub_i, skeleton_j) >= min_dist for all i, j (inequality constraints)
       (excluding the link-link pair which is handled by equality constraint)
    
    Objective
    ---------
    Maximize: sum of all pairwise distances between sub atoms and skeleton atoms
    
    Raises
    ------
    IndexError
        If skeleton_link_index or sub_link_index is out of range
    """
    n_skeleton = skeleton_coord.shape[0]
    n_sub = sub_coord.shape[0]
    
    if skeleton_link_index < 0 or skeleton_link_index >= n_skeleton:
        raise IndexError(f"skeleton_link_index {skeleton_link_index} out of range [0, {n_skeleton-1}]")
    if sub_link_index < 0 or sub_link_index >= n_sub:
        raise IndexError(f"sub_link_index {sub_link_index} out of range [0, {n_sub-1}]")
    
    # Calculate relative coordinates of sub atoms with respect to sub_link
    # This gives us the fixed relative positions that must be maintained
    sub_relative = calc_relative_coords(sub_coord, sub_link_index)
    relative_offsets = sub_relative[:, 1:4]  # shape (m, 3), offsets from sub_link
    
    # Get element indices
    sub_elements = sub_coord[:, 0].astype(int)
    skeleton_elements = skeleton_coord[:, 0].astype(int)
    
    # Get skeleton coordinates
    skeleton_coords = skeleton_coord[:, 1:4]  # shape (n, 3)
    skeleton_link_coords = skeleton_coords[skeleton_link_index]
    
    # Calculate bond distance for equality constraint (sub_link to skeleton_link)
    sub_link_element = sub_elements[sub_link_index]
    skeleton_link_element = skeleton_elements[skeleton_link_index]
    bond_dist = covalent_radii[sub_link_element] + covalent_radii[skeleton_link_element]
    
    # Calculate minimum distance matrix for inequality constraints
    # min_dist[i, j] = covalent_radii[sub_i] + covalent_radii[skeleton_j] + buffer
    sub_radii = np.array([covalent_radii[z] for z in sub_elements])  # (m,)
    skeleton_radii = np.array([covalent_radii[z] for z in skeleton_elements])  # (n,)
    min_dist_matrix = sub_radii[:, np.newaxis] + skeleton_radii[np.newaxis, :] + buffer  # (m, n)
    
    # Mask to exclude the link-link pair from inequality constraints
    # because it is already constrained by the equality constraint
    ineq_mask = np.ones((n_sub, n_skeleton), dtype=bool)
    ineq_mask[sub_link_index, skeleton_link_index] = False
    
    # =========================================================================
    # Lagrange Multiplier Formulation:
    # 
    # Decision variables: 
    #   x[:3] = position of sub_link (3 values)
    #   x[3:] = rotation angles (Euler angles: alpha, beta, gamma) (3 values)
    # 
    # Derived positions: 
    #   R = Rotation(x[3:])
    #   sub_i_position = x[:3] + (relative_offsets[i] @ R.T)
    # 
    # Objective: max sum_{i,j} ||sub_i - skeleton_j||
    #          = min -sum_{i,j} ||sub_i - skeleton_j||
    # 
    # Constraints:
    #   Equality:   ||x[:3] - skeleton_link||^2 = bond_dist^2
    #   Inequality: ||sub_i - skeleton_j||^2 >= min_dist[i,j]^2 for all i, j (except link pair)
    # =========================================================================
    
    def get_sub_positions(x):
        """Helper to calculate sub positions from optimization variables"""
        pos_link = x[:3]
        angles = x[3:]
        
        # Create rotation matrix
        r = Rotation.from_euler('xyz', angles)
        R_matrix = r.as_matrix()  # (3, 3)
        
        # Rotate offsets: (m, 3) @ (3, 3).T -> (m, 3)
        rotated_offsets = relative_offsets @ R_matrix.T
        
        # Translate
        sub_positions = pos_link + rotated_offsets
        return sub_positions

    def objective(x):
        """
        Negative sum of all pairwise distances (to maximize via minimization)
        """
        sub_positions = get_sub_positions(x)
        
        # Compute all pairwise distances between sub and skeleton
        # sub_positions: (m, 3) -> (m, 1, 3)
        # skeleton_coords: (n, 3) -> (1, n, 3)
        diff = sub_positions[:, np.newaxis, :] - skeleton_coords[np.newaxis, :, :]  # (m, n, 3)
        distances = np.linalg.norm(diff, axis=2)  # (m, n)
        
        return -np.sum(distances)
    
    # Note: We use numerical gradients (finite difference) for simplicity and robustness
    # with rotation variables.
    
    def eq_constraint(x):
        """Equality constraint: ||x[:3] - skeleton_link||^2 - bond_dist^2 = 0"""
        pos_link = x[:3]
        diff = pos_link - skeleton_link_coords
        return np.dot(diff, diff) - bond_dist**2
    
    def ineq_constraint(x):
        """
        Inequality constraints: ||sub_i - skeleton_j||^2 - min_dist[i,j]^2 >= 0
        Returns a 1D array of shape (m * n - 1,)
        """
        sub_positions = get_sub_positions(x)
        
        # diff[i, j, :] = sub_i - skeleton_j
        diff = sub_positions[:, np.newaxis, :] - skeleton_coords[np.newaxis, :, :]  # (m, n, 3)
        dist_sq = np.sum(diff**2, axis=2)  # (m, n)
        
        # Apply mask to exclude link-link pair
        constraints = (dist_sq - min_dist_matrix**2)[ineq_mask]
        return constraints
    
    # Initial guess
    # 1. Position: project sub_link onto the constraint sphere
    x_init_pos = sub_coord[sub_link_index, 1:4].copy()
    init_diff = x_init_pos - skeleton_link_coords
    init_dist = np.linalg.norm(init_diff)
    
    if init_dist < 1e-6:
        # Default direction: opposite to the center of mass of skeleton
        skeleton_com = np.mean(skeleton_coords, axis=0)
        direction = skeleton_link_coords - skeleton_com
        dir_norm = np.linalg.norm(direction)
        if dir_norm < 1e-6:
            direction = np.array([1.0, 0.0, 0.0])
        else:
            direction = direction / dir_norm
        x_init_pos = skeleton_link_coords + direction * bond_dist
    else:
        # Project onto the constraint sphere
        x_init_pos = skeleton_link_coords + (init_diff / init_dist) * bond_dist
        
    # 2. Rotation: start with 0 angles (original orientation)
    # Ideally we could try multiple random rotations if this fails
    x_init_rot = np.zeros(3)
    
    x_init = np.concatenate([x_init_pos, x_init_rot])
    
    # Define constraints for scipy
    constraints = [
        {
            'type': 'eq',
            'fun': eq_constraint
        },
        {
            'type': 'ineq',
            'fun': ineq_constraint
        }
    ]
    
    # Solve using SLSQP
    # We remove 'jac' to let scipy estimate gradients numerically
    result = minimize(
        objective,
        x_init,
        method='SLSQP',
        constraints=constraints,
        options={'ftol': 1e-6, 'maxiter': 1000}
    )
    
    # Extract optimal position
    optimal_x = result.x
    optimal_sub_positions = get_sub_positions(optimal_x)
    
    # Build result array
    optimal_sub = np.zeros_like(sub_coord)
    optimal_sub[:, 0] = sub_coord[:, 0]  # Keep atomic numbers
    optimal_sub[:, 1:4] = optimal_sub_positions
    
    return optimal_sub
    # Compute optimal positions for all sub atoms
    optimal_sub_positions = optimal_sub_link + relative_offsets  # (m, 3)
    
    # Build result array
    optimal_sub = np.zeros_like(sub_coord)
    optimal_sub[:, 0] = sub_coord[:, 0]  # Keep atomic numbers
    optimal_sub[:, 1:4] = optimal_sub_positions
    
    return optimal_sub


def read_xyz(filepath: str) -> np.ndarray:
    """
    Read XYZ format file and return a numpy array containing atomic information.
    
    Parameters
    ----------
    filepath : str
        Path to the XYZ file
    
    Returns
    -------
    np.ndarray
        Array of shape (n_atoms, 4), each row is [atomic_number, x, y, z]
        atomic_number is an integer, x, y, z are float coordinates (unit: Angstrom)
    
    Example
    -------
    >>> atoms = read_xyz('molecule.xyz')
    >>> print(atoms)
    [[ 6.   0.000  0.000  0.000]   # C atom
     [ 1.   1.089  0.000  0.000]   # H atom
     ...]
    """
    atoms_data = []
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # First line is the number of atoms
    n_atoms = int(lines[0].strip())
    
    # Second line is a comment, skip it
    # Atom data starts from the third line
    for i in range(2, 2 + n_atoms):
        line = lines[i].strip()
        if not line:
            continue
        
        # Split line and handle possible comments
        parts = line.split('#')[0].split()  # Remove comments after #
        
        element = parts[0]
        x = float(parts[1])
        y = float(parts[2])
        z = float(parts[3])
        
        # Convert element symbol to atomic number
        # Handle case sensitivity (e.g., 'c' -> 'C', 'cl' -> 'Cl')
        element = element.capitalize()
        atomic_number = ELEMENT_TO_ATOMIC_NUMBER.get(element, 0)
        
        if atomic_number == 0:
            raise ValueError(f"Unknown element symbol: {element}")
        
        atoms_data.append([atomic_number, x, y, z])
    
    return np.array(atoms_data, dtype=np.float64)


def get_covalent_radius(element):
    """
    Returns the covalent radius of an element in Å.

    Args:
        element (str): Atomic symbol (e.g., "C", "O", "H").

    Returns:
        float: Covalent radius in Å, or None if not found.
    """
    element = element.capitalize()  # Ensure correct capitalization
    atomic_number = atomic_numbers.get(element)
    if atomic_number is None:
        raise ValueError(f"Unknown element: {element}")
    return covalent_radii[atomic_number]


def get_bond_cutoff(element1, element2, buffer=DEFAULT_BUFFER):
    """
    Calculates bond cutoff distance based on covalent radii and buffer.
    A good bond cutoff distance for molecular graphs depends on the type of
    chemical bonds and the elements involved. Here are some guidelines:

        General Bond Cutoffs (Approximate Values)
            Bond Type	Cutoff Distance (Å)
            Covalent Bonds	1.0 – 1.6 Å
            Hydrogen Bonds	2.5 – 3.5 Å
            Van der Waals	3.0 – 4.5 Å

        Element-Specific Covalent Cutoffs
        A common way to estimate a covalent bond cutoff is:
            𝑅_cutoff = 𝑅_𝐴 + 𝑅_𝐵 + tolerance
        where:
        R_A and R_B are the covalent radii of atoms A and B,
        A tolerance (~0.2-0.4 Å) accounts for bond flexibility.
        For example:
            C–C bond: 1.54 Å
            C–O bond: ~1.43 Å
            C–H bond: ~1.1 Å
    Recommended Bond Cutoffs for Molecular Graphs
        For Covalent Bond Networks → 1.6 Å – 2.0 Å
        For Molecular Clusters (H-bonding included) → 2.5 Å – 3.5 Å
        For Van der Waals interactions (e.g., protein-ligand) → 3.5 Å – 4.5 Å

    Args:
        element1 (str): Atomic symbol of first element.
        element2 (str): Atomic symbol of second element.
        buffer (float): Additional buffer for bond length (default: 0.3 Å).

    Returns:
        float: Bond cutoff distance in Å.
    """
    r1 = get_covalent_radius(element1)
    r2 = get_covalent_radius(element2)
    return r1 + r2 + buffer


def calc_pairwise_bond_cutoffs(mol_a: np.ndarray, mol_b: np.ndarray, buffer: float = DEFAULT_BUFFER) -> np.ndarray:
    """
    Calculate pairwise bond cutoff distances between all atoms in two molecules.
    
    Parameters
    ----------
    mol_a : np.ndarray
        First molecule, shape (n, 4), each row is [atomic_number, x, y, z]
    mol_b : np.ndarray
        Second molecule, shape (m, 4), each row is [atomic_number, x, y, z]
    buffer : float
        Additional buffer for bond length (default: 0.3 Å)
    
    Returns
    -------
    np.ndarray
        Bond cutoff matrix of shape (n, m), where element [i, j] is the bond cutoff
        distance between atom i in mol_a and atom j in mol_b
    """
    n = mol_a.shape[0]
    m = mol_b.shape[0]
    
    # Get atomic numbers
    atomic_nums_a = mol_a[:, 0].astype(int)
    atomic_nums_b = mol_b[:, 0].astype(int)
    
    # Get covalent radii for all atoms
    radii_a = np.array([covalent_radii[z] for z in atomic_nums_a])
    radii_b = np.array([covalent_radii[z] for z in atomic_nums_b])
    
    # Use broadcasting to compute all pairwise bond cutoffs
    # radii_a: (n,) -> (n, 1)
    # radii_b: (m,) -> (1, m)
    cutoffs = radii_a[:, np.newaxis] + radii_b[np.newaxis, :] + buffer
    
    return cutoffs


def write_xyz(filepath: str, mol: np.ndarray, comment: str = "") -> None:
    """
    Write molecular data to an XYZ format file.
    
    Parameters
    ----------
    filepath : str
        Path to the output XYZ file
    mol : np.ndarray
        Molecule array, shape (n, 4), each row is [atomic_number, x, y, z]
    comment : str, optional
        Comment line (second line in XYZ file), default is empty string
    
    Returns
    -------
    None
    
    Example
    -------
    >>> mol = np.array([[6, 0.0, 0.0, 0.0], [1, 1.089, 0.0, 0.0]])
    >>> write_xyz('output.xyz', mol, comment='My molecule')
    """
    n_atoms = mol.shape[0]
    
    with open(filepath, 'w') as f:
        # First line: number of atoms
        f.write(f"{n_atoms}\n")
        
        # Second line: comment
        f.write(f"{comment}\n")
        
        # Atom data
        for row in mol:
            atomic_num = int(row[0])
            x, y, z = row[1], row[2], row[3]
            
            # Convert atomic number to element symbol
            element = ATOMIC_NUMBER_TO_ELEMENT.get(atomic_num, 'X')
            
            # Write atom line with 6 decimal places
            f.write(f"{element:<2} {x:>14.6f} {y:>14.6f} {z:>14.6f}\n")


if __name__ == "__main__":
    # Test
    test_path = '/Users/wanglewen/Desktop/Project/test/ase/test.xyz'
    atoms = read_xyz(test_path)
    print(f"Read {len(atoms)} atoms")
    print(f"Array shape: {atoms.shape}")
    print("\nFirst 5 atoms:")
    print(f"{'Atomic #':<10}{'X':<12}{'Y':<12}{'Z':<12}")
    for row in atoms[:5]:
        print(f"{int(row[0]):<10}{row[1]:<12.6f}{row[2]:<12.6f}{row[3]:<12.6f}")


