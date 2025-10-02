"""
Data for molecules
"""

from ase.data import atomic_numbers, covalent_radii

# Default buffer for bond cutoff calculations in Angstroms
DEFAULT_BUFFER = 0.3


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
