"""
Data for molecules
"""

from ase.data import atomic_numbers, covalent_radii

# Common chemical abbreviations and their SMILES representations
# Used for expanding abbreviations in chemical structure images
CHEMICAL_ABBREVIATIONS = {
    "Ad": "C1C2CC3CC1CC(C2)C3",  # Adamantyl (1-adamantyl)
    "Ph": "c1ccccc1",  # Phenyl
    "Me": "C",  # Methyl
    "Et": "CC",  # Ethyl
    "nPr": "CCC",  # n-Propyl
    "iPr": "C(C)C",  # Isopropyl
    "Bu": "CCCC",  # Butyl
    "iBu": "CC(C)C",  # Isobutyl
    "sBu": "C(C)CC",  # sec-Butyl
    "tBu": "C(C)(C)C",  # tert-Butyl
    "Bn": "Cc1ccccc1",  # Benzyl
    "Ac": "C(=O)C",  # Acetyl
    "Bz": "C(=O)c1ccccc1",  # Benzoyl
    "Ts": "S(=O)(=O)c1ccc(C)cc1",  # Tosyl
    "Ms": "S(=O)(=O)C",  # Mesyl
    "Tf": "S(=O)(=O)C(F)(F)F",  # Triflyl
    "Cy": "C1CCCCC1",  # Cyclohexyl
    "OMe": "OC",  # Methoxy
    "OEt": "OCC",  # Ethoxy
    "NMe2": "N(C)C",  # Dimethylamino
    "CF3": "C(F)(F)F",  # Trifluoromethyl
    "NO2": "N(=O)=O",  # Nitro
    "Boc": "C(=O)OC(C)(C)C",  # tert-Butoxycarbonyl
    "Cbz": "C(=O)OCc1ccccc1",  # Benzyloxycarbonyl
    "Fmoc": "C(=O)OCC1c2ccccc2-c2ccccc21",  # Fluorenylmethyloxycarbonyl
}

# Dash characters commonly used in chemical drawings
DASH_CHARACTERS = ["-", "–", "—"]  # hyphen, en-dash, em-dash

# Minimum length for a valid SMILES string
MIN_VALID_SMILES_LENGTH = 3

# Mapping of functional group text to SMILES atom
SUBSTITUENT_MAPPING = {
    "SH": "S",  # Thiol
    "OH": "O",  # Hydroxyl
    "NH2": "N",  # Amine
    "NH₂": "N",  # Amine (with subscript)
}


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
