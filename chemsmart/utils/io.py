import re

import numpy as np

from chemsmart.io.molecules.structure import Molecule


def create_molecule_list(
    orientations,
    orientations_pbc,
    energies,
    forces,
    symbols,
    charge,
    multiplicity,
    frozen_atoms,
    pbc_conditions,
    num_structures=None,
):
    """Helper function to create Molecule objects."""
    num_structures = num_structures or len(orientations)
    print(energies, num_structures)

    return [
        Molecule(
            symbols=symbols,
            positions=orientations[i],
            translation_vectors=orientations_pbc[i],
            charge=charge,
            multiplicity=multiplicity,
            frozen_atoms=frozen_atoms,
            pbc_conditions=pbc_conditions,
            energy=energies[i] if energies else None,
            forces=forces[i] if forces else None,
        )
        for i in range(num_structures)
    ]


def clean_duplicate_structure(orientations):
    """Remove the last structure if it's a duplicate of the previous one."""
    if orientations and len(orientations) > 1:
        if np.allclose(orientations[-1], orientations[-2], rtol=1e-5):
            orientations.pop(-1)


def increment_numbers(s, increment=1):
    # Use re.sub with a lambda to increment each number
    return re.sub(r"\d+", lambda m: str(int(m.group()) + increment), s)


def remove_keyword(text, keyword):
    return re.sub(
        r"\b" + re.escape(keyword) + r"\b", "", text, flags=re.IGNORECASE
    )
