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

    return [
        Molecule(
            symbols=symbols,
            positions=orientations[i],
            translation_vectors=orientations_pbc[i],
            charge=charge,
            multiplicity=multiplicity,
            frozen_atoms=frozen_atoms,
            pbc_conditions=pbc_conditions,
            energy=energies[i],
            forces=forces[i],
        )
        for i in range(num_structures)
    ]


def clean_duplicate_structure(orientations):
    """Remove the last structure if it's a duplicate of the previous one."""
    if orientations and len(orientations) > 1:
        if np.allclose(orientations[-1], orientations[-2], rtol=1e-5):
            orientations.pop(-1)
