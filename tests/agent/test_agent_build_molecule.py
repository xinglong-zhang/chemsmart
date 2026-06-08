import numpy as np

from chemsmart.agent.tools import build_molecule
from chemsmart.io.molecules.structure import Molecule


class TestBuildMolecule:
    def test_returns_molecule_from_xyz_fixture(self, single_molecule_xyz_file):
        molecule = build_molecule(single_molecule_xyz_file)

        assert isinstance(molecule, Molecule)
        assert molecule.empirical_formula == "C37H25Cl3N3O3"

    def test_preserves_1_based_indexing(self, multiple_molecules_xyz_file):
        molecule = build_molecule(multiple_molecules_xyz_file, index="1")

        assert isinstance(molecule, Molecule)
        assert np.allclose(
            molecule.positions[0],
            np.array([-1.0440166707, -2.3921211654, -1.1765767093]),
            rtol=1e-5,
        )
