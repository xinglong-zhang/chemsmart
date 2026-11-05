"""In-memory molecule fixtures."""

import os

import numpy as np
import pytest
import rdkit.Chem.rdDistGeom as rdDistGeom
from ase import Atoms
from rdkit import Chem

from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.xyz.xyzfile import XYZFile


## pytest fixtures for molecules
@pytest.fixture()
def methanol_molecule():
    symbols = ["C", "O", "H", "H", "H", "H"]
    coords = np.array(
        [
            [0.000000, 0.000000, 0.000000],  # C
            [1.430000, 0.000000, 0.000000],  # O
            [1.109545, 0.904936, 0.000000],  # H (hydroxyl)
            [-0.363849, 1.027479, 0.000000],  # H
            [-0.363849, -0.513740, 0.889823],  # H
            [-0.363849, -0.513740, -0.889823],  # H
        ],
        dtype=float,
    )
    methanol = Molecule(symbols=symbols, positions=coords)
    return methanol


@pytest.fixture()
def ethanol_molecule():
    symbols = ["O", "C", "C", "H", "H", "H", "H", "H", "H"]
    coords = np.array(
        [
            [-1.1712, 0.2997, 0.0000],
            [-0.0463, -0.5665, 0.0000],
            [1.2175, 0.2668, 0.0000],
            [-0.0958, -1.2120, 0.8819],
            [-0.0952, -1.1938, -0.8946],
            [2.1050, -0.3720, -0.0177],
            [1.2426, 0.9307, -0.8704],
            [1.2616, 0.9052, 0.8886],
            [-1.1291, 0.8364, 0.8099],
        ],
        dtype=float,
    )
    ethanol = Molecule(symbols=symbols, positions=coords)
    return ethanol


@pytest.fixture()
def methanol_molecules(methanol_molecule):
    methanol = methanol_molecule

    # rotated methanol
    ase_atoms = methanol.to_ase()
    ase_atoms.rotate(90, [0, 0, 1])
    methanol_rot1 = Molecule.from_ase_atoms(ase_atoms)

    ase_atoms = methanol.to_ase()
    ase_atoms.rotate(20, [1, 1, 1])
    methanol_rot2 = Molecule.from_ase_atoms(ase_atoms)

    methanol_molecules = [methanol, methanol_rot1, methanol_rot2]

    return methanol_molecules


@pytest.fixture()
def methyl3hexane_molecule():
    symbols = [
        "C",
        "C",
        "C",
        "C",
        "C",
        "C",
        "C",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
    ]
    coords = np.array(
        [
            [0.828, -0.5939, -0.4105],
            [-0.6074, -0.1292, -0.7341],
            [1.6188, 0.3738, 0.4896],
            [-1.5335, 0.0189, 0.4761],
            [0.8379, -1.9975, 0.2043],
            [1.7539, 1.7749, -0.0888],
            [-2.8977, 0.553, 0.0634],
            [1.3601, -0.6621, -1.3694],
            [-1.0561, -0.849, -1.4317],
            [-0.5655, 0.8214, -1.2791],
            [2.6287, -0.0297, 0.6378],
            [1.1655, 0.4353, 1.4856],
            [-1.0999, 0.707, 1.2088],
            [-1.6692, -0.9475, 0.973],
            [1.8635, -2.3744, 0.2836],
            [0.4049, -2.004, 1.2096],
            [0.2737, -2.7005, -0.4174],
            [2.1704, 1.7451, -1.1005],
            [2.4243, 2.3747, 0.5353],
            [0.7895, 2.2899, -0.1235],
            [-3.384, -0.1214, -0.6486],
            [-2.8085, 1.5396, -0.4023],
            [-3.5477, 0.6486, 0.9389],
        ]
    )
    methyl3hexane = Molecule(symbols=symbols, positions=coords)
    return methyl3hexane


@pytest.fixture()
def tetrahydrofuran_molecule():
    symbols = ["O", "C", "C", "C", "C", "H", "H", "H", "H", "H", "H", "H", "H"]
    coords = np.array(
        [
            [1.2328, -0.0005, 0.0000],
            [-1.0107, -0.7202, -0.2205],
            [-1.0102, 0.7210, 0.2205],
            [0.3936, -1.1560, 0.1374],
            [0.3946, 1.1557, -0.1375],
            [-1.7823, -1.3279, 0.2593],
            [-1.1544, -0.7757, -1.3060],
            [-1.7812, 1.3292, -0.2593],
            [-1.1537, 0.7766, 1.3061],
            [0.4518, -1.4889, 1.1792],
            [0.7622, -1.9589, -0.5071],
            [0.4532, 1.4885, -1.1793],
            [0.7639, 1.9583, 0.5070],
        ]
    )
    thf_molecule = Molecule(symbols=symbols, positions=coords)
    return thf_molecule


@pytest.fixture()
def constrained_atoms():
    """Fixture to create a simple Ar2 dimer with constraints."""
    from ase import Atoms
    from ase.calculators.lj import LennardJones
    from ase.constraints import FixAtoms, FixBondLength

    # Simple Ar2 dimer with a reasonable separation
    r0 = 3.5  # Å
    atoms = Atoms(
        "Ar2", positions=[(0.0, 0.0, 0.0), (r0, 0.0, 0.0)], pbc=False
    )

    # Light-weight calculator for tests
    atoms.calc = LennardJones()  # defaults are fine for unit tests

    # Constraints:
    #  - Fix the first atom in space
    #  - Keep the Ar–Ar bond length fixed at its initial value
    constraints = [
        FixAtoms(indices=[0]),
        FixBondLength(0, 1),
    ]
    # set the constraints on the Atoms object
    atoms.set_constraint(constraint=constraints)

    # set velocity
    atoms.set_velocities([[0, 0, 0], [0, 0, 0]])  # Set zero velocities

    return atoms


@pytest.fixture()
def methanol_and_ethanol(methanol_molecule, ethanol_molecule):
    # molecules for testing
    # methanol
    methanol = methanol_molecule

    # ethanol
    ethanol = ethanol_molecule

    methanol_and_ethanol = [methanol, ethanol]
    return methanol_and_ethanol


@pytest.fixture()
def conformers_from_rdkit():
    """Generate multiple conformers for a complex molecule using RDKit."""
    smiles = "O=C([O-])CCn1c(=O)c(=O)[nH]c2cc([N+](=O)[O-])c(-n3ccc(C=NOCc4ccccc4)c3)cc21"
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # Generate 3D conformers
    ps = rdDistGeom.ETKDGv3()
    ps.randomSeed = 0xD06F00D
    ps.numThreads = 10
    conf_ids = rdDistGeom.EmbedMultipleConfs(mol, numConfs=300, params=ps)

    # Ensure each conformer is extracted into its own unique RDKit Mol object
    conformers = []
    for conf_id in conf_ids:
        single_conf_mol = Chem.Mol(mol)  # Copy molecule structure
        single_conf_mol.RemoveAllConformers()  # Remove all existing conformers
        single_conf_mol.AddConformer(
            mol.GetConformer(conf_id), assignId=True
        )  # Add only this conformer
        conformers.append(single_conf_mol)

    # Verify that each conformer contains exactly one conformer
    for conf in conformers:
        assert (
            conf.GetNumConformers() == 1
        ), "Each conformer should contain exactly one conformer."

    # Convert to Molecule instances
    conformers_from_rdkit = [
        Molecule.from_rdkit_mol(conf) for conf in conformers
    ]

    return conformers_from_rdkit


@pytest.fixture()
def water_molecule():
    """H₂O with realistic geometry."""
    return Molecule(
        symbols=["O", "H", "H"],
        positions=np.array(
            [
                [0.0000, 0.0000, 0.1173],
                [0.0000, 0.7572, -0.4692],
                [0.0000, -0.7572, -0.4692],
            ]
        ),
    )


@pytest.fixture()
def methane_molecule():
    """CH₄ tetrahedral."""
    return Molecule(
        symbols=["C", "H", "H", "H", "H"],
        positions=np.array(
            [
                [0.0000, 0.0000, 0.0000],
                [0.6276, 0.6276, 0.6276],
                [0.6276, -0.6276, -0.6276],
                [-0.6276, 0.6276, -0.6276],
                [-0.6276, -0.6276, 0.6276],
            ]
        ),
    )


@pytest.fixture()
def ethylene_molecule():
    """C₂H₄ – contains a C=C double bond."""
    return Molecule(
        symbols=["C", "C", "H", "H", "H", "H"],
        positions=np.array(
            [
                [0.0000, 0.0000, 0.6695],
                [0.0000, 0.0000, -0.6695],
                [0.0000, 0.9289, 1.2321],
                [0.0000, -0.9289, 1.2321],
                [0.0000, 0.9289, -1.2321],
                [0.0000, -0.9289, -1.2321],
            ]
        ),
    )


@pytest.fixture()
def acetylene_molecule():
    """C₂H₂ – contains a C≡C triple bond."""
    return Molecule(
        symbols=["C", "C", "H", "H"],
        positions=np.array(
            [
                [0.0000, 0.0000, 0.6013],
                [0.0000, 0.0000, -0.6013],
                [0.0000, 0.0000, 1.6644],
                [0.0000, 0.0000, -1.6644],
            ]
        ),
    )


@pytest.fixture()
def hydrogen_molecule():
    """H₂ – minimal molecule with H-H bond."""
    return Molecule(
        symbols=["H", "H"],
        positions=np.array(
            [
                [0.0, 0.0, 0.0],
                [0.74, 0.0, 0.0],
            ]
        ),
    )


@pytest.fixture()
def co2_molecule():
    """CO₂ – linear molecule with two C=O double bonds."""
    return Molecule(
        symbols=["C", "O", "O"],
        positions=np.array(
            [
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 1.16],
                [0.0, 0.0, -1.16],
            ]
        ),
    )


@pytest.fixture()
def single_atom_molecule():
    """Single argon atom – edge case for bond detection."""
    return Molecule(
        symbols=["Ar"],
        positions=np.array([[0.0, 0.0, 0.0]]),
    )


@pytest.fixture()
def water_with_metadata_molecule():
    """Water molecule with PDB residue metadata."""
    mol = Molecule(
        symbols=["O", "H", "H"],
        positions=np.array(
            [
                [0.0, 0.0, 0.0],
                [0.96, 0.0, 0.0],
                [-0.24, 0.93, 0.0],
            ]
        ),
        info={
            "record_type": ["HETATM", "HETATM", "HETATM"],
            "atom_name": ["O", "H1", "H2"],
            "residue_name": ["HOH", "HOH", "HOH"],
            "residue_number": [1, 1, 1],
            "chain_id": ["A", "A", "A"],
        },
    )
    mol.atom_names = ["O", "H1", "H2"]
    mol.residue_names = ["HOH", "HOH", "HOH"]
    mol.residue_numbers = [1, 1, 1]
    mol.chain_ids = ["A", "A", "A"]
    return mol


@pytest.fixture()
def chiral_molecule():
    """Molecule with a chiral center (C with 4 different substituents)."""
    return Molecule(
        symbols=["C", "Cl", "F", "Br", "I"],
        positions=np.array(
            [
                [0.0, 0.0, 0.0],
                [1.2, 0.0, -0.5],
                [-0.6, 1.0, 0.5],
                [-0.6, -1.0, 0.5],
                [0.0, 0.0, 1.3],
            ]
        ),
    )


@pytest.fixture()
def single_model_pdb_file(tmpdir):
    """PDB file with one implicit model (no MODEL/ENDMDL records)."""
    pdb_content = (
        "HETATM    1  O   HOH A   7       0.000   0.000   0.000  1.00  0.00           O\n"
        "HETATM    2  H1  HOH A   7       0.960   0.000   0.000  1.00  0.00           H\n"
        "HETATM    3  H2  HOH A   7      -0.240   0.930   0.000  1.00  0.00           H\n"
        "END\n"
    )
    filepath = os.path.join(str(tmpdir), "water.pdb")
    with open(filepath, "w") as f:
        f.write(pdb_content)
    return filepath


@pytest.fixture()
def multi_model_pdb_file(tmpdir):
    """PDB file with two explicit MODEL/ENDMDL blocks."""
    pdb_content = (
        "MODEL        1\n"
        "ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM      2  H1  HOH A   1       0.960   0.000   0.000  1.00  0.00           H\n"
        "ENDMDL\n"
        "MODEL        2\n"
        "ATOM      1  O   HOH B   2       1.500   2.500   3.500  1.00  0.00           O\n"
        "ATOM      2  H1  HOH B   2       2.460   2.500   3.500  1.00  0.00           H\n"
        "ENDMDL\n"
        "END\n"
    )
    filepath = os.path.join(str(tmpdir), "multi_model.pdb")
    with open(filepath, "w") as f:
        f.write(pdb_content)
    return filepath


@pytest.fixture()
def blank_element_pdb_file(tmpdir):
    """PDB file where element columns (77-78) are blank, requiring inference."""
    pdb_content = (
        "HETATM    1 FE   HEM A   1       0.000   0.000   0.000  1.00  0.00\n"
        "HETATM    2 ZN   ZN  A   2       1.000   0.000   0.000  1.00  0.00\n"
        "HETATM    3 CL   CL  A   3       2.000   0.000   0.000  1.00  0.00\n"
        "ATOM      4  CA  ALA A   4       3.000   0.000   0.000  1.00  0.00\n"
        "END\n"
    )
    filepath = os.path.join(str(tmpdir), "blank_elements.pdb")
    with open(filepath, "w") as f:
        f.write(pdb_content)
    return filepath


@pytest.fixture()
def empty_pdb_file(tmpdir):
    """PDB file with no ATOM/HETATM records."""
    filepath = os.path.join(str(tmpdir), "empty.pdb")
    with open(filepath, "w") as f:
        f.write("REMARK  This PDB has no atoms.\nEND\n")
    return filepath


@pytest.fixture()
def water_pdb_file(structure_test_directory):
    """Single-molecule SDF (V2000 molfile) for water."""
    return os.path.join(structure_test_directory, "pdb", "water.pdb")


@pytest.fixture()
def water_pdb_file_content(water_pdb_file):
    """Returns the content of the water.sdf file as a string."""
    with open(water_pdb_file, "r") as f:
        content = f.read()
    return content


############ Molecule and Structure Fixtures ##################
@pytest.fixture()
def a_molecule(single_molecule_xyz_file):
    return Molecule.from_filepath(single_molecule_xyz_file)


@pytest.fixture()
def ethanol_conformers(ethanol_molecule):
    return [ethanol_molecule.copy() for _ in range(4)]


@pytest.fixture()
def ethanol_irc_molecules(ethanol_molecule):
    # A short synthetic "IRC trajectory": 5 copies of the same
    # 9-atom ethanol geometry, standing in for points along a path.
    return [ethanol_molecule.copy() for _ in range(5)]


@pytest.fixture()
def ethanol_traj_molecules(ethanol_molecule):
    mols = []
    for i in range(20):
        mol = ethanol_molecule.copy()
        mol.energy = float(i)
        mols.append(mol)
    return mols


@pytest.fixture()
def simple_ase_atoms():
    return Atoms("Ar2", positions=[(0.0, 0.0, 0.0), (3.5, 0.0, 0.0)])


@pytest.fixture()
def water_with_vib_modes():
    """Water molecule with a single synthetic mass-weighted mode."""
    mol = Molecule(
        symbols=["O", "H", "H"],
        positions=np.array(
            [
                [0.0000, 0.0000, 0.1173],
                [0.0000, 0.7572, -0.4692],
                [0.0000, -0.7572, -0.4692],
            ]
        ),
        vibrational_modes=[
            np.array(
                [
                    [0.0, 0.0, 0.1],
                    [0.0, 0.2, -0.1],
                    [0.0, -0.2, -0.1],
                ]
            )
        ],
    )
    return mol


@pytest.fixture()
def crest_molecules(xyz_directory):
    xyz_file = XYZFile(
        filename=os.path.join(xyz_directory, "crest_conformers.xyz")
    )
    return xyz_file.get_molecules(index=":", return_list=True)


@pytest.fixture
def minimal_struct_dict():
    return {
        "chemical_symbols": ["H", "H"],
        "positions": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]],
        "charge": 0,
        "multiplicity": 1,
        "energy": -1.0,
    }
