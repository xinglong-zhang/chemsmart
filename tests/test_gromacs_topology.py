from rdkit import Chem
from rdkit.Chem import AllChem

from chemsmart.io.gromacs import build_gromacs_topology
from chemsmart.io.molecules.structure import Molecule


def _build_rdkit_mol(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)
    return mol


def test_gromacs_topology_from_molecule():
    rdkit_mol = _build_rdkit_mol("CCO")
    molecule = Molecule.from_rdkit_mol(rdkit_mol)

    topo = build_gromacs_topology(molecule, resname="MOL")

    assert "[ atoms ]" in topo.itp
    assert "[ bonds ]" in topo.itp
    assert "[ system ]" in topo.top

    lines = topo.itp.splitlines()
    atoms_start = lines.index("[ atoms ]")
    next_section = next(
        idx
        for idx in range(atoms_start + 1, len(lines))
        if lines[idx].startswith("[")
    )
    atom_lines = [
        line
        for line in lines[atoms_start + 1 : next_section]
        if line.strip() and line.strip()[0].isdigit()
    ]
    assert len(atom_lines) == molecule.num_atoms
