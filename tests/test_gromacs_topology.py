from rdkit import Chem
from rdkit.Chem import AllChem

from chemsmart.io.gromacs import build_gromacs_topology
from chemsmart.io.molecules.structure import Molecule


def _build_rdkit_mol(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDG()
    params.randomSeed = 0
    embed_status = AllChem.EmbedMolecule(mol, params)
    assert embed_status == 0, f"RDKit embedding failed for {smiles!r} with code {embed_status}"
    optimize_status = AllChem.UFFOptimizeMolecule(mol)
    assert optimize_status == 0, f"RDKit UFF optimization failed for {smiles!r} with code {optimize_status}"
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
