"""Tiny runner for the GROMACS topology writer."""

from rdkit import Chem
from rdkit.Chem import AllChem

from chemsmart.io.gromacs import write_gromacs_itp_top
from chemsmart.io.molecules.structure import Molecule


def main() -> None:
    mol = Chem.MolFromSmiles("CCO")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)

    molecule = Molecule.from_rdkit_mol(mol)
    write_gromacs_itp_top(
        molecule, "molecule.itp", "molecule.top", resname="MOL"
    )


if __name__ == "__main__":
    main()
