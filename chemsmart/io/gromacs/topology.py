"""Lightweight GROMACS topology writer for Molecule.

This module provides a minimal, extensible skeleton for atom typing and
.itp/.top output. It is intentionally conservative and leaves full parameter
coverage to user-supplied force-field tables.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

from rdkit import Chem

from chemsmart.io.molecules.structure import Molecule

AtomTypeKey = Tuple[str, str, bool]


# Simple starter atom-typing rules (extend per force field)
ATOM_TYPE_RULES: Dict[AtomTypeKey, str] = {
    ("C", "SP3", False): "CT",
    ("C", "SP2", False): "C",
    ("C", "SP2", True): "CA",
    ("N", "SP3", False): "NT",
    ("N", "SP2", False): "N",
    ("O", "SP3", False): "OS",
    ("O", "SP2", False): "O",
    ("S", "SP3", False): "S",
    ("H", "SP3", False): "HC",
}

ELEMENT_FALLBACK: Dict[str, str] = {
    "C": "CT",
    "N": "N",
    "O": "O",
    "S": "S",
    "H": "HC",
}


@dataclass
class GromacsTopology:
    itp: str
    top: str


def _assign_gasteiger_charges(mol: Chem.Mol) -> List[float]:
    Chem.rdPartialCharges.ComputeGasteigerCharges(mol)
    charges: List[float] = []
    for atom in mol.GetAtoms():
        try:
            value = float(atom.GetProp("_GasteigerCharge"))
        except Exception:
            value = 0.0
        if value != value:  # NaN check without extra deps
            value = 0.0
        charges.append(value)
    return charges


def _infer_atom_type(atom: Chem.Atom) -> str:
    element = atom.GetSymbol()
    hybrid = str(atom.GetHybridization()).replace("SP", "SP")
    aromatic = atom.GetIsAromatic()
    key = (element, hybrid, aromatic)
    return ATOM_TYPE_RULES.get(key, ELEMENT_FALLBACK.get(element, element))


def _iter_bonds(graph) -> Iterable[Tuple[int, int]]:
    for i, j in graph.edges:
        if i < j:
            yield i, j


def _iter_angles(graph) -> Iterable[Tuple[int, int, int]]:
    angles = set()
    for j in graph.nodes:
        neighbors = list(graph.neighbors(j))
        for idx_i, i in enumerate(neighbors):
            for k in neighbors[idx_i + 1 :]:
                angles.add((i, j, k))
    for i, j, k in sorted(angles):
        yield i, j, k


def _iter_dihedrals(graph) -> Iterable[Tuple[int, int, int, int]]:
    dihedrals = set()
    for i, j in graph.edges:
        for k in graph.neighbors(j):
            if k == i:
                continue
            for idx_l in graph.neighbors(i):
                if idx_l == j:
                    continue
                dihedrals.add((k, j, i, idx_l))
    for k, j, i, idx_l in sorted(dihedrals):
        yield k, j, i, idx_l


def build_gromacs_topology(
    molecule: Molecule,
    resname: str = "MOL",
    charges: Optional[List[float]] = None,
    add_hs: bool = False,
) -> GromacsTopology:
    """Build minimal .itp/.top content for a Molecule.

    Parameters
    ----------
    molecule
        Molecule object with symbols/positions and valid connectivity.
    resname
        Residue name used in output sections.
    charges
        Optional per-atom charges. If None, Gasteiger charges are assigned
        using RDKit.
    add_hs
        If True, add explicit hydrogens before charge assignment.
    """

    rdkit_mol = molecule.to_rdkit(add_bonds=True)
    if add_hs:
        rdkit_mol = Chem.AddHs(rdkit_mol, addCoords=True)
    graph = molecule.to_graph()

    if charges is None:
        charges = _assign_gasteiger_charges(rdkit_mol)
    if len(charges) != molecule.num_atoms:
        raise ValueError("Charges length must match number of atoms")

    atom_types: List[str] = []
    for atom in rdkit_mol.GetAtoms():
        atom_types.append(_infer_atom_type(atom))

    itp_lines: List[str] = []
    itp_lines.append("[ moleculetype ]")
    itp_lines.append(f"{resname}  3")
    itp_lines.append("")

    itp_lines.append("[ atoms ]")
    itp_lines.append("; nr  type  resnr  residue  atom  cgnr  charge  mass")
    for idx, (atype, charge, mass) in enumerate(
        zip(atom_types, charges, molecule.masses), start=1
    ):
        atom_name = f"{molecule.symbols[idx - 1]}{idx}"
        itp_lines.append(
            f"{idx:5d}  {atype:4s}  1  {resname:4s}  {atom_name:4s}"
            f"  {idx:5d}  {charge:8.4f}  {mass:8.3f}"
        )

    itp_lines.append("")
    itp_lines.append("[ bonds ]")
    itp_lines.append("; i  j  funct")
    for i, j in _iter_bonds(graph):
        itp_lines.append(f"{i + 1:5d} {j + 1:5d}  1")

    itp_lines.append("")
    itp_lines.append("[ angles ]")
    itp_lines.append("; i  j  k  funct")
    for i, j, k in _iter_angles(graph):
        itp_lines.append(f"{i + 1:5d} {j + 1:5d} {k + 1:5d}  1")

    itp_lines.append("")
    itp_lines.append("[ dihedrals ]")
    itp_lines.append("; i  j  k  l  funct")
    for i, j, k, idx_l in _iter_dihedrals(graph):
        itp_lines.append(f"{i + 1:5d} {j + 1:5d} {k + 1:5d} {idx_l + 1:5d}  3")

    top_lines = [
        '#include "oplsaa.ff/forcefield.itp"',
        f'#include "{resname.lower()}.itp"',
        "",
        "[ system ]",
        f"{resname} system",
        "",
        "[ molecules ]",
        f"{resname}  1",
    ]

    return GromacsTopology(itp="\n".join(itp_lines), top="\n".join(top_lines))


def write_gromacs_itp_top(
    molecule: Molecule,
    itp_path: str,
    top_path: str,
    resname: str = "MOL",
    charges: Optional[List[float]] = None,
    add_hs: bool = False,
) -> GromacsTopology:
    """Write .itp/.top files and return the in-memory content."""

    topo = build_gromacs_topology(
        molecule, resname=resname, charges=charges, add_hs=add_hs
    )
    with open(itp_path, "w") as itp_file:
        itp_file.write(topo.itp)
    with open(top_path, "w") as top_file:
        top_file.write(topo.top)
    return topo
