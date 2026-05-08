# GROMACS Topology Writer

This folder contains a lightweight, RDKit-backed topology writer that builds
minimal `.itp` and `.top` content directly from a `Molecule` instance.

## Usage

```python
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.gromacs import build_gromacs_topology

molecule = Molecule.from_filepath("input.xyz")
topo = build_gromacs_topology(molecule, resname="MOL")
print(topo.itp)
print(topo.top)
```

## Notes

- Atom typing is rule-based and intentionally small. Extend
  `ATOM_TYPE_RULES` for full coverage.
- Charges default to RDKit Gasteiger charges if none are supplied.

