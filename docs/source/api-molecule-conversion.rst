Format Conversion
=================

CHEMSMART provides bidirectional conversion with popular chemistry libraries:

.. note::

   Always start a format-conversion workflow from a **calculation output file**
   (``xTB .out``, ``Gaussian .log``, ``ORCA .out``, etc.). Plain ``.xyz``
   files contain only atomic coordinates — they carry *no* charge or spin
   multiplicity information, so ``charge`` and ``multiplicity`` will be
   ``None`` and downstream conversions to libraries such as **pymatgen** or
   **ASE** will fail with a ``TypeError`` (e.g.
   ``unsupported operand type(s) for -=: 'float' and 'NoneType'``).

.. code:: python

   # xTB output file — contains charge, multiplicity, coordinates, energy, forces
   mol = Molecule.from_filepath(
       "tests/data/XTBTests/outputs/co2_ohess/co2_ohess.out"
   )
   print(mol.chemical_formula)   # CO2
   print(mol.charge, mol.multiplicity)   # 0 1

ASE
---

.. code:: python

   # To ASE Atoms (energy converted from Hartree to eV, forces to eV/Å)
   atoms = mol.to_ase()
   print(type(atoms))   # <class 'ase.Atoms'>

:meth:`~chemsmart.io.molecules.structure.Molecule.to_ase()`
Convert the molecule to an ASE :class:`~ase.Atoms` object. Energy is converted from Hartree to eV and forces from Hartree/Bohr to eV/Å. The result is an :class:`AtomsChargeMultiplicity` instance that preserves ``charge`` and ``multiplicity``.

- Returns: :class:`ase.Atoms` — Wrapped as :class:`AtomsChargeMultiplicity` for round-trip fidelity.

pymatgen
--------

.. code:: python

   # To pymatgen Molecule
   pymatgen_mol = mol.to_pymatgen()
   print(type(pymatgen_mol))   # <class 'pymatgen.core.Molecule'>

:meth:`~chemsmart.io.molecules.structure.Molecule.to_pymatgen()`
Convert to a pymatgen :class:`Molecule` via the ASE adapter. Requires the ``pymatgen`` package.

- Returns: :class:`pymatgen.core.structure.Molecule`

RDKit
-----

.. code:: python

   # To RDKit Mol (with stereochemistry from 3D)
   rdkit_mol = mol.to_rdkit()
   print(type(rdkit_mol))   # <class 'rdkit.Chem.rdchem.Mol'>

:meth:`~chemsmart.io.molecules.structure.Molecule.to_rdkit(add_bonds=True, bond_cutoff_buffer=0.05, adjust_H=True)`
Build an RDKit :class:`Mol` from the 3D geometry. Bonds and bond orders are inferred from interatomic distances; stereochemistry is assigned from 3D. Requires the ``rdkit`` package.

- Returns: :class:`rdkit.Chem.Mol`

SMILES
------

.. code:: python

   # To SMILES string
   smiles = mol.to_smiles()
   print(smiles)   # e.g. "O=C=O"

:meth:`~chemsmart.io.molecules.structure.Molecule.to_smiles()`
Generate a SMILES string by way of RDKit. Requires RDKit.

- Returns: ``str`` — Canonical SMILES, or ``None`` if conversion fails.

PDB
---

.. code:: python

   # To PDB string
   pdb_string = mol.to_pdb()
   print(pdb_string[:80])   # first line

:meth:`~chemsmart.io.molecules.structure.Molecule.to_pdb()`
Return a PDB-format string representation of the molecule.

- Returns: ``str`` — PDB block string.

COSMORS-XYZ
-----------

:meth:`~chemsmart.io.molecules.structure.Molecule.to_cosmorsxyz()`
Return a COSMORS-XYZ format string. Used by COSMO-RS solvation calculations.

- Returns: ``str`` — COSMORS-XYZ block.

NetworkX
--------

.. code:: python

   # To NetworkX graph (nodes = atoms, edges = bonds with bond_order)
   graph = mol.to_graph()
   print(graph.number_of_nodes(), graph.number_of_edges())   # 3 2

:meth:`~chemsmart.io.molecules.structure.Molecule.to_graph(bond_cutoff_buffer=0.3)`
Convert to a NetworkX :class:`Graph` where nodes are atoms and edges are bonds with a ``bond_order`` attribute.

- Returns: :class:`networkx.Graph`

ML Feature Vector
-----------------

.. code:: python

   # To ML feature vector
   X = mol.to_X_data()
   print(len(X))   # feature vector length

:meth:`~chemsmart.io.molecules.structure.Molecule.to_X_data(wbo=False)`
Build a feature vector for ML/descriptor workflows.

- ``wbo`` (bool) — If ``True``, include Wiberg bond-order features.
- Returns: ``list`` — Fixed-length feature vector.

.. note::

   The :meth:`to_ase` method converts energy from Hartree to eV and forces from Hartree/Bohr to eV/Å, following ASE
   conventions. The :meth:`to_rdkit` method infers bonds from interatomic distances and assigns stereochemistry from
   the 3D geometry.