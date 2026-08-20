########################
 Working with Molecules
########################

*********************
 Indexing and Copies
*********************

CHEMSMART uses **1-based indexing** for atom indices throughout the API, matching the convention used by Gaussian and
ORCA input files.

.. code:: python

   import numpy as np
   from chemsmart.io.molecules.structure import Molecule

   # Create a water molecule for indexing examples
   mol = Molecule(
       symbols=["O", "H", "H"],
       positions=np.array([
           [0.0,  0.0,  0.119],
           [0.0,  0.76, -0.477],
           [0.0, -0.76, -0.477],
       ]),
   )

   # Number of atoms
   print(len(mol))                   # 3

   # Subset of atoms (1-based indices)
   subset = mol[[1, 2]]             # atoms 1 and 2
   print(subset.chemical_formula)   # 'HO' (Hill notation: alphabetical when no C)

   # Deep copy
   mol_copy = mol.copy()
   print("Copy formula:", mol_copy.chemical_formula)

   # Delete atoms by 1-based index
   # Delete atom 3 (the second H), keeping OH
   oh = mol.delete_atoms_by_indices(atom_indices=3)
   print("After deleting atom 3:", oh.chemical_formula)   # HO
   # Delete atoms 2 and 3 (both H atoms), keeping just O
   o_only = mol.delete_atoms_by_indices(atom_indices=[2, 3])
   print("After deleting atoms 2,3:", o_only.chemical_formula)   # O
   # Delete atom 0 (0-based, with one_based=False) — removes O, keeping H2
   h2 = mol.delete_atoms_by_indices(atom_indices=0, one_based=False)
   print("After deleting atom 0 (0-based):", h2.chemical_formula)   # H2

.. warning::

   CHEMSMART uses **1-based indexing** to match most molecular visualization software, unlike Python's 0-based indexing.

:meth:`~chemsmart.io.molecules.structure.Molecule.copy()` Returns a deep copy of the molecule.

-  Returns: :class:`Molecule` — A new molecule with independent copies of all arrays and metadata.

:meth:`~chemsmart.io.molecules.structure.Molecule.delete_atoms_by_indices(atom_indices, one_based=True)` Remove atoms by
index and return a new molecule.

-  ``atom_indices`` (int or iterable of int) — Atom indices to remove. Default is **1-based** (matching Gaussian/Orca
   conventions).
-  ``one_based`` (bool) — If ``True`` (default), indices are 1-based. Set to ``False`` for 0-based Python-style
   indexing.
-  Returns: :class:`Molecule` — New molecule with the specified atoms removed.
-  Raises: ``ValueError`` if any index is out of range or removing all atoms would leave an empty molecule.

.. warning::

   CHEMSMART uses **1-based indexing** to match most molecular visualization software, unlike Python's 0-based indexing.

****************************
 Writing Molecules to Files
****************************

Use :meth:`Molecule.write` to export molecules to various formats:

:meth:`~chemsmart.io.molecules.structure.Molecule.write(filename, format="xyz", mode="w", **kwargs)` Write the molecule
to a file in one of the supported formats.

.. code:: python

   import numpy as np
   from chemsmart.io.molecules.structure import Molecule

   # Create a simple molecule for writing examples
   mol = Molecule(
       symbols=["C", "O", "O"],
       positions=np.array([[0, 0, 0], [0, 0, 1.16], [0, 0, -1.16]]),
       charge=0,
       multiplicity=1,
   )

   # Write to XYZ
   mol.write("output.xyz")
   print("Wrote output.xyz")

   # Write to extended XYZ (includes energy and forces)
   mol.write("output.extxyz", format="extxyz")
   print("Wrote output.extxyz")

   # Write to Gaussian input file
   mol.write("input.com", format="com")
   print("Wrote input.com")

   # Write to PDB
   mol.write("output.pdb", format="pdb")
   print("Wrote output.pdb")

Supported formats:

+--------------------+------------------------------------------+
| Format             | Description                              |
+====================+==========================================+
| ``"xyz"``          | Standard XYZ format (default)            |
+--------------------+------------------------------------------+
| ``"extxyz"``       | Extended XYZ with energy, forces,        |
|                    | metadata                                 |
+--------------------+------------------------------------------+
| ``"com"``          | Gaussian input file with route section   |
+--------------------+------------------------------------------+
| ``"pdb"``          | Protein Data Bank format                 |
+--------------------+------------------------------------------+

For formats not listed above (e.g. COSMORS-XYZ, MOL), use the dedicated writer methods directly:

.. code:: python

   # COSMORS-XYZ (direct method, not via write())
   mol.write_cosmorsxyz("output.cosmorsxyz")
   print("Wrote output.cosmorsxyz")

   # RDKit's MOL block (via RDKit conversion)
   from rdkit import Chem
   rdkit_mol = mol.to_rdkit()
   molblock = Chem.MolToMolBlock(rdkit_mol)
   with open("output.mol", "w") as f:
       f.write(molblock)
   print("Wrote output.mol")

-  ``filename`` (str) — Output file path.
-  ``format`` (str) — Output format. Defaults to ``"xyz"``.
-  ``mode`` (str) — File write mode (default: ``"w"``).
-  Raises: ``ValueError`` if the requested format is not supported.

**Dedicated writer methods:**

-  ``write_xyz(filename, mode="w")`` — Write standard XYZ.
-  ``write_extxyz(filename, mode="w")`` — Write extended XYZ.
-  ``write_com(filename)`` — Write Gaussian input (``.com``).
-  ``write_pdb(filename, mode="w")`` — Write PDB format.
-  ``write_cosmorsxyz(filename, mode="w")`` — Write COSMORS-XYZ format.
-  ``write_coordinates(filename, mode="w", atom_indices=None)`` — Write a text coordinate block.
-  ``write_pdb_pybabel(filename)`` — Write PDB via PyBabel.

*******************
 Format Conversion
*******************

CHEMSMART provides bidirectional conversion with popular chemistry libraries (ASE, pymatgen, RDKit, NetworkX, etc.) via
a set of :meth:`to_*` methods.

.. note::

   Always start a format-conversion workflow from a **calculation output file** (``xTB .out``, ``Gaussian .log``, ``ORCA
   .out``, etc.). Plain ``.xyz`` files carry no charge or spin multiplicity, so downstream conversions may fail.

.. code:: python

   # Start from a calculation output (has charge, multiplicity, energy, forces)
   mol = Molecule.from_filepath(
       "tests/data/XTBTests/outputs/co2_ohess/co2_ohess.out"
   )
   print(mol.chemical_formula)   # CO2

   # ASE Atoms (energy → eV, forces → eV/Å)
   atoms = mol.to_ase()
   print(type(atoms).__name__)          # AtomsChargeMultiplicity

   # pymatgen Molecule
   pymatgen_mol = mol.to_pymatgen()
   print(type(pymatgen_mol).__name__)   # Molecule

   # RDKit Mol (bonds & stereo from 3D)
   rdkit_mol = mol.to_rdkit()
   print(type(rdkit_mol).__name__)      # Mol

   # SMILES string
   print(mol.to_smiles())               # O=C=O

   # PDB string
   print(mol.to_pdb()[:50])             # first line of PDB block

   # NetworkX graph (nodes = atoms, edges = bonds)
   graph = mol.to_graph()
   print(graph.number_of_nodes(), graph.number_of_edges())   # 3 2

   # ML feature vector
   X = mol.to_X_data()
   print(len(X))                        # feature vector length

**Conversion methods:** :meth:`to_ase`, :meth:`to_pymatgen`, :meth:`to_rdkit`, :meth:`to_smiles`, :meth:`to_pdb`,
:meth:`to_cosmorsxyz`, :meth:`to_graph`, :meth:`to_X_data`.
