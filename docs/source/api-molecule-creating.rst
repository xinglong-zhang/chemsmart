#####################
 Creating a Molecule
#####################

****************************
 From Symbols and Positions
****************************

The simplest way to create a :class:`Molecule` is to specify atomic symbols and Cartesian positions directly:

.. code:: python

   from chemsmart.io.molecules.structure import Molecule
   import numpy as np

   mol = Molecule(
       symbols=["O", "H", "H"],
       positions=np.array([
           [0.0,  0.0,  0.119],
           [0.0,  0.76, -0.477],
           [0.0, -0.76, -0.477],
       ]),
       charge=0,
       multiplicity=1,
   )

   print(mol.chemical_formula)     # H2O
   print(mol.num_atoms)            # 3
   print(mol.mass)                 # ~18.015 amu
   print(mol.get_distance(2, 3))   # distance between H2 and H3 in Å

************
 From Files
************

CHEMSMART supports a wide range of file formats via :meth:`Molecule.from_filepath`:

The following examples use real test files shipped with the CHEMSMART repository (under ``tests/data/``):

.. code:: python

   from chemsmart.io.molecules.structure import Molecule

   # From an XYZ file (methane, 5 atoms)
   mol = Molecule.from_filepath("tests/data/StructuresTests/canonical/methane.xyz")
   print(mol.chemical_formula)   # CH4

   # From a Gaussian output file (benzene, last structure by default)
   mol = Molecule.from_filepath("tests/data/GaussianTests/outputs/benzene.log")
   print(mol.chemical_formula)   # C6H6

   # Use a specific structure from a multi-step optimization (1-based index)
   # Note: index must be a string ("5", "-1", or ":" for all)
   mol = Molecule.from_filepath(
       "tests/data/GaussianTests/outputs/collidine_opt.log", index="5"
   )
   print(f"Structure 5: E = {mol.energy:.6f} Ha")   # -365.455947

   # Compare with the first and last structures
   mol_first = Molecule.from_filepath(
       "tests/data/GaussianTests/outputs/collidine_opt.log", index="1"
   )
   mol_last = Molecule.from_filepath(
       "tests/data/GaussianTests/outputs/collidine_opt.log", index="-1"
   )
   print(f"Structure 1:  E = {mol_first.energy:.6f} Ha")   # -365.446140
   print(f"Structure 13: E = {mol_last.energy:.6f} Ha")   # -365.456178

   # Get all conformers as a list (18 structures in this file)
   mols = Molecule.from_filepath(
       "tests/data/StructuresTests/xyz/crest_conformers.xyz",
       index=":",
       return_list=True,
   )
   print(len(mols))   # 18

.. tip::

   The ``index`` parameter uses **1-based** indexing to match Gaussian / ORCA conventions. Use ``"-1"`` for the last
   structure (default), ``":"`` for all structures, or a specific integer.

**************
 From PubChem
**************

Query PubChem directly by name, CID, or SMILES string:

.. code:: python

   # By CID (Compound ID) — unambiguous, no fallback errors
   mol = Molecule.from_pubchem("2244")
   print(mol.chemical_formula)   # C9H8O4

   # By SMILES string — unambiguous, no fallback errors
   mol = Molecule.from_pubchem("CC(=O)OC1=CC=CC=C1C(=O)O")
   print(mol.chemical_formula)   # C9H8O4

   # By compound name — CHEMSMART tries SMILES first (which fails for a
   # plain name) and then falls back to name search. The intermediate
   # "400 PUGREST.BadRequest" message for the SMILES attempt is normal
   # and the result is still returned correctly:
   mol = Molecule.from_pubchem("aspirin")
   print(mol.chemical_formula)   # C9H8O4

.. note::

   PubChem queries require an **internet connection** and the ``tenacity`` Python package (``pip install tenacity``). If
   a 3D conformer is unavailable, CHEMSMART automatically falls back to the 2D structure and generates 3D coordinates
   using RDKit (``pip install rdkit``).

.. tip::

   To avoid the harmless-but-confusing "400 BadRequest" logger message for name-based queries, prefer CID or SMILES
   identifiers when running in a notebook or CI.

***************************
 From Other Python Objects
***************************

CHEMSMART provides seamless conversion from popular chemistry libraries:

.. code:: python

   from chemsmart.io.molecules.structure import Molecule

   # From ASE Atoms
   from ase import Atoms
   atoms = Atoms("H2O", positions=[[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]])
   mol = Molecule.from_ase_atoms(atoms)
   print(mol.chemical_formula)   # H2O

   # Deep copy (works for any Molecule)
   mol_copy = mol.copy()
   print(mol_copy.chemical_formula)   # H2O

   # From RDKit Mol
   from rdkit import Chem
   from rdkit.Chem import AllChem
   rdkit_mol = Chem.MolFromSmiles("CCO")
   AllChem.EmbedMolecule(rdkit_mol)
   mol = Molecule.from_rdkit_mol(rdkit_mol)
   print(mol.chemical_formula)   # C2H6O

   # From a coordinate block text (Gaussian-style)
   mol = Molecule.from_coordinate_block_text(
       "C  0.0  0.0  0.0\n"
       "O  1.2  0.0  0.0\n"
       "O -1.2  0.0  0.0"
   )
   print(mol.chemical_formula)   # CO2

.. important::

   :doc:`molecule-input-formats` for a full overview of supported input formats and CLI usage.

************************
 Constructor Parameters
************************

The constructor accepts the following parameters. Parameters typed in **bold** are required; the rest are optional and
are typically populated automatically when reading calculation output files (Gaussian ``.log``, ORCA ``.out``, xTB,
etc.).

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   -  -  Parameter
      -  Type
      -  Description

   -  -  **symbols**
      -  list
      -  List of atomic symbols (e.g. ``["C", "O", "O"]``).

   -  -  **positions**
      -  ndarray
      -  ``(N, 3)`` array of Cartesian coordinates in Ångströms.

   -  -  ``charge``
      -  int
      -  Molecular charge

   -  -  ``multiplicity``
      -  int
      -  Spin multiplicity (2S + 1)

   -  -  ``frozen_atoms``
      -  list
      -  Per-atom constraint flags: ``-1`` = frozen, ``0`` = relaxed (Gaussian convention)

   -  -  ``pbc_conditions``
      -  list
      -  Periodic boundary conditions ``[x, y, z]``

   -  -  ``translation_vectors``
      -  list
      -  Cell / translation vectors for periodic systems

   -  -  ``energy``
      -  float
      -  Total energy in Hartree

   -  -  ``forces``
      -  ndarray
      -  ``(N, 3)`` array of forces in Hartree / Bohr

   -  -  ``velocities``
      -  ndarray
      -  ``(N, 3)`` array of atomic velocities

   -  -  ``vibrational_frequencies``
      -  list
      -  Vibrational frequencies in cm⁻¹

   -  -  ``vibrational_modes``
      -  list
      -  Each entry is an ``(N, 3)`` mass-weighted normal-mode displacement

   -  -  ``dipole_moment``
      -  ndarray
      -  Dipole moment vector ``[X, Y, Z]`` in Debye

   -  -  ``point_group``
      -  str
      -  Point group label (e.g. ``"C2V"``, ``"CS"``)

   -  -  ``info``
      -  dict
      -  Arbitrary extra metadata stored with the molecule
