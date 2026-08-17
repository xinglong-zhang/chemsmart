############################
 The Molecule Object (API)
############################

The :class:`~chemsmart.io.molecules.structure.Molecule` object is the central data structure in CHEMSMART for
representing a molecular structure. It stores atomic symbols, 3D coordinates, electronic properties (charge,
multiplicity), calculated properties (energy, forces, vibrations), and provides a rich set of methods for geometry
analysis, format conversion, and interop with popular chemistry libraries (RDKit, ASE, pymatgen, Open Babel).

This page is intended for Python users who want to use CHEMSMART programmatically via its API. For the command-line
interface, see :doc:`cli-overview` and :doc:`molecule-input-formats`.

.. note::

   For more information about CHEMSMART's design and capabilities, see the preprint: https://arxiv.org/abs/2508.20042

************************
 Creating a Molecule
************************

From Symbols and Positions
==========================

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

From Files
==========

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
   # Note: ``index`` must be a string (``"5"``, ``"-1"``, or ``":"`` for all)
   mol = Molecule.from_filepath(
       "tests/data/GaussianTests/outputs/collidine_opt.log", index="5"
   )
   print(mol.chemical_formula)   # C8H11N

   # Get all conformers as a list (18 structures in this file)
   mols = Molecule.from_filepath(
       "tests/data/StructuresTests/xyz/crest_conformers.xyz",
       index=":",
       return_list=True,
   )
   print(len(mols))   # 18

.. list-table::
   :header-rows: 1
   :widths: 25 75

   -  -  Extension
      -  Description

   -  -  ``.xyz``
      -  Standard XYZ format

   -  -  ``.sdf``
      -  MDL MOL / SDF format

   -  -  ``.pdb``
      -  Protein Data Bank format

   -  -  ``.com``, ``.gjf``
      -  Gaussian input files

   -  -  ``.log``
      -  Gaussian output files (multi-step trajectories supported)

   -  -  ``.inp``
      -  ORCA input files

   -  -  ``.out``
      -  ORCA / xTB / Gaussian output (auto-detected by file header)

   -  -  ``.cdx``, ``.cdxml``
      -  ChemDraw files (binary or XML)

   -  -  ``.gro``, ``.trr``
      -  GROMACS files

   -  -  ``.db``
      -  CHEMSMART database or ASE database (auto-detected)

   -  -  *other*
      -  Falls back to ASE :func:`ase.io.read`

.. tip::

   The ``index`` parameter uses **1-based** indexing to match Gaussian / ORCA conventions. Use ``"-1"`` for the last
   structure (default), ``":"`` for all structures, or a specific integer.

From a Calculation Directory
============================

For xTB calculations, use :meth:`Molecule.from_directorypath`. The example below uses a real
xTB Hessian calculation directory shipped with the CHEMSMART test suite:

.. code:: python

   # CO2 molecule from a real xTB calculation directory
   mol = Molecule.from_directorypath(
       "tests/data/XTBTests/outputs/co2_ohess/", program="xtb"
   )
   print(mol.chemical_formula)   # CO2
   print(mol.num_atoms)          # 3

From PubChem
============

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

   PubChem queries require an **internet connection** and the ``tenacity``
   Python package (``pip install tenacity``). If a 3D conformer is
   unavailable, CHEMSMART automatically falls back to the 2D structure and
   generates 3D coordinates using RDKit (``pip install rdkit``).

.. tip::

   To avoid the harmless-but-confusing "400 BadRequest" logger message for
   name-based queries, prefer CID or SMILES identifiers when running in a
   notebook or CI.

From Other Python Objects
=========================

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

.. seealso::

   :doc:`molecule-input-formats` for a full overview of supported input formats and CLI usage.

******************************
 Constructor Parameters
******************************

The constructor accepts the following parameters. Parameters typed in
**bold** are required; the rest are optional and are typically populated
automatically when reading calculation output files (Gaussian ``.log``,
ORCA ``.out``, xTB, etc.).

.. automethod:: chemsmart.io.molecules.structure.Molecule.__init__
   :noindex:

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

***********************
 Indexing and Copies
***********************

CHEMSMART uses **1-based indexing** for atom indices throughout the API, matching the convention used by Gaussian and
ORCA input files.

.. code:: python

   # Number of atoms
   len(mol)                         # 3

   # Subset of atoms (1-based indices)
   subset = mol[[1, 2]]             # atoms 1 and 2
   subset.chemical_formula          # 'HO' (Hill notation: alphabetical when no C)

   # Deep copy
   mol_copy = mol.copy()

   # Delete atoms by 1-based index
   phenoxide = mol.delete_atoms_by_indices(atom_indices=13)
   phenoxide = mol.delete_atoms_by_indices(atom_indices=[13])
   phenoxide = mol.delete_atoms_by_indices(atom_indices=12, one_based=False)

.. warning::

   CHEMSMART uses **1-based indexing** to match most molecular visualization software, unlike Python's 0-based indexing.

.. automethod:: chemsmart.io.molecules.structure.Molecule.copy

.. automethod:: chemsmart.io.molecules.structure.Molecule.delete_atoms_by_indices

***************************
 Writing Molecules to Files
***************************

Use :meth:`Molecule.write` to export molecules to various formats:

.. code:: python

   # Write to XYZ
   mol.write("output.xyz")

   # Write to extended XYZ (includes energy and forces)
   mol.write("output.extxyz", format="extxyz")

   # Write to Gaussian input file
   mol.write("input.com", format="com")

   # Write to PDB
   mol.write("output.pdb", format="pdb")

   # Write to COSMORS-XYZ
   mol.write("output.cosmorsxyz", format="cosmorsxyz")

.. list-table::
   :header-rows: 1
   :widths: 20 80

   -  -  Format
      -  Description

   -  -  ``"xyz"``
      -  Standard XYZ format

   -  -  ``"extxyz"``
      -  Extended XYZ with energy, forces, and metadata

   -  -  ``"com"``
      -  Gaussian input file (with route section)

   -  -  ``"pdb"``
      -  Protein Data Bank format with bond connectivity

.. automethod:: chemsmart.io.molecules.structure.Molecule.write

.. automethod:: chemsmart.io.molecules.structure.Molecule.write_xyz

.. automethod:: chemsmart.io.molecules.structure.Molecule.write_extxyz

.. automethod:: chemsmart.io.molecules.structure.Molecule.write_com

.. automethod:: chemsmart.io.molecules.structure.Molecule.write_pdb

.. automethod:: chemsmart.io.molecules.structure.Molecule.write_pdb_pybabel

.. automethod:: chemsmart.io.molecules.structure.Molecule.write_cosmorsxyz

.. automethod:: chemsmart.io.molecules.structure.Molecule.write_coordinates

************************
 Format Conversion
************************

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

   # To ASE Atoms (energy converted from Hartree to eV, forces to eV/Å)
   atoms = mol.to_ase()
   print(type(atoms))   # <class 'ase.Atoms'>

   # To pymatgen Molecule
   pymatgen_mol = mol.to_pymatgen()
   print(type(pymatgen_mol))   # <class 'pymatgen.core.Molecule'>

   # To RDKit Mol (with stereochemistry from 3D)
   rdkit_mol = mol.to_rdkit()
   print(type(rdkit_mol))   # <class 'rdkit.Chem.rdchem.Mol'>

   # To SMILES string
   smiles = mol.to_smiles()
   print(smiles)   # e.g. "O=C=O"

   # To PDB string
   pdb_string = mol.to_pdb()
   print(pdb_string[:80])   # first line

   # To NetworkX graph (nodes = atoms, edges = bonds with bond_order)
   graph = mol.to_graph()
   print(graph.number_of_nodes(), graph.number_of_edges())   # 3 2

   # To ML feature vector
   X = mol.to_X_data()
   print(len(X))   # feature vector length

.. automethod:: chemsmart.io.molecules.structure.Molecule.to_ase

.. automethod:: chemsmart.io.molecules.structure.Molecule.to_pymatgen

.. automethod:: chemsmart.io.molecules.structure.Molecule.to_rdkit

.. automethod:: chemsmart.io.molecules.structure.Molecule.to_smiles

.. automethod:: chemsmart.io.molecules.structure.Molecule.to_pdb

.. automethod:: chemsmart.io.molecules.structure.Molecule.to_cosmorsxyz

.. automethod:: chemsmart.io.molecules.structure.Molecule.to_graph

.. automethod:: chemsmart.io.molecules.structure.Molecule.to_X_data

.. note::

   The :meth:`to_ase` method converts energy from Hartree to eV and forces from Hartree/Bohr to eV/Å, following ASE
   conventions. The :meth:`to_rdkit` method infers bonds from interatomic distances and assigns stereochemistry from
   the 3D geometry.

***********************
 Molecular Identifiers
***********************

CHEMSMART provides several identifiers for molecules and structures:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   -  -  Property
      -  Description

   -  -  :attr:`~Molecule.chemical_formula`
      -  Chemical formula in Hill notation (e.g. ``"H2O"``)

   -  -  :attr:`~Molecule.empirical_formula`
      -  Empirical formula in Hill notation

   -  -  :attr:`~Molecule.smiles`
      -  SMILES string via RDKit

   -  -  :attr:`~Molecule.cxsmiles`
      -  CXSMILES string (with 3D coordinates and stereo)

   -  -  :attr:`~Molecule.inchi`
      -  Full InChI string (requires Open Babel)

   -  -  :attr:`~Molecule.inchikey`
      -  InChIKey string (requires Open Babel)

   -  -  :attr:`~Molecule.molecule_id`
      -  Unique molecular identifier (= InChIKey)

   -  -  :attr:`~Molecule.structure_id`
      -  SHA-256 hash of canonical geometry + charge + multiplicity

   -  -  :attr:`~Molecule.structure_label`
      -  Human-readable label: ``"str-H2O-a1b2c3d4e5f6"``

   -  -  :attr:`~Molecule.molecule_label`
      -  Human-readable label: ``"mol-H2O-UHOVQNZJYSORNB-UHFFFAOYSA-N"``

.. code:: python

   print(mol.chemical_formula)       # H2O
   print(mol.smiles)                 # O
   print(mol.inchikey)               # UHOVQNZJYSORNB-UHFFFAOYSA-N
   print(mol.structure_id)           # a1b2c3d4e5f6...
   print(mol.structure_label)        # str-H2O-a1b2c3d4e5f6

**Attributes:**

-  ``smiles`` (str or None) — SMILES string (convenience, equivalent to ``to_smiles()``).
-  ``cxsmiles`` (str or None) — CXSMILES string (extended SMILES with 3D/stereo info, via RDKit).
-  ``inchi`` (str or None) — InChI string (via RDKit).
-  ``inchikey`` (str or None) — InChIKey string (27 characters, via Open Babel).
-  ``molecule_id`` (str or None) — Unique molecular identifier (InChIKey).
-  ``molecule_label`` (str or None) — Human-readable label, e.g. ``"mol-C6H6-UHOVQNZJYSORNB-UHFFFAOYSA-N"``.
-  ``structure_id`` (str) — SHA-256 hex digest of canonical geometry + charge + multiplicity.
-  ``structure_label`` (str) — Human-readable label, e.g. ``"str-C6H6-a1b2c3d4e5f6"``.
-  ``chemical_formula`` (str) — Chemical formula in Hill notation.
-  ``empirical_formula`` (str) — Empirical formula in Hill notation.
-  ``elements`` (list) — Sorted unique element symbols, e.g. ``["C", "O"]``.
-  ``element_counts`` (dict) — Per-element counts, e.g. ``{"C": 1, "O": 2}``.

.. automethod:: chemsmart.io.molecules.structure.Molecule.get_chemical_formula

***********************************
 Atomic and Molecular Properties
***********************************

Counts and Symbols
===================

.. code:: python

   # Use a real calculation output so all metadata is populated
   mol = Molecule.from_filepath(
       "tests/data/XTBTests/outputs/co2_ohess/co2_ohess.out"
   )

   print(mol.num_atoms)            # 3
   print(list(mol.symbols))        # ['C', 'O', 'O']
   print(mol.chemical_symbols)     # ['C', 'O', 'O']
   print(mol.positions.shape)      # (3, 3)
   print(mol.positions)            # Cartesian coordinates in Å

**Attributes:**

-  ``num_atoms`` (int) — Number of atoms in the molecule.
-  ``symbols`` (Symbols) — A ``Symbols`` object wrapping the chemical symbols.
-  ``positions`` (numpy.ndarray) — ``(N, 3)`` array of Cartesian coordinates in Ångströms.
-  ``chemical_symbols`` (list) — List of chemical symbol strings (e.g. ``["C", "O", "O"]``).

Masses
=======

CHEMSMART supports three mass conventions:

-  **Standard atomic masses** (:attr:`mass`, :attr:`masses`)
-  **Natural abundance weighted masses** (:attr:`natural_abundance_weighted_mass`, :attr:`natural_abundance_weighted_masses`)
-  **Most abundant isotope masses** (:attr:`most_abundant_mass`, :attr:`most_abundant_masses`)

.. code:: python

   print(mol.mass)                              # ~18.015 amu
   print(mol.natural_abundance_weighted_mass)  # weighted by isotope abundance
   print(mol.most_abundant_mass)                # using most abundant isotopes
   print(mol.masses)                            # list of per-atom masses
   print(mol.most_abundant_masses)              # per-atom most-abundant isotope masses

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.mass

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.masses

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.natural_abundance_weighted_mass

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.natural_abundance_weighted_masses

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.most_abundant_mass

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.most_abundant_masses

Energy, Forces, and Velocities
================================

These are populated from calculation output files (Gaussian ``.log``,
ORCA ``.out``, xTB, etc.).

.. code:: python

   # Total energy (Hartree)
   print(mol.energy)                  # e.g. -4.17657788 Hartree for CO2 at xTB

   # Forces on each atom (Hartree/Bohr), shape (N, 3)
   print(mol.forces.shape)            # (3, 3) for CO2
   print(mol.forces)                  # force components per atom

   # Atomic velocities (if available, N, 3)
   if mol.velocities is not None:
       print(mol.velocities.shape)

**Attributes:**

-  ``energy`` (float or None) — Total molecular energy in Hartree.
-  ``forces`` (numpy.ndarray or None) — ``(N, 3)`` array of forces in Hartree/Bohr.
-  ``velocities`` (numpy.ndarray or None) — ``(N, 3)`` array of atomic velocities.

.. note::

   Energy is stored in **Hartree** and forces in **Hartree/Bohr**, following quantum chemistry conventions. Use
   :meth:`to_ase` to convert to eV and eV/Å respectively.

Electronic and Thermochemical Properties
=========================================

These properties are typically populated when reading calculation output files.
The example below uses a real Gaussian output file:

.. code:: python

   mol = Molecule.from_filepath("tests/data/GaussianTests/outputs/co2.log")

   # Basic electronic properties
   print(mol.charge, mol.multiplicity)          # 0 1

   # Dipole moment
   if mol.dipole_moment is not None:
       print("Dipole (Debye):", mol.dipole_moment)         # [X, Y, Z] components
       print("|μ|:", mol.dipole_moment_magnitude)           # total magnitude

   # Point group / rotational constants / symmetry number
   print("Point group:", mol.point_group)                    # e.g. "D∞H" or "CS"
   print("Rot. const. (Hz):", mol.rotational_constants)     # [A, B, C]
   print("Rot. symm. #:", mol.rotational_symmetry_number)   # e.g. 2

   # Per-atom Mulliken charges and spin densities
   if mol.mulliken_atomic_charges is not None:
       print("Mulliken charges:", mol.mulliken_atomic_charges)  # {"C1": 0.84, "O2": -0.42, ...}
   if mol.mulliken_spin_densities is not None:
       print("Mulliken spins:", mol.mulliken_spin_densities)

   # Optimization flag
   print("Is optimized?", mol.is_optimized_structure)       # True / False / None

**Attributes:**

-  ``charge`` (int or None) — Molecular charge.
-  ``multiplicity`` (int or None) — Spin multiplicity (2S + 1).
-  ``dipole_moment`` (numpy.ndarray or None) — Dipole moment ``[X, Y, Z]`` in Debye.
-  ``dipole_moment_magnitude`` (float or None) — Total dipole moment magnitude in Debye.
-  ``rotational_constants`` (numpy.ndarray or None) — Rotational constants ``[A, B, C]`` in Hz.
-  ``point_group`` (str or None) — Molecular point group (e.g. ``"CS"``, ``"C2V"``).
-  ``rotational_symmetry_number`` (int or None) — Rotational symmetry number.
-  ``mulliken_atomic_charges`` (dict or None) — Per-atom Mulliken charges keyed like ``"O1"``, ``"C2"``.
-  ``mulliken_spin_densities`` (dict or None) — Per-atom Mulliken spin densities.
-  ``is_optimized_structure`` (bool or None) — Whether the structure is an optimized / final geometry.

Atomic Radii
=============

.. code:: python

   # Covalent radii (Å)
   print(mol.atomic_radii_list)      # e.g. [0.77, 0.73, 0.73] for CO2
   # Van der Waals radii (Å)
   print(mol.vdw_radii_list)         # e.g. [1.70, 1.52, 1.52] for CO2

-  ``atomic_radii_list`` (list) — Covalent radii (Å) for each atom, from ASE's ``covalent_radii`` table.
-  ``vdw_radii_list`` (list) — Van der Waals radii (Å) for each atom, from ASE's ``vdW_radii`` table.

**************************
 Structure Classification
**************************

CHEMSMART provides several boolean properties for classifying molecular structure. These are based on RDKit's
perception of the 3D geometry (bonds, rings, stereochemistry).

.. list-table::
   :header-rows: 1
   :widths: 25 75

   -  -  Property
      -  Description

   -  -  :attr:`~Molecule.is_aromatic`
      -  ``True`` if any atom is aromatic and belongs to a ring

   -  -  :attr:`~Molecule.is_ring`
      -  ``True`` if the molecule contains any ring

   -  -  :attr:`~Molecule.is_monoatomic`
      -  ``True`` for single-atom molecules

   -  -  :attr:`~Molecule.is_diatomic`
      -  ``True`` for two-atom molecules

   -  -  :attr:`~Molecule.is_linear`
      -  ``True`` if all atoms are collinear

   -  -  :attr:`~Molecule.is_chiral`
      -  ``True`` if the molecule has chiral centers

   -  -  :attr:`~Molecule.is_multicomponent`
      -  ``True`` if the molecule has multiple fragments (salts, complexes)

.. code:: python

   print(mol.is_aromatic)         # False
   print(mol.is_linear)           # True
   print(mol.is_chiral)           # False
   print(mol.chiral_centers)      # {}
   print(mol.num_components)      # 1

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.is_aromatic

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.is_ring

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.is_monoatomic

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.is_diatomic

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.is_linear

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.is_chiral

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.chiral_centers

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.is_multicomponent

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.num_components

.. note::

   Bond orders are inferred from 3D geometry, not from an electronic structure calculation, so aromaticity detection is
   **model-dependent**. For borderline or unusual systems, the result should be treated as a heuristic estimate.

Periodic Boundary Conditions
=============================

.. code:: python

   import numpy as np
   from chemsmart.io.molecules.structure import Molecule

   # Build a periodic molecule with an orthorhombic cell
   symbols = ["Si", "Si"]
   positions = np.array([[0.0, 0.0, 0.0], [1.36, 1.36, 1.36]])
   cell = np.array([[2.72, 0.0, 0.0], [0.0, 2.72, 0.0], [0.0, 0.0, 2.72]])

   mol = Molecule(
       symbols=symbols,
       positions=positions,
       pbc_conditions=[True, True, True],
       translation_vectors=cell.tolist(),
   )

   print(mol.pbc)                    # [True, True, True]
   print(mol.pbc_conditions)          # [True, True, True]
   print(mol.translation_vectors)     # 3x3 cell matrix

**Attributes:**

-  ``pbc`` (bool) — ``True`` if the molecule has no periodic boundary conditions (all entries in ``pbc_conditions`` are 0).
-  ``pbc_conditions`` (list) — Periodic boundary conditions ``[x, y, z]``.
-  ``translation_vectors`` (list) — Cell / translation vectors for periodic systems.

Frozen Atoms
=============

.. code:: python

   from chemsmart.io.molecules.structure import Molecule
   import numpy as np

   # Build a CO2 molecule with the two oxygen atoms frozen
   symbols = ["C", "O", "O"]
   positions = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.16], [0.0, 0.0, -1.16]])

   mol = Molecule(
       symbols=symbols,
       positions=positions,
       frozen_atoms=[0, -1, -1],  # 0 = relax C, -1 = freeze O atoms
   )

   print(mol.frozen_atoms)            # [0, -1, -1]
   frozen_indices = [i+1 for i, f in enumerate(mol.frozen_atoms) if f == -1]
   print("Frozen atoms:", frozen_indices)  # [2, 3]

-  ``frozen_atoms`` (list) — Per-atom constraint flags following the Gaussian convention: ``-1`` = frozen, ``0`` = relaxed.

.. note::

   Frozen atom flags follow the Gaussian convention: ``-1`` denotes a frozen atom and ``0`` denotes a relaxed atom.

Miscellaneous
=============

.. code:: python

   mol = Molecule.from_filepath(
       "tests/data/XTBTests/outputs/co2_ohess/co2_ohess.out"
   )

   # User-defined metadata dict
   if mol.info is not None:
       print("Info keys:", list(mol.info.keys()))

   # 1-based position of this structure within the source file
   print("Structure index in file:", mol.structure_index_in_file)  # e.g. 1

-  ``info`` (dict) — Arbitrary extra metadata stored with the molecule.
-  ``structure_index_in_file`` (int or None) — 1-based position of this structure within the source file.

************************
 Geometry Analysis
************************

All geometry methods use **1-based** atom indices, matching Gaussian / ORCA conventions.

Distances, Angles, and Dihedrals
=================================

.. code:: python

   # Bond distance between atoms 1 and 2 (in Å)
   dist = mol.get_distance(1, 2)
   print(f"Distance 1-2: {dist:.3f} Å")

   # Bond angle H2–O1–H3 (in degrees)
   angle = mol.get_angle(2, 1, 3)
   print(f"Angle 2-1-3: {angle:.1f}°")

   # Dihedral angle about bond 2–3 (in degrees)
   dihedral = mol.get_dihedral(1, 2, 3, 4)
   print(f"Dihedral 1-2-3-4: {dihedral:.1f}°")

   # All pairwise distances
   all_distances = mol.get_all_distances()
   print(f"All distances: {len(all_distances)} pairs")

   # Pairwise distance matrix (N × N)
   dist_matrix = mol.distance_matrix
   print(f"Distance matrix shape: {dist_matrix.shape}")

.. automethod:: chemsmart.io.molecules.structure.Molecule.get_distance

.. automethod:: chemsmart.io.molecules.structure.Molecule.get_angle

.. automethod:: chemsmart.io.molecules.structure.Molecule.get_angle_from_positions

.. automethod:: chemsmart.io.molecules.structure.Molecule.get_dihedral

.. automethod:: chemsmart.io.molecules.structure.Molecule.get_dihedral_from_positions

.. automethod:: chemsmart.io.molecules.structure.Molecule.get_all_distances

.. automethod:: chemsmart.io.molecules.structure.Molecule.bond_lengths

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.distance_matrix

Center of Mass
==============

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.center_of_mass

Canonical Geometry
===================

These properties yield representations invariant to translation, rotation, and atom-index permutation:

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.canonical_positions

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.canonical_geometry

.. note::

   Canonical positions are translated to the centre of mass and rotated into the principal-axes frame of the
   moment-of-inertia tensor. The :attr:`structure_id` is a SHA-256 hash of the canonical geometry, charge, and
   multiplicity.

Moments of Inertia and Rotational Temperatures
===============================================

Moments of inertia are in units of amu·Å². Rotational temperatures are in Kelvin.

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.moments_of_inertia_tensor

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.moments_of_inertia

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.moments_of_inertia_weighted_mass

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.moments_of_inertia_most_abundant_mass

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.moments_of_inertia_principal_axes

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.rotational_temperatures

.. tip::

   For linear molecules, the moment of inertia along the molecular axis is zero, so the axial rotational temperature
   is infinite. CHEMSMART returns only the finite perpendicular-axis value for linear molecules.

************************
 Molecular Volume
************************

CHEMSMART provides several volume calculation methods, ranging from fast approximations to grid-based numerical
integration. All volumes are in cubic Ångströms (Å³).

.. list-table::
   :header-rows: 1
   :widths: 35 65

   -  -  Property
      -  Description

   -  -  :attr:`~Molecule.crude_volume_by_atomic_radii`
      -  Crude volume using covalent radii (no overlap correction)

   -  -  :attr:`~Molecule.crude_volume_by_vdw_radii`
      -  Crude volume using van der Waals radii (no overlap correction)

   -  -  :attr:`~Molecule.vdw_volume`
      -  VDW volume with pairwise overlap correction (fastest)

   -  -  :attr:`~Molecule.grid_vdw_volume`
      -  Grid-based VDW volume (recommended for complex molecules)

   -  -  :attr:`~Molecule.vdw_volume_from_rdkit`
      -  RDKit's native grid-based volume (requires RDKit)

   -  -  :attr:`~Molecule.voronoi_dirichlet_occupied_volume`
      -  Voronoi-Dirichlet tessellation volume

   -  -  :attr:`~Molecule.voronoi_dirichlet_polyhedra_occupied_volume`
      -  Voronoi-Dirichlet Polyhedra (VDP) volume

.. code:: python

   mol = Molecule.from_filepath(
       "tests/data/XTBTests/outputs/co2_ohess/co2_ohess.out"
   )

   # Crude estimates (no overlap correction)
   print(mol.crude_volume_by_atomic_radii)   # via covalent radii (Å³)
   print(mol.crude_volume_by_vdw_radii)      # via van der Waals radii (Å³)

   # Better estimates
   print(mol.vdw_volume)                     # pairwise overlap-corrected (Å³)
   print(mol.grid_vdw_volume)                # grid-based VDW, recommended (Å³)
   print(mol.vdw_volume_from_rdkit)          # RDKit grid-based volume (Å³)

   # Advanced
   print(mol.voronoi_dirichlet_occupied_volume)   # V-D tessellation (Å³)
   print(mol.voronoi_dirichlet_polyhedra_occupied_volume)  # VDP (Å³)
   print(mol.estimated_dispersion)                 # dispersion energy estimate

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.crude_volume_by_atomic_radii

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.crude_volume_by_vdw_radii

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.vdw_volume

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.grid_vdw_volume

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.vdw_volume_from_rdkit

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.voronoi_dirichlet_occupied_volume

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.voronoi_dirichlet_polyhedra_occupied_volume

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.estimated_dispersion

************************
 Bonds and Bond Orders
************************

Bond orders are inferred from interatomic distances using covalent radii and a tolerance buffer.

.. code:: python

   mol = Molecule.from_filepath(
       "tests/data/XTBTests/outputs/co2_ohess/co2_ohess.out"
   )
   print(mol.chemical_formula)   # CO2

   # --- All bonds (heuristic from 3D geometry) ---
   print("Bond orders:", mol.bond_orders)          # e.g. [2.0, 2.0] (C=O double bonds)

   # --- Bond orders extracted from RDKit connectivity ---
   bo_graph = mol.get_bond_orders_from_graph()
   print("Bond orders (from graph):", bo_graph)    # same structure as bond_orders

   # --- RDKit-based bond orders (requires RDKit) ---
   bo_rdkit = mol.get_bond_orders_from_rdkit_mol()
   if bo_rdkit is not None:
       print("Bond orders (RDKit):", bo_rdkit)

   # --- Query a single bond pair (1-based atom indices) ---
   # Between atom 1 (C) and atom 2 (O)
   single_bo = mol.determine_bond_order_one_bond(1, 2)
   print("Bond order 1–2:", single_bo)             # 2.0 for C=O

   # --- Full bond-order list (lower-level) ---
   all_bos = mol.determine_bond_order()
   print("All bond orders list:", all_bos)         # all detected bonds

   # --- RDKit Morgan fingerprint (for ML / similarity tasks) ---
   fp = mol.rdkit_fingerprints
   if fp is not None:
       print("Fingerprint shape:", fp.shape)        # e.g. (1, 2048)

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.bond_orders

.. automethod:: chemsmart.io.molecules.structure.Molecule.get_bond_orders_from_graph

.. automethod:: chemsmart.io.molecules.structure.Molecule.get_bond_orders_from_rdkit_mol

.. automethod:: chemsmart.io.molecules.structure.Molecule.determine_bond_order_one_bond

.. automethod:: chemsmart.io.molecules.structure.Molecule.determine_bond_order

.. autoattribute:: chemsmart.io.molecules.structure.Molecule.rdkit_fingerprints

************************
 Vibrational Data
************************

Vibrational properties are typically populated when reading frequency calculation output files.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   -  -  Property
      -  Description

   -  -  :attr:`~Molecule.vibrational_frequencies`
      -  Vibrational frequencies in cm⁻¹

   -  -  :attr:`~Molecule.vibrational_reduced_masses`
      -  Reduced masses per mode in amu

   -  -  :attr:`~Molecule.vibrational_force_constants`
      -  Force constants per mode in mDyne/Å

   -  -  :attr:`~Molecule.vibrational_ir_intensities`
      -  IR intensities in km/mol

   -  -  :attr:`~Molecule.vibrational_mode_symmetries`
      -  Irreducible representation labels (e.g. ``"A1"``)

   -  -  :attr:`~Molecule.vibrational_modes`
      -  Mass-weighted normal-mode displacements, each ``(N, 3)``

.. code:: python

   if mol.has_vibrations:
       print(mol.num_vib_frequencies)        # e.g. 3
       print(mol.vibrational_frequencies)     # e.g. [1595, 3657, 3756]

       # Generate a geometry displaced along mode 1
       displaced = mol.vibrationally_displaced(mode_idx=1, amp=0.1)

       # Generate a 20-frame trajectory for animation
       frames = mol.vibrationally_displaced(mode_idx=1, amp=0.1, nframes=20)

**Attributes:**

-  ``vibrational_frequencies`` (list) — Vibrational frequencies in cm⁻¹.
-  ``vibrational_reduced_masses`` (list) — Reduced masses per mode in amu.
-  ``vibrational_force_constants`` (list) — Force constants per mode in mDyne/Å.
-  ``vibrational_ir_intensities`` (list) — IR intensities in km/mol.
-  ``vibrational_mode_symmetries`` (list) — Irreducible representation labels (e.g. ``"A1"``).
-  ``vibrational_modes`` (list) — Each entry is an ``(N, 3)`` mass-weighted normal-mode displacement.
-  ``has_vibrations`` (bool) — ``True`` if any vibrational frequencies are present.
-  ``num_vib_frequencies`` (int) — Number of vibrational frequencies available.
-  ``num_vib_modes`` (int) — Number of normal modes available.

.. automethod:: chemsmart.io.molecules.structure.Molecule.vibrationally_displaced

************************
 Molecule Subclasses
************************

PKaMolecule
===========

A :class:`Molecule` subclass for pKa calculations that attaches a resolved ``proton_index`` (1-based) identifying the
acidic proton to be removed during deprotonation.

.. code:: python

   from chemsmart.io.molecules.structure import Molecule, PKaMolecule

   mol = Molecule(symbols=["O", "H", "H"],
                  positions=np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]]))
   pka_mol = PKaMolecule(molecule=mol, proton_index=2)
   print(pka_mol.proton_index)     # 2
   print(pka_mol.chemical_formula)  # H2O

.. autoclass:: chemsmart.io.molecules.structure.PKaMolecule
   :members:
   :inherited-members:
   :show-inheritance:

QMMMMolecule
============

A :class:`Molecule` subclass for QM/MM (ONIOM) calculations that stores partition level information (high / medium /
low level atoms), link bonds, and scale factors.

.. code:: python

   from chemsmart.io.molecules.structure import Molecule, QMMMMolecule

   # Load a molecule with at least 8 atoms
   mol = Molecule.from_filepath(
       "tests/data/GaussianTests/outputs/collidine_opt.log"
   )
   print(mol.chemical_formula, mol.num_atoms)   # C8H11N 20

   # When low_level_atoms is not specified, all atoms not in high/medium
   # are automatically assigned to low level
   qmmm = QMMMMolecule(
       molecule=mol,
       high_level_atoms=[1, 2, 3],
       bonded_atoms=[(3, 4)],
   )
   print(qmmm.partition_level_strings[:8])
   # ['H', 'H', 'H', 'L', 'L', 'L', 'L', 'L']

.. autoclass:: chemsmart.io.molecules.structure.QMMMMolecule
   :members:
   :inherited-members:
   :show-inheritance:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   -  -  Attribute
      -  Description

   -  -  ``high_level_atoms``
      -  Atom indices (1-based) treated at the high level of theory

   -  -  ``medium_level_atoms``
      -  Atom indices treated at the medium level (ONIOM 3-layer)

   -  -  ``low_level_atoms``
      -  Atom indices treated at the low level (auto-assigned if only high is given)

   -  -  ``bonded_atoms``
      -  List of ``(atom1, atom2)`` tuples describing QM/MM link bonds

   -  -  ``scale_factors``
      -  Dict mapping ``(atom1, atom2)`` → ``[low, medium, high]`` scale factors

   -  -  ``real_charge``, ``real_multiplicity``
      -  Charge / multiplicity of the full (real) QM/MM system

CoordinateBlock
================

A helper class that parses a Gaussian-style coordinate block text into a :class:`Molecule`:

.. code:: python

   from chemsmart.io.molecules.structure import CoordinateBlock, Molecule

   # Example: a Gaussian-style Z-matrix / Cartesian block with constraints
   block_text = """
   C  0  0  0.000000  0.000000  0.000000
   O  0  0  0.000000  0.000000  1.160000
   O -1  0  0.000000  0.000000 -1.160000
   """

   cb = CoordinateBlock(block_text)

   # Parsed chemical symbols
   print("Symbols:", cb.chemical_symbols)                 # ['C', 'O', 'O']

   # Per-atom constraints (-1 = frozen, 0 = relaxed)
   print("Constraints:", cb.constrained_atoms)           # [0, 0, -1]

   # Cartesian positions
   print("Positions (Å):", cb.positions)

   # Atom partition levels (for QM/MM ONIOM-like partitions)
   print("Partitions:", cb.partitions)                   # [0, 0, 0]

   # Directly convert to a Molecule
   mol = cb.molecule
   mol_or = cb.convert_coordinate_block_list_to_molecule(block_text.splitlines())
   print("Formula via Molecule:", mol.chemical_formula)  # CO2

.. autoclass:: chemsmart.io.molecules.structure.CoordinateBlock
   :members:
   :inherited-members:
   :show-inheritance:

AtomsChargeMultiplicity
=======================

A thin wrapper around ASE :class:`ase.Atoms` that additionally stores ``charge``, ``multiplicity``, ``frozen_atoms``,
``energy``, and ``forces``. Used for round-tripping between CHEMSMART :class:`Molecule` and ASE:

.. code:: python

   from chemsmart.io.molecules.structure import Molecule
   from chemsmart.io.molecules.atoms import AtomsChargeMultiplicity

   # 1. Build a CHEMSMART Molecule
   mol = Molecule(symbols=["C", "O", "O"],
                  positions=[[0, 0, 0], [0, 0, 1.16], [0, 0, -1.16]],
                  charge=0,
                  multiplicity=1)

   # 2. Molecule → ASE Atoms (wrapped as AtomsChargeMultiplicity, preserves charge/mult)
   atoms = mol.to_ase()
   print(type(atoms).__name__)                # AtomsChargeMultiplicity
   print("Atoms charge / mult:", atoms.charge, atoms.multiplicity)   # 0 1

   # 3. ASE Atoms → Molecule (round-trip back)
   mol2 = Molecule.from_ase_atoms(atoms)
   print("Round-trip formula:", mol2.chemical_formula)               # CO2
   print(mol2.charge, mol2.multiplicity)                              # 0 1

   # 4. Direct AtomsChargeMultiplicity construction
   import numpy as np
   acm = AtomsChargeMultiplicity(
       symbols=["H", "H", "O"],
       positions=np.array([[0.96, 0.0, 0.0], [-0.96, 0.0, 0.0], [0.0, 0.0, 0.0]]),
       charge=0,
       multiplicity=1,
   )
   print("ACM symbols:", list(acm.symbols))     # ['H', 'H', 'O']

.. autoclass:: chemsmart.io.molecules.atoms.AtomsChargeMultiplicity
   :members:
   :inherited-members:
   :show-inheritance:

************************
 Helper Functions
************************

Bond Cutoffs and Covalent Radii
=================================

Helper functions for computing covalent radii and bond cutoff distances:

.. code:: python

   from chemsmart.io.molecules import get_covalent_radius, get_bond_cutoff

   # Get covalent radius of an element (in Ångströms)
   r_C = get_covalent_radius("C")   # ~0.77 Å
   r_O = get_covalent_radius("O")   # ~0.73 Å

   # Compute a bond cutoff for a pair of elements (includes a buffer)
   cutoff_CO = get_bond_cutoff("C", "O")   # r_C + r_O + 0.3 (default buffer)
   cutoff_CH = get_bond_cutoff("C", "H", buffer=0.2)   # custom buffer

.. autofunction:: chemsmart.io.molecules.get_covalent_radius

.. autofunction:: chemsmart.io.molecules.get_bond_cutoff

PubChem Search
==============

Programmatic search of the PubChem database:

.. code:: python

   from chemsmart.io.molecules.pubchem import pubchem_search

   # Search by compound name
   result = pubchem_search("aspirin", fail_silently=True)

   # Search by CID
   result = pubchem_search(cid="2244", fail_silently=True)

   # Search by SMILES
   result = pubchem_search(smiles="CC(=O)OC1=CC=CC=C1C(=O)O", fail_silently=True)

   # Low-level raw search
   from chemsmart.io.molecules.pubchem import search_pubchem_raw
   raw = search_pubchem_raw("aspirin", "name", suffix="3d")

.. autofunction:: chemsmart.io.molecules.pubchem.pubchem_search

.. autofunction:: chemsmart.io.molecules.pubchem.search_pubchem_raw

**********
 See Also
**********

-  :doc:`molecule-input-formats` — Overview of supported input formats and CLI usage
-  :doc:`cli-overview` — Command-line interface reference
-  :doc:`thermochemistry-analysis` — Thermochemistry analysis using the Python API
-  :doc:`modules` — Full auto-generated module reference

For more technical details on the implementation, see the CHEMSMART preprint: https://arxiv.org/abs/2508.20042
