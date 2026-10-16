######################
 Molecular Properties
######################

***********************
 Molecular Identifiers
***********************

CHEMSMART provides several identifiers for molecules and structures. The four most commonly used ones are:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   -  -  Property
      -  Description
   -  -  :attr:`~Molecule.smiles`
      -  SMILES string via RDKit
   -  -  :attr:`~Molecule.cxsmiles`
      -  CXSMILES string (with 3D coordinates and stereo)
   -  -  :attr:`~Molecule.inchi`
      -  Full InChI string (requires Open Babel)
   -  -  :attr:`~Molecule.inchikey`
      -  InChIKey string (requires Open Babel)

All other identifiers — :attr:`~Molecule.chemical_formula`, :attr:`~Molecule.empirical_formula`,
:attr:`~Molecule.molecule_id`, :attr:`~Molecule.structure_id`, etc. — are also available as attributes.

.. code:: python

   mol = Molecule(
       symbols=["O", "H", "H"],
       positions=np.array([
           [0.0,  0.0,  0.119],
           [0.0,  0.76, -0.477],
           [0.0, -0.76, -0.477],
       ]),
   )

   print(mol.smiles)       # O
   print(mol.cxsmiles)     # O
   print(mol.inchi)        # InChI=1S/H2O/h1H2O
   print(mol.inchikey)     # UHOVQNZJYSORNB-UHFFFAOYSA-N

*********************************
 Atomic and Molecular Properties
*********************************

Counts and Symbols
==================

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
======

CHEMSMART supports three mass conventions:

-  **Standard atomic masses** (:attr:`mass`, :attr:`masses`)
-  **Natural abundance weighted masses** (:attr:`natural_abundance_weighted_mass`,
   :attr:`natural_abundance_weighted_masses`)
-  **Most abundant isotope masses** (:attr:`most_abundant_mass`, :attr:`most_abundant_masses`)

.. code:: python

   # Use water molecule for mass examples
   mol = Molecule(
       symbols=["O", "H", "H"],
       positions=np.array([
           [0.0,  0.0,  0.119],
           [0.0,  0.76, -0.477],
           [0.0, -0.76, -0.477],
       ]),
   )

   print(mol.mass)                              # ~18.015 amu
   print(mol.natural_abundance_weighted_mass)  # weighted by isotope abundance
   print(mol.most_abundant_mass)                # using most abundant isotopes
   print(mol.masses)                            # list of per-atom masses
   print(mol.most_abundant_masses)              # per-atom most-abundant isotope masses

**Attributes:**

-  ``mass`` (float) — Total molecular mass in amu (standard atomic masses).
-  ``masses`` (list) — Per-atom standard atomic masses in amu.
-  ``natural_abundance_weighted_mass`` (float) — Total mass using natural isotope-abundance weighted atomic masses.
-  ``natural_abundance_weighted_masses`` (list) — Per-atom abundance-weighted masses.
-  ``most_abundant_mass`` (float) — Total mass using the most abundant isotope of each element.
-  ``most_abundant_masses`` (list) — Per-atom most-abundant isotope masses.

Energy, Forces, and Velocities
==============================

These are populated from calculation output files (Gaussian ``.log``, ORCA ``.out``, xTB, etc.).

.. code:: python

   # Use CO2 molecule from xTB calculation (has energy and forces)
   mol = Molecule.from_filepath(
       "tests/data/XTBTests/outputs/co2_ohess/co2_ohess.out"
   )

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
========================================

These properties are typically populated when reading calculation output files. The example below uses a real Gaussian
output file:

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
============

.. code:: python

   # Use CO2 molecule
   mol = Molecule.from_filepath("tests/data/GaussianTests/outputs/co2.log")

   # Covalent radii (Å)
   print(mol.atomic_radii_list)      # e.g. [0.77, 0.73, 0.73] for CO2
   # Van der Waals radii (Å)
   print(mol.vdw_radii_list)         # e.g. [1.70, 1.52, 1.52] for CO2

-  ``atomic_radii_list`` (list) — Covalent radii (Å) for each atom, from ASE's ``covalent_radii`` table.
-  ``vdw_radii_list`` (list) — Van der Waals radii (Å) for each atom, from ASE's ``vdW_radii`` table.

**************************
 Structure Classification
**************************

CHEMSMART provides several boolean properties for classifying molecular structure. These are based on RDKit's perception
of the 3D geometry (bonds, rings, stereochemistry).

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

   # Use CO2 molecule (linear, non-aromatic)
   mol = Molecule.from_filepath("tests/data/GaussianTests/outputs/co2.log")

   print(mol.is_aromatic)         # False
   print(mol.is_linear)           # True
   print(mol.is_chiral)           # False
   print(mol.chiral_centers)      # {}
   print(mol.num_components)      # 1

**Attributes:**

-  ``is_aromatic`` (bool) — ``True`` if any atom is aromatic and part of a ring (requires RDKit connectivity).
-  ``is_ring`` (bool) — ``True`` if the molecular graph contains at least one cycle.
-  ``is_monoatomic`` (bool) — ``True`` when ``num_atoms == 1``.
-  ``is_diatomic`` (bool) — ``True`` when ``num_atoms == 2``.
-  ``is_linear`` (bool) — ``True`` if all atoms are collinear (within tolerance).
-  ``is_chiral`` (bool) — ``True`` if any tetrahedral centre has no improper-rotation symmetry.
-  ``chiral_centers`` (dict) — Map of 1-based atom index → stereodescriptor (``"R"``/``"S"``) for detected chiral
   centres.
-  ``is_multicomponent`` (bool) — ``True`` if the connectivity graph has more than one connected component (salts,
   solvates, etc.).
-  ``num_components`` (int) — Number of connected components in the molecular graph.

.. note::

   Bond orders are inferred from 3D geometry, not from an electronic structure calculation, so aromaticity detection is
   **model-dependent**. For borderline or unusual systems, the result should be treated as a heuristic estimate.

Periodic Boundary Conditions
============================

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

-  ``pbc`` (bool) — ``True`` if the molecule has no periodic boundary conditions (all entries in ``pbc_conditions`` are
   0).
-  ``pbc_conditions`` (list) — Periodic boundary conditions ``[x, y, z]``.
-  ``translation_vectors`` (list) — Cell / translation vectors for periodic systems.

Frozen Atoms
============

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

-  ``frozen_atoms`` (list) — Per-atom constraint flags following the Gaussian convention: ``-1`` = frozen, ``0`` =
   relaxed.

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
