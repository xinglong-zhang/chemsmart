#####################
 API Reference Index
#####################

This page provides quick navigation targets for all ``Molecule`` methods, attributes, helper classes and helper
functions mentioned in the API guide. Each entry is a manually declared Sphinx signature — no source-code docstrings are
pulled in (so there are no ``>>>`` doctest blocks here).

****************
 Molecule Class
****************

.. py:class:: Molecule
   :module: chemsmart.io.molecules.structure

*************
 Constructor
*************

.. py:method:: Molecule.__init__(symbols, positions, charge=None, multiplicity=None, frozen_atoms=None, pbc_conditions=None, translation_vectors=None, energy=None, forces=None, velocities=None, vibrational_frequencies=None, vibrational_modes=None, dipole_moment=None, point_group=None, info=None)
   :module: chemsmart.io.molecules.structure

Constructor Parameters
======================

**symbols** (*list*) — List of atomic symbols (e.g. ``["C", "O", "O"]``).

**positions** (*ndarray*) — ``(N, 3)`` array of Cartesian coordinates in Ångströms.

``charge`` (*int*) — Molecular charge

``multiplicity`` (*int*) — Spin multiplicity (2S + 1)

``frozen_atoms`` (*list*) — Per-atom constraint flags: ``-1`` = frozen, ``0`` = relaxed (Gaussian convention)

``pbc_conditions`` (*list*) — Periodic boundary conditions ``[x, y, z]``

``translation_vectors`` (*list*) — Cell / translation vectors for periodic systems

``energy`` (*float*) — Total energy in Hartree

``forces`` (*ndarray*) — ``(N, 3)`` array of forces in Hartree / Bohr

``velocities`` (*ndarray*) — ``(N, 3)`` array of atomic velocities

``vibrational_frequencies`` (*list*) — Vibrational frequencies in cm⁻¹

``vibrational_modes`` (*list*) — Each entry is an ``(N, 3)`` mass-weighted normal-mode displacement

``dipole_moment`` (*ndarray*) — Dipole moment vector ``[X, Y, Z]`` in Debye

``point_group`` (*str*) — Point group label (e.g. ``"C2V"``, ``"CS"``)

``info`` (*dict*) — Arbitrary extra metadata stored with the molecule

Constructor Examples
--------------------

.. code:: python

   # From symbols and positions
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

*****************
 Reading Methods
*****************

.. py:method:: Molecule.from_filepath(filepath, charge=None, multiplicity=None, **kwargs)
   :module: chemsmart.io.molecules.structure

.. py:method:: Molecule.from_directorypath(directorypath, index=-1, return_list=False, program=None)
   :module: chemsmart.io.molecules.structure

.. py:method:: Molecule.from_coordinate_block_text(coordinate_block_text, charge=None, multiplicity=None, **kwargs)
   :module: chemsmart.io.molecules.structure

.. py:method:: Molecule.from_ase_atoms(atoms, charge=None, multiplicity=None)
   :module: chemsmart.io.molecules.structure

.. py:method:: Molecule.from_rdkit_mol(rdkit_mol, charge=None, multiplicity=None)
   :module: chemsmart.io.molecules.structure

.. py:method:: Molecule.from_pubchem(query, search_type="name")
   :module: chemsmart.io.molecules.structure

**********************
 Manipulation Methods
**********************

.. py:method:: Molecule.copy()
   :module: chemsmart.io.molecules.structure

.. py:method:: Molecule.delete_atoms_by_indices(atom_indices, one_based=True)
   :module: chemsmart.io.molecules.structure

*****************
 Writing Methods
*****************

.. py:method:: Molecule.write(filename, format="xyz", mode="w", **kwargs)
   :module: chemsmart.io.molecules.structure

.. py:method:: Molecule.write_cosmorsxyz(filename)
   :module: chemsmart.io.molecules.structure

********************
 Conversion Methods
********************

.. py:method:: Molecule.to_ase()
   :module: chemsmart.io.molecules.structure

.. py:method:: Molecule.to_pymatgen()
   :module: chemsmart.io.molecules.structure

.. py:method:: Molecule.to_rdkit(add_bonds=True, bond_cutoff_buffer=0.05, adjust_H=True)
   :module: chemsmart.io.molecules.structure

.. py:method:: Molecule.to_smiles()
   :module: chemsmart.io.molecules.structure

.. py:method:: Molecule.to_pdb()
   :module: chemsmart.io.molecules.structure

.. py:method:: Molecule.to_cosmorsxyz()
   :module: chemsmart.io.molecules.structure

.. py:method:: Molecule.to_graph(bond_cutoff_buffer=0.3)
   :module: chemsmart.io.molecules.structure

.. py:method:: Molecule.to_X_data(wbo=False)
   :module: chemsmart.io.molecules.structure

********************
 Identifier Methods
********************

.. py:method:: Molecule.get_chemical_formula(hill=True)
   :module: chemsmart.io.molecules.structure

******************
 Geometry Methods
******************

.. py:method:: Molecule.get_distance(idx1, idx2)
   :module: chemsmart.io.molecules.structure

.. py:method:: Molecule.get_angle(idx1, idx2, idx3)
   :module: chemsmart.io.molecules.structure

.. py:method:: Molecule.get_dihedral(idx1, idx2, idx3, idx4)
   :module: chemsmart.io.molecules.structure

.. py:method:: Molecule.get_all_distances()
   :module: chemsmart.io.molecules.structure

.. py:method:: Molecule.bond_lengths()
   :module: chemsmart.io.molecules.structure

.. py:method:: Molecule.get_angle_from_positions(position1, position2, position3)
   :module: chemsmart.io.molecules.structure

.. py:method:: Molecule.get_dihedral_from_positions(position1, position2, position3, position4)
   :module: chemsmart.io.molecules.structure

***********************
 Bond Analysis Methods
***********************

.. py:method:: Molecule.get_bond_orders_from_graph(bond_cutoff_buffer=0.3)
   :module: chemsmart.io.molecules.structure

.. py:method:: Molecule.get_bond_orders_from_rdkit_mol(bond_cutoff_buffer=0.3, **kwargs)
   :module: chemsmart.io.molecules.structure

.. py:method:: Molecule.determine_bond_order_one_bond(bond_length, bond_cutoff)
   :module: chemsmart.io.molecules.structure

.. py:method:: Molecule.determine_bond_order(bond_lengths, bond_cutoff)
   :module: chemsmart.io.molecules.structure

*******************
 Vibration Methods
*******************

.. py:method:: Molecule.vibrationally_displaced(mode_idx, amp=0.1, nframes=None)
   :module: chemsmart.io.molecules.structure

************
 Attributes
************

**Identifiers**

.. py:attribute:: Molecule.chemical_formula
   :module: chemsmart.io.molecules.structure

.. py:attribute:: Molecule.smiles
   :module: chemsmart.io.molecules.structure

.. py:attribute:: Molecule.inchikey
   :module: chemsmart.io.molecules.structure

.. py:attribute:: Molecule.structure_id
   :module: chemsmart.io.molecules.structure

.. py:attribute:: Molecule.structure_label
   :module: chemsmart.io.molecules.structure

**Mass Properties**

.. py:attribute:: Molecule.mass
   :module: chemsmart.io.molecules.structure

.. py:attribute:: Molecule.natural_abundance_weighted_mass
   :module: chemsmart.io.molecules.structure

.. py:attribute:: Molecule.most_abundant_mass
   :module: chemsmart.io.molecules.structure

**Structure Classification**

.. py:attribute:: Molecule.is_aromatic
   :module: chemsmart.io.molecules.structure

.. py:attribute:: Molecule.is_ring
   :module: chemsmart.io.molecules.structure

.. py:attribute:: Molecule.is_linear
   :module: chemsmart.io.molecules.structure

.. py:attribute:: Molecule.is_chiral
   :module: chemsmart.io.molecules.structure

**Spatial Properties**

.. py:attribute:: Molecule.center_of_mass
   :module: chemsmart.io.molecules.structure

.. py:attribute:: Molecule.canonical_positions
   :module: chemsmart.io.molecules.structure

.. py:attribute:: Molecule.moments_of_inertia
   :module: chemsmart.io.molecules.structure

.. py:attribute:: Molecule.rotational_temperatures
   :module: chemsmart.io.molecules.structure

.. py:attribute:: Molecule.distance_matrix
   :module: chemsmart.io.molecules.structure

**Volume Properties**

.. py:attribute:: Molecule.crude_volume_by_atomic_radii
   :module: chemsmart.io.molecules.structure

.. py:attribute:: Molecule.vdw_volume
   :module: chemsmart.io.molecules.structure

.. py:attribute:: Molecule.estimated_dispersion
   :module: chemsmart.io.molecules.structure

**Electronic Properties**

.. py:attribute:: Molecule.rdkit_fingerprints
   :module: chemsmart.io.molecules.structure

**Bond Properties**

.. py:attribute:: Molecule.bond_orders
   :module: chemsmart.io.molecules.structure

****************
 Helper Classes
****************

.. py:class:: CoordinateBlock
   :module: chemsmart.io.molecules.structure

Attributes
==========

.. py:attribute:: CoordinateBlock.atom_labels
   :module: chemsmart.io.molecules.structure

.. py:attribute:: CoordinateBlock.coordinates
   :module: chemsmart.io.molecules.structure

Constructor
===========

.. py:method:: CoordinateBlock.__init__(atom_labels, coordinates)
   :module: chemsmart.io.molecules.structure

Methods
=======

.. py:method:: CoordinateBlock.convert_to_molecule()
   :module: chemsmart.io.molecules.structure

.. py:method:: CoordinateBlock.convert_coordinate_block_list_to_molecule()
   :module: chemsmart.io.molecules.structure

.. py:class:: AtomsChargeMultiplicity
   :module: chemsmart.io.molecules.atoms

Constructor
===========

.. py:method:: AtomsChargeMultiplicity.__init__(symbols, positions, charge, multiplicity, frozen_atoms, energy=None, forces=None, velocities=None, vibrational_frequencies=None, vibrational_modes=None, dipole_moment=None, point_group=None, pbc_conditions=None, translation_vectors=None, info=None)
   :module: chemsmart.io.molecules.atoms

******************
 Helper Functions
******************

Bond Cutoffs and Covalent Radii
===============================

.. py:function:: get_covalent_radius(element_symbol)
   :module: chemsmart.io.molecules

.. py:function:: get_bond_cutoff(el1, el2, buffer=0.3)
   :module: chemsmart.io.molecules

PubChem Search
==============

.. py:function:: pubchem_search(name=None, cid=None, smiles=None, fail_silently=False, timeout=10)
   :module: chemsmart.io.molecules.pubchem

.. py:function:: search_pubchem_raw(query, namespace, suffix="2d", timeout=10)
   :module: chemsmart.io.molecules.pubchem
