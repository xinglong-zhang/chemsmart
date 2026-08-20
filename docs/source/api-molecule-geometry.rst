###################
 Geometry Analysis
###################

*******************
 Geometry Analysis
*******************

All geometry methods use **1-based** atom indices, matching Gaussian / ORCA conventions.

Distances, Angles, and Dihedrals
================================

.. code:: python

   # Use H2O molecule for geometry examples (3 atoms)
   from chemsmart.io.molecules.structure import Molecule
   import numpy as np

   mol = Molecule(
       symbols=["O", "H", "H"],
       positions=np.array([
           [0.0,  0.0,  0.119],
           [0.0,  0.76, -0.477],
           [0.0, -0.76, -0.477],
       ]),
   )

   # Bond distance between atoms 1 and 2 (H2O: O–H)
   dist = mol.get_distance(1, 2)
   print(f"Distance 1-2: {dist:.3f} Å")

   # Bond angle H2–O1–H3 (in degrees)
   angle = mol.get_angle(2, 1, 3)
   print(f"Angle 2-1-3: {angle:.1f}°")

   # Dihedral angle requires at least 4 atoms. Create a 4-atom molecule:
   butane = Molecule(
       symbols=["C", "C", "C", "C"],
       positions=np.array([
           [-1.26, 0.0, 0.0],
           [0.0, 0.0, 0.0],
           [1.26, 0.35, 0.0],
           [2.52, 0.0, 0.0],
       ]),
   )
   dihedral = butane.get_dihedral(1, 2, 3, 4)
   print(f"Dihedral 1-2-3-4: {dihedral:.1f}°")

   # All pairwise distances (H2O: 3 atoms → 3 pairs)
   all_distances = mol.get_all_distances()
   print(f"All distances: {len(all_distances)} pairs")

   # Pairwise distance matrix (N × N)
   dist_matrix = mol.distance_matrix
   print(f"Distance matrix shape: {dist_matrix.shape}")

:meth:`~chemsmart.io.molecules.structure.Molecule.get_distance(idx1, idx2)` Calculate the distance between two atoms
(1-based indices).

-  ``idx1`` (int) — First atom index (1-based).
-  ``idx2`` (int) — Second atom index (1-based).
-  Returns: ``float`` — Distance in Ångströms.

:meth:`~chemsmart.io.molecules.structure.Molecule.get_angle(idx1, idx2, idx3)` Calculate the angle between three atoms
(1-based indices), measured at atom ``idx2``.

-  ``idx1``, ``idx2``, ``idx3`` (int) — Atom indices forming the angle (1-based).
-  Returns: ``float`` — Angle in degrees.

:meth:`~chemsmart.io.molecules.structure.Molecule.get_dihedral(idx1, idx2, idx3, idx4)` Calculate the dihedral (torsion)
angle between four atoms (1-based indices).

-  ``idx1``, ``idx2``, ``idx3``, ``idx4`` (int) — Atom indices defining the dihedral (1-based).
-  Returns: ``float`` — Dihedral angle in degrees.

:meth:`~chemsmart.io.molecules.structure.Molecule.get_all_distances()` Compute all pairwise interatomic distances.

-  Returns: ``list`` — List of distances in Ångströms.

:attr:`~chemsmart.io.molecules.structure.Molecule.distance_matrix` Pairwise distance matrix of shape ``(N, N)``.

-  Returns: ``numpy.ndarray`` — Distance matrix in Ångströms.

:meth:`~chemsmart.io.molecules.structure.Molecule.bond_lengths()` Compute bond lengths for all detected bonds.

-  Returns: ``list`` — List of bond lengths in Ångströms.

:meth:`~chemsmart.io.molecules.structure.Molecule.get_angle_from_positions(position1, position2, position3)` Calculate
an angle from three positions directly (no molecule needed).

:meth:`~chemsmart.io.molecules.structure.Molecule.get_dihedral_from_positions(position1, position2, position3,
position4)` Calculate a dihedral angle from four positions directly (no molecule needed).

Center of Mass
==============

-  ``center_of_mass`` (numpy.ndarray) — Shape ``(3,)`` vector of the molecular centre of mass in Å (standard atomic
   masses).

Canonical Geometry
==================

These properties yield representations invariant to translation, rotation, and atom-index permutation:

-  ``canonical_positions`` (numpy.ndarray) — Coordinates translated to centre-of-mass and rotated into principal axes
   frame.
-  ``canonical_geometry`` (tuple) — A hashable, sorted descriptor of the canonical geometry (used for ``structure_id`` /
   ``structure_label``).

.. note::

   Canonical positions are translated to the centre of mass and rotated into the principal-axes frame of the
   moment-of-inertia tensor. The :attr:`structure_id` is a SHA-256 hash of the canonical geometry, charge, and
   multiplicity.

Moments of Inertia and Rotational Temperatures
==============================================

Moments of inertia are in units of amu·Å². Rotational temperatures are in Kelvin.

-  ``moments_of_inertia_tensor`` (numpy.ndarray) — ``(3, 3)`` symmetric inertia tensor (standard masses).
-  ``moments_of_inertia`` (numpy.ndarray) — ``(3,)`` sorted principal moments (ascending, standard masses).
-  ``moments_of_inertia_weighted_mass`` (numpy.ndarray) — Principal moments using natural-abundance weighted masses.
-  ``moments_of_inertia_most_abundant_mass`` (numpy.ndarray) — Principal moments using most-abundant isotope masses.
-  ``moments_of_inertia_principal_axes`` (numpy.ndarray) — ``(3, 3)`` matrix with principal axes as rows.
-  ``rotational_temperatures`` (numpy.ndarray) — Rotational temperatures for each principal moment (K); linear molecules
   return only the finite perpendicular value.

.. tip::

   For linear molecules, the moment of inertia along the molecular axis is zero, so the axial rotational temperature is
   infinite. CHEMSMART returns only the finite perpendicular-axis value for linear molecules.

******************
 Molecular Volume
******************

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

**Volume Attributes:**

-  ``crude_volume_by_atomic_radii`` (float) — Sum of spheres volumes from covalent radii (no overlap correction).
-  ``crude_volume_by_vdw_radii`` (float) — Sum of spheres volumes from van der Waals radii (no overlap correction).
-  ``vdw_volume`` (float) — VDW volume with pairwise overlap subtraction (closed-form).
-  ``grid_vdw_volume`` (float) — VDW volume evaluated on a Cartesian grid (most reliable for complex shapes).
-  ``vdw_volume_from_rdkit`` (float or None) — RDKit grid-based volume (requires RDKit).
-  ``voronoi_dirichlet_occupied_volume`` (float) — Occupied volume from Voronoi-Dirichlet tessellation.
-  ``voronoi_dirichlet_polyhedra_occupied_volume`` (float) — VDP polyhedra-based occupied volume.
-  ``estimated_dispersion`` (float) — London dispersion energy estimate derived from atomic radii and volume.

***********************
 Bonds and Bond Orders
***********************

Bond orders are inferred from interatomic distances using covalent radii and a tolerance buffer.

.. code:: python

   import numpy as np
   from chemsmart.io.molecules import get_covalent_radius

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

   # --- Query a single bond (low-level: needs bond length + cutoff) ---
   # Get bond length between atoms 1 and 2, plus the C–O cutoff
   r_C = get_covalent_radius("C")
   r_O = get_covalent_radius("O")
   cutoff_CO = r_C + r_O + 0.3
   bond_len = mol.get_distance(1, 2)
   single_bo = mol.determine_bond_order_one_bond(bond_len, cutoff_CO)
   print("Bond order 1–2:", single_bo)             # 2.0 for C=O

   # --- Full bond-order list (lower-level, vectorized) ---
   all_bos = mol.determine_bond_order(np.array(mol.get_all_distances()), cutoff_CO)
   print("All bond orders list:", all_bos)         # all detected bonds

   # --- RDKit Morgan fingerprint (for ML / similarity tasks) ---
   fp = mol.rdkit_fingerprints
   if fp is not None:
       print("Fingerprint bits:", fp.GetNumBits())   # 2048 for RDK fingerprint
       print("Fingerprint as array:", list(fp.GetOnBits()))   # list of active bit indices

**Attributes:**

-  ``bond_orders`` (list) — Bond orders for all detected bonds.

:meth:`~chemsmart.io.molecules.structure.Molecule.get_bond_orders_from_graph(bond_cutoff_buffer=0.3)` Infer bond orders
from the molecular graph using covalent radii cutoffs.

-  ``bond_cutoff_buffer`` (float) — Extra tolerance for bond detection (default: 0.3 Å).
-  Returns: ``list`` — Bond order values.

:meth:`~chemsmart.io.molecules.structure.Molecule.get_bond_orders_from_rdkit_mol(bond_cutoff_buffer=0.3, **kwargs)` Get
bond orders from RDKit's molecular graph perception.

-  Returns: ``list`` — Bond order values (e.g. ``[1.0]``, ``[2.0]``).

:meth:`~chemsmart.io.molecules.structure.Molecule.determine_bond_order_one_bond(bond_length, bond_cutoff)` Determine the
bond order of a single bond given its length and cutoff.

-  ``bond_length`` (float) — Length of the bond in Ångströms.
-  ``bond_cutoff`` (float) — Cutoff distance for bond detection in Ångströms.
-  Returns: ``float`` — Bond order (0 = no bond, 1 = single, 1.5 = aromatic, 2 = double, 3 = triple).

:meth:`~chemsmart.io.molecules.structure.Molecule.determine_bond_order(bond_lengths, bond_cutoff)` Vectorized bond-order
determination for multiple bond lengths.

-  ``bond_lengths`` (array-like) — Bond lengths in Ångströms.
-  ``bond_cutoff`` (float) — Cutoff distance in Ångströms.
-  Returns: ``numpy.ndarray`` — Bond orders for each input length.

:attr:`~chemsmart.io.molecules.structure.Molecule.rdkit_fingerprints` Morgan fingerprint from RDKit for ML / similarity
tasks.

-  Returns: ``numpy.ndarray`` or ``None`` — Fingerprint vector (e.g. shape ``(1, 2048)``).

:meth:`~chemsmart.io.molecules.structure.Molecule.get_all_distances()` Compute all pairwise interatomic distances.

-  Returns: ``list`` — List of distances in Ångströms.
