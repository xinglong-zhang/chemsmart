##################
 SOAP Descriptors
##################

CHEMSMART can compute Smooth Overlap of Atomic Positions (SOAP) descriptors from a
:class:`~chemsmart.io.molecules.structure.Molecule` geometry. SOAP is a local atomic environment fingerprint based on
nuclear coordinates and elemental species. No bonding / connectivity information is required.

**************
 Installation
**************

SOAP support depends on the optional `DScribe <https://singroup.github.io/dscribe/>`_ package and is **not** installed
with the default CHEMSMART dependencies:

.. code:: bash

   pip install 'chemsmart[soap]'

Without DScribe, importing ``Molecule`` still works; calling SOAP methods raises an actionable ``ImportError``.

*************
 Basic usage
*************

.. code:: python

   from chemsmart.io.molecules.structure import Molecule

   mol = Molecule.from_filepath("water.xyz")

   # Local (per-atom) SOAP vectors: shape (n_atoms, n_features)
   local = mol.calculate_soap(r_cut=6.0, n_max=8, l_max=6, sigma=1.0)

   # Global fingerprint by outer averaging / summing local power spectra
   # over the selected centers (DScribe average="outer" semantics).
   # Prefer "mean" for size-comparable fingerprints; "sum" is extensive.
   global_mean = mol.calculate_soap(aggregation="mean")
   global_sum = mol.calculate_soap(aggregation="sum")

   # SOAP centered on selected atoms (1-based indices)
   metal_centered = mol.calculate_soap(centers=[1])

*****************
 Hyperparameters
*****************

+------------------+----------------------------------------------------------+
| Parameter        | Meaning                                                  |
+==================+==========================================================+
| ``r_cut``        | Cutoff radius of the local atomic environment (Å). Must  |
|                  | be greater than 1 Å for the default GTO basis. Typical   |
|                  | values are about 3–8 Å.                                  |
+------------------+----------------------------------------------------------+
| ``n_max``        | Number of radial basis functions.                        |
+------------------+----------------------------------------------------------+
| ``l_max``        | Maximum angular momentum (spherical harmonics degree).   |
+------------------+----------------------------------------------------------+
| ``sigma``        | Gaussian width / standard deviation in Å.                |
+------------------+----------------------------------------------------------+
| ``species``      | Elemental basis of the SOAP feature space.               |
+------------------+----------------------------------------------------------+
| ``centers``      | Optional 1-based atom indices used as SOAP centers.      |
+------------------+----------------------------------------------------------+
| ``aggregation``  | ``None`` (local), ``"mean"``, or ``"sum"`` (outer        |
|                  | average / sum of local power spectra).                   |
+------------------+----------------------------------------------------------+

****************************
 Species basis for datasets
****************************

When ``species`` is omitted, CHEMSMART infers a sorted unique list of elements present in the **current** molecule. That
is convenient for a single structure, but feature length and meaning then depend on the molecule's composition.

For comparable descriptors across a dataset (for example, machine-learning models), always pass an explicit shared
species list that covers every element present in the dataset:

.. code:: python

   species = ["H", "C", "N", "O", "Fe"]
   features = mol.calculate_soap(species=species, aggregation="mean")

Only standard chemical element symbols recognized by ASE are accepted (for example ``"H"``, not deuterium ``"D"``).

***********************************
 Supported systems and limitations
***********************************

SOAP from geometry works well for molecular conformers, transition states, organometallic complexes, and trajectory
frames that are treated as finite molecules.

This CHEMSMART interface currently supports **finite (non-periodic)** molecules only. Structures with active periodic
boundary conditions are rejected.

Aggregation modes ``"mean"`` and ``"sum"`` implement **outer** averaging of local power-spectrum vectors over the
selected centers. They are **not** DScribe's ``average="inner"`` global fingerprint (average expansion coefficients
before forming the power spectrum). Inner averaging is not exposed by this API.

Standard geometric SOAP does **not** distinguish structures that share the same nuclear coordinates but differ
electronically, for example different:

-  spin states (high-spin vs low-spin Fe)
-  oxidation states with nearly identical geometries
-  charge states
-  electronic configurations

Conventional power-spectrum SOAP is also reflection-invariant and does not encode molecular chirality / enantiomer
identity.

For catalytic chemistry, SOAP is therefore often augmented with electronic or topological descriptors such as formal
oxidation state, partial charges, spin density, frontier-orbital energies, NBO descriptors, or local coordination
numbers. Those augmentations are outside the scope of the current geometry-only API.

***************
 API reference
***************

The low-level entry point is :func:`chemsmart.analysis.soap.calculate_soap`. The convenience method
:meth:`~chemsmart.io.molecules.structure.Molecule.calculate_soap` delegates to it.
