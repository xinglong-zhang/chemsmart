##################
 SOAP Descriptors
##################

CHEMSMART can compute Smooth Overlap of Atomic Positions (SOAP) descriptors from a
:class:`~chemsmart.io.molecules.structure.Molecule` geometry. SOAP is a local atomic environment fingerprint based on
nuclear coordinates and elemental species. No bonding / connectivity information is required.

SOAP is implemented in-tree with NumPy and SciPy (no optional dependencies). The GTO formulation matches DScribe 2.1.2
with ``rbf="gto"``, ``average="off"``, and default ``compression="off"``. Matching that release, neighbor atoms within
``r_cut + sigma * sqrt(-2 * ln(1e-3))`` (about ``r_cut + 3.72 * sigma``) contribute through Gaussian density tails;
there is no hard exclusion at exactly ``r_cut``.

*************
 Basic usage
*************

.. code:: python

   from chemsmart.io.molecules.structure import Molecule

   mol = Molecule.from_filepath("water.xyz")

   # Local (per-atom) SOAP vectors: shape (n_atoms, n_features)
   local = mol.calculate_soap(r_cut=6.0, n_max=8, l_max=6, sigma=1.0)

   # Global fingerprint by post-hoc outer averaging / summing of local
   # power spectra (local per-center vectors, then mean/sum over centers).
   # Prefer "mean" for size-comparable fingerprints; "sum" is extensive.
   global_mean = mol.calculate_soap(aggregation="mean")
   global_sum = mol.calculate_soap(aggregation="sum")

   # SOAP centered on selected atoms (1-based indices; duplicates kept)
   metal_centered = mol.calculate_soap(centers=[1])

*****************
 Hyperparameters
*****************

+------------------+----------------------------------------------------------+
| Parameter        | Meaning                                                  |
+==================+==========================================================+
| ``r_cut``        | Cutoff radius of the local atomic environment (Å). Must  |
|                  | be greater than 1 Å for the default GTO basis. Typical   |
|                  | values are about 3–8 Å. Matching DScribe 2.1.2, atoms    |
|                  | out to ``r_cut + ~3.72·sigma`` also contribute.          |
+------------------+----------------------------------------------------------+
| ``n_max``        | Number of radial basis functions.                        |
+------------------+----------------------------------------------------------+
| ``l_max``        | Maximum angular momentum (spherical harmonics degree).   |
|                  | Must be between 0 and 20 inclusive.                      |
+------------------+----------------------------------------------------------+
| ``sigma``        | Gaussian width / standard deviation in Å.                |
+------------------+----------------------------------------------------------+
| ``species``      | Elemental basis of the SOAP feature space.               |
+------------------+----------------------------------------------------------+
| ``centers``      | Optional 1-based atom indices used as SOAP centers.      |
|                  | Order is preserved; duplicate indices are kept and       |
|                  | overweight ``"mean"`` / ``"sum"`` aggregations.          |
+------------------+----------------------------------------------------------+
| ``aggregation``  | ``None`` (local), ``"mean"``, or ``"sum"``. Mean/sum are |
|                  | post-hoc over local power spectra (outer-average /       |
|                  | extensive-sum equivalent).                               |
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

Only standard chemical element symbols are accepted (for example ``"H"``, not deuterium ``"D"``). Symbols are
case-normalized (``"fe"`` → ``"Fe"``). Internally, species channels are ordered by atomic number (matching DScribe's
convention).

***********************************
 Supported systems and limitations
***********************************

SOAP from geometry works well for molecular conformers, transition states, organometallic complexes, and trajectory
frames that are treated as finite molecules.

This CHEMSMART interface currently supports **finite (non-periodic)** molecules only. Structures with active periodic
boundary conditions or non-empty translation vectors (cell) are rejected.

Aggregation modes ``"mean"`` and ``"sum"`` are computed **post-hoc** over local power-spectrum vectors from each
selected center. That is equivalent to outer averaging / summing over the selected centers. They are **not** inner
averaging (average expansion coefficients before forming the power spectrum). Inner averaging is not exposed by this
API.

Duplicate ``centers`` entries are allowed and preserved; they overweight ``"mean"`` and ``"sum"`` aggregations.

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
