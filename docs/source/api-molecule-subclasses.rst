Subclasses, Helpers & Reference
================================

Molecule Subclasses
-------------------

PKaMolecule
^^^^^^^^^^^

A :class:`Molecule` subclass for pKa calculations that attaches a resolved ``proton_index`` (1-based) identifying the
acidic proton to be removed during deprotonation.

.. code:: python

   import numpy as np
   from chemsmart.io.molecules.structure import Molecule, PKaMolecule

   mol = Molecule(symbols=["O", "H", "H"],
                  positions=np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]]))
   pka_mol = PKaMolecule(molecule=mol, proton_index=2)
   print(pka_mol.proton_index)     # 2
   print(pka_mol.chemical_formula)  # H2O

**Constructor:**

- ``molecule`` (:class:`Molecule`) — The parent molecule whose data is inherited.
- ``proton_index`` (int) — 1-based index of the acidic proton in the molecule.

**Key attributes:**

- ``proton_index`` (int) — The 1-based proton index.
- All other :class:`Molecule` attributes are inherited.

QMMMMolecule
^^^^^^^^^^^^

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

**Constructor:**

- ``molecule`` (:class:`Molecule`) — The parent molecule.
- ``high_level_atoms`` (list) — 1-based atom indices in the high-level QM region.
- ``medium_level_atoms`` (list) — 1-based atom indices in the medium-level region.
- ``low_level_atoms`` (list) — 1-based atom indices in the low-level (MM) region.
- ``bonded_atoms`` (list of tuples) — Link bonds between QM and MM regions.
- ``scale_factors`` (dict) — Scaling factors for electrostatic embedding.

**Key attributes:**

- ``partition_level_strings`` (list) — Per-atom partition labels (``"H"``, ``"M"``, ``"L"``).
- ``link_bonds`` (list) — Link bond tuples.
- All other :class:`Molecule` attributes are inherited.

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
^^^^^^^^^^^^^^^^

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
   mol_or = cb.convert_coordinate_block_list_to_molecule()
   print("Formula via Molecule:", mol.chemical_formula)  # CO2

**CoordinateBlock Attributes and Methods:**

- ``chemical_symbols`` (list) — Parsed element symbol strings.
- ``constrained_atoms`` (list) — Per-atom constraint flags (-1 = frozen, 0 = relaxed).
- ``positions`` (numpy.ndarray) — Parsed Cartesian coordinates (N, 3).
- ``partitions`` (list) — Per-atom ONIOM-style partition levels.
- ``freeze_flags`` (list) — Legacy alias for ``constrained_atoms``.
- ``molecule`` (Molecule) — Convert directly into a CHEMSMART Molecule (charge/multiplicity inferred from partition column).
- ``convert_coordinate_block_list_to_molecule()`` — Alternate constructor-path conversion into a Molecule.
- ``from_coordinate_block_block_string(block_string)`` — Classmethod: build a Molecule directly from a coordinate block string.

AtomsChargeMultiplicity
^^^^^^^^^^^^^^^^^^^^^^^^

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
       charge=0,
       multiplicity=1,
       frozen_atoms=None,
       energy=None,
       forces=None,
       symbols=["H", "H", "O"],
       positions=np.array([[0.96, 0.0, 0.0], [-0.96, 0.0, 0.0], [0.0, 0.0, 0.0]]),
   )
   print("ACM symbols:", list(acm.symbols))     # ['H', 'H', 'O']

**AtomsChargeMultiplicity Constructor:**

:class:`~chemsmart.io.molecules.atoms.AtomsChargeMultiplicity(charge, multiplicity, frozen_atoms, energy, forces, **kwargs)`
- ``charge`` (int) — Molecular charge.
- ``multiplicity`` (int) — Spin multiplicity (``2S + 1``).
- ``frozen_atoms`` (list or None) — 1-based indices of frozen atoms.
- ``energy`` (float or None) — Total energy in eV (ASE convention).
- ``forces`` (numpy.ndarray or None) — Per-atom forces in eV/Å.
- ``**kwargs`` — Additional ASE Atoms keyword arguments (``symbols``, ``positions``, ``cell``, ``pbc``, etc.).

**Attributes:**

- ``charge`` (int) — Get/set molecular charge.
- ``multiplicity`` (int) — Get/set spin multiplicity.
- ``frozen_atoms``, ``energy``, ``forces`` — Direct access to the corresponding fields.
- Also inherits every method and property from :class:`ase.Atoms`.

Helper Functions
----------------

Bond Cutoffs and Covalent Radii
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Helper functions for computing covalent radii and bond cutoff distances:

.. code:: python

   from chemsmart.io.molecules import get_covalent_radius, get_bond_cutoff

   # Get covalent radius of an element (in Ångströms)
   r_C = get_covalent_radius("C")   # ~0.77 Å
   r_O = get_covalent_radius("O")   # ~0.73 Å
   print("r_C:", r_C, "r_O:", r_O)

   # Compute a bond cutoff for a pair of elements (includes a buffer)
   cutoff_CO = get_bond_cutoff("C", "O")   # r_C + r_O + 0.3 (default buffer)
   cutoff_CH = get_bond_cutoff("C", "H", buffer=0.2)   # custom buffer
   print("cutoff_C–O:", cutoff_CO, "cutoff_C–H:", cutoff_CH)

**Functions:**

:func:`~chemsmart.io.molecules.get_covalent_radius(element_symbol)`
Return the covalent radius of a given element. Source: ASE's standard covalent-radii lookup table.

- ``element_symbol`` (str) — Chemical symbol (case-insensitive).
- Returns: ``float`` — Covalent radius in Ångströms.

:func:`~chemsmart.io.molecules.get_bond_cutoff(el1, el2, buffer=0.3)`
Compute a bond-detection cutoff distance between two elements as ``r(el1) + r(el2) + buffer``.

- ``el1``, ``el2`` (str) — Element symbols for the two atoms.
- ``buffer`` (float) — Extra distance tolerance added to the sum of covalent radii (default 0.3 Å).
- Returns: ``float`` — Bond cutoff in Ångströms.

PubChem Search
^^^^^^^^^^^^^^

Programmatic search of the PubChem database:

.. code:: python

   from chemsmart.io.molecules.pubchem import pubchem_search

   # Search by compound name
   result = pubchem_search("aspirin", fail_silently=True)
   if result is not None:
       print("Aspirin formula:", result.chemical_formula)   # C9H8O4

   # Search by CID
   result = pubchem_search(cid="2244", fail_silently=True)
   if result is not None:
       print("CID 2244 formula:", result.chemical_formula)

   # Search by SMILES
   result = pubchem_search(smiles="CC(=O)OC1=CC=CC=C1C(=O)O", fail_silently=True)
   if result is not None:
       print("SMILES formula:", result.chemical_formula)

   # Low-level raw search
   from chemsmart.io.molecules.pubchem import search_pubchem_raw
   raw = search_pubchem_raw("aspirin", "name", suffix="3d")
   print("Raw result (CID):", raw.get("CID") if isinstance(raw, dict) else "skipped")

**Functions:**

:func:`~chemsmart.io.molecules.pubchem.pubchem_search(name=None, cid=None, smiles=None, fail_silently=False, timeout=10)`
High-level PubChem search that returns a CHEMSMART :class:`Molecule` on success. Exactly one of ``name``, ``cid``, or ``smiles`` must be provided.

- ``name`` (str) — Search by common or IUPAC compound name.
- ``cid`` (str) — Search by PubChem CID (string or integer).
- ``smiles`` (str) — Search by SMILES string.
- ``fail_silently`` (bool) — If ``True`` return ``None`` on error instead of raising.
- ``timeout`` (int) — Request timeout in seconds.
- Returns: :class:`Molecule` or ``None``.

:func:`~chemsmart.io.molecules.pubchem.search_pubchem_raw(query, namespace, suffix="2d", timeout=10)`
Low-level raw JSON search against PubChem PUG REST. Useful if you need metadata beyond the parsed Molecule.

- ``query`` (str) — Search term (name / CID / SMILES / InChIKey, etc.).
- ``namespace`` (str) — PubChem namespace, e.g. ``"name"``, ``"cid"``, ``"smiles"``, ``"inchikey"``.
- ``suffix`` (str) — Representation flavour: ``"2d"`` or ``"3d"``.
- ``timeout`` (int) — Request timeout in seconds.
- Returns: ``dict`` — Raw JSON data from PubChem, or ``None`` if the request fails.

See Also
--------

-  :doc:`molecule-input-formats` — Overview of supported input formats and CLI usage
-  :doc:`cli-overview` — Command-line interface reference
-  :doc:`thermochemistry-analysis` — Thermochemistry analysis using the Python API
-  :doc:`modules` — Full auto-generated module reference

For more technical details on the implementation, see the CHEMSMART preprint: https://arxiv.org/abs/2508.20042