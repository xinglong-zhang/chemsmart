########################
 Molecule Input Formats
########################

Chemsmart provides flexible molecule input capabilities, supporting multiple file formats and molecular representations.
This page describes the various ways you can create molecules for use in quantum chemistry calculations.

.. note::

   For more information about CHEMSMART's design and capabilities, see the preprint: https://arxiv.org/abs/2508.20042

*******************
 Workflow Overview
*******************

The following diagram illustrates how different input formats are processed to create molecules for quantum chemistry
calculations:

.. code:: text

   ┌──────────────────┐
   │  Input Sources   │
   └──────────────────┘
           │
           ├─── ASE Atoms
           ├─── Pymatgen Molecule
           ├─── RDKit Molecule
           ├─── File Formats (.xyz, .sdf, .com/.log, .inp/.out, .cdx/.cdxml, etc.)
           └─── PubChem queries (by name, CID, or SMILES)
           │
           ▼
   ┌──────────────────────────────┐
   │  Input Processing            │  ──────►  Creates Molecule Object
   │  (chemsmart.io)              │           (parses all supported input formats)
   └──────────────────────────────┘
           │
           ▼
   ┌──────────────────┐
   │    Molecule      │  ──────►  Central molecular representation
   │     Object       │           (symbols, positions, charge, multiplicity)
   └──────────────────┘
           │
           ├─────────────────────────────────┐
           ▼                                 ▼
   ┌──────────────────┐            ┌──────────────────┐
   │ chemsmart.jobs   │            │ chemsmart.jobs   │
   │ .gaussian.writer │            │   .orca.writer   │
   └──────────────────┘            └──────────────────┘
           │                                 │
           ▼                                 ▼
   ┌──────────────────┐            ┌──────────────────┐
   │  .com/.gjf files │            │    .inp files    │
   │  (Gaussian, Inc.)│            │      (ORCA)      │
   └──────────────────┘            └──────────────────┘

************************
 Supported File Formats
************************

Coordinate Files
================

XYZ Files
---------

Standard XYZ format for atomic coordinates:

.. code:: bash

   # Single structure
   chemsmart sub -s server gaussian -p project -f molecule.xyz -c 0 -m 1 opt

   # Multi-structure file (select 5th structure)
   chemsmart sub -s server gaussian -p project -f molecules.xyz -i 5 -c 0 -m 1 opt

.. note::

   XYZ files do not contain charge or multiplicity information. You must specify these using ``-c`` and ``-m`` flags.

SDF Files
---------

Structure-Data Format files with 2D or 3D coordinates:

.. code:: bash

   chemsmart sub -s server gaussian -p project -f molecule.sdf -c 0 -m 1 opt

Gaussian Files
==============

Gaussian Input Files (.com, .gjf)
---------------------------------

Gaussian input files contain route, charge, multiplicity, and coordinates:

.. code:: bash

   chemsmart sub -s server gaussian -p project -f input.com opt

.. tip::

   Chemsmart reads the existing charge and multiplicity from ``.com`` files. Override with ``-c`` and ``-m`` if needed.

Gaussian Output Files (.log, .out)
----------------------------------

Use optimized structures from Gaussian output:

.. code:: bash

   # Use last structure from optimization
   chemsmart sub -s server gaussian -p project -f opt.log sp

   # Use specific structure (e.g., 10th structure from scan)
   chemsmart sub -s server gaussian -p project -f scan.log -i 10 sp

ORCA Files
==========

ORCA Input Files (.inp)
-----------------------

ORCA input files with coordinates and calculation settings:

.. code:: bash

   chemsmart sub -s server gaussian -p project -f input.inp opt

ORCA Output Files (.out)
------------------------

Extract structures from ORCA calculations:

.. code:: bash

   chemsmart sub -s server gaussian -p project -f orca_opt.out sp

ChemDraw Files
==============

ChemDraw XML (.cdxml) and Binary (.cdx) Files
---------------------------------------------

Chemsmart supports reading molecular structures directly from ChemDraw files, including **organometallic complexes**
with aromatic ligands such as Cp, Cp\*, and benzene rings.

.. code:: bash

   # Organic molecule
   chemsmart sub -s server gaussian -p project -f molecule.cdxml -c 0 -m 1 opt

   # Organometallic complex (charge and multiplicity must be specified explicitly)
   chemsmart sub -s server gaussian -p project -f ferrocene.cdxml -c 0 -m 1 opt

.. tip::

   You can submit a Gaussian optimization directly from a ChemDraw file in a single command:

   .. code:: bash

      chemsmart sub -s server gaussian -p project -f molecule.cdxml -c 0 -m 1 opt

   This will:

   #. Read the molecular structure from the ChemDraw file
   #. Generate 3D coordinates using RDKit's EmbedMolecule
   #. Create a Gaussian input file with specified charge and multiplicity
   #. Submit the geometry optimization job to the HPC cluster

.. note::

   -  Both binary (``.cdx``) and XML-based (``.cdxml``) ChemDraw formats are supported.
   -  RDKit is used internally to parse ChemDraw files and generate 3D coordinates.
   -  For multi-molecule ChemDraw files, use ``-i`` to select a specific molecule.
   -  3D coordinates are automatically generated from 2D structures.
   -  Reading binary ``.cdx`` files requires Open Babel (``obabel``) to be installed. If Open Babel is not available,
      save the file as ``.cdxml`` instead.
   -  Charge and multiplicity of organometallic complexes are **not** inferred from the ChemDraw file – always specify
      ``-c`` and ``-m`` explicitly.

**Example: Multi-molecule ChemDraw file**

.. code:: bash

   # Select the 2nd molecule from a ChemDraw file with multiple structures
   chemsmart sub -s server gaussian -p project -f molecules.cdxml -i 2 -c 0 -m 1 opt

For full details on organometallic complex support and its restrictions, see :doc:`chemdraw-organometallic`.

Image Files
===========

Image Files (.png, .jpg, .jpeg, .tif, .tiff)
---------------------------------------------

Chemsmart can read 2D molecular drawings from image files and convert them to a
``Molecule`` object.

.. code:: bash

   chemsmart sub -s server gaussian -p project -f molecule.png -c 0 -m 1 opt

Processing pipeline
~~~~~~~~~~~~~~~~~~~

When ``Molecule.from_filepath()`` receives an image path, chemsmart:

#. Loads the image in grayscale and applies preprocessing (resize + Otsu threshold).
#. Runs DECIMER to predict a SMILES string from the structure drawing.
#. Optionally runs OCR (``pytesseract``) to detect text abbreviations in the image.
#. If DECIMER output is missing/too short and recognized abbreviations are detected,
   uses abbreviation-based fallback assembly for simple dash-separated labels.
#. Converts the final SMILES to a chemsmart ``Molecule``.

Python dependencies
~~~~~~~~~~~~~~~~~~~

Required for image parsing:

- ``decimer``
- ``opencv-python``
- ``Pillow`` (PIL backend)

Optional but recommended for abbreviation detection:

- ``pytesseract`` (Python wrapper)

Install commands:

.. code:: bash

   # Core image parsing stack
   pip install decimer opencv-python Pillow

   # Optional OCR support for abbreviation expansion
   pip install chemsmart[image]
   # or
   pip install pytesseract

.. note::

   ``pytesseract`` also requires the Tesseract OCR binary to be available on the
   system PATH.

Recognized abbreviations
~~~~~~~~~~~~~~~~~~~~~~~~

The image parser currently recognizes these built-in abbreviations from
``chemsmart.io.molecules.CHEMICAL_ABBREVIATIONS``:

- ``Ad, Ph, Me, Et, nPr, iPr, Bu, nBu, iBu, sBu, tBu, nPent, nHex, Bn, Ac, Bz``
- ``Ts, Ms, Tf, Cy, OMe, OEt, NMe2, CF3, NO2, CHO, CO2H, CO2Me, CO2Et, CN, N3``
- ``Vinyl, Allyl, Propargyl, Piv, OMs, OTs, OTf``
- ``Boc, Cbz, Alloc, Fmoc, MOM, TMS, TBDMS, TBDPS, Trt``

For dash-separated shorthand, the fallback parser also recognizes these terminal
substituent labels from ``SUBSTITUENT_MAPPING`` (both ``NH2`` and Unicode
``NH₂`` are supported):

- ``SH, OH, NH2, NH₂``

Limitations
~~~~~~~~~~~

- Image input currently supports only a **single** molecule per file.
- For image files, ``index`` and ``return_list`` arguments are ignored.
- Abbreviation fallback is heuristic and limited to currently defined labels.
- OCR quality and image quality (resolution, contrast, font clarity) strongly
  affect recognition accuracy.
- Without ``pytesseract``, abbreviation-aware fallback is disabled.
- Incorrect DECIMER/OCR interpretations can still produce wrong structures and
  should be manually verified for critical workflows.

*********************
 Molecular Databases
*********************

PubChem Integration
===================

Query PubChem directly by name, CID, or SMILES:

.. code:: bash

   # By compound name
   chemsmart sub -s server gaussian -p project -P water -c 0 -m 1 -l water opt

   # By CID (Compound ID)
   chemsmart sub -s server gaussian -p project -P 962 -c 0 -m 1 -l water opt

   # By SMILES string
   chemsmart sub -s server gaussian -p project -P "O=Cc1ccccc1" -c 0 -m 1 -l benzaldehyde opt

.. note::

   -  When using ``-P`` (PubChem), the ``-l`` flag specifies the output filename.
   -  PubChem queries support compound names, CIDs (Compound IDs), and SMILES strings.
   -  SMILES strings are processed by PubChem to retrieve 3D structures and generate coordinates.

ASE Database Files (.db, .traj)
===============================

Use structures from ASE database or trajectory files:

.. code:: bash

   # From ASE database
   chemsmart sub -s server gaussian -p project -f molecules.db -i 5 -c 0 -m 1 opt

   # From ASE trajectory
   chemsmart sub -s server gaussian -p project -f md.traj -i 10 -c 0 -m 1 opt

***************************
 Python Object Integration
***************************

Chemsmart's ``Molecule`` class provides seamless integration with popular Python chemistry libraries:

From ASE Atoms
==============

.. code:: python

   from ase import Atoms
   from chemsmart.io.molecules.structure import Molecule

   # Create ASE Atoms object
   atoms = Atoms('H2O', positions=[[0, 0, 0], [0, 0, 1], [0, 1, 0]])

   # Convert to Chemsmart Molecule
   molecule = Molecule.from_ase_atoms(atoms)

From RDKit Mol
==============

.. code:: python

   from rdkit import Chem
   from rdkit.Chem import AllChem
   from chemsmart.io.molecules.structure import Molecule

   # Create RDKit molecule from SMILES
   rdkit_mol = Chem.MolFromSmiles('CCO')
   AllChem.EmbedMolecule(rdkit_mol)

   # Convert to Chemsmart Molecule
   molecule = Molecule.from_rdkit_mol(rdkit_mol)

Aromaticity Detection (``is_aromatic``)
---------------------------------------

The ``is_aromatic`` property detects aromaticity by converting the molecule to an RDKit representation and checking
whether any atom both carries the aromatic flag **and** belongs to a ring. This guards against false positives in
acyclic molecules (e.g. H₂O, MgI₂) that can arise when the geometry-based bond-order heuristic assigns a bond order of
1.5 to short single bonds.

.. note::

   **Known limitations of aromaticity detection:**

   -  Bond orders are inferred purely from 3D geometry (interatomic distances), not from an electronic structure
      calculation. This means the detection is **model-dependent** and may not match formal aromaticity criteria in all
      cases.

   -  **Edge cases** such as the cyclopropenyl cation (aromatic) versus the cyclopropenyl radical (non-aromatic) may not
      be distinguished correctly, because the outcome depends on how bond orders and electron counts are assigned from
      the geometry alone.

   -  For borderline or unusual systems (strained rings, metal-organic frameworks, non-Kekulé structures, etc.) the
      result should be treated as a heuristic estimate rather than a definitive answer.

   -  If precise aromaticity information is required, consider constructing the RDKit molecule directly from a SMILES
      string or from an output file that encodes explicit bond orders.

From Pymatgen
=============

Chemsmart molecules can be converted to and from Pymatgen format:

.. code:: python

   from chemsmart.io.molecules.structure import Molecule

   # Convert Chemsmart Molecule to Pymatgen
   molecule = Molecule.from_filepath('input.xyz')
   pymatgen_mol = molecule.to_pymatgen()

.. note::

   For converting Pymatgen molecules to Chemsmart, you can use the ASE Atoms adaptor as an intermediate format.

****************
 Best Practices
****************

Charge and Multiplicity
=======================

Always specify charge and multiplicity for:

-  XYZ files
-  SDF files
-  ASE database/trajectory files
-  PubChem queries
-  ChemDraw files

.. code:: bash

   # Correct usage with charge and multiplicity
   chemsmart sub -s server gaussian -p project -f molecule.xyz -c -1 -m 2 opt

Structure Selection
===================

For multi-structure files, use ``-i`` to select a specific structure:

.. code:: bash

   # Use 1-based indexing
   chemsmart sub -s server gaussian -p project -f scan.log -i 15 sp

.. warning::

   Chemsmart uses **1-based indexing** to match most molecular visualization software, unlike Python's 0-based indexing.

File Format Auto-Detection
==========================

Chemsmart automatically detects file formats based on extensions:

-  ``.xyz`` → XYZ format
-  ``.sdf`` → SDF format
-  ``.com``, ``.gjf`` → Gaussian input
-  ``.log`` → Gaussian output
-  ``.inp`` → ORCA input
-  ``.out`` → ORCA/Gaussian output (auto-detected by reading file header)
-  ``.cdx``, ``.cdxml`` → ChemDraw format
-  ``.db``, ``.traj`` → ASE database/trajectory

.. note::

   For ``.out`` files, Chemsmart automatically detects whether the file is from ORCA or Gaussian by examining the file
   header. If detection fails, an error will be raised indicating the unsupported format.

For unsupported extensions, Chemsmart falls back to ASE's file reading capabilities.

**********
 See Also
**********

-  :doc:`gaussian-cli-options`
-  :doc:`orca-cli-options`
-  :doc:`cli-overview`

For more technical details on the implementation, see the CHEMSMART preprint: https://arxiv.org/abs/2508.20042
