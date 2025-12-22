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

**NEW**: Chemsmart now supports reading molecular structures directly from ChemDraw files! This enables seamless
integration with chemical drawing tools.

.. code:: bash

   # Basic usage with ChemDraw file
   chemsmart sub -s server gaussian -p project -f molecule.cdxml -c 0 -m 1 opt

   # Complete single-command workflow
   chemsmart sub -s server gaussian -p project -f benzene.cdxml -c 0 -m 1 opt

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

   -  Both binary (``.cdx``) and XML-based (``.cdxml``) ChemDraw formats are supported
   -  RDKit is used internally to parse ChemDraw files and generate 3D coordinates
   -  For multi-molecule ChemDraw files, use ``-i`` to select a specific molecule
   -  3D coordinates are automatically generated from 2D structures

**Example: Multi-molecule ChemDraw file**

.. code:: bash

   # Select the 2nd molecule from a ChemDraw file with multiple structures
   chemsmart sub -s server gaussian -p project -f molecules.cdxml -i 2 -c 0 -m 1 opt

CIF Files (Crystallographic Information Files)
==============================================

**NEW**: Chemsmart now supports reading molecular structures from CIF (Crystallographic Information File) format,
commonly used for crystal structures from databases like the Cambridge Crystallographic Data Centre (CCDC).

Reading Local CIF Files
-----------------------

CIF files contain crystallographic data including atomic coordinates and unit cell parameters:

.. code:: bash

   # Gaussian job with CIF file
   chemsmart sub -s server gaussian -p project -f structure.cif -c 0 -m 1 opt

   # ORCA job with CIF file
   chemsmart sub -s server orca -p project -f structure.cif -c 0 -m 1 opt

**Python API:**

.. code:: python

   from chemsmart.io.molecules.structure import Molecule

   # Read a single structure from a CIF file
   mol = Molecule.from_cif_file('structure.cif')

   # Read with return_list option
   mols = Molecule.from_cif_file('structure.cif', return_list=True)

   # Also works with from_filepath (auto-detects .cif extension)
   mol = Molecule.from_filepath('structure.cif')

.. note::

   -  CIF files are commonly used for crystallographic data
   -  Chemsmart uses ASE (Atomic Simulation Environment) to parse CIF files
   -  For multi-structure CIF files, use ``-i`` to select a specific structure
   -  Charge and multiplicity must be specified with ``-c`` and ``-m`` flags

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

CCDC Database Integration
=========================

**NEW**: Retrieve crystal structures directly from the Cambridge Crystallographic Data Centre (CCDC) database by
deposition number.

.. code:: bash

   # Gaussian optimization with CCDC structure
   chemsmart sub -s server gaussian -p project -C 1428476 -c 0 -m 1 opt

   # ORCA single point with CCDC structure
   chemsmart sub -s server orca -p project -C 1428476 -c 0 -m 1 sp

   # With custom label
   chemsmart sub -s server gaussian -p project -C 1428476 -l complex -c 0 -m 1 opt

**Python API:**

.. code:: python

   from chemsmart.io.molecules.structure import Molecule

   # Retrieve structure from CCDC by deposition number
   mol = Molecule.from_ccdc(1428476)

   # With return_list option
   mols = Molecule.from_ccdc(1428476, return_list=True)

.. tip::

   The CCDC option ``-C`` works similarly to the PubChem option ``-P``:

   -  Use ``-l`` to specify a custom output filename
   -  Downloaded structures are automatically cached to avoid repeated downloads
   -  Cache location: ``/tmp/chemsmart_ccdc_cache/``

**Configuring CCDC URL:**

The CCDC base URL can be configured via environment variable if you need to use a different endpoint:

.. code:: bash

   export CCDC_BASE_URL="https://your-ccdc-endpoint.com"
   chemsmart sub -s server gaussian -p project -C 1428476 -c 0 -m 1 opt

.. warning::

   The default CCDC URL is a placeholder. Users should verify and configure the actual CCDC API endpoint based on their
   institutional access and licensing requirements.

**Caching:**

CCDC structures are automatically cached to improve performance:

-  Default cache location: ``/tmp/chemsmart_ccdc_cache/``
-  Cache files are named by deposition number (e.g., ``1428476.cif``)
-  Cached files are reused across sessions to minimize network requests

**Error Handling:**

Clear error messages are provided for common issues:

-  Invalid deposition numbers
-  Network failures (with automatic retry)
-  Missing or restricted structures
-  CIF parsing errors

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
-  CIF files
-  ASE database/trajectory files
-  PubChem queries
-  CCDC queries
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
-  ``.cif`` → CIF (Crystallographic Information File) format
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
