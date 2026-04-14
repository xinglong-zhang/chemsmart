.. _chemdraw-organometallic:

#######################################
 ChemDraw Organometallic Complex Files
#######################################

Chemsmart can read organometallic complexes drawn in ChemDraw (``.cdx`` and ``.cdxml``) and generate 3D structures
suitable for quantum chemistry calculations with Gaussian or ORCA.

.. note::

   Organometallic complex support requires RDKit and, for binary ``.cdx`` files, Open Babel. See :ref:`installation
   requirements <chemdraw-organometallic-requirements>` below.

*****
 Why
*****

Organometallic complexes present several challenges when read from ChemDraw files:

-  RDKit raises ``Can't kekulize mol`` errors for aromatic ligands coordinated to a metal centre.
-  RDKit raises ``UFFTYPER: Unrecognized atom type`` errors for transition metals.
-  ChemDraw can store aromatic ligands (e.g. Cp, benzene) as **separate fragments** that need to be combined with the
   metal fragment before 3D coordinates can be generated.

Chemsmart handles all of these cases automatically.

*********************************
 Supported Organometallic Inputs
*********************************

The following types of organometallic complexes drawn in ChemDraw are supported:

-  Transition-metal complexes with **η5-cyclopentadienyl (Cp)** ligands, including **Cp\***
   (pentamethylcyclopentadienyl)
-  Transition-metal complexes with **η6-arene** (e.g. benzene) ligands
-  General transition-metal complexes with ancillary phosphine, amine, carbonyl, or halide ligands
-  Mixed complexes combining aromatic and non-aromatic ligands

Example usage:

.. code:: bash

   # Submit a Gaussian optimization for a ferrocene-like complex
   chemsmart sub -s server gaussian -p project -f ferrocene.cdxml -c 0 -m 1 opt

   # Binary CDX format (requires Open Babel)
   chemsmart sub -s server gaussian -p project -f complex.cdx -c 0 -m 1 opt

   # Multi-molecule file: select the second molecule
   chemsmart sub -s server gaussian -p project -f complexes.cdxml -i 2 -c 2 -m 3 opt

.. _chemdraw-organometallic-requirements:

**************
 Requirements
**************

+------------------+------------------------------------------+-----------------------------------+
| Dependency       | Purpose                                  | Required for                      |
+==================+==========================================+===================================+
| RDKit            | Parse ``.cdxml``; sanitize molecules     | Both ``.cdx`` and ``.cdxml``      |
+------------------+------------------------------------------+-----------------------------------+
| Open Babel CLI   | Convert binary ``.cdx`` to SDF for RDKit | ``.cdx`` files only               |
| (``obabel``)     |                                          |                                   |
+------------------+------------------------------------------+-----------------------------------+

If ``obabel`` is not installed and a ``.cdx`` file is provided, Chemsmart raises a ``ValueError`` with instructions to
install Open Babel or re-save the file as ``.cdxml``.

**************
 How It Works
**************

Chemsmart applies the following pipeline when reading ChemDraw files:

#. **Parse without sanitization** – ``sanitize=False`` is passed to RDKit (``MolsFromCDXMLFile``) to avoid kekulization
   errors during initial parsing. For ``.cdx`` files, Open Babel converts the binary format to SDF first.

#. **Update property cache** – ``mol.UpdatePropertyCache(strict=False)`` is called on every molecule before processing
   to avoid *pre-condition violation* errors.

#. **Combine metal and ligand fragments** – ChemDraw sometimes stores the metal centre and its aromatic ligands as
   separate molecules in the same file. Chemsmart detects a small metal-containing fragment followed by aromatic ring
   fragments and merges them into a single molecule automatically.

#. **Normalize metal bonds** – Aromatic bond flags on any bond involving a metal atom are removed (set to single bonds).
   RDKit does not support aromatic bonds to metal centres.

#. **Fix cyclopentadienyl aromaticity** – Neutral aromatic 5-membered carbon rings (``c1cccc1``) that cannot be
   kekulized are converted to the cyclopentadienyl anion (``[cH-]``) so that RDKit can process them correctly.

#. **Add η5 coordination bond** – One single bond is added from the metal to an anchor carbon in each remaining aromatic
   5-membered ring to represent η5 coordination. The ring is de-aromatized to an alternating single/double pattern
   before the bond is added.

#. **Selective sanitization** – Kekulization is skipped for molecules that contain metals (to avoid ``Can't kekulize
   mol`` errors) while all other sanitization steps are applied normally.

#. **Add hydrogens and generate 3D coordinates** – Explicit hydrogens are added with ``AddHs``, then 3D coordinates are
   generated with ``EmbedMolecule`` (ETKDG) and optionally refined with MMFF. If MMFF fails (common for exotic metal
   types), the ETKDG geometry is kept.

**********************
 Current Restrictions
**********************

.. warning::

   The following restrictions apply to the current organometallic complex support.

Coordination Geometry
   3D coordinate generation uses the RDKit ETKDG embedding algorithm, which was designed for organic molecules. Metal
   coordination geometries (octahedral, square-planar, tetrahedral at metal centres) are approximated rather than
   explicitly enforced. The resulting 3D structure should be used as a **starting geometry for a DFT optimization**, not
   as a final structure.

η5 Coordination Representation
   Cp and Cp\* η5 coordination is represented with a **single metal–carbon σ-bond** to one ring carbon. This is a
   structural approximation that allows RDKit to build a valid connectivity graph and generate 3D coordinates. The bond
   order is **not** meaningful for electronic structure purposes.

η6 Arene Coordination
   Full η6 metal–arene coordination (e.g. benzene complexes) is represented as a disconnected fragment combined with the
   metal centre. No explicit metal–carbon bonds are added for 6-membered aromatic rings; the metal and the arene ring
   are present in the molecular formula but their interaction is only implied.

Multi-Hapto Ligands Beyond Cp/Benzene
   Higher-order hapticity ligands (η7-cycloheptatrienyl, η8-cyclooctatetraene, etc.) and non-carbon η-donors are not
   explicitly handled and may produce errors or incomplete structures.

Force-Field Optimization
   The MMFF94 force field does not have parameters for most transition metals. Geometry optimization with MMFF is
   attempted but silently skipped on failure. The ETKDG embedded geometry is used in that case and may require further
   DFT optimization.

Charge and Multiplicity
   Charge and multiplicity of organometallic complexes are **not** inferred from the ChemDraw file. You must always
   specify them explicitly with ``-c`` and ``-m``:

   .. code:: bash

      # Iron(II) complex with overall charge 2+ and singlet multiplicity
      chemsmart sub -s server gaussian -p project -f fecpx.cdxml -c 2 -m 1 opt

   Incorrect charge or multiplicity will lead to a failed quantum chemistry calculation.

Fragment Combination Heuristic
   The fragment-combination step uses the heuristic that a small (< 10 heavy atoms) metal-containing fragment followed
   immediately by aromatic ring fragments should be merged. This may not work correctly for very unusual ChemDraw
   layouts or complex multi-metal systems.

**********
 Examples
**********

Ferrocene Derivative (Cp₂Fe)
============================

Draw the complex in ChemDraw (or use an existing ``.cdxml`` file) and run:

.. code:: bash

   chemsmart sub -s server gaussian -p project -f ferrocene.cdxml -c 0 -m 1 opt B3LYP/def2-SVP

Half-Sandwich Complex
=====================

For a half-sandwich complex such as [CpFe(CO)₂Cl]:

.. code:: bash

   chemsmart sub -s server gaussian -p project -f half_sandwich.cdxml -c 0 -m 2 opt B3LYP/def2-SVP

.. tip::

   For open-shell transition-metal complexes, always verify the multiplicity. A d⁵ Fe(III) centre in a weak-field
   environment typically has multiplicity 6 (high-spin), whereas in a strong-field environment it may be 2 (low-spin).

**********
 See Also
**********

-  :doc:`molecule-input-formats` – all supported input formats
-  :doc:`gaussian-cli-options` – available Gaussian calculation options
-  :doc:`orca-cli-options` – available ORCA calculation options
