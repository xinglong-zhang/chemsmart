#####################
 Database Management
#####################

CHEMSMART provides a unified database subsystem for managing, querying, and reusing quantum chemistry calculation data.
The ``chemsmart run database`` command group converts heterogeneous computational output files and calculation folders
into a structured, identity-aware, and queryable database for downstream analysis, machine learning, and CHEMSMART
workflow automation.

************
 Motivation
************

Computational chemistry projects typically generate hundreds to thousands of output files spread across different
directories, programs, and levels of theory. Without a unified data model, it can be difficult to manage, compare,
query, and report these results in a consistent and reproducible way.

The CHEMSMART database module addresses this challenge by converting fragmented calculation outputs into normalized,
reusable database records. It is designed to:

-  parse computational output files into schema-driven records;
-  separate molecular identity from calculation context;
-  deduplicate chemical species and 3D conformers across calculations;
-  store everything in a portable, single-file CHEMSMART database (``.db``);
-  provide CLI tools for inspecting, querying, filtering, and exporting data;
-  enable stored records and structures to be reused directly in downstream CHEMSMART workflows.

Together, these features allow computational results to be treated not as isolated output files, but as structured,
traceable, and reusable scientific data.

************
 Data Model
************

The CHEMSMART database separates calculation data into three core entities: ``record``, ``molecule``, and ``structure``.
These entities represent *what was computed*, *what chemical species was studied*, and *which 3D geometry was used*.

.. list-table::
   :header-rows: 1
   :widths: 18 32 50

   -  -  Object
      -  Represents
      -  Description

   -  -  ``record``
      -  One calculation
      -  Identified by ``record_id``. Stores the calculation context, parsed results, runtime information, and
         provenance for one calculation.

   -  -  ``molecule``
      -  One chemical species
      -  Identified by ``molecule_id``. Stores connectivity-level molecular identity, such as formula, SMILES, InChI,
         composition, and basic molecular descriptors.

   -  -  ``structure``
      -  One 3D geometry
      -  Identified by ``structure_id``. Stores geometry-level information, including coordinates, charge, multiplicity,
         and simple geometric descriptors.

This organization allows CHEMSMART to recognize when the same molecule or conformer appears in multiple calculations,
while preserving calculation-specific information needed for analysis and reproducibility.

Identifiers
===========

.. list-table::
   :header-rows: 1
   :widths: 22 18 60

   -  -  Identifier
      -  CLI option
      -  Notes

   -  -  ``record_id``
      -  ``--rid``
      -  Hash-derived identifier for a calculation record. Prefixes of at least 12 characters are accepted.

   -  -  ``record_index``
      -  ``--ri``
      -  1-based index assigned in insertion order. Useful for quick interactive inspection.

   -  -  ``molecule_id``
      -  ``--mid``
      -  InChIKey-derived identifier for a molecular identity. Prefixes of at least 16 characters are accepted.

   -  -  ``structure_id``
      -  ``--sid``
      -  Hash-derived identifier based on geometry, charge, and multiplicity. Prefixes of at least 12 characters are
         accepted.

.. note::

   Partial IDs are resolved automatically. If a prefix uniquely matches one entry, the command proceeds; otherwise an
   explicit error is raised with the matching candidates.

**********
 Workflow
**********

A typical CHEMSMART database workflow is:

#. **Assemble** output files from a project directory into a single ``.db`` file using ``chemsmart run database
   assemble`` (see :doc:`database-assemble`).
#. **Query** database entries, with or without structured filters, using ``chemsmart run database query`` (see
   :doc:`database-query`).
#. **Inspect** selected entries in detail using ``chemsmart run database inspect`` (see :doc:`database-inspect`).
#. **Export** selected content as JSON, CSV, XYZ, or extended XYZ (with per-atom forces) for downstream analysis using
   ``chemsmart run database export`` (see :doc:`database-export`).
#. **Reuse** database entries in downstream chemsmart workflows (see :doc:`database-workflow`).

*****************
 Command Summary
*****************

.. list-table::
   :header-rows: 1
   :widths: 20 80

   -  -  Subcommand
      -  Purpose
   -  -  ``assemble``
      -  Parse computational chemistry output files and build a CHEMSMART database.
   -  -  ``inspect``
      -  Print detailed information for databases, records, molecules, and structures.
   -  -  ``query``
      -  Filter records, molecules, or structures using a structured query expression.
   -  -  ``export``
      -  Export the database or selected entries to JSON, CSV, XYZ, or extended XYZ.

*************
 Usage Notes
*************

-  **Supported programs.** Gaussian and ORCA output files and xTB calculation folders are currently parsed.
   Files/folders from other programs are skipped during ``assemble``. Note that xTB results are organized as calculation
   folders rather than single output files.

-  **Idempotent assembly.** Re-running ``assemble`` on the same directory will not duplicate records. Entries with the
   same ``record_id`` are replaced, and the number of replaced duplicates is reported.

-  **Custom basis identity.** Gaussian GEN/GENECP output stores the printed orbital functions in
   ``heavy_elements_basis_functions`` and any effective core potentials in ``heavy_elements_ecp``.
   ``heavy_elements_basis`` stores a common inferred name, ``mixed`` for different identified basis sets, or null when
   identification is incomplete. ``heavy_elements_basis_by_element`` records each element's inferred name; unidentified
   names remain null without discarding the printed functions or ECPs. For unchanged parsed functions and ECPs, record
   IDs remain compatible with the earlier functions-only metadata format; choosing a different name for the same
   functions does not create a new calculation record.

-  **Gaussian basis reuse.** For verbose ``#p`` Gen/GenECP output, basis names are inferred from Basis Set Exchange only
   when both the orbital functions and the ECP match, including the absence of an ECP. Orbital-only matches remain
   available in ``heavy_elements_basis_candidates``. If no complete match exists, supply an explicit basis name or a
   ``gen_genecp_file`` containing the intended basis and ECP before preparing another calculation.

-  **Database validation.** Every database operation first verifies that the input is a valid CHEMSMART database with
   the expected schema. Foreign ``.db`` files are rejected with a clear error.

-  **Portability.** A CHEMSMART database is stored as a single portable ``.db`` file. It can be copied between machines
   and shared alongside the source outputs for reproducibility, since each record carries provenance metadata.

***********************
 Important Assumptions
***********************

-  Molecular identity depends on the connectivity information that can be inferred from the parsed output. For unusual
   bonding situations, organometallic complexes, ion pairs, or transition states, users should inspect the assigned
   molecule and structure identifiers carefully.

-  Structure identity is geometry-based and includes charge and multiplicity. The same Cartesian geometry with different
   charge or spin state is treated as a different structure.

-  Molecular descriptors such as ``is_aromatic``, ``is_multicomponent``, ``smiles``, and ``inchi`` are derived from 3D
   geometry using RDKit heuristics, not from first-principles electronic structure. Bond orders are inferred from
   interatomic distances, which may give incorrect results for unusual bonding, metal-ligand interactions, or strained
   systems. These tags should be treated as approximate.

ORCA custom basis metadata
==========================

ORCA assembly stores named element-specific NewGTO assignments and expanded NewGTO functions printed in the output.
Unprinted functions remain unavailable. Printed ECP summaries retain the name, replaced core-electron count and lmax;
``parameters_available: false`` indicates that potential coefficients were not parsed and the summary is not a complete
reconstructable ECP.

Custom definitions echoed in the output are retained as ``custom_basis_input`` even when expanded functions are absent.
This includes complete tokenized ``%basis`` blocks, with comments and echo line numbers removed. For atom-specific or
other overrides outside those blocks, the tokenized input is retained to preserve assignment context. These definitions
participate in record identity so different unexpanded custom bases cannot overwrite each other. Incomplete printed
function blocks are not presented as complete functions. Standard route-basis expansions remain metadata and do not
change record identity across ORCA versions.

A named ORCA basis augmented with extra shells (inside NewGTO or via AddGTO) has an unidentified final name rather than
retaining the original library nickname. Both element symbols and atomic numbers are accepted for element assignments.
All echoed ``%basis`` definitions, including global Basis/ECP settings and GTOName file references, participate in
custom basis identity. A file reference preserves its filename, not the contents of an external file absent from the
output.
