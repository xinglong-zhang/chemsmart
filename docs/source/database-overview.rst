#####################
 Database Management
#####################

Chemsmart provides a unified database subsystem for managing, querying, and reusing quantum chemistry calculation data.
The ``chemsmart run database`` command group converts heterogeneous computational output files into a structured,
identity-aware, and queryable database for downstream analysis, machine learning, and chemsmart workflow automation.

************
 Motivation
************

Computational chemistry projects typically generate hundreds to thousands of output files spread across different
directories, programs, and levels of theory. Without a unified data model, it can be difficult to manage, compare,
query, and report these results in a consistent and reproducible way.

The chemsmart database module addresses this challenge by converting fragmented calculation outputs into normalized,
reusable database records. It is designed to:

-  parse computational output files into schema-driven records;
-  separate molecular identity from calculation context;
-  deduplicate chemical species and 3D conformers across calculations;
-  store everything in a portable, single-file chemsmart database (``.db``);
-  provide CLI tools for inspecting, querying, filtering, and exporting data;
-  enable stored records and structures to be reused directly in downstream chemsmart workflows.

Together, these features allow computational results to be treated not as isolated output files, but as structured,
traceable, and reusable scientific data.

************
 Data Model
************

The chemsmart database separates calculation data into three core entities: ``record``, ``molecule``, and ``structure``.
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
         provenance for one output file.

   -  -  ``molecule``
      -  One chemical species
      -  Identified by ``molecule_id``. Stores connectivity-level molecular identity, such as formula, SMILES, InChI,
         composition, and basic molecular descriptors.

   -  -  ``structure``
      -  One 3D geometry
      -  Identified by ``structure_id``. Stores geometry-level information, including coordinates, charge, multiplicity,
         and simple geometric descriptors.

This organization allows chemsmart to recognize when the same molecule or conformer appears in multiple calculations,
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

A typical chemsmart database workflow is:

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
      -  Parse computational chemistry output files and build a chemsmart database.
   -  -  ``inspect``
      -  Print detailed information for databases, records, molecules, and structures.
   -  -  ``query``
      -  Filter records, molecules, or structures using a structured query expression.
   -  -  ``export``
      -  Export the database or selected entries to JSON, CSV, XYZ, or extended XYZ.

*************
 Usage Notes
*************

-  **Supported programs.** Only Gaussian and ORCA output files are currently parsed. Files from other programs are
   skipped during ``assemble``.
-  **Idempotent assembly.** Re-running ``assemble`` on the same directory will not duplicate records. Entries with the
   same ``record_id`` are replaced, and the number of replaced duplicates is reported.
-  **Database validation.** Every database operation first verifies that the input is a valid chemsmart database with
   the expected schema. Foreign ``.db`` files are rejected with a clear error.
-  **Portability.** A chemsmart database is stored as a single portable ``.db`` file. It can be copied between machines
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
