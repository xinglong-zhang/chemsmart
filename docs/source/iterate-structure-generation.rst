######################
 Structure Generation
######################

This section describes how to use the ``iterate`` command to systematically generate molecular libraries.

**********
 Overview
**********

The ``iterate`` tool automates the process of attaching substituents to skeleton molecules. It is designed to generate a
library of structures by combinatorially attaching a list of substituents to specific sites on a list of skeletons.

The tool handles: - Alignment of the substituent to the specified attachment vector. - Rotation optimization to minimize
steric clashes (using algorithms like Lagrange multipliers). - Batch processing with parallel execution.

********************
 Configuration File
********************

The core of an iterate job is the TOML configuration file, which defines the input molecules and attachment points.

.. code:: toml

   # Example 1: Basic skeleton definition
   [[skeletons]]
   file_path = "core_simple.xyz"
   label = "Core_Simple"
   link_index = "1, 5"          # Atom indices (1-based) on skeleton to attach to

   # Example 2: Skeleton with specific efficient region defined
   [[skeletons]]
   file_path = "core_complex.xyz"
   label = "Core_Complex"
   link_index = [10]
   skeleton_indices = "1-15"    # Only atoms 1-15 are treated as the rigid part
   # Note: The link_index (10) MUST be included within the skeleton_indices (1-15).
   # Otherwise, the program will raise an error.

   [[substituents]]
   file_path = "group.xyz"
   label = "GroupA"
   link_index = 1               # Atom index (1-based) on substituent to attach

   [[substituents]]
   file_path = "group2.gjf"
   label = "GroupB"
   link_index = [2]             # List format is also supported

   [[substituents]]
   file_path = "group3.xyz"
   label = "GroupC"
   link_index = "4-7, 17"       # Range format "start-end" is supported

Explanation of Keys
===================

skeletons
---------

A list of tables defining the core structures.

-  **file_path**: Path to the molecule file (XYZ, GJF, etc.).

-  **label**: Unique identifier for the molecule.

-  **link_index**: One or more atom indices (1-based) where the attachment should occur. This defines *where*
   substituents will be attached.

-  **skeleton_indices** (Optional): Indices defining the rigid atoms of the skeleton. If not provided, the program will
   **automatically detect** and remove the original group (e.g., a hydrogen atom) attached at the ``link_index``,
   treating the remaining structure as the skeleton.

substituents
------------

A list of tables defining the functional groups to attach.

-  **file_path**: Path to the molecule file.
-  **label**: Unique identifier.
-  **link_index**: The atom index (1-based) on the substituent that will form the bond with the skeleton.

Supported Index Formats
-----------------------

For ``link_index`` (and ``skeleton_indices``), the following formats are supported:

-  **Integer**: ``1``
-  **List of integers**: ``[1, 5]``
-  **String with comma-separated values**: ``"1, 5"``
-  **String with ranges**: ``"1-5, 8"`` (meaning 1, 2, 3, 4, 5, and 8)

****************
 Usage Examples
****************

Generate Template
=================

To quickly get started, you can generate a template configuration file:

.. code:: bash

   chemsmart run iterate -g my_config.toml

Run Generation (Merged Output)
==============================

To run the job and save all generated molecules into a single XYZ file (useful for visualization):

.. code:: bash

   chemsmart run iterate -f my_config.toml -o results

This will create ``results.xyz``.

Run Generation (Separate Outputs)
=================================

To generate a separate file for each valid combination (useful for subsequent distinct calculations):

.. code:: bash

   chemsmart run iterate -f my_config.toml --separate-outputs -d ./output_dir -np 4

This will: 1. Use 4 parallel processes to speed up generation. 2. Save each generated molecule as an individual file in
``./output_dir``.
