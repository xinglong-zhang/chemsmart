######################
 Structure Generation
######################

The ``iterate`` command builds molecular libraries by attaching substituents to selected sites on one or more skeleton
molecules. It enumerates the requested substitution patterns, optimizes each attached structure with RDKit ETKDGv3 or
Joint Lagrange Geometry Optimization (JLGO), and writes the successful structures as XYZ files.

For the complete command-line option reference, see :doc:`Iterate CLI Options <iterate-cli-options>`.

*************
 Quick Start
*************

Generate an annotated YAML template, edit its molecule paths and attachment sites, and then run it:

.. code:: bash

   # Generate iterate_template.yaml
   chemsmart run iterate yaml -g

   # Generate structures with the default ETKDG algorithm
   chemsmart run iterate yaml -f iterate_template.yaml -o results

   # Write one XYZ file per structure, using four worker processes
   chemsmart run iterate yaml -f iterate_template.yaml \
       --separate-outputs -d ./library -np 4

The command has an input-format layer (currently ``yaml``) and an optional algorithm subcommand (``etkdg`` or ``jlgo``):

.. code:: text

   chemsmart run iterate yaml -f CONFIG.yaml [YAML_OPTIONS] \
       [ALGORITHM [ALGORITHM_OPTIONS]]

If the algorithm subcommand is omitted, ``iterate`` uses the algorithm in the YAML ``algorithm`` block. If neither is
provided, ETKDG is used.

.. note::

   The ``cdxml`` input subcommand is reserved for future support and is not implemented yet.

********************
 Configuration File
********************

The ``-f`` file is a YAML manifest containing molecule paths and attachment rules. Relative molecule paths are resolved
relative to the YAML file, not the directory from which the command is run. For supported molecule file formats, see
:doc:`Molecule Input Formats <molecule-input-formats>`.

The following example allows ``Me`` and ``OH`` to be placed at either of two sites on one skeleton:

.. code:: yaml

   skeletons:
     - file_path: "benzene.xyz"
       label: "benzene"
       skeleton_indices: "1-6"
       slots:
         - group: 1
           link_indices: "1,3"

   substituents:
     - file_path: "methane.xyz"
       label: "Me"
       link_index: 1
       groups: [1]
     - file_path: "water.xyz"
       label: "OH"
       link_index: 1
       groups: [1]

   algorithm:
     name: jlgo

Each of the two sites has three choices: keep the original branch, attach ``Me``, or attach ``OH``. The all-original
choice is excluded, so this example generates ``3² - 1 = 8`` combinations: four single substitutions and four double
substitutions.

Skeleton entries
================

Each entry under ``skeletons`` accepts:

``file_path``
   Path to the skeleton molecule. This field is required.

``label``
   Identifier used in combination labels and separate-output filenames. Labels are optional, but explicit, unique labels
   are recommended. A label must contain only characters accepted for safe filenames.

``link_index`` or ``slots``
   The skeleton attachment sites. Every skeleton must define exactly one of these forms. Atom indices are 1-based.

``skeleton_indices``
   Optional 1-based indices identifying the skeleton core during branch preprocessing. They help distinguish core atoms
   from branches to remove at the selected attachment sites; they do not mean that all other atoms are unconditionally
   discarded. Every attachment atom must be included.

Substituent entries
===================

Each entry under ``substituents`` accepts:

``file_path``
   Path to the substituent molecule. This field is required.

``label``
   Identifier used in combination labels and separate-output filenames. The same filename-safety and uniqueness rules as
   skeleton labels apply.

``link_index``
   The single 1-based atom that bonds to the skeleton. A substituent needs one attachment atom to participate in
   generation.

``groups``
   One or more group numbers that the substituent may fill. This field is required for YAML input. A substituent can be
   made available to several groups, for example ``groups: [1, 2]``.

Index formats
=============

``link_index``, ``link_indices`` and ``skeleton_indices`` accept the following forms:

-  an integer: ``1``;
-  a list of integers: ``[1, 5]``;
-  a comma-separated string: ``"1,5"``;
-  a range string: ``"1-5,8"``.

All indices are 1-based. Duplicate attachment sites and sites shared by multiple slots are rejected during configuration
validation.

*********************************
 Attachment Groups and Expansion
*********************************

Groups determine which substituents are candidates for which skeleton sites. Group numbers are global across the entire
YAML file, contiguous, and start at 1:

-  a skeleton using ``link_index`` occupies one implicit group;
-  a skeleton using ``slots`` occupies one group per slot.

For example, if the first skeleton uses ``link_index`` and the second has two slots, their group numbers are 1 and then
2--3 respectively. Every substituent's ``groups`` values must refer to groups defined by the skeletons.

``link_index`` shorthand
========================

``link_index`` is convenient for a single-site screen, or for screening several sites independently:

.. code:: yaml

   skeletons:
     - file_path: "core.xyz"
       label: "core"
       link_index: "1,3"

   substituents:
     - file_path: "methyl.xyz"
       label: "Me"
       link_index: 1
       groups: [1]
     - file_path: "hydroxyl.xyz"
       label: "OH"
       link_index: 1
       groups: [1]

In the current implementation, each shorthand site is expanded separately. This example therefore generates four
single-site structures (``1Me``, ``1OH``, ``3Me`` and ``3OH``), not double-substituted structures. The
``--combination-mode`` option does not change shorthand expansion.

Use explicit ``slots`` whenever several sites must occur in the same Cartesian product.

Explicit slots
==============

A slot contains one group number and one or more skeleton positions:

.. code:: yaml

   skeletons:
     - file_path: "core.xyz"
       label: "core"
       slots:
         - group: 1
           link_indices: "1,3"
         - group: 2
           link_indices: "8"

All positions within one slot are expanded together. In the first slot above, positions 1 and 3 therefore participate in
the same Cartesian product. The second slot can use a different set of substituents through group 2.

Combination modes
=================

``-cm/--combination-mode`` controls whether *different explicit slots* are cross-combined. Each position always has an
implicit "keep the original branch" choice, and the all-original result is excluded.

Suppose slot ``R1`` has candidates ``A`` and ``B``, and slot ``R2`` has candidate ``C``:

``independent`` (default)
   Expand each slot separately and take the union. This gives three structures: ``R1=A``, ``R1=B`` and ``R2=C``. If a
   single slot contains several positions, those positions are still expanded together.

``global``
   Combine all slots in one Cartesian product. This gives five structures: the three single substitutions plus
   ``R1=A/R2=C`` and ``R1=B/R2=C``.

.. important::

   ``--combination-mode`` applies to explicit slots. It does not turn a multi-value ``link_index`` shorthand into
   multi-site combinations.

*********************
 Algorithm Selection
*********************

ETKDG (default)
===============

The ``etkdg`` algorithm uses RDKit ETKDGv3 distance-geometry embedding. In the default local mode, the processed
skeleton atoms are fixed while the attached substituent atoms are embedded. Global mode re-embeds the complete combined
molecule. ETKDG is used when neither the YAML file nor the CLI selects an algorithm explicitly.

Algorithm options can be declared in YAML:

.. code:: yaml

   algorithm:
     name: etkdg
     options:
       use_global_optimization: false
       num_conformers: 10
       random_seed: 42

They can also be supplied through an algorithm subcommand:

.. code:: bash

   chemsmart run iterate yaml -f config.yaml etkdg \
       --num-conformers 10 --random-seed 42

Explicit CLI algorithm options take precedence over YAML options. See :doc:`Iterate CLI Options <iterate-cli-options>`
for all names, defaults and validation rules.

JLGO
====

Joint Lagrange Geometry Optimization (``jlgo``) places all substituents in a combination through a joint 6K-dimensional
optimization, where K is the number of attached substituents. It keeps each linking atom on its bond sphere while
minimizing steric clashes. Single- and multi-substituent structures use the same joint optimization path.

Because ETKDG is the built-in default, select JLGO explicitly on the CLI:

.. code:: bash

   chemsmart run iterate yaml -f config.yaml jlgo

JLGO can also be selected in the YAML file:

.. code:: yaml

   algorithm:
     name: jlgo
     options:
       use_adaptive_sampling: true
       max_starts: 8

*************************
 Execution and Debugging
*************************

``-np/--nprocs`` sets the maximum number of worker processes. ``-t/--timeout`` sets the timeout in seconds for each
combination, not for the entire run.

In normal interactive use, the terminal shows a progress bar followed by the total, successful and failed combination
counts and the report/output paths. When stdout is redirected, ``iterate`` prints the run header but suppresses the
dynamic bar to avoid carriage-return output.

For development diagnostics, enable the existing top-level debug flag:

.. code:: bash

   chemsmart run --debug iterate yaml -f config.yaml

Debug mode disables the progress bar and exposes detailed Iterate worker logs and RDKit warnings. Failure details are
also recorded in the run report in normal mode.

********************
 Outputs and Report
********************

Structure output
================

By default, successful structures are written to one multi-structure XYZ file named by ``-o/--outputfile``
(``iterate_out.xyz`` by default). Each structure's XYZ comment line contains its combination label, for example
``benzene_1Me_3OH``.

With ``--separate-outputs``, each successful combination is written as ``<combination_label>.xyz`` in the directory
selected by ``-d/--directory``. Existing files with the same names may be overwritten.

The generated atom ordering consists of the processed skeleton followed by the attached substituent blocks. XYZ stores
elements, coordinates and the combination label; it does not encode molecular charge or multiplicity.

If all combinations fail, ``iterate`` does not create an empty merged XYZ file. In a partially successful completed run,
successful structures are kept and written even though the process exits with an error code.

Run report
==========

After YAML parsing and validation succeeds and the runner starts, ``iterate`` attempts to write a Gaussian-log-style
plain-text report named ``<config_stem>_iterate.out``. For example, ``123.yaml`` produces ``123_iterate.out``. The
report is placed beside the merged XYZ, in the separate-output directory, or in the current directory when no output
directory is specified.

The report records the input entries, resolved algorithm options, combination counts, per-combination results,
input/execution/timeout/write failures, output paths and final statistics. It ends with one of:

.. code:: text

   Normal termination of CHEMSMART Iterate at <time>.
   Error termination of CHEMSMART Iterate at <time>.

Configuration and usage errors (exit code 2) occur before the runner starts and do not produce a run report. Reports for
an unexpected internal error or SIGINT are best effort, and report writing can itself fail.

****************************
 Exit Codes and Error Codes
****************************

A run terminates normally only when every requested combination and output write succeeds. A partially successful run is
an error termination, but its successfully delivered structures are retained.

.. list-table::
   :header-rows: 1
   :widths: 15 85

   -  -  Exit code
      -  Meaning

   -  -  ``0``
      -  Normal termination with no errors.

   -  -  ``1``
      -  Runtime error termination, including input load/preprocessing, execution, timeout, structure write, internal,
         or report-write failure. Partial structure output may exist.

   -  -  ``2``
      -  CLI usage or YAML configuration error before molecule generation. No Iterate run report is created.

   -  -  ``130``
      -  Interrupted by the user (SIGINT). The interrupt report is best effort; partial in-memory results may not have
         been written.

The report can list several error codes when several failure types occur:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   -  -  Error code
      -  Meaning
   -  -  ``ITR-INPUT-001``
      -  Input molecule load error before combination execution.
   -  -  ``ITR-EXEC-001``
      -  Combination preprocessing or algorithm execution failure.
   -  -  ``ITR-TIMEOUT-001``
      -  Combination execution timeout.
   -  -  ``ITR-WRITE-001``
      -  Structure output write failure.
   -  -  ``ITR-INTERNAL-001``
      -  Unexpected internal error.
   -  -  ``ITR-INTERRUPTED-001``
      -  Run interrupted by the user.

**********
 Examples
**********

.. code:: bash

   # Merged output with the default ETKDG algorithm
   chemsmart run iterate yaml -f config.yaml -o results

   # Cross-combine explicit slots and write one file per structure
   chemsmart run iterate yaml -f config.yaml -cm global \
       --separate-outputs -d ./library -np 4

   # ETKDG local embedding with a fixed random seed
   chemsmart run iterate yaml -f config.yaml etkdg --random-seed 42

   # ETKDG global embedding
   chemsmart run iterate yaml -f config.yaml etkdg --global

   # Select JLGO instead of the default ETKDG algorithm
   chemsmart run iterate yaml -f config.yaml jlgo

For all common and algorithm-specific options, see :doc:`Iterate CLI Options <iterate-cli-options>`.
