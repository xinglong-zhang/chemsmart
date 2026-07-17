#####################
 Iterate CLI Options
#####################

This page documents the CLI options available for the ``iterate`` command.

*************************
 Basic Command Structure
*************************

The ``iterate`` command is organized into two independent layers:

#. **Input format layer** — how the skeletons/substituents are provided.
   This is selected with a subcommand of ``iterate`` (currently ``yaml``;
   ``cdxml`` is coming soon).

#. **Algorithm layer** — how substituent positions are optimized. This is
   selected with an *optional* algorithm subcommand of the input command
   (``lagrange``, ``etkdg``).

.. code:: bash

   chemsmart run iterate yaml -f <CONFIG_FILE> [YAML_OPTIONS] \
       [ALGORITHM [ALGORITHM_OPTIONS]]

If no algorithm subcommand is given, the algorithm declared in the YAML
``algorithm`` block is used, falling back to the built-in default
(``lagrange_multipliers``).

********************************
 Role of the Configuration File
********************************

The configuration file (specified by ``-f``) acts as a **manifest** or
**catalog** rather than a standard structure input. Instead of containing raw
molecular data, it lists the **file paths** for skeleton and substituent
molecules, along with the **link indices** that define their connection
points. Optionally, it may also declare the optimization ``algorithm`` and
its options.

You can generate a template configuration file using the ``-g`` flag:

.. code:: bash

   chemsmart run iterate yaml -g my_template.yaml

For a detailed explanation of each section within the configuration file,
please refer to the :doc:`Structure Generation <iterate-structure-generation>`
page.

***********************
 Common (yaml) Options
***********************

These options belong to the ``yaml`` command and are **algorithm-agnostic**.

.. list-table::
   :header-rows: 1
   :widths: 35 15 50

   -  -  Option
      -  Type
      -  Description

   -  -  ``-f, --filename``
      -  string
      -  Path to the YAML configuration file (required)

   -  -  ``-np, --nprocs``
      -  int
      -  Number of parallel processes to use (default: 1)

   -  -  ``-t, --timeout``
      -  int
      -  Timeout in seconds for each molecule to be generated (default: 120)

   -  -  ``-cm, --combination-mode``
      -  choice
      -  Combination strategy for skeleton slots: ``independent`` (default)
         or ``global``

   -  -  ``-g, --generate-template``
      -  string
      -  Generate a template configuration file (in YAML format) and exit.
         Defaults to ``iterate_template.yaml`` if not specified.

Output Options
==============

.. list-table::
   :header-rows: 1
   :widths: 35 15 50

   -  -  Option
      -  Type
      -  Description

   -  -  ``--separate-outputs`` / ``--no-separate-outputs``
      -  flag
      -  Whether to save each structure as a separate XYZ file or merge into
         one file (default: ``--no-separate-outputs``)

   -  -  ``-o, --outputfile``
      -  string
      -  Output filename for merged output (default: ``iterate_out``). Only
         valid with ``--no-separate-outputs``.

   -  -  ``-d, --directory``
      -  string
      -  Directory to save separate output files. Only valid with
         ``--separate-outputs``.

.. note::

   The legacy options ``-m/--method``,
   ``-s/--sphere-direction-samples-number`` and
   ``-a/--axial-rotations-sample-number`` have been **removed**. Algorithm
   parameters now live exclusively on the algorithm subcommands (see below).

***********************
 Algorithm Subcommands
***********************

Algorithm parameters are isolated per algorithm; each algorithm exposes only
its own options.

``lagrange`` (Joint Lagrange)
=============================

Attaches one or more substituents in a single joint (6K-dimensional)
optimization; the same path handles single- and multi-substituent
combinations.

.. list-table::
   :header-rows: 1
   :widths: 40 12 48

   -  -  Option
      -  Type
      -  Description

   -  -  ``--adaptive-sampling`` / ``--no-adaptive-sampling``
      -  flag
      -  Run a fixed coarse sampling stage first (default:
         ``--adaptive-sampling``).

   -  -  ``--link-sphere-samples``
      -  int
      -  Full-stage number of linking-atom bond-sphere position samples
         (default: 48)

   -  -  ``--orientation-sphere-samples``
      -  int
      -  Full-stage number of substituent principal-axis direction samples
         (default: 24)

   -  -  ``--axial-samples``
      -  int
      -  Number of axial rotations per orientation direction (default: 4)

   -  -  ``--candidate-pool-size``
      -  int
      -  Per-substituent candidate pool size kept after region exclusion
         (default: 20)

   -  -  ``--preselect``
      -  int
      -  Top joint combinations fed into greedy start selection (default: 48)

   -  -  ``--beam-width``
      -  int
      -  Beam width retained per layer during feasible-domain pruning
         (default: 4096)

   -  -  ``--max-starts``
      -  int
      -  Maximum number of 6K-dimensional joint starts handed to SLSQP
         (default: 8)

   -  -  ``--slsqp-maxiter``
      -  int
      -  Maximum SLSQP iterations per start (default: 200)

.. note::

   With adaptive sampling enabled (the default), a fixed coarse stage runs
   first and the six full-stage sampling/pruning options
   (``--link-sphere-samples``, ``--orientation-sphere-samples``,
   ``--axial-samples``, ``--candidate-pool-size``, ``--preselect``,
   ``--beam-width``) only take effect when the coarse stage does not yield
   an acceptable optimized structure; ``--max-starts`` and
   ``--slsqp-maxiter`` always apply. Use ``--no-adaptive-sampling`` to always
   apply the full sampling parameters. Increasing the sampling numbers
   increases the time required for molecule generation.

``etkdg`` (RDKit ETKDGv3)
=========================

Generates the attached-substituent geometry with RDKit's ETKDGv3 distance
geometry. In the default **local** mode the skeleton atoms are held fixed and
only the substituent is re-embedded; ``--global`` re-embeds every atom.

.. list-table::
   :header-rows: 1
   :widths: 40 12 48

   -  -  Option
      -  Type
      -  Description

   -  -  ``--global`` / ``--local``
      -  flag
      -  Embedding mode. ``--local`` (default) keeps the skeleton fixed;
         ``--global`` re-embeds the whole molecule.

   -  -  ``--num-conformers``
      -  int
      -  Number of ETKDG conformers to try per attachment; the lowest-energy
         one is kept (default: 10)

   -  -  ``--random-seed``
      -  int
      -  Base RDKit random seed; ``-1`` for a non-reproducible seed
         (default: 42)

   -  -  ``--max-iterations``
      -  int
      -  Maximum ETKDG embedding iterations, ``0`` uses the RDKit default
         (default: 2000)

   -  -  ``--random-coords`` / ``--no-random-coords``
      -  flag
      -  Start embedding from random coordinates (default: ``--random-coords``)

   -  -  ``--enforce-chirality`` / ``--no-enforce-chirality``
      -  flag
      -  Enforce the input chirality during embedding
         (default: ``--no-enforce-chirality``)

   -  -  ``--force-field``
      -  choice
      -  Optional force-field post-optimization: ``none`` (default), ``uff``,
         ``mmff94`` or ``mmff94s``

.. note::

   Every ETKDG option is available both as a CLI flag (above) and as a key in
   the YAML ``algorithm`` block, using the same names (the YAML keys are
   ``use_global_optimization``, ``num_conformers``, ``random_seed``,
   ``max_iterations``, ``use_random_coordinates``, ``enforce_chirality`` and
   ``force_field``).

********************************
 Algorithm Selection & Priority
********************************

The effective algorithm is resolved from three sources, in increasing order
of priority:

#. **Built-in default** — ``lagrange_multipliers`` with its default options.
#. **YAML ``algorithm`` block** — overrides the default name and/or options.
#. **CLI algorithm subcommand** — the subcommand name overrides the YAML
   algorithm name; explicitly-passed options override the matching YAML
   options.

Additional rules:

-  CLI options that are **not** explicitly passed do **not** override YAML
   values (Click defaults never clobber the YAML configuration).
-  If the CLI selects a **different** algorithm than the YAML block, the YAML
   options (which belong to the other algorithm) are **discarded**, not mixed
   into the newly selected algorithm.

YAML ``algorithm`` block (ETKDG example):

.. code:: yaml

   algorithm:
     name: etkdg
     options:
       use_global_optimization: false   # true re-embeds every atom
       num_conformers: 10
       random_seed: 42
       max_iterations: 2000
       force_field: none                # or uff / mmff94 / mmff94s

   skeletons:
     ...

   substituents:
     ...

**********
 Examples
**********

.. code:: bash

   # Use the YAML algorithm block (or the default if none is declared)
   chemsmart run iterate yaml -f config.yaml

   # Explicitly select Joint Lagrange (equivalent to the default)
   chemsmart run iterate yaml -f config.yaml lagrange

   # Override Lagrange sampling parameters from the CLI
   chemsmart run iterate yaml -f config.yaml lagrange \
       --no-adaptive-sampling \
       --link-sphere-samples 48 \
       --axial-samples 4

   # ETKDG, local mode (skeleton fixed) is the default
   chemsmart run iterate yaml -f config.yaml etkdg

   # ETKDG, global mode (re-embed the whole molecule)
   chemsmart run iterate yaml -f config.yaml etkdg --global

   # ETKDG with more conformers and a fixed seed
   chemsmart run iterate yaml -f config.yaml etkdg \
       --num-conformers 50 --random-seed 1

   # Generate a template YAML configuration file
   chemsmart run iterate yaml -g my_config.yaml
