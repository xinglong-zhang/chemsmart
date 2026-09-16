.. _chain-jobs:

#################
 Chain Workflows
#################

Chain projects tie together CREST, xTB, Gaussian, and ORCA calculations. A chain YAML file holds **project aliases**
only. Combined YAML is used only by ``chain`` jobs.

``chemsmart run/sub chain`` has two shapes:

-  A **custom pipeline** from repeatable ``-s/--steps`` (or the ``custom`` subcommand).
-  **Workflow subcommands** ``pka``, ``fukui``, ``redox``, and ``reaction``. These are not ``-s`` pipeline steps.

.. contents:: Table of Contents
   :local:
   :depth: 2

********************
 Chain Project YAML
********************

Chain project files live in ``~/.chemsmart/chain/``. The ``-p`` / ``--project`` option takes the file stem (without
``.yaml``), the same convention as ``gaussian -p test`` → ``~/.chemsmart/gaussian/test.yaml``.

Example (``~/.chemsmart/chain/combined.yaml``):

.. code:: yaml

   crest: crest_project1
   xtb: xtb_project1
   gaussian: gaussian_project2
   orca: orca_project3

Top-level keys ``crest``, ``xtb``, ``gaussian``, and ``orca`` are aliases to existing per-program project names. Any
other top-level key is rejected. A legacy ``steps`` key is also rejected — specify the pipeline with ``-s/--steps``.

Each pipeline step's program, and each workflow ``--program``, must have an alias in the file. The target per-program
YAML must exist. ``chemsmart sub gaussian -p test`` loads ``~/.chemsmart/gaussian/test.yaml`` only — it does not read
chain YAML.

Packaged templates are copied to ``~/.chemsmart/chain/`` by ``chemsmart config`` and ``chemsmart update projects``. See
:doc:`configuration-project-settings`.

********************
 CLI Pipeline Steps
********************

Repeat ``-s`` on the chain command for each phase. A step is ``PROGRAM:JOB`` or a quoted ``PROGRAM: [OPTIONS] JOB``
string that matches the program CLI after the program name:

.. code:: bash

   chemsmart run chain -p combined -f mol.xyz -c 0 -m 1 \
     -s "gaussian: -o maxstep=8,maxsize=12 opt" \
     -s gaussian:sp

   chemsmart sub chain -p combined -f mol.xyz -s gaussian:opt custom

The following ``(program, job)`` pairs are supported:

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   -  -  Program
      -  Job
      -  Description

   -  -  ``crest``
      -  ``conformers``
      -  CREST conformer search; best conformer is passed to the next step

   -  -  ``xtb``
      -  ``opt``
      -  xTB geometry optimization (with frequency by default)

   -  -  ``xtb``
      -  ``sp``
      -  xTB single-point energy

   -  -  ``xtb``
      -  ``hess``
      -  xTB Hessian / frequency calculation

   -  -  ``gaussian``
      -  ``opt``, ``ts``, ``sp``
      -  Gaussian optimization, transition state, or single point

   -  -  ``orca``
      -  ``opt``, ``ts``, ``sp``
      -  ORCA optimization, transition state, or single point

**pKa**, **Fukui**, **redox**, and **reaction** are **not** ``-s/--steps`` jobs. ``chain -s gaussian:pka`` is rejected.
Run them as :ref:`chain-workflow-subcommands`, or on the program CLI with that program's own project YAML.

Geometry handoff: optimization, TS, and SP use the previous job's optimized structure; CREST passes the best conformer
(``crest_best.xyz``). If a phase is incomplete, the chain stops so HPC resubmission with ``--skip-completed`` can
continue later.

Child labels follow ``{chain_label}_{index:02d}_{program}_{job}``. Charge and multiplicity from ``-c`` / ``-m`` are
merged into each child. Do not combine ``-s/--steps`` with ``pka``, ``fukui``, ``redox``, or ``reaction``.

The whole pipeline is **one** ``ChainJob`` submission. ``chemsmart sub`` writes one scheduler script for all phases.
``--fake`` uses the chain fake runner for the parent; each child gets the matching program fake runner when its type
differs.

.. _chain-workflow-subcommands:

**********************
 Workflow Subcommands
**********************

Submit requires ``-p``, ``-f``, and ``--program {gaussian,orca}``. Theory and solvent come from that program's alias in
the chain YAML. Analyze (pKa, Fukui, redox) does not need ``-p``, ``-f``, or ``--program``. Reaction is submit/batch
only.

.. code:: bash

   chemsmart sub chain -p combined -f mol.xyz -c 0 -m 1 \
     pka --program gaussian [pka options] [submit]
   chemsmart run chain pka analyze [pka analyze options]

Program CLI remains ``chemsmart sub gaussian|orca -p <program-project> … pka|fukui|redox``. Reaction is ``chemsmart sub
chain … reaction --program {gaussian,orca}`` only. See :doc:`pka-calculations`, :ref:`fukui-jobs`,
:doc:`redox-calculations`, and :doc:`reaction`.

***************
 CLI Reference
***************

See :doc:`chain-cli-options` for chain-specific options.
