###################
 Chain CLI Options
###################

This page documents the CLI options for the ``chain`` command group. Use ``chemsmart sub chain --help`` for the complete
list, including ``custom`` and the workflow subcommands ``pka``, ``fukui``, ``redox``, and ``reaction``.

*************************
 Basic Command Structure
*************************

**Custom pipeline** (no nested subcommand, or ``custom``):

.. code:: bash

   chemsmart sub [OPTIONS] chain [CHAIN_OPTIONS] -f FILE -s STEP [-s ...] [custom]

**Workflow submit** (``pka``, ``fukui``, ``redox``, ``reaction``):

.. code:: bash

   chemsmart sub [OPTIONS] chain -p PROJECT -f FILE -c CHARGE -m MULT \
     pka|fukui|redox|reaction --program {gaussian,orca} [WORKFLOW_OPTIONS]

**Workflow analyze** (pKa, Fukui, and redox only; no ``-p``, ``-f``, or ``--program``). Redox infers ``-r`` from
Ref_ox/Ref_red formulas; pass ``-r`` to override:

.. code:: bash

   chemsmart run chain pka|fukui analyze [ANALYZE_OPTIONS]
   chemsmart run chain redox analyze [ANALYZE_OPTIONS]

***************
 Chain Options
***************

Project and File Options
========================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-p, --project``
      -  string
      -  Chain project settings from ``~/.chemsmart/chain/*.yaml``. Required for the ``-s/--steps`` pipeline and for
         ``pka`` / ``fukui`` / ``redox`` / ``reaction`` **submit**. Not required for analyze.

   -  -  ``-f, --filename``
      -  string
      -  Input structure file (required for the CLI pipeline and for workflow submit)

   -  -  ``-s, --steps``

      -  string (repeatable)

      -  One pipeline step as ``PROGRAM:JOB`` or quoted ``PROGRAM: [OPTIONS] JOB`` (e.g. ``gaussian:opt`` or
         ``"gaussian: -o maxstep=8,maxsize=12 opt"``). Required for the CLI pipeline; repeat for each phase. Not the
         same as ``run``/``sub`` ``-s/--server``. Cannot be combined with ``pka``, ``fukui``, ``redox``, or
         ``reaction``.

   -  -  ``-l, --label``
      -  string
      -  Custom job label (without extension)

   -  -  ``-a, --append-label``
      -  string
      -  String appended to the base filename for the label

   -  -  ``-i, --index``
      -  string
      -  Structure index (1-based; default: last structure)

   -  -  ``-c, --charge``
      -  int
      -  Molecular charge (pipeline and workflow submit; merged into each child job)

   -  -  ``-m, --multiplicity``
      -  int
      -  Spin multiplicity (pipeline and workflow submit; merged into each child job)

Workflow Submit Options
=======================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``--program``
      -  ``gaussian`` or ``orca``
      -  QC program for workflow **submit**. Required. Theory and solvent come from that program's alias in the chain
         YAML. Not used for analyze.

.. note::

   -  ``-p`` uses the chain project name without the ``.yaml`` extension.

   -  For the CLI pipeline, ``-f`` and at least one ``-s/--steps`` are required. The ``custom`` subcommand runs the same
      pipeline.

   -  For workflow submit, ``-p``, ``-f``, and ``--program`` are required. Analyze does not need ``-p``, ``-f``, or
      ``--program``.

   -  ``-l`` and ``-a`` are mutually exclusive.

   -  Workflow option tables (proton index, Fukui mode, redox reference files) match the corresponding program CLI.
      Reaction options (``--product``, …) are on ``chain … reaction`` only. See :doc:`pka-calculations`,
      :ref:`fukui-jobs`, :doc:`redox-calculations`, and :doc:`reaction`.

Execution Control
=================

The chain group accepts the standard job options documented in :doc:`cli-overview`:

-  ``-S/-R, --skip-completed/--no-skip-completed`` — applies to the ``-s/--steps`` pipeline and to workflow **submit**
   (pKa, Fukui, redox, reaction). ``-R`` on the chain group is enough; repeating it after the workflow name also works.
   Analyze is unchanged.

-  ``--fake/--no-fake``

-  ``--scratch/--no-scratch``

**********
 Examples
**********

Pipeline (optional ``custom``):

.. code:: bash

   chemsmart sub -s mycluster chain -p combined -f ligand.xyz -c 0 -m 1 \
     -s crest:conformers -s xtb:opt -s gaussian:opt -s orca:sp
   chemsmart run chain -p combined -f mol.xyz -c 0 -m 1 \
     -s "gaussian: -o maxstep=8,maxsize=12 opt" \
     -s gaussian:sp custom

Workflow submit:

.. code:: bash

   chemsmart sub chain -p combined -f ligand.xyz -c 0 -m 1 \
     pka --program gaussian -pi 10 submit

Program CLI uses that program's own project YAML, not chain YAML:

.. code:: bash

   chemsmart sub gaussian -p gaussian_project2 -f ligand.xyz -c 0 -m 1 pka submit

************
 Next Steps
************

-  :doc:`chain-jobs` — chain YAML aliases, ``-s/--steps`` pipeline, and workflow subcommands
-  :doc:`configuration-project-settings` — chain project file location and alias rules
-  :doc:`pka-calculations` — pKa submit and analyze
-  :doc:`redox-calculations` — exchange redox submit and analyze
-  :doc:`reaction` — reaction submit (no analyze)
