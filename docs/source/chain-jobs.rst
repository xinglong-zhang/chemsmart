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

Each ``crest``, ``xtb``, ``gaussian``, and ``orca`` step, and each workflow ``--program``, must have an alias in the
file. A ``grouper`` step does not. The target per-program YAML must exist. ``chemsmart sub gaussian -p test`` loads
``~/.chemsmart/gaussian/test.yaml`` only — it does not read chain YAML.

Packaged templates are copied to ``~/.chemsmart/chain/`` by ``chemsmart config`` and ``chemsmart update projects``. See
:doc:`configuration-project-settings`.

********************
 CLI Pipeline Steps
********************

Repeat ``-s`` on the chain command for each phase. A step is ``PROGRAM:JOB`` or a quoted ``PROGRAM: [OPTIONS] JOB``
string. Options before the job name are that program's options. Options after the job name are that command's options:

.. code:: bash

   chemsmart run chain -p combined -f mol.xyz -c 0 -m 1 \
     -s "gaussian: -o maxstep=8,maxsize=12 opt" \
     -s gaussian:sp

   chemsmart sub chain -p combined -f mol.xyz -c 0 -m 1 \
     -s "crest:conformers -nc 10" -s gaussian:opt

   chemsmart sub chain -p combined -f mol.xyz -c 0 -m 1 \
     -s "crest:conformers -nc 10" \
     -s "gaussian: -l structure_from_lowest -c 0 -m 1 crest -j opt -nc 10" \
     -s gaussian:sp

   chemsmart sub chain -p combined -f ligand.xyz \
     -s "crest:conformers -nc 10" \
     -s "grouper: -T 0.3 irmsd --inversion on" \
     -s gaussian:opt

A step for ``crest``, ``xtb``, ``gaussian``, or ``orca`` is any subcommand of that program except ``pka``, ``fukui``,
``redox``, and ``reaction``. Those four stay on :ref:`chain-workflow-subcommands`. A ``grouper`` step is a grouper
strategy, such as ``irmsd``. Grouper has no alias in the chain YAML. ``gaussian:crest`` and ``grouper`` must immediately
follow ``crest:conformers``.

``-nc N`` after ``crest:conformers`` selects the lowest N conformers from ``crest_conformers.xyz`` for the next step and
is not passed to CREST. ``N`` is an integer of at least 1. Omitting ``-nc`` passes ``crest_best.xyz``. There is no
chain-level ``-N``.

**Geometry handoff.** Each phase reads the previous completed child. If a phase is incomplete, the chain stops so HPC
resubmission with ``--skip-completed`` can continue later.

-  An optimization, transition-state, or single-point step passes its optimized structure, or the output molecule when
   the job did not optimize.

-  A Gaussian, ORCA, or xTB step passes its own output file (``.log`` or ``.out``) as ``-f`` with ``-i -1``, for example
   ``gaussian:ts`` followed by ``gaussian:qrc``. The first step, CREST frames, and grouper representatives are written
   to ``.xyz`` files.

-  ``crest:conformers -nc N`` passes the lowest N frames of ``crest_conformers.xyz``. Frames are ordered by energy when
   every frame has an energy; otherwise the file order is kept. CREST writes that file lowest-energy first. Fewer frames
   than ``N`` uses every available frame. Omitting ``-nc`` passes the structure in ``crest_best.xyz``. A missing file
   leaves the next phase empty until the CREST step finishes.

-  ``gaussian:crest`` and ``grouper`` read that same selection. ``-nc N`` on the CREST step gives them
   ``crest_lowest_N.xyz`` written beside ``crest_conformers.xyz``. Omitting ``-nc`` gives them ``crest_best.xyz``.
   ``-nc`` on ``gaussian:crest`` is that command's own option and chooses how many of its input structures to optimize.

-  The step after Gaussian crest receives that job's lowest-energy optimized conformer.

-  The step after ``gaussian:irc`` receives the final structure of each directional log. Both directions hand those two
   endpoint jobs to the next step as a list, so that step builds one child per log. A forward-only or reverse-only IRC
   hands off that one log.

-  The step after grouper receives one representative, the first member of each group. Fewer than two CREST structures
   skips grouping and passes the available structure on.

Child labels start from ``{chain_label}_{index:02d}_{program}_{job}``. More than one group representative appends
``_g1``, ``_g2``, and so on. More than one other structure appends ``_c1``, ``_c2``, and so on. A single structure has
no suffix. The program command may extend the label it is given. Gaussian crest keeps the labels that command already
uses (``{label}_{jobtype}_c1``, …). Chain ``-c`` / ``-m`` go to the first program step. A later Gaussian, ORCA, or xTB
step reads charge and multiplicity from the previous output file. A step fed by an xyz file, including CREST frames and
grouper representatives, gets the charge carried with that structure. A step's own ``-c`` / ``-m`` wins for that step
and is written into its output, so the following step sees it. Grouper does not take those flags: on the grouper command
``-c`` is concentration and ``-m`` is the matrix file format. Do not combine ``-s/--steps`` with ``pka``, ``fukui``,
``redox``, or ``reaction``.

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
