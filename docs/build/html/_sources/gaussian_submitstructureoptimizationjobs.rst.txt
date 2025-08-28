Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to
contact our team if you have questions or feedback.

####################################
 Submit Structure Optimization Jobs
####################################

ChemSmart provides comprehensive structure optimization capabilities using Gaussian. This section covers basic geometry
optimization, constrained optimization with frozen atoms, modified redundant coordinate optimization, and CREST
optimization workflows.

***********************
 Geometry Optimization
***********************

Run basic geometry optimization to find the minimum energy structure using the ``opt`` command.

.. code:: console

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] opt

Run constrained optimization with specific atoms frozen in place using the ``opt`` command with frozen atom options
(such as https://www.researchgate.net/post/Freezing-atoms-in-gaussian-how-to-do-it).

.. code:: console

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] opt [SUBCMD_OPTIONS]

Opt-Specific OPTIONS
====================

.. list-table:: Frozen Atom Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-f, --freeze-atoms``
      -  string
      -  Indices of atoms to freeze for constrained optimization (1-based indexing)

Basic Usage
===========

**Basic geometry optimization**

   .. code:: console

      chemsmart sub gaussian -p project_name -f input.gjf -c 0 -m 1 opt

**Constrained optimization with frozen atoms**

-  To submit the geometry optimization job with atoms frozen, one can do

   .. code:: console

      chemsmart sub -s shared gaussian -p frozen_opt -f input.com opt -f 1-10

   and

   .. code:: console

      chemsmart sub -s shared gaussian -p frozen_opt -f input.com opt -f 1-3,5,7

Examples
========

**optimize structure directly from a Gaussian optimization output file with different charge and multiplicity**

-  After a Gaussian optimization job is done, one can use the output file from the previous optimization job to run a
   new optimization with different charge and multiplicity, e.g.:

   .. code:: console

      chemsmart sub -s SLURM gaussian -p project1 -f k_atom_opt.log -c 1 -m 1 -l k_cation_opt opt

   ``-c`` and ``-m`` will override the charge and multiplicity in the original Gaussian log file (from K_atom to
   K_cation).

*************************
 CREST Optimization Jobs
*************************

Run CREST-based optimization on multiple conformers using the ``crestopt`` command.

.. code:: console

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] crestopt [SUBCMD_OPTIONS]

CRESTOPT-Specific OPTIONS
=========================

.. list-table:: CREST Optimization Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-n, --num-confs-to-opt``
      -  int
      -  Number of conformers to optimize from the CREST ensemble

Basic Usage
===========

**CREST Optimization for all conformers**

   .. code:: console

      chemsmart sub gaussian -p crest_optimization -f conformers.xyz crestopt

**with specific number of conformers**

   .. code:: console

      chemsmart sub gaussian -p crest_opt -f molecule.xyz crestopt -n 10
