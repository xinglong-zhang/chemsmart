###################################
 Structure Optimization (Gaussian)
###################################

This page covers geometry optimization workflows using Gaussian.

***********************
 Geometry Optimization
***********************

Run basic geometry optimization to find the minimum energy structure:

.. code:: bash

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] opt [SUBCMD_OPTIONS]

Optimization Options
====================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-f, --freeze-atoms``
      -  string
      -  Atom indices to freeze (1-based indexing)

Basic Usage
===========

Standard geometry optimization:

.. code:: bash

   chemsmart sub gaussian -p project -f input.gjf -c 0 -m 1 opt

Constrained optimization with frozen atoms:

.. code:: bash

   # Freeze atoms 1-10
   chemsmart sub gaussian -p project -f input.com opt -f 1-10

   # Freeze atoms 1, 2, 3, 5, and 7
   chemsmart sub gaussian -p project -f input.com opt -f 1-3,5,7

Examples
========

**Optimization from a previous output file:**

Use an optimized structure with different charge/multiplicity:

.. code:: bash

   chemsmart sub -s SLURM gaussian -p project -f k_atom_opt.log -c 1 -m 1 -l k_cation_opt opt

This takes the optimized K atom structure and runs a cation optimization.
