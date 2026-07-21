##############################
 Structure Optimization (xTB)
##############################

This page covers geometry optimization, single point, and Hessian/frequency workflows using xTB.

***********************
 Geometry Optimization
***********************

Find a local minimum structure with GFN-xTB / GFN-FF.

.. code:: bash

   chemsmart run [OPTIONS] xtb [XTB_OPTIONS] opt [SUBCMD_OPTIONS]

Optimization Options
====================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``--optimization-level``
      -  choice
      -  Convergence level: ``crude``, ``sloppy``, ``loose``, ``lax``, ``normal``, ``tight``, ``vtight``, ``extreme``

Basic Usage
===========

Standard geometry optimization:

.. code:: bash

   chemsmart run xtb -p project -f input.xyz -c 0 -m 1 opt

Override the project optimization level:

.. code:: bash

   chemsmart run xtb -p project -f molecule.xyz opt --optimization-level loose

HPC job submission:

.. code:: bash

   chemsmart sub -s server xtb -p project -f molecule.xyz opt

Solvated Optimization
=====================

Solvent options are specified at the xTB group level (before ``opt``). See :doc:`xtb-cli-options` for the full solvent
option reference.

.. code:: bash

   # ALPB(toluene)
   chemsmart run xtb -p project -f molecule.xyz -sm alpb -si toluene -a solv opt

   # Gas phase when the project YAML already sets a solvent
   chemsmart run xtb -p solv_project -f molecule.xyz --remove-solvent -a gas opt

Examples
========

**Optimization with GFN1 and gradients:**

.. code:: bash

   chemsmart run xtb -p project -f molecule.xyz -g gfn1 --grad opt

**Optimization from a Gaussian or ORCA output:**

.. code:: bash

   chemsmart sub -s server xtb -p project -f water_opt.out -c 0 -m 1 -l water_xtbopt opt

***************************
 Single Point Calculations
***************************

Compute energy (and optional properties) at a fixed geometry.

.. code:: bash

   chemsmart run [OPTIONS] xtb [XTB_OPTIONS] sp

Basic Usage
===========

.. code:: bash

   chemsmart run xtb -p project -f input.xyz sp

Single point with gradient output and extra xTB flags:

.. code:: bash

   chemsmart run xtb -p project -f molecule.xyz --grad -r "--json" sp

Solvated single point:

.. code:: bash

   chemsmart run xtb -p project -f molecule.xyz -sm alpb -si water sp

**********************************
 Hessian / Frequency Calculations
**********************************

Compute the numerical Hessian and vibrational frequencies on the input geometry.

.. code:: bash

   chemsmart run [OPTIONS] xtb [XTB_OPTIONS] hess

Basic Usage
===========

.. code:: bash

   chemsmart run xtb -p project -f molecule.xyz hess

Typical workflow: optimize first, then run frequencies on the optimized xTB output:

.. code:: bash

   chemsmart run xtb -p project -f molecule.xyz -l mol_opt opt
   chemsmart run xtb -p project -f mol_opt/mol_opt.out -l mol_hess hess

***************
 Output Layout
***************

Each xTB job runs in a dedicated folder named after the job label. CHEMSMART writes ``{label}.xyz`` and redirects
stdout/stderr to ``{label}.out`` / ``{label}.err``. Additional xTB artifacts (``charges``, ``xtbopt.xyz``, ``g98.out``,
``hessian``, etc.) appear in the same folder depending on the job type and flags.

The resulting files and folders can later be:

-  visualized or reused as geometry input (``-f path/to/{label}.out``)
-  analyzed with :doc:`thermochemistry-analysis` individually (``-f path/to/{label}.out``) or in batch (``-d
   path/to/parent_folder`` + ``-p xtb``)
-  assembled in batch into a database with :doc:`database-assemble` (``-d path/to/parent_folder`` + ``-p xtb``)

************
 Next Steps
************

-  :doc:`xtb-cli-options` — full xTB CLI option reference
-  :doc:`molecule-input-formats` — using xTB ``.out`` as Gaussian/ORCA geometry input
-  :doc:`thermochemistry-analysis`
-  :doc:`database-assemble`
