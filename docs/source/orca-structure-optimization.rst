####################################
 Structure Optimization (ORCA)
####################################

This page covers geometry optimization and single point calculations using ORCA.

************************
 Geometry Optimization
************************

Find the minimum energy structure by adjusting atomic positions.

.. code-block:: bash

   chemsmart sub [OPTIONS] orca [ORCA_OPTIONS] opt [SUBCMD_OPTIONS]

Optimization Options
====================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-f, --freeze-atoms``
     - string
     - Atom indices to freeze (1-based indexing)
   * - ``-i, --invert-constraints/--no-invert-constraints``
     - bool
     - Invert frozen atom constraints (default: disabled)

Basic Usage
===========

Standard geometry optimization:

.. code-block:: bash

   chemsmart sub orca -p project -f input.xyz opt

Optimization with frozen atoms:

.. code-block:: bash

   chemsmart sub orca -p project -f molecule.xyz opt -f 1,2,3

Optimization with inverted constraints:

.. code-block:: bash

   chemsmart sub orca -p project -f molecule.xyz opt -f 1,2,3 -i

***************************
 Single Point Calculations
***************************

Compute energy and properties at a fixed geometry.

.. code-block:: bash

   chemsmart sub [OPTIONS] orca [ORCA_OPTIONS] sp

Basic Usage
===========

Standard single point calculation:

.. code-block:: bash

   chemsmart sub orca -p project -f input.xyz sp
