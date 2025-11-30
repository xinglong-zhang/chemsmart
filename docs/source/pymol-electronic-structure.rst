########################################
 Electronic Structure Analysis (PyMOL)
########################################

This page covers electronic structure visualization using PyMOL, including
molecular orbitals and spin density plots.

*****************************
 Molecular Orbital (MO) Jobs
*****************************

Generate molecular orbital visualizations for frontier orbitals and other
electronic states.

.. code-block:: bash

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] mo [SUBCMD_OPTIONS]

MO Options
==========

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-n, --number``
     - int
     - Specific MO number to visualize
   * - ``-h, --homo``
     - bool
     - Plot HOMO (default: disabled)
   * - ``-l, --lumo``
     - bool
     - Plot LUMO (default: disabled)

.. note::

   MO jobs inherit all visualization options including styling, ray tracing,
   and surface rendering.

Basic Usage
===========

Standard MO visualization:

.. code-block:: bash

   chemsmart run mol -f calculation.log mo

HOMO visualization:

.. code-block:: bash

   chemsmart run mol -f molecule.log mo -h

LUMO visualization:

.. code-block:: bash

   chemsmart run mol -f molecule.log mo -l

Specific orbital:

.. code-block:: bash

   chemsmart run mol -f molecule.log mo -n 5

*******************
 Spin Density Jobs
*******************

Generate spin density visualizations for open-shell systems.

.. code-block:: bash

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] spin

.. note::

   Requires both ``.log`` and ``.chk`` files in the same folder.

Spin density jobs inherit all visualization options.

Basic Usage
===========

Standard spin density:

.. code-block:: bash

   chemsmart run mol -f radical.log spin

With ray tracing:

.. code-block:: bash

   chemsmart run mol -f radical.log spin -t
