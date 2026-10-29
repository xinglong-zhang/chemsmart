###########################
 Reaction Analysis (PyMOL)
###########################

This page covers reaction pathway visualization using PyMOL, including IRC trajectories and molecular dynamics during
chemical reactions.

**********
 IRC Jobs
**********

Generate IRC movies and trajectory visualizations.

.. code:: bash

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] irc [SUBCMD_OPTIONS]

IRC Options
===========

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-r, --reactant``
      -  string
      -  IRC file for reactant side

   -  -  ``-p, --product``
      -  string
      -  IRC file for product side

   -  -  ``-a, --all``
      -  string
      -  File containing complete IRC trajectory

.. note::

   IRC jobs inherit all visualization options.

Basic Usage
===========

With separate reactant/product files:

.. code:: bash

   chemsmart run mol irc -r irc_reactant_side.log -p irc_product_side.log

From complete trajectory:

.. code:: bash

   chemsmart run mol irc -a full_irc_trajectory.log

For reactant side pathway movie:

.. code:: bash

   chemsmart run mol irc -r irc_reactant_side.log

For product side pathway movie:

.. code:: bash

   chemsmart run mol irc -p irc_product_side.log
