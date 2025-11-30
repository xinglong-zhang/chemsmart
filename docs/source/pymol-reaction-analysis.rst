############################
 Reaction Analysis (PyMOL)
############################

This page covers reaction pathway visualization using PyMOL, including IRC
trajectories and molecular dynamics during chemical reactions.

**********
 IRC Jobs
**********

Generate IRC movies and trajectory visualizations.

.. code-block:: bash

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] irc [SUBCMD_OPTIONS]

IRC Options
===========

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-r, --reactant``
     - string
     - IRC file for reactant side
   * - ``-p, --product``
     - string
     - IRC file for product side
   * - ``-a, --all``
     - string
     - File containing complete IRC trajectory

.. note::

   IRC jobs inherit all visualization options.

Basic Usage
===========

Standard IRC visualization:

.. code-block:: bash

   chemsmart run mol -f irc_output.log irc

With separate reactant/product files:

.. code-block:: bash

   chemsmart run mol -f ts_structure.xyz irc -r reactant.log

From complete trajectory:

.. code-block:: bash

   chemsmart run mol -f ts.xyz irc -a full_irc_trajectory.log
