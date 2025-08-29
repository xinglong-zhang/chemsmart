Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to
contact our team if you have questions or feedback.

############################
 Run Reaction Analysis Jobs
############################

ChemSmart provides powerful reaction analysis capabilities using PyMOL for visualizing IRC trajectories, reaction
pathways, and molecular dynamics during chemical reactions.

**********
 IRC Jobs
**********

Generate automatic PyMOL IRC movies and trajectory visualizations for reaction pathway analysis.

.. code:: console

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] irc [SUBCMD_OPTIONS]

IRC-Specific OPTIONS
====================

.. list-table:: IRC Job Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-r, --reactant``
      -  string
      -  IRC file leading to the reactant side (default=None)

   -  -  ``-p, --product``
      -  string
      -  IRC file leading to the product side (default=None)

   -  -  ``-a, --all``
      -  string
      -  File containing all structures in the IRC, from full IRC run (default=None)

.. note::

   IRC jobs also inherit all options from visualization jobs.

IRC Basic Usage
===============

**Basic IRC trajectory visualization**

   .. code:: console

      chemsmart run mol -f irc_output.log irc

**IRC with separate reactant and product files**

   .. code:: console

      chemsmart run mol -f ts_structure.xyz irc -r reactant.log

**IRC from complete trajectory file**

   .. code:: console

      chemsmart run mol -f transition_state.xyz irc -a full_irc_trajectory.log
