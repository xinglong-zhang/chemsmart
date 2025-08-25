Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to contact our team if you have questions or feedback.

Run Reaction Analysis Jobs
==========================

ChemSmart provides powerful reaction analysis capabilities using PyMOL for visualizing IRC trajectories, reaction pathways, and molecular dynamics during chemical reactions.

IRC Jobs
--------

Generate automatic PyMOL IRC movies and trajectory visualizations for reaction pathway analysis.

IRC Specific Options
^^^^^^^^^^^^^^^^^^^^

.. list-table:: IRC Job Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-r, --reactant <string>``
     - string
     - IRC file leading to the reactant side (default=None)
   * - ``-p, --product <string>``
     - string
     - IRC file leading to the product side (default=None)
   * - ``-a, --all <string>``
     - string
     - File containing all structures in the IRC, from full IRC run (default=None)

**Note**: IRC jobs also inherit all options from visualization jobs.

Basic Usage
^^^^^^^^^^^

* **Basic IRC trajectory visualization**:

    .. code-block:: console

        chemsmart run mol -f irc_output.log irc

* **IRC with separate reactant and product files**:

    .. code-block:: console

        chemsmart run mol -f ts_structure.xyz irc -r reactant.log -p product.log

* **IRC from complete trajectory file**:

    .. code-block:: console

        chemsmart run mol -f transition_state.xyz irc -a full_irc_trajectory.log

* **High-quality IRC movie with ray tracing**:

    .. code-block:: console

        chemsmart run mol -f irc_data.log irc -t

* **IRC visualization with custom style**:

    .. code-block:: console

        chemsmart run mol -f reaction_path.log irc -s cylview -v

Usage Examples
--------------

Complete Reaction Analysis Workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

?
