Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to contact our team if you have questions or feedback.

Submit Structure Optimization Jobs
==================================

ChemSmart provides comprehensive tools for structure optimization calculations using ORCA. This section covers geometry optimization workflows including constrained optimizations.

Structure Optimization
----------------------

Geometry optimization is used to find the minimum energy structure of a molecule by adjusting atomic positions until forces are minimized.

.. code-block:: console

    chemsmart sub [OPTIONS] orca [ORCA_OPTIONS] opt [SUBCMD_OPTIONS]


Opt-Specific OPTIONS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: Structure Optimization Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-f, --freeze-atoms``
     - string
     - Indices of atoms to freeze for constrained optimization. 1-indexed (default=None)
   * - ``-i, --invert-constraints/--no-invert-constraints``
     - bool
     - Invert the constraints for frozen atoms in optimization (default=False)

Basic Usage
^^^^^^^^^^^

**Basic geometry optimization**:

    .. code-block:: console

        chemsmart sub orca -p project_name -f input.xyz opt

**Optimization with inverted constraints**:

    .. code-block:: console

        chemsmart sub orca -p inverted_opt -f molecule.xyz opt -f 1,2,3 -i


Single Point Calculations
-------------------------

Single point calculations compute the energy and properties of a molecule at a fixed geometry without optimization.

.. code-block:: console

    chemsmart sub [OPTIONS] orca [ORCA_OPTIONS] sp [SUBCMD_OPTIONS]

Basic Usage
^^^^^^^^^^^

**Basic single point calculation**:

    .. code-block:: console

        chemsmart sub orca -p project_name -f input.xyz sp


