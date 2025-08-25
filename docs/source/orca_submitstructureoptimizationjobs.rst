Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to contact our team if you have questions or feedback.

Submit Structure Optimization Jobs
==================================

ChemSmart provides comprehensive tools for structure optimization calculations using ORCA. This section covers geometry optimization workflows including constrained optimizations.

Structure Optimization
----------------------

Geometry optimization is used to find the minimum energy structure of a molecule by adjusting atomic positions until forces are minimized.

Run "chemsmart sub [GENERAL_OPTIONS] orca [ORCA_OPTIONS] opt [OPTIONS]" to perform structure optimization.

Optimization-Specific Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: Structure Optimization Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-f, --freeze-atoms <string>``
     - string
     - Indices of atoms to freeze for constrained optimization. 1-indexed (default=None)
   * - ``-i, --invert-constraints/--no-invert-constraints``
     - flag
     - Invert the constraints for frozen atoms in optimization (default=False)

Basic Usage
^^^^^^^^^^^

* **Basic geometry optimization**:

    .. code-block:: console

        chemsmart sub orca -p project_name -f input.xyz opt

* **Optimization with frozen atoms**:

    .. code-block:: console

        chemsmart sub orca -p constrained_opt -f molecule.xyz opt -f "1,2,3"

* **Optimization with inverted constraints**:

    .. code-block:: console

        chemsmart sub orca -p inverted_opt -f molecule.xyz opt -f "1,2,3" -i

* **Optimization with specific method and basis set**:

    .. code-block:: console

        chemsmart sub orca -p dft_opt -f molecule.xyz -x B3LYP -b def2-TZVP opt


Modred jobs
----------------------

Jobtype mission
Whats the Difference
