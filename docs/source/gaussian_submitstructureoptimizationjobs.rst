Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to contact our team if you have questions or feedback.

Submit Structure Optimization Jobs
===================================

ChemSmart provides comprehensive structure optimization capabilities using Gaussian. This section covers basic geometry optimization, constrained optimization with frozen atoms, and CREST optimization workflows.

Basic Geometry Optimization
----------------------------

Run basic geometry optimization to find the minimum energy structure.

Basic Usage
^^^^^^^^^^^

* **Basic geometry optimization**:

    .. code-block:: console

        chemsmart sub gaussian -p project_name -f input.xyz opt

* **Example with specific charge and multiplicity**:

    .. code-block:: console

        chemsmart sub gaussian -p optimization -f molecule.xyz -c 0 -m 1 opt

Frozen Atom Optimization
-------------------------

Run constrained optimization with specific atoms frozen in place.

Frozen Atom Options
^^^^^^^^^^^^^^^^^^^

.. list-table:: Frozen Atom Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-f, --freeze-atoms <string>``
     - string
     - Indices of atoms to freeze for constrained optimization

Basic Usage
^^^^^^^^^^^

* **Optimization with frozen atoms**:

    .. code-block:: console

        chemsmart sub gaussian -p constrained_opt -f molecule.xyz opt -f "1 2 3"

* **Example freezing specific atom indices**:

    .. code-block:: console

        chemsmart sub gaussian -p frozen_opt -f complex.xyz opt -f "5-10 15 20"

CREST Optimization Jobs
-----------------------

Run CREST-based optimization on multiple conformers.

CREST Optimization Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: CREST Optimization Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-n, --num-confs-to-opt <int>``
     - int
     - Number of conformers to optimize

Basic Usage
^^^^^^^^^^^

* **Basic CREST optimization**:

    .. code-block:: console

        chemsmart sub gaussian -p crest_optimization -f conformers.xyz crestopt

* **CREST optimization with specific number of conformers**:

    .. code-block:: console

        chemsmart sub gaussian -p crest_opt -f molecule.xyz crestopt -n 10

Modified Redundant Coordinate (Modred) Jobs
--------------------------------------------

Run modified redundant coordinate optimization for constrained geometry optimization and transition state searches.

Modred-Specific Options
^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: Modred Job Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-c, --coordinates <string>``
     - string
     - List of coordinates to be fixed for modred job. 1-indexed

Basic Usage
^^^^^^^^^^^


