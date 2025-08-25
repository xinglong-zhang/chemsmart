Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to contact our team if you have questions or feedback.

Submit Transition State Search Jobs
===================================

ChemSmart provides comprehensive transition state search capabilities using Gaussian. This section covers transition state optimization, IRC calculations, and potential energy surface scanning.

Transition State Optimization
-----------------------------

Run transition state optimization to find saddle points on the potential energy surface.

Basic Usage
^^^^^^^^^^^

* **Basic transition state optimization**:

    .. code-block:: console

        chemsmart sub gaussian -p ts_search -f guess_structure.xyz ts

* **Transition state with specific charge and multiplicity**:

    .. code-block:: console

        chemsmart sub gaussian -p ts_opt -f ts_guess.xyz -c 0 -m 1 ts

Transition State Search with Frozen Atoms
------------------------------------------

Run transition state optimization with specific atoms frozen to constrain the search.

Frozen Atom Options for TS
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: TS Frozen Atom Options
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

* **TS optimization with frozen atoms**:

    .. code-block:: console

        chemsmart sub gaussian -p constrained_ts -f ts_guess.xyz ts -f "1 2 3"

* **Example freezing specific atom indices during TS search**:

    .. code-block:: console

        chemsmart sub gaussian -p frozen_ts -f complex.xyz ts -f "5-10 15 20"

modred ts
---------


Intrinsic Reaction Coordinate (IRC) Calculations
------------------------------------------------

Run IRC calculations to follow the reaction path from a transition state.

IRC-Specific Options
^^^^^^^^^^^^^^^^^^^^

.. list-table:: IRC Job Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-fl, --flat-irc/--no-flat-irc``
     - bool
     - Whether to run flat IRC or not (default=False)
   * - ``-pt, --predictor``
     - string
     - Type of predictors used for IRC. Options: LQA, HPC, EulerPC, DVV, Euler (default=none)
   * - ``-rc, --recorrect``
     - string
     - Recorrection step of HPC and EulerPC IRCs. Options: Never, Always, Test (default=none)
   * - ``-rs, --recalc-step``
     - int
     - Compute the Hessian analytically every N predictor steps or every |N| corrector steps if N<0 (default=6)
   * - ``-p, --maxpoints``
     - int
     - Number of points along reaction path to examine (default=512)
   * - ``-c, --maxcycles``
     - int
     - Maximum number of steps along IRC to run (default=128)
   * - ``-s, --stepsize``
     - int
     - Step size along reaction path, in units of 0.01 Bohr (default=20)

Basic Usage
^^^^^^^^^^^


Potential Energy Surface Scanning
----------------------------------

Run coordinate scanning to explore potential energy surfaces and locate transition states.

Scanning coordinates, step size and number of steps of scan required!

Scan-Specific Options
^^^^^^^^^^^^^^^^^^^^^

.. list-table:: Scan Job Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-c, --coordinates``
     - string
     - List of coordinates to be fixed for scan job. 1-indexed (default=None)
   * - ``-s, --step-size``
     - float
     - Step size of coordinates to scan (default=None)
   * - ``-n, --num-steps``
     - int
     - Number of steps to scan (default=None)

Basic Usage
^^^^^^^^^^^

* **Basic coordinate scan**:

    .. code-block:: console

        chemsmart sub gaussian -p pes_scan -f molecule.xyz scan -c "[[2,3]]" -s 0.1 -n 15

* **Multi-coordinate scan**:

    .. code-block:: console

        chemsmart sub gaussian -p multi_scan -f complex.xyz scan -c "[[2,3],[6,7]]" -s 0.1 -n 10

* **Bond distance scan example**:

    .. code-block:: console

        chemsmart sub gaussian -p bond_scan -f reactant.xyz scan -c "[[1,2]]" -s 0.05 -n 20
