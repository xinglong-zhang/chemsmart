Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to contact our team if you have questions or feedback.

Run Basic Visualization Jobs
=============================

ChemSmart provides powerful molecular visualization capabilities using PyMOL for creating high-quality molecular graphics, movies, and interactive visualizations.

Visualization Jobs
------------------

Create static PyMOL visualizations and interactive session files.

Visualization Specific Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: Visualization Job Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-f, --file <string>``
     - string
     - PyMOL file script or style. If not specified, defaults to use zhang_group_pymol_style.py (default=None)
   * - ``-s, --style``
     - string
     - PyMOL render style. Options: pymol, cylview (default=None)
   * - ``-t, --trace/--no-trace``
     - bool
     - PyMOL options to ray trace or not (default=True)
   * - ``-v, --vdw``
     - bool
     - Add Van der Waals surface (default=False)
   * - ``-q, --quiet/--no-quiet``
     - bool
     - Run PyMOL in quiet mode (default=False)
   * - ``--command-line-only/--no-command-line-only``
     - bool
     - Run PyMOL in command line only (default=True)
   * - ``-c, --coordinates <string>``
     - string
     - List of coordinates (bonds, angles and dihedrals) for labelling. 1-indexed (default=None)


Basic Usage
^^^^^^^^^^^

* **Basic molecular visualization**:

    .. code-block:: console

        chemsmart run mol -f molecule.xyz visualize

* **Quiet mode visualization**:

    .. code-block:: console

        chemsmart run mol -f calculation.log visualize -q

* **Visualization with coordinate labeling**:

    .. code-block:: console

        chemsmart run mol -f structure.xyz visualize -c "1,2,3"

* **Command line only mode**:

    .. code-block:: console

        chemsmart run mol -f molecule.log visualize --command-line-only

Movie Jobs
----------

Generate rotating movie animations of molecular structures.

Movie Specific Options
^^^^^^^^^^^^^^^^^^^^^^

Movie jobs inherit all options from visualization jobs and use the same parameters.

tips: Specific OPTIONS are inherited from visualization jobs

Basic Usage
^^^^^^^^^^^

* **Basic rotating movie**:

    .. code-block:: console

        chemsmart run mol -f molecule.xyz movie

* **Quiet mode movie generation**:

    .. code-block:: console

        chemsmart run mol -f optimized.log movie -q

* **Movie with coordinate labels**:

    .. code-block:: console

        chemsmart run mol -f structure.log movie -c "1,2,5,10"

