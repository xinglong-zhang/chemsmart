Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to contact our team if you have questions or feedback.

Run Electronic Structure Analysis Jobs
======================================

ChemSmart provides powerful electronic structure visualization capabilities using PyMOL for generating molecular orbitals, spin density plots, and other electronic property visualizations.

Molecular Orbital (MO) Jobs
---------------------------

Generate molecular orbital visualizations for frontier orbitals and other electronic states.

.. code-block:: console

    chemsmart run [OPTIONS] mol [MOL_OPTIONS] mo [SUBCMD_OPTIONS]

MO-Specific OPTIONS
^^^^^^^^^^^^^^^^^^^

.. list-table:: MO Job Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-n, --number``
     - int
     - Molecular Orbital number to be visualized (e.g., 31 will visualize MO #31) (default=None)
   * - ``-h, --homo``
     - bool
     - Plot the highest occupied molecular orbital (HOMO) (default=False)
   * - ``-l, --lumo``
     - bool
     - Plot the lowest unoccupied molecular orbital (LUMO) (default=False)

.. note::

    MO jobs also inherit all options from visualization jobs including styling, ray tracing, and surface rendering options.

Basic Usage
^^^^^^^^^^^

**Basic molecular orbital visualization**

    .. code-block:: console

        chemsmart run mol -f calculation.log mo

**HOMO visualization**

    .. code-block:: console

        chemsmart run mol -f molecule.log mo -h

**LUMO visualization**

    .. code-block:: console

        chemsmart run mol -f molecule.log mo -l

**Specific orbital visualization**

    .. code-block:: console

        chemsmart run mol -f molecule.log mo -n 5


Spin Density Jobs
-----------------

Generate spin density visualizations for open-shell systems and radical species.

.. code-block:: console

    chemsmart run [OPTIONS] mol [MOL_OPTIONS] mo [SUBCMD_OPTIONS]

.. note::

    2 files are needed in one folder (log and chk files)

Spin-Specific OPTIONS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Spin density jobs inherit all options from visualization jobs including styling, ray tracing, and surface rendering options.

Basic Usage
^^^^^^^^^^^

* **Basic spin density visualization**

    .. code-block:: console

        chemsmart run mol -f radical_calculation.log spin

* **High-quality spin density with ray tracing**

    .. code-block:: console

        chemsmart run mol -f radical.log spin -t



