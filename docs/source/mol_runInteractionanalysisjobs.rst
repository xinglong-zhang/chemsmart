Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to contact our team if you have questions or feedback.

Run Interaction Analysis Jobs
=============================

ChemSmart provides powerful interaction analysis capabilities using PyMOL for visualizing non-covalent interactions (NCI) and intermolecular interaction patterns.

NCI Jobs
--------

Generate non-covalent interaction (NCI) plots for analyzing intermolecular and intramolecular interactions.

.. code-block:: console

    chemsmart run [OPTIONS] mol [MOL_OPTIONS] nci [SUBCMD_OPTIONS]

NCI-Specific OPTIONS
^^^^^^^^^^^^^^^^^^^^

.. list-table:: NCI Job Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-i, --isosurface``
     - float
     - Isosurface value for NCI plot (default=0.5)
   * - ``-r, --color-range``
     - float
     - Ramp range for NCI plot (default=1.0)
   * - ``-b, --binary``
     - bool
     - Plot NCI plots with two colors only (default=False)
   * - ``--intermediate``
     - bool
     - Plot NCI plots with intermediate range colors (default=False)

.. note::

    NCI jobs also inherit all options from visualization jobs including styling, ray tracing, and surface rendering options.

Basic Usage
^^^^^^^^^^^

**Basic NCI plot**

    .. code-block:: console

        chemsmart run mol -f complex.xyz nci

**NCI plot with custom color range**

    .. code-block:: console

        chemsmart run mol -f complex.xyz nci -r 1.5

**Binary color NCI plot**

    .. code-block:: console

        chemsmart run mol -f hydrogen_bond.xyz nci -b

**High-quality NCI plot with ray tracing**

    .. code-block:: console

        chemsmart run mol -f complex.xyz nci -i 0.4 -t -s pymol
