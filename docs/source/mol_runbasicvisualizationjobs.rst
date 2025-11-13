Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to
contact our team if you have questions or feedback.

##############################
 Run Basic Visualization Jobs
##############################

ChemSmart provides powerful molecular visualization capabilities using PyMOL for creating high-quality molecular
graphics, movies, and interactive visualizations.

********************
 Visualization Jobs
********************

Create static PyMOL visualizations and interactive session files.

.. code:: console

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] visualize [SUBCMD_OPTIONS]

Visualization-Specific OPTIONS
==============================

.. list-table:: Visualization Job Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-f, --file``
      -  string
      -  PyMOL file script or style. If not specified, defaults to use zhang_group_pymol_style.py (default=None)

   -  -  ``-s, --style``
      -  string
      -  PyMOL render style. Options: pymol, cylview (default=None)

   -  -  ``-t, --trace/--no-trace``
      -  bool
      -  PyMOL options to ray trace or not (default=True)

   -  -  ``-v, --vdw``
      -  bool
      -  Add Van der Waals surface (default=False)

   -  -  ``-q, --quiet/--no-quiet``
      -  bool
      -  Run PyMOL in quiet mode (default=False)

   -  -  ``--command-line-only/--no-command-line-only``
      -  bool
      -  Run PyMOL in command line only (default=True)

   -  -  ``-c, --coordinates <string>``
      -  string
      -  List of coordinates (bonds, angles and dihedrals) for labelling. 1-indexed (default=None)

Visualization Basic Usage
=========================

**Basic molecular visualization**

   .. code:: console

      chemsmart run mol -f molecule.xyz visualize

**Quiet mode visualization**

   .. code:: console

      chemsmart run mol -f calculation.log visualize -q

**Visualization with coordinate labeling**

   .. code:: console

      chemsmart run mol -f structure.xyz visualize -c [[1,2,3]]

**Visualization using custom style or script file**

   .. code:: console

      chemsmart run mol -f molecule.log visualize -f custom_style.py

Visualization Examples
======================

************
 Movie Jobs
************

Generate rotating movie animations of molecular structures.

.. code:: console

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] movie [SUBCMD_OPTIONS]

Movie-Specific OPTIONS
======================

Movie jobs inherit all options from visualization jobs and use the same parameters.

Movie Basic Usage
=================

**Basic rotating movie**

   .. code:: console

      chemsmart run mol -f molecule.xyz movie

Movie Examples
==============
##############################
 Run Hybrid Visualization Jobs
##############################

Apart from basic visualization jobs, Chemsmart can also creating graphics with different substrates displayed
in different style (what 'hybrid' refers to).

********************
 Hybrid Visualization Jobs
********************

Create static PyMOL visualizations and interactive session files.

.. code:: console

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] visualize [SUBCMD_OPTIONS] --hybrid [SUBCMD_OPTIONS]

Visualization-Specific OPTIONS
==============================

.. list-table:: Hybrid Visualization Job Options
   :header-rows: 1
   :widths: 25 15 60

   -  -  Option
      -  Type
      -  Description

   -  -  ``--hybrid``
      -  bool
      -  Enable hybrid visualization mode. Allows drawing different groups in different styles (default=False).

   -  -  ``-g, --group <string>``
      -  string, multiple
      -  Indexes of atoms to select for a group. Repeatable for multiple groups, e.g., -g '1-5' -g '6,7,8'.

   -  -  ``-C, --color <string>``
      -  string, multiple
      -  Color for each group. Repeatable to match -g options. Example color schemes: ['cbap', 'cbac', 'cbay', 'cbag', ...].

   -  -  ``-sc, --surface-color <string>``
      -  string
      -  Customize the surface color of the molecule (default=None).

   -  -  ``-st, --surface-transparency <string>``
      -  string
      -  Customize the surface transparency of the molecule (default=None).

   -  -  ``-nc, --new-color-carbon <string>``
      -  string
      -  Set a custom color for carbon atoms. Example: -nc [0.8, 0.8, 0.9] (default=None).

   -  -  ``-nn, --new-color-nitrogen <string>``
      -  string
      -  Set a custom color for nitrogen atoms. Example: -nn [0.6, 0.8, 1.0] (default=None).

   -  -  ``-no, --new-color-oxygen <string>``
      -  string
      -  Set a custom color for oxygen atoms. Example: -no [1.0, 0.7, 0.7] (default=None).

   -  -  ``-ns, --new-color-sulfur <string>``
      -  string
      -  Set a custom color for sulfur atoms. Example: -ns [0.8, 0.8, 0.9] (default=None).

   -  -  ``-np, --new-color-phosphorus <string>``
      -  string
      -  Set a custom color for phosphorus atoms. Example: -np [1.0, 0.85, 0.6] (default=None).

Visualization Basic Usage
=========================

**Basic hybrid visualization**

   .. code:: console

      chemsmart run mol -f molecule.xyz visualize --hybrid -g [1,2,3]
This is the minimal command needed for hybrid visualization job. Only the highlighted groups is specified, others will be assigned by default.

**Hybrid visualization with customized color scheme**

   .. code:: console

      chemsmart run mol -f calculation.log visualize --hybrid -g [1,2,3] -g [4,5,6] -c cbay -c cbak

**Visualization with coordinate labeling**

   .. code:: console

      chemsmart run mol -f structure.xyz visualize -c [[1,2,3]]

**Visualization using custom style or script file**

   .. code:: console

      chemsmart run mol -f molecule.log visualize -f custom_style.py

Visualization Examples
======================

************
 Movie Jobs
************

Generate rotating movie animations of molecular structures.

.. code:: console

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] movie [SUBCMD_OPTIONS]

Movie-Specific OPTIONS
======================

Movie jobs inherit all options from visualization jobs and use the same parameters.

Movie Basic Usage
=================

**Basic rotating movie**

   .. code:: console

      chemsmart run mol -f molecule.xyz movie

Movie Examples
==============
