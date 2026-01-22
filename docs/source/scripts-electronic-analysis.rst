###############################
 Electronic Structure Analysis
###############################

This page documents scripts for analyzing electronic structure properties from Gaussian and ORCA output files.

*********************
 FMO Analysis Script
*********************

The ``fmo.py`` script performs Frontier Molecular Orbital (FMO) analysis, extracting and calculating:

-  HOMO and LUMO energies
-  HOMO-LUMO energy gap
-  Chemical potential
-  Chemical hardness
-  Electrophilicity index

The script supports both closed-shell (restricted) and open-shell (unrestricted) calculations. For open-shell systems,
it additionally provides:

-  Alpha and beta spin channels (α-HOMO, α-LUMO, β-HOMO, β-LUMO)
-  Singly Occupied Molecular Orbitals (SOMOs)
-  Spin-specific reactivity descriptors
-  Overall FMO gap (highest SOMO to lowest LUMO)

Usage
=====

.. code:: bash

   fmo.py [-f filename] [-u eV|kcal/mol]

Options
=======

.. list-table::
   :header-rows: 1
   :widths: 20 10 70

   -  -  Option
      -  Type
      -  Description

   -  -  ``-f, --filename``
      -  string
      -  Output file to analyze (required)

   -  -  ``-u, --unit``
      -  string
      -  Energy units: eV or kcal/mol (default: eV)

Examples
========

Closed-Shell System
-------------------

.. code:: text

   fmo.py -f fe2_singlet.out -u kcal/mol

   System multiplicity: 1

   ============================================================
   CLOSED-SHELL SYSTEM FMO ANALYSIS
   ============================================================

   HOMO energy: -463.5399 kcal/mol
   LUMO energy: -253.2063 kcal/mol
   HOMO-LUMO gap: 210.3336 kcal/mol
   Chemical potential, μ: -358.3731 kcal/mol
   Chemical hardness, η: 105.1668 kcal/mol
   Electrophilicity index, ω: 610.6075 kcal/mol

Open-Shell System
-----------------

.. code:: text

   fmo.py -f fe2_triplet.out

   System multiplicity: 3

   ============================================================
   OPEN-SHELL SYSTEM FMO ANALYSIS
   ============================================================

   --- Alpha Spin Channel ---
   α-HOMO energy: -20.7301 eV
   α-LUMO energy: -9.4291 eV
   α-HOMO-LUMO gap: 11.3010 eV
   Chemical potential, μ_α: -15.0796 eV
   Chemical hardness, η_α: 5.6505 eV
   Electrophilicity index, ω_α: 20.1216 eV

   --- Beta Spin Channel ---
   β-HOMO energy: -19.6200 eV
   β-LUMO energy: -10.3698 eV
   β-HOMO-LUMO gap: 9.2502 eV
   Chemical potential, μ_β: -14.9949 eV
   Chemical hardness, η_β: 4.6251 eV
   Electrophilicity index, ω_β: 24.3072 eV

   --- Singly Occupied Molecular Orbitals (SOMOs) ---
   Number of unpaired electrons: 2
   SOMO-1 energy: -20.8980 eV
   SOMO-2 energy: -20.7301 eV
   Lowest SOMO energy: -20.8980 eV
   Highest SOMO energy: -20.7301 eV

   --- Overall FMO Gap ---
   FMO gap (highest SOMO to lowest LUMO): 10.3603 eV

*************************************
 Mulliken Population Analysis Script
*************************************

The ``mulliken.py`` script extracts Mulliken atomic charges and spin densities. Supports Gaussian and ORCA output files.

Usage
=====

.. code:: bash

   mulliken.py [-f filename] [-n atom_number]

Options
=======

.. list-table::
   :header-rows: 1
   :widths: 20 10 70

   -  -  Option
      -  Type
      -  Description

   -  -  ``-f, --filename``
      -  string
      -  Output file to analyze (required)

   -  -  ``-n, --numbers``
      -  int
      -  Atom number(s) to report (1-indexed)

Examples
========

.. code:: bash

   mulliken.py -f molecule.log

   Mulliken Charges:
   C1      :     0.062
   C2      :    -0.008
   ...

.. code:: bash

   mulliken.py -f triplet.log -n 3

   Mulliken Charges:
   ...
   Mulliken Spin densities:
   ...
   Mulliken Charge at P3 is 0.481.
   Mulliken Spin density at P3 is -0.058.

**************************************
 Hirshfeld Population Analysis Script
**************************************

The ``hirshfeld.py`` script extracts Hirshfeld charges and spin densities.

Usage
=====

.. code:: bash

   hirshfeld.py [-f filename] [-n atom_number]

Options
=======

.. list-table::
   :header-rows: 1
   :widths: 20 10 70

   -  -  Option
      -  Type
      -  Description

   -  -  ``-f, --filename``
      -  string
      -  Output file with Hirshfeld analysis (required)

   -  -  ``-n, --numbers``
      -  int
      -  Atom number(s) to report (1-indexed)

Example
=======

.. code:: bash

   hirshfeld.py -f calculation_hirshfeld.out

   Hirshfeld Charges:
   O1      :    -0.177
   O2      :    -0.319
   ...

   Hirshfeld Spins:
   O1      :     0.000
   ...

***********************************
 Wiberg Bond Index Analysis Script
***********************************

The ``wbi_analysis.py`` script extracts Natural Population Analysis (NPA) data from Wiberg Bond Index calculations. See
:ref:`wbi-jobs` for running WBI calculations.

Usage
=====

.. code:: bash

   wbi_analysis.py [-f filename] [-n atom_number]

Options
=======

.. list-table::
   :header-rows: 1
   :widths: 20 10 70

   -  -  Option
      -  Type
      -  Description

   -  -  ``-f, --filename``
      -  string
      -  WBI output file (required)

   -  -  ``-n, --numbers``
      -  int
      -  Atom number(s) to report (1-indexed)

Example
=======

.. code:: bash

   wbi_analysis.py -f ts_wbi.log -n 2 -n 5

   Natural Charges:
   Ni1     :     0.528
   P2      :     0.930
   ...

   Natural Charge at atom 2 is 0.930.
   Natural Charge at atom 5 is -0.562.
