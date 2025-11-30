###############################
 Electronic Structure Analysis
###############################

This page documents scripts for analyzing electronic structure properties
from Gaussian and ORCA output files.

*********************
 FMO Analysis Script
*********************

The ``fmo.py`` script performs Frontier Molecular Orbital (FMO) analysis,
extracting:

- HOMO and LUMO energies
- HOMO-LUMO energy gap
- Chemical potential
- Chemical hardness
- Electrophilicity index

Currently supports restricted (closed-shell) calculations only.

Usage
=====

.. code-block:: bash

   fmo.py [-f filename] [-u eV|kcal/mol]

Options
=======

.. list-table::
   :header-rows: 1
   :widths: 20 10 70

   * - Option
     - Type
     - Description
   * - ``-f, --filename``
     - string
     - Output file to analyze (required)
   * - ``-u, --unit``
     - string
     - Energy units: eV or kcal/mol (default: eV)

Examples
========

.. code-block:: bash

   fmo.py -f water_mp2.log

   HOMO energy: -13.8816 eV
   LUMO energy: 0.7992 eV
   HOMO-LUMO gap: 14.6808 eV
   Chemical potential, mu: -6.541 eV
   Chemical hardness, eta: 7.34 eV
   Electrophilicity Index, omega: 2.915 eV

.. code-block:: bash

   fmo.py -f water_mp2.log -u kcal/mol

   HOMO energy: -320.1176 kcal/mol
   LUMO energy: 18.4299 kcal/mol
   HOMO-LUMO gap: 338.5475 kcal/mol
   ...

*************************************
 Mulliken Population Analysis Script
*************************************

The ``mulliken.py`` script extracts Mulliken atomic charges and spin densities.
Supports Gaussian and ORCA output files.

Usage
=====

.. code-block:: bash

   mulliken.py [-f filename] [-n atom_number]

Options
=======

.. list-table::
   :header-rows: 1
   :widths: 20 10 70

   * - Option
     - Type
     - Description
   * - ``-f, --filename``
     - string
     - Output file to analyze (required)
   * - ``-n, --numbers``
     - int
     - Atom number(s) to report (1-indexed)

Examples
========

.. code-block:: bash

   mulliken.py -f molecule.log

   Mulliken Charges:
   C1      :     0.062
   C2      :    -0.008
   ...

.. code-block:: bash

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

.. code-block:: bash

   hirshfeld.py [-f filename] [-n atom_number]

Options
=======

.. list-table::
   :header-rows: 1
   :widths: 20 10 70

   * - Option
     - Type
     - Description
   * - ``-f, --filename``
     - string
     - Output file with Hirshfeld analysis (required)
   * - ``-n, --numbers``
     - int
     - Atom number(s) to report (1-indexed)

Example
=======

.. code-block:: bash

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

The ``wbi_analysis.py`` script extracts Natural Population Analysis (NPA) data
from Wiberg Bond Index calculations. See :ref:`wbi-jobs` for running WBI
calculations.

Usage
=====

.. code-block:: bash

   wbi_analysis.py [-f filename] [-n atom_number]

Options
=======

.. list-table::
   :header-rows: 1
   :widths: 20 10 70

   * - Option
     - Type
     - Description
   * - ``-f, --filename``
     - string
     - WBI output file (required)
   * - ``-n, --numbers``
     - int
     - Atom number(s) to report (1-indexed)

Example
=======

.. code-block:: bash

   wbi_analysis.py -f ts_wbi.log -n 2 -n 5

   Natural Charges:
   Ni1     :     0.528
   P2      :     0.930
   ...

   Natural Charge at atom 2 is 0.930.
   Natural Charge at atom 5 is -0.562.
