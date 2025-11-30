##########################
 Thermochemistry Analysis
##########################

Chemsmart provides thermochemistry analysis capabilities for computing
thermodynamic properties from Gaussian and ORCA output files.

*******************************
 Thermochemistry Analysis Jobs
*******************************

The ``thermochemistry`` command parses output files and calculates
thermochemical properties including enthalpy, entropy, and Gibbs free
energy.

Usage
=====

.. code:: bash

   chemsmart run thermochemistry [-d path/to/directory] [-t log|out] [-f filename(s)]
                                 [-csg s_freq_cutoff] [-cst s_freq_cutoff]
                                 [-ch h_freq_cutoff] [-c concentration] [-p pressure] [-w]
                                 [-T temperature] [-a alpha] [-u hartree|eV|kcal/mol|kJ/mol]
                                 [-o outfile.dat] [-O] [-i] [-S|-R]

Options
=======

**Input Options:**

.. list-table::
   :header-rows: 1
   :widths: 15 10 75

   -  -  Option
      -  Type
      -  Description

   -  -  ``-d, --directory``
      -  string
      -  Directory for batch processing (mutually exclusive with -f)

   -  -  ``-t, --filetype``
      -  string
      -  File type: log (Gaussian) or out (ORCA)

   -  -  ``-f, --filenames``
      -  string
      -  Specific file(s) to analyze (repeatable, mutually exclusive
         with -d)

**Quasi-RRHO Corrections:**

.. list-table::
   :header-rows: 1
   :widths: 15 10 75

   -  -  Option
      -  Type
      -  Description

   -  -  ``-csg, --cutoff-entropy-grimme``
      -  float
      -  Grimme's quasi-RRHO entropy cutoff frequency (cm⁻¹)

   -  -  ``-cst, --cutoff-entropy-truhlar``
      -  float
      -  Truhlar's quasi-RRHO entropy cutoff frequency (cm⁻¹)

   -  -  ``-ch, --cutoff-enthalpy``
      -  float
      -  Head-Gordon's quasi-RRHO enthalpy cutoff frequency (cm⁻¹)

.. note::

   ``-csg`` and ``-cst`` are mutually exclusive. If neither is
   specified, no quasi-RRHO entropy correction is applied.

**Thermodynamic Conditions:**

.. list-table::
   :header-rows: 1
   :widths: 15 10 75

   -  -  Option
      -  Type
      -  Description

   -  -  ``-c, --concentration``
      -  float
      -  Solution concentration in mol/L (mutually exclusive with -p)

   -  -  ``-p, --pressure``
      -  float
      -  Gas-phase pressure in atm (default: 1.0)

   -  -  ``-w, --weighted``
      -  bool
      -  Use isotopically weighted masses

   -  -  ``-T, --temperature``
      -  float
      -  Temperature in Kelvin (required)

   -  -  ``-a, --alpha``
      -  int
      -  Interpolator exponent (default: 4)

**Output Options:**

.. list-table::
   :header-rows: 1
   :widths: 15 10 75

   -  -  Option
      -  Type
      -  Description

   -  -  ``-u, --energy-units``
      -  string
      -  Units: hartree, eV, kcal/mol, kJ/mol (default: hartree)

   -  -  ``-o, --outputfile``
      -  string
      -  Output filename (default: <basename>.dat)

   -  -  ``-O, --overwrite``
      -  bool
      -  Overwrite existing output files

   -  -  ``-i, --check-imaginary-frequencies``
      -  bool
      -  Check for imaginary frequencies

   -  -  ``-S/-R, --skip-completed/--no-skip-completed``
      -  bool
      -  Skip or rerun completed jobs

Examples
========

**Single file analysis:**

.. code:: bash

   chemsmart run thermochemistry -T 298.15 -f water_opt.out

   Structure                        E        ZPE          H       T.S        G(T)
   ================================================================================
   water_opt               -76.323311   0.021581  -76.297951   0.021430  -76.319381

**Multiple files:**

.. code:: bash

   chemsmart run thermochemistry -T 298.15 -f he_gaussian.log -f he_orca.out

**Batch processing:**

.. code:: bash

   chemsmart run thermochemistry -T 298.15 -d . -t log -o thermo.dat

**Solution phase:**

.. code:: bash

   chemsmart run thermochemistry -T 298.15 -f co2.log -c 1.0

**With quasi-RRHO corrections:**

.. code:: bash

   # Grimme's entropy correction
   chemsmart run thermochemistry -T 298.15 -f water_mp2.log -csg 100

   # Head-Gordon's enthalpy correction
   chemsmart run thermochemistry -T 298.15 -f water_mp2.log -ch 200

   # Both corrections
   chemsmart run thermochemistry -T 298.15 -f water_mp2.log -cst 40 -ch 100

***********************************
 Boltzmann Weighted Averaging Jobs
***********************************

The ``boltzmann`` subcommand performs Boltzmann-weighted averaging of
thermochemical results across multiple conformers.

Usage
=====

.. code:: bash

   chemsmart run thermochemistry [THERMO_OPTIONS] boltzmann [-w gibbs|electronic] [-S|-R]

Options
=======

.. list-table::
   :header-rows: 1
   :widths: 15 10 75

   -  -  Option
      -  Type
      -  Description

   -  -  ``-w, --energy-type-for-weighting``
      -  string
      -  Weighting scheme: gibbs or electronic (default: gibbs)

   -  -  ``-S/-R, --skip-completed/--no-skip-completed``
      -  bool
      -  Skip or rerun completed jobs

Examples
========

**Boltzmann averaging of conformers:**

.. code:: bash

   chemsmart run thermochemistry -T 298.15 -f conf1.log -f conf2.log boltzmann -w electronic

   Structure                                        E        ZPE          H       T.S        G(T)
   ================================================================================================
   conformer_boltzmann_avg_by_electronic   -2189.631938   0.288732  -2189.312517   0.097016  -2189.409533

**Batch directory averaging:**

.. code:: bash

   chemsmart run thermochemistry -T 298.15 -d . -t log boltzmann

Results are saved to ``thermochemistry_job_boltzmann.dat``.
