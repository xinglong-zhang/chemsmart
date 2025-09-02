Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to
contact our team if you have questions or feedback.

###################################
 Run Thermochemistry Analysis Jobs
###################################

*******************************
 Thermochemistry Analysis Jobs
*******************************

In CHEMSMART, the ``thermochemistry`` command currently supports parsing output files from **Gaussian** and **ORCA**,
enabling users to recalculate and apply corrections to thermochemical properties such as enthalpy, entropy, and Gibbs
free energy based on quantum chemical results.

USAGE
=====

.. code:: console

   chemsmart run thermochemistry [-d path/to/directory] [-t log|out] [-f filename(s)]
                                 [-csg s_freq_cutoff] [-cst s_freq_cutoff]
                                 [-ch h_freq_cutoff] [-c concentration] [-p pressure] [-w]
                                 [-T temperature] [-a alpha] [-u hartree|eV|kcal/mol|kJ/mol]
                                 [-o outfile.dat] [-O] [-i] [-S|-R]

THERMO_OPTIONS
==============

.. list-table::
   :header-rows: 1
   :widths: 15 10 75

   -  -  Option
      -  Type
      -  Description

   -  -  ``-d, --directory``
      -  string
      -  |  Specifies a directory containing files to be processed in batch.
         |  This option is mutually exclusive with ``-f``.

   -  -  ``-t, --filetype``
      -  string
      -  |  Specifies the file type in the directory, either ``log`` (Gaussian) or ``out`` (ORCA).
         |  This option is only used together with ``-d``.

   -  -  ``-f, --filenames``

      -  string

      -  |  Specifies one or more individual files to analyze.
         |  This option can be used multiple times and is mutually exclusive with ``-d``
         |  and ``-t``.

   -  -  ``-csg, --cutoff-entropy-grimme``

      -  float

      -  |  Specifies the cut-off frequency (in wavenumbers) for applying Grimme’s
         |  quasi-rigid-rotor-harmonic-oscillator (quasi-RRHO) correction to entropy.
         |  For example, ``-csg 100`` applies cut-off frequency of 100 :math:`\rm cm^{-1}` to entropy
         |  when computing corrected entropy (*qh-S*).

   -  -  ``-cst, --cutoff-entropy-truhlar``

      -  float

      -  |  Specifies the cut-off frequency (in wavenumbers) for applying Truhlar’s
         |  quasi-RRHO approach to entropy.
         |  Vibrational modes with frequencies below the cut-off are raised to the
         |  cut-off value instead of being treated by the standard RRHO approximation.

.. note::

   The ``-csg`` and ``-cst`` options are mutually exclusive, and only one of them can be specified at most. If neither
   is specified, the quasi-RRHO entropy correction is not applied and only the standard entropy (*S*) is calculated.

.. list-table::
   :header-rows: 0
   :widths: 15 10 75

   -  -  ``-ch, --cutoff-enthalpy``

      -  float

      -  |  Specifies the cut-off frequency (in wavenumbers) for applying
         |  Head-Gordon’s quasi-RRHO correction to enthalpy.
         |  For example, ``-ch 50`` applies a cut-off frequency of 50 :math:`\rm cm^{-1}` to enthalpy
         |  when calculating corrected enthalpy (*qh-H*).
         |  If not specified, no quasi-RRHO correction is applied to enthalpy and only
         |  the standard enthalpy (*H*) is calculated.

   -  -  ``-c, --concentration``

      -  float

      -  |  Specifies the solution concentration (in mol/L).
         |  In this case, the translational partition function is evaluated in the
         |  concentration-dependent form, instead of the pressure-dependent form
         |  used in the gas phase.
         |  This correction follows directly from the Sackur-Tetrode equation of
         |  translational entropy.
         |  If ``-c`` is not specified, the calculation defaults to the gas phase.

   -  -  ``-p, --pressure``
      -  float
      -  |  Specifies the gas-phase pressure (in atm), with a default value of 1.0 atm.

.. note::

   The ``-c`` and ``-p`` options are mutually exclusive. If both are provided, ``-c`` takes precedence and ``-p`` is
   ignored.

.. list-table::
   :header-rows: 0
   :widths: 15 10 75

   -  -  ``-w, --weighted``
      -  bool
      -  |  Uses masses weighted by natural isotopic abundance.
         |  By default, the most abundant isotope masses are used.

   -  -  ``-T, --temperature``
      -  float
      -  |  **Required**. Specifies the temperature (in Kelvin).

   -  -  ``-a, --alpha``
      -  int
      -  |  Specifies the interpolator exponent (alpha) in the quasi-RRHO scheme.
         |  The default value is 4.

   -  -  ``-u, --energy-units``

      -  string

      -  |  Specifies the units for energetic values.
         |  Options include ``hartree``, ``eV``, ``kcal/mol``, and ``kJ/mol``.
         |  The default unit is ``hartree``.

   -  -  ``-o, --outputfile``
      -  string
      -  |  Specifies the output file to save the thermochemistry results.
         |  By default, results are saved to a file named ``file_basename.dat``.

   -  -  ``-O, --overwrite``
      -  bool
      -  |  Overwrites existing output files if they already exist.
         |  By default, the results are appended to the existing file.

   -  -  ``-i, --check-imaginary-frequencies``
      -  bool
      -  |  Checks for imaginary frequencies in the calculations.
         |  By default, this check is performed.

   -  -  ``-S, --skip-completed / -R, --no-skip-completed``

      -  bool

      -  |  Controls how completed jobs are handled.
         |  ``-R`` forces rerun completed Thermochemistry jobs, ``-S`` skips completed jobs.
         |  By default, ``-S`` is enabled.

EXAMPLES
========

**Thermochemical calculation for a single file**

.. code:: console

   chemsmart run thermochemistry -T 298.15 -f water_opt.out

   Structure                                           E        ZPE             H        T.S          G(T)
   =======================================================================================================
   water_opt                                  -76.323311   0.021581    -76.297951   0.021430    -76.319381

Calculate standard thermochemical properties (*E*, *ZPE*, *H*, *T.S*, *G(T)*) at 298.15 K for ``water_opt.out``. Results
are saved as ``water_opt.dat``.

**Thermochemical calculation for multiple files**

.. code:: console

   chemsmart run thermochemistry -T 298.15 -f he_gaussian.log -f he_orca.out

   Structure                                           E        ZPE             H        T.S          G(T)
   =======================================================================================================
   he_gaussian                                 -2.915130   0.000000     -2.912769   0.014313     -2.927083

   Structure                                           E        ZPE             H        T.S          G(T)
   =======================================================================================================
   he_orca                                     -2.899161   0.000000     -2.896800   0.014313     -2.911114

Compute standard thermochemical properties for both ``he_gaussian.log`` and ``he_orca.out``.

Results are saved as ``he_gaussian.dat`` and ``he_orca.dat``, respectively.

**Batch calculation for all files of a specified type**

.. code:: console

   chemsmart run thermochemistry -T 298.15 -d </path/to/directory> -t log -o thermo.dat

Process all ``.log`` files in the specified directory and summarize results in ``thermo.dat``.

**Solution vs gas-phase calculation**

.. code:: console

   chemsmart run thermochemistry -T 298.15 -f co2.log -c 1.0

   Structure                                           E        ZPE             H        T.S          G(T)
   =======================================================================================================
   co2                                       -188.444680   0.011776   -188.429325   0.021262   -188.450587

Calculate standard thermochemical properties in the solution at a concentration of 1.0 mol/L.

.. code:: console

   chemsmart run thermochemistry -T 298.15 -f co2.log -p 1.0

   Structure                                           E        ZPE             H        T.S          G(T)
   =======================================================================================================
   co2                                       -188.444680   0.011776   -188.429325   0.024281   -188.453606

Calculate standard thermochemical properties in the gas phase at a pressure of 1.0 atm.

**Thermochemical calculation with quasi-RRHO corrections**

.. code:: console

   chemsmart run thermochemistry -T 298.15 -f water_mp2.log -csg 100

   Structure                                           E        ZPE             H        T.S     T.qh-S          G(T)       qh-G(T)
   ================================================================================================================================
   water_mp2                                  -76.328992   0.021410    -76.303803   0.021424   0.021424    -76.325227    -76.325227

Apply Grimme’s quasi-RRHO correction with a frequency cut-off of 100 :math:`\rm cm^{-1}` for entropy to obtain the
corrected entropy (*qh-S*).

Output includes *E*, *ZPE*, *H*, *T.S*, *T.qh-S*, *G(T)*, and *qh-G(T)*, where *qh-G(T)* is calculated using *qh-S*
together with *H*.

.. code:: console

   chemsmart run thermochemistry -T 298.15 -f water_mp2.log -ch 200

   Structure                                           E        ZPE             H          qh-H        T.S          G(T)       qh-G(T)
   ===================================================================================================================================
   water_mp2                                  -76.328992   0.021410    -76.303803    -76.303804   0.021424    -76.325227    -76.325228

Apply Head-Gordon’s quasi-RRHO correction with a frequency cut-off of 200 :math:`\rm cm^{-1}` for enthalpy to calculate
the corrected enthalpy (*qh-H*).

Output includes *E*, *ZPE*, *H*, *qh-H*, *T.S*, *G(T)*, and *qh-G(T)*, where *qh-G(T)* is calculated using *qh-H*
together with *S*.

.. code:: console

   chemsmart run thermochemistry -T 298.15 -f water_mp2.log -cst 40 -ch 100

   Structure                                           E        ZPE             H          qh-H        T.S     T.qh-S          G(T)       qh-G(T)
   ==============================================================================================================================================
   water_mp2                                  -76.328992   0.021410    -76.303803    -76.303803   0.021424   0.021424    -76.325227    -76.325227

Perform thermochemical analysis by applying Truhlar’s quasi-RRHO method with a frequency cut-off of 40 :math:`\rm
cm^{-1}` to entropy and Head-Gordon’s method with a frequency cut-off of 100 :math:`\rm cm^{-1}` to enthalpy.

Output includes *E*, *ZPE*, *H*, *qh-H*, *T.S*, *T.qh-S*, *G(T)*, and *qh-G(T)*, where *qh-G(T)* is calculated using
both *qh-H* and *qh-S*.

***********************************
 Boltzmann Weighted Averaging Jobs
***********************************

The ``boltzmann`` subcommand performs Boltzmann-weighted averaging of thermochemical results across multiple conformers
to obtain the overall thermochemical properties of the system.

USAGE
=====

.. code:: console

   chemsmart run thermochemistry [THERMO_OPTIONS]
                 boltzmann [-w gibbs|electronic] [-S|-R]

BOLTZMANN_OPTIONS
=================

.. list-table::
   :header-rows: 1
   :widths: 15 10 75

   -  -  Option
      -  Type
      -  Description

   -  -  ``-w, --energy-type-for-weighting``

      -  string

      -  |  Specifies the weighting scheme for Boltzmann averaging based on the type of
         |  energy.
         |  Options include ``gibbs`` (default) and ``electronic``.

   -  -  ``-S, --skip-completed / -R, --no-skip-completed``

      -  bool

      -  |  Controls how completed jobs are handled.
         |  ``-R`` forces rerun completed Boltzmann jobs, ``-S`` skips completed jobs.
         |  By default, ``-S`` is enabled.

EXAMPLES
========

**Boltzmann-weighted averaging for two conformers**

.. code:: console

   chemsmart run thermochemistry -T 298.15 -f udc3_mCF3_monomer_c1.log -f udc3_mCF3_monomer_c4.log boltzmann -w electronic

   Structure                                           E        ZPE             H        T.S          G(T)
   =======================================================================================================
   udc3_mCF3_monomer_c_boltzmann_avg_by_electronic  -2189.631938   0.288732  -2189.312517   0.097016  -2189.409533

Perform thermochemical analysis for ``udc3_mCF3_monomer_c1.log`` and ``udc3_mCF3_monomer_c4.log``, followed by
Boltzmann-weighted averaging using electronic energies.

Results will be saved to ``thermochemistry_job_boltzmann.dat``.

**Boltzmann-weighted averaging for batch directory calculation**

.. code:: console

   chemsmart run thermochemistry -T 298.15 -d . -t log boltzmann

   Structure                                           E        ZPE             H        T.S          G(T)
   =======================================================================================================
   udc3_mCF3_monomer_c_boltzmann_avg_by_gibbs  -2189.631914   0.288696  -2189.312512   0.097155  -2189.409667

Perform thermochemical calculations for all ``.log`` files in current directory and apply Boltzmann-weighted average
using Gibbs free energies.

Results will be saved to ``thermochemistry_job_boltzmann.dat``.
