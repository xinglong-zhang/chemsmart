Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to
contact our team if you have questions or feedback.

###################################
 Run Thermochemistry Analysis Jobs
###################################

*******************************
 Thermochemistry Analysis Jobs
*******************************

In CHEMSMART, the ``thermochemistry`` command currently supports parsing output files from **Gaussian** and **ORCA**,
enabling users to recalculate and apply corrections to thermochemical properties such as enthalpy, entropy, and Gibbs free energy
based on quantum chemical results.

USAGE
=====

.. code-block:: console

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

   * - Option
     - Type
     - Description
   * - ``-d, --directory``
     - string
     - | Specifies a directory containing files to be processed in batch.
       | This option is mutually exclusive with ``-f``.
   * - ``-t, --filetype``
     - string
     - | Specifies the file type in the directory, either ``log`` (Gaussian) or ``out`` (ORCA).
       | This option is only used together with ``-d``.
   * - ``-f, --filenames``
     - string
     - | Specifies one or more individual files to analyze.
       | This option can be used multiple times and is mutually exclusive with ``-d``
       | and ``-t``.
   * - ``-csg, --cutoff-entropy-grimme``
     - float
     - | Specifies the cut-off frequency (in wavenumbers) for applying Grimme’s
       | quasi-rigid-rotor-harmonic-oscillator (quasi-RRHO) correction to entropy.
       | For example, ``-csg 100`` applies cut-off frequency of 100 :math:`\rm cm^{-1}` to entropy
       | when computing corrected entropy (*qh-S*).
   * - ``-cst, --cutoff-entropy-truhlar``
     - float
     - | Specifies the cut-off frequency (in wavenumbers) for applying Truhlar’s
       | quasi-RRHO approach to entropy.
       | Vibrational modes with frequencies below the cut-off are raised to the
       | cut-off value instead of being treated by the standard RRHO approximation.
.. note::
       The ``-csg`` and ``-cst`` options are mutually exclusive, and only one of them can be specified at most.
       If neither is specified, the quasi-RRHO entropy correction is not applied and only the standard entropy (*S*) is calculated.
.. list-table::
   :header-rows: 0
   :widths: 15 10 75

   * - ``-ch, --cutoff-enthalpy``
     - float
     - | Specifies the cut-off frequency (in wavenumbers) for applying
       | Head-Gordon’s quasi-RRHO correction to enthalpy.
       | For example, ``-ch 50`` applies a cut-off frequency of 50 :math:`\rm cm^{-1}` to enthalpy
       | when calculating corrected enthalpy (*qh-H*).
       | If not specified, no quasi-RRHO correction is applied to enthalpy and only
       | the standard enthalpy (*H*) is calculated.
   * - ``-c, --concentration``
     - float
     - | Specifies the solution concentration (in mol/L).
       | In this case, the translational partition function is evaluated in the
       | concentration-dependent form, instead of the pressure-dependent form
       | used in the gas phase.
       | This correction follows directly from the Sackur-Tetrode equation of
       | translational entropy.
       | If ``-c`` is not specified, the calculation defaults to the gas phase.
   * - ``-p, --pressure``
     - float
     - | Specifies the gas-phase pressure (in atm), with a default value of 1.0 atm.
.. note::
       The ``-c`` and ``-p`` options are mutually exclusive. If both are provided, ``-c`` takes precedence and ``-p`` is ignored.
.. list-table::
   :header-rows: 0
   :widths: 15 10 75

   * - ``-w, --weighted``
     - bool
     - | Uses masses weighted by natural isotopic abundance.
       | By default, the most abundant isotope masses are used.
   * - ``-T, --temperature``
     - float
     - | **Required**. Specifies the temperature (in Kelvin).
   * - ``-a, --alpha``
     - int
     - | Specifies the interpolator exponent (alpha) in the quasi-RRHO scheme.
       | The default value is 4.
   * - ``-u, --energy-units``
     - string
     - | Specifies the units for energetic values.
       | Options include ``hartree``, ``eV``, ``kcal/mol``, and ``kJ/mol``.
       | The default unit is ``hartree``.
   * - ``-o, --outputfile``
     - string
     - | Specifies the output file to save the thermochemistry results.
       | By default, results are saved to a file named ``file_basename.dat``.
   * - ``-O, --overwrite``
     - bool
     - | Overwrites existing output files if they already exist.
       | By default, the results are appended to the existing file.
   * - ``-i, --check-imaginary-frequencies``
     - bool
     - | Checks for imaginary frequencies in the calculations.
       | By default, this check is performed.
   * - ``-S, --skip-completed / -R, --no-skip-completed``
     - bool
     - | Controls how completed jobs are handled.
       | ``-R`` forces rerun completed Thermochemistry jobs, ``-S`` skips completed jobs.
       | By default, ``-S`` is enabled.

EXAMPLES
========

**Thermochemical calculation for a single file**

.. code-block:: console

   chemsmart run thermochemistry -T 298.15 -f file.out

Calculate standard thermochemical properties (*E*, *ZPE*, *H*, *T.S*, *G(T)*) at 298.15 K for ``file.out``. Results are saved as ``file.dat``.

**Thermochemical calculation for multiple files**

.. code-block:: console

   chemsmart run thermochemistry -T 298.15 -f file1.log -f file2.out

Compute standard thermochemical properties for both ``file1.log`` and ``file2.out``.

Results are saved as ``file1.dat`` and ``file2.dat``, respectively.

**Batch calculation for all files of a specified type**

.. code-block:: console

   chemsmart run thermochemistry -T 298.15 -d </path/to/directory> -t log -o thermo.dat

Process all ``.log`` files in the specified directory and summarize results in ``thermo.dat``.

**Solution vs gas-phase calculation**

.. code-block:: console

   chemsmart run thermochemistry -T 298.15 -f file.log -c 0.5

Calculate standard thermochemical properties in the solution at a concentration of 0.5 mol/L.

.. code-block:: console

   chemsmart run thermochemistry -T 298.15 -f file.log -p 2.0

Calculate standard thermochemical properties in the gas phase at a pressure of 2.0 atm.

**Thermochemical calculation with quasi-RRHO corrections**

.. code-block:: console

   chemsmart run thermochemistry -T 298.15 -f file.out -csg 100

Apply Grimme’s quasi-RRHO correction with a frequency cut-off of 100 :math:`\rm cm^{-1}` for entropy to obtain the corrected entropy (*qh-S*).

Output includes *E*, *ZPE*, *H*, *T.S*, *T.qh-S*, *G(T)*, and *qh-G(T)*, where *qh-G(T)* is calculated using *qh-S* together with *H*.

.. code-block:: console

   chemsmart run thermochemistry -T 298.15 -f file.out -ch 200

Apply Head-Gordon’s quasi-RRHO correction with a frequency cut-off of 200 :math:`\rm cm^{-1}` for enthalpy to calculate the corrected enthalpy (*qh-H*).

Output includes *E*, *ZPE*, *H*, *qh-H*, *T.S*, *G(T)*, and *qh-G(T)*, where *qh-G(T)* is calculated using *qh-H* together with *S*.

.. code-block:: console

   chemsmart run thermochemistry -T 298.15 -f file.out -cst 40 -ch 100

Perform thermochemical analysis by applying Truhlar’s quasi-RRHO method with a frequency cut-off of 40 :math:`\rm cm^{-1}` to entropy and Head-Gordon’s method with a frequency cut-off of 100 :math:`\rm cm^{-1}` to enthalpy.

Output includes *E*, *ZPE*, *H*, *qh-H*, *T.S*, *T.qh-S*, *G(T)*, and *qh-G(T)*, where *qh-G(T)* is calculated using both *qh-H* and *qh-S*.


***********************************
 Boltzmann Weighted Averaging Jobs
***********************************

The ``boltzmann`` subcommand performs Boltzmann-weighted averaging of thermochemical results across multiple conformers to obtain the overall thermochemical properties of the system.

USAGE
=====

.. code-block:: console

    chemsmart run thermochemistry [THERMO_OPTIONS]
                  boltzmann [-w gibbs|electronic] [-S|-R]

BOLTZMANN_OPTIONS
=================

.. list-table::
   :header-rows: 1
   :widths: 15 10 75

   * - Option
     - Type
     - Description
   * - ``-w, --energy-type-for-weighting``
     - string
     - | Specifies the weighting scheme for Boltzmann averaging based on the type of energy.
       | Options include ``gibbs`` (default) and ``electronic``.
   * - ``-S, --skip-completed / -R, --no-skip-completed``
     - bool
     - | Controls how completed jobs are handled.
       | ``-R`` forces rerun completed Boltzmann jobs, ``-S`` skips completed jobs.
       | By default, ``-S`` is enabled.

EXAMPLES
========

**Boltzmann-weighted averaging for two conformers**

.. code-block:: console

   chemsmart run thermochemistry -T 298.15 -f conformer1.log -f conformer2.log boltzmann -w electronic

Perform thermochemical analysis for ``conformer1.log`` and ``conformer2.log``, followed by Boltzmann-weighted averaging using electronic energies.

Results will be saved to ``thermochemistry_job_boltzmann.dat``.

**Boltzmann-weighted averaging for batch directory calculation**

.. code-block:: console

   chemsmart run thermochemistry -T 298.15 -d </path/to/directory> -t out boltzmann

Perform thermochemical calculations for all ``.out`` files in the directory and apply Boltzmann-weighted average using Gibbs free energies.

Results will be saved to ``thermochemistry_job_boltzmann.dat``.
