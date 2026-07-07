##########################
 Thermochemistry Analysis
##########################

CHEMSMART provides thermochemistry analysis capabilities for computing thermodynamic properties from Gaussian and ORCA
output files and xTB calculation folders.

*******************************
 Thermochemistry Analysis Jobs
*******************************

The ``thermochemistry`` command parses output files and calculates thermochemical properties including enthalpy,
entropy, and Gibbs free energy.

When processing a directory, use ``-p/--program`` when chemsmart needs to identify Gaussian versus ORCA output from the
file content, and use ``-t/--filetype`` when you only want to filter by filename suffix such as ``.log`` or ``.out``.
For xTB, use ``-d`` with ``-p xtb`` (calculation directories, not individual files).

Usage
=====

.. code:: text

   chemsmart run thermochemistry [-d path/to/directory] [-p gaussian|orca|xtb] [-t filetype] [-f filename(s)]
                                 [-csg s_freq_cutoff] [-cst s_freq_cutoff]
                                 [-ch h_freq_cutoff] [-c concentration] [-P pressure] [--weighted | --no-weighted]
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
      -  Directory for batch processing (mutually exclusive with ``-f``)

   -  -  ``-p, --program``
      -  string
      -  Program that produced the output files: ``gaussian``, ``orca``, and ``xtb``. Use this when parsing depends on
         program identity, for example to distinguish Gaussian and ORCA outputs that may share extensions.

   -  -  ``-t, --filetype``
      -  string
      -  File extension to filter when using ``-d`` (e.g. ``log``, ``out``). Use this when selection is purely
         suffix-based, regardless of which program generated the files.

   -  -  ``-f, --filenames``
      -  string
      -  Specific Gaussian or ORCA file(s) to analyze (repeatable, mutually exclusive with ``-d``)

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

   ``-csg`` and ``-cst`` are mutually exclusive. If neither is specified, no quasi-RRHO entropy correction is applied.

**Thermodynamic Conditions:**

.. list-table::
   :header-rows: 1
   :widths: 15 10 75

   -  -  Option
      -  Type
      -  Description

   -  -  ``-c, --concentration``
      -  float
      -  Solution concentration in mol/L (mutually exclusive with -P)

   -  -  ``-P, --pressure``
      -  float
      -  Gas-phase pressure in atm (default: 1.0)

   -  -  ``-w, --weighted, --no-weighted``
      -  bool
      -  Toggle between natural abundance weighted masses (``--weighted``) and most abundant isotope masses
         (``--no-weighted``). Default: ``--weighted``.

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

Analyze a single Gaussian/ORCA output file.

.. code:: bash

   chemsmart run thermochemistry -T 298.15 -f water_opt.out

Output:

.. code:: text

   Structure                        E        ZPE          H       T.S        G(T)
   ================================================================================
   water_opt               -76.323311   0.021581  -76.297951   0.021430  -76.319381

**Single xTB calculation directory:**

Analyze a single xTB calculation directory containing all relevant xTB output files.

.. code:: bash

   chemsmart run thermochemistry -T 298.15 -d water_ohess/ -p xtb

**Multiple files:**

Analyze multiple Gaussian/ORCA output files in a single run.

.. code:: bash

   chemsmart run thermochemistry -T 298.15 -f he_gaussian.log -f he_orca.out

**Batch processing of Gaussian/ORCA output files:**

Automatically discover and analyze Gaussian or ORCA output files within a directory.

.. code:: bash

   chemsmart run thermochemistry -T 298.15 -d . -p gaussian -o thermo.dat

Select files by extension only:

.. code:: bash

   chemsmart run thermochemistry -T 298.15 -d . -t log -o thermo_logs.dat

**Batch processing of xTB calculation directories:**

Analyze multiple xTB calculation directories located within a parent directory.

Example directory layout:

.. code::

   molecules/
   ├── molecule1_ohess/
   ├── molecule2_hess/
   ├── molecule3_opt/
   └── molecule4_sp/

.. code:: bash

   chemsmart run thermochemistry -T 298.15 -d molecules/ -p xtb -o thermo_xtb.dat

The specified directory is scanned for xTB calculation folders. If the given directory itself is not an xTB calculation
directory, only its immediate subdirectories (one level deep) are inspected. Thermochemistry analysis is performed for
each identified calculation directory.

**Batch processing with custom pressure:**

.. code:: bash

   chemsmart run thermochemistry -T 298.15 -P 2.5 -d . -p gaussian -o thermo.dat

.. note::

   Both ``-p`` (program) and ``-P`` (pressure) are case-sensitive options that can be used together. ``-p`` specifies
   the quantum chemistry program that generated the files, while ``-P`` specifies the pressure in atm.

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

*******************************************
 Linear and Quasi-Linear Molecule Handling
*******************************************

CHEMSMART correctly handles linear and quasi-linear molecules when computing rotational thermochemistry. The
``Thermochemistry`` Python class exposes a ``rotational_mode`` parameter that controls how rotational constants are
treated.

Why This Matters
================

Gaussian may print an overflow token (``*************``) for the axial rotational constant of a linear or near-linear
molecule because the field width in the output is too narrow. CHEMSMART converts such tokens to ``numpy.inf`` during
parsing and can then recover physically meaningful rotational constants from the molecular geometry.

Additionally, a molecule whose bond angle is very close to 180° (e.g. a quasi-linear triatomic with a 179.99° angle) may
still print three finite rotational constants in the Gaussian output; the axial constant will simply be many orders of
magnitude larger than the perpendicular constant. CHEMSMART detects this quasi-linear case and handles it correctly in
``"physical"`` mode.

When a quasi-linear species is treated as a linear rotor, its vibrational degrees of freedom should be 3N-5, but
Gaussian and ORCA typically report only 3N-6 frequencies (e.g. three modes for a triatomic). In ``"physical"`` mode,
CHEMSMART pads ``cleaned_frequencies`` by duplicating the lowest positive frequency to restore the missing degenerate
bending mode.

Rotational Modes
================

Two modes are available via the ``rotational_mode`` argument of ``Thermochemistry``:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   -  -  Mode
      -  Description

   -  -  ``"physical"`` *(default)*

      -  Derives rotational constants from the molecular geometry. Collapses effectively linear or quasi-linear ``[A, B,
         C]`` triples to the single perpendicular constant ``[B_perp]``, treating the molecule as a linear rotor.
         Nonlinear molecules are kept unchanged. For quasi-linear species, also pads vibrational frequencies from 3N-6
         to 3N-5 when needed.

   -  -  ``"gaussian"``
      -  Preserves the rotational constants exactly as printed in the Gaussian output. If Gaussian overflowed one or
         more values (``*************``), CHEMSMART falls back to geometry-derived constants automatically.

Python API Examples
===================

**Physical mode (default) — linear rotor treatment:**

.. code:: python

   from chemsmart.analysis.thermochemistry import Thermochemistry

   thermo = Thermochemistry(
       filename="koh_opt.log",
       temperature=298.15,
       rotational_mode="physical",   # default
   )
   print(thermo.rotational_partition_function)

The molecule is analysed geometrically; any near-linear or truly linear structure is automatically collapsed to a
linear-rotor model.

**Gaussian mode — reproduce Gaussian output file values:**

.. code:: python

   thermo = Thermochemistry(
       filename="koh_opt.log",
       temperature=298.15,
       rotational_mode="gaussian",
   )
   print(thermo.rotational_partition_function)

This mode matches the partition function and thermochemical corrections that Gaussian itself would report. When Gaussian
printed overflow tokens for the axial constant, CHEMSMART silently falls back to geometry-derived constants to avoid a
division-by-zero error while keeping the nonlinear-rotor model.

.. note::

   ``rotational_mode`` is a Python-API parameter of ``Thermochemistry`` and is not yet exposed as a CLI option. The
   ``chemsmart run thermochemistry`` command always uses the ``"physical"`` mode.

***********************************
 Boltzmann Weighted Averaging Jobs
***********************************

The ``boltzmann`` subcommand performs Boltzmann-weighted averaging of thermochemical results across multiple conformers.

Usage
=====

.. code:: text

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

Output:

.. code:: text

   Structure                                        E        ZPE          H       T.S        G(T)
   ================================================================================================
   conformer_boltzmann_avg_by_electronic   -2189.631938   0.288732  -2189.312517   0.097016  -2189.409533

**Batch directory averaging:**

.. code:: bash

   chemsmart run thermochemistry -T 298.15 -d . -p gaussian boltzmann

Results are saved to ``thermochemistry_job_boltzmann.dat``.
