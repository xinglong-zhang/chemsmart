Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to
contact our team if you have questions or feedback.

###################################
 Run Thermochemistry Analysis Jobs
###################################

ChemSmart provides comprehensive thermochemistry analysis capabilities for calculating thermodynamic properties,
Boltzmann weighting, and statistical mechanics data from quantum chemistry calculations.

*************************
 Basic Command Structure
*************************

The basic command structure for thermochemistry jobs is:

.. code:: console

   chemsmart run [OPTIONS] thermochemistry [THERMO_OPTIONS]

.. code:: console

   chemsmart run thermochemistry [THERMO_OPTIONS] boltzmann [SUB_OPTIONS]”

*******************************
 Thermochemistry Analysis Jobs
*******************************

Run comprehensive thermochemistry analysis on frequency calculation outputs.

THERMO_OPTIONS
==============

Works for all thermochemistry jobs

.. list-table:: Thermochemistry Job Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-o, --outputfile``
      -  string
      -  Output file to save the thermochemistry results. Defaults to file_basename.dat (default=None)

   -  -  ``-O, --overwrite``
      -  bool
      -  Overwrite existing output files if they already exist (default=False)

   -  -  ``-i, --check-imaginary-frequencies``
      -  bool
      -  Check for imaginary frequencies in the calculations (default=True)

.. list-table:: Temperature and Pressure Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-t, --temperature``
      -  float
      -  Temperature in Kelvin (default=None)

   -  -  ``-p, --pressure``
      -  float
      -  Pressure in atm (default=None)

   -  -  ``-c, --concentration``
      -  float
      -  Concentration in mol/L (default=None)

.. list-table:: Frequency Cutoff and Correction Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-cs, --cutoff-entropy``
      -  float
      -  Cutoff frequency for entropy in wavenumbers (default=None)

   -  -  ``-ch, --cutoff-enthalpy``
      -  float
      -  Cutoff frequency for enthalpy in wavenumbers (default=None)

   -  -  ``-a, --alpha``
      -  int
      -  Interpolator exponent used in the quasi-RRHO approximation (default=4)

.. list-table:: Mass and Energy Unit Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-w, --weighted``
      -  bool
      -  Use natural abundance weighted masses (True) or use most abundant masses (False) (default=False)

   -  -  ``-u, --energy-units``
      -  string
      -  Units of energetic values. Options: hartree, eV, kcal/mol, kJ/mol (default="hartree")

.. list-table:: Output and Validation Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-o, --outputfile``
      -  string
      -  Output file to save the thermochemistry results. Defaults to file_basename.dat (default=None)

   -  -  ``-O, --overwrite``
      -  bool
      -  Overwrite existing output files if they already exist (default=False)

   -  -  ``-i, --check-imaginary-frequencies``
      -  bool
      -  Check for imaginary frequencies in the calculations (default=True)

Basic Usage
===========

**Basic thermochemistry analysis**

   .. code:: console

      chemsmart run thermochemistry -p thermo_analysis -f freq_output.log analysis

**Thermochemistry with custom output file**

   .. code:: console

      chemsmart run thermochemistry -p thermo_custom -f molecule_freq.log analysis -o custom_results.dat

**Thermochemistry with overwrite option**

   .. code:: console

      chemsmart run thermochemistry -p thermo_overwrite -f calculation.log analysis -O

Examples
========

***********************************
 Boltzmann Weighted Averaging Jobs
***********************************

Run Boltzmann weighted averaging for thermochemistry jobs with multiple conformers.

.. list-table:: Boltzmann Weighting Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-w, --energy-type-for-weighting``
      -  string
      -  Type of energy to use for Boltzmann weighting. Options: gibbs, electronic (default="gibbs")

Boltzmann Basic Usage
=====================

Boltzmann Examples
==================
