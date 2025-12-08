Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to
contact our team if you have questions or feedback.

#######################################
 General CLI Options for All ORCA Jobs
#######################################

ChemSmart provides comprehensive command-line options for ORCA quantum chemistry calculations. Use ``chemsmart sub orca
--help`` for detailed information about all available options.

*************************
 Basic Command Structure
*************************

The basic command structure for ORCA jobs is:

.. code:: console

   chemsmart sub [OPTIONS] orca [ORCA_OPTIONS] <SUBCMD> [SUBCMD_OPTIONS]

.. tip::

   Apart from some differences in certain parameters, the logic of using ORCA is largely similar to that of Gaussian.

**************
 ORCA_OPTIONS
**************

Works for all ORCA jobs

.. list-table:: Project and File Management Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-p, --project <string>``
      -  string
      -  Project settings (required)

   -  -  ``-f, --filename <string>``
      -  string
      -  Filename from which new ORCA input is prepared (required)

   -  -  ``-l, --label <string>``
      -  string
      -  Write user input filename for the job (without extension). Will overwrite your original filename (default=None)

   -  -  ``-a, --append-label <string>``
      -  string
      -  Name to be appended to file for the job. Will append name to current filename (default=None)

   -  -  ``-t, --title <string>``
      -  string
      -  ORCA job title (default=None)

   -  -  ``-i, --index <string>``
      -  string
      -  Index of molecule to use; default to the last molecule structure (default=None)

   -  -  ``-P, --pubchem <string>``
      -  string
      -  Get molecule structure from PubChem database using identifier (default=None)

.. note::

   |  ``-p`` is followed by the one of the projects specified in ``~/.chemsmart/project/*.yaml`` files, without
      ``.yaml`` extension.
   |  ``-f`` is followed by an input file the user wishes to run job on. Note that this input file can be any format,
      such as ``.xyz``, Gaussian ``.com``, ``.gjf`` or ``.log`` file or ORCA ``.inp`` or ``.out`` file.

**Specify Name of the File**

-  Users can specify the name of the file to be created for the job, without file extension, they want to run by using
   the option ``-l``, e.g.:

   .. code:: console

      chemsmart sub -s shared gaussian -p test -f test.com -l custom_job_name opt

   will create input file named ``custom_job_name.com`` instead of the default ``test_opt.com``.

**Append String to Input File Name**

-  Users can also simply append a string to the base name of the filename supplied, e.g.:

   .. code:: console

      chemsmart sub -s shared gaussian -p test -f test.com -a append_string ts

   will create input file named ``test_append_string.com`` instead of the default ``test_ts.com``.

**Select the Particular Structure in file**

-  If one has more than one structure in the supplied file for input preparation, one can select the particular
   structure to perform job on by using the ``-i/--index`` option, e.g.:

   .. code:: console

      chemsmart sub -s shared gaussian -p test -f small.db -i 5 -c 0 -m 1 opt

   will take the 5th structure (1-indexed, as in chemsmart) from ase database file, ``small.db``, to create the input
   file for geometry optimization.

.. Warning::

   1-indexed numbers are used, instead of 0-indexed numbers in Python language, since most visualization softwares for
   molecules are are 1-indexed.

**Get molecule from pubchem directly**

-  Users can also get the molecule structure directly from PubChem using name, smiles, cid and conformer information,
   e.g.:

   .. code:: console

      chemsmart sub -s shared gaussian -p test -P 962 -c 0 -m 1 -l water opt

   will create input file named ``water.com`` for optimization calculation of water.

.. list-table:: Molecular Properties Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-c, --charge <int>``
      -  int
      -  Charge of the molecule (default=None)

   -  -  ``-m, --multiplicity <int>``
      -  int
      -  Multiplicity of the molecule (default=None)

.. list-table:: Method and Basis Set Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-A, --ab-initio <string>``
      -  string
      -  Ab initio method to be used (default=None)

   -  -  ``-x, --functional <string>``
      -  string
      -  New functional to run (default=None)

   -  -  ``-D, --dispersion <string>``
      -  string
      -  Dispersion for DFT functional (default=None)

   -  -  ``-b, --basis <string>``
      -  string
      -  New basis set to run (default=None)

   -  -  ``-a, --aux-basis <string>``
      -  string
      -  Auxiliary basis set (default=None)

   -  -  ``-e, --extrapolation-basis <string>``
      -  string
      -  Extrapolation basis set (default=None)

.. list-table:: SCF and Grid Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-d, --defgrid``
      -  choice
      -  Grid for numerical integration. Options: defgrid1, defgrid2, defgrid3 (default="defgrid2")

   -  -  ``--scf-tol <choice>``
      -  choice
      -  SCF convergence tolerance. Options: NormalSCF, LooseSCF, SloppySCF, StrongSCF, TightSCF, VeryTightSCF,
         ExtremeSCF (default=None)

   -  -  ``--scf-algorithm <choice>``
      -  choice
      -  SCF algorithm to use. Options: GDIIS, DIIS, SOSCF, AutoTRAH (default=None)

   -  -  ``--scf-maxiter <int>``
      -  int
      -  Maximum number of SCF iterations (default=None)

   -  -  ``--scf-convergence <float>``
      -  float
      -  SCF convergence criterion (default=None)

.. list-table:: Property Calculation Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``--dipole/--no-dipole``
      -  bool
      -  Dipole moment calculation (default=None)

   -  -  ``--quadrupole/--no-quadrupole``
      -  bool
      -  Quadrupole moment calculation (default=None)

   -  -  ``--forces/--no-forces``
      -  bool
      -  Forces calculation (default=False)

.. list-table:: MDCI Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``--mdci-cutoff``
      -  choice
      -  MDCI cutoff. Options: loose, normal, tight (default=None)

   -  -  ``--mdci-density``
      -  choice
      -  MDCI density. Options: none, unrelaxed, relaxed (default=None)

.. list-table:: Additional Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-r, --additional-route-parameters <string>``
      -  string
      -  Additional route parameters (default=None)

********************************
 SUBCMD for Different ORCA Jobs
********************************

.. list-table:: Structure Optimization and Singlet Point Jobs
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``opt``
      -  CLI for optimization calculation for ORCA
   -  -  ``sp``
      -  CLI for single point calculation for ORCA

.. list-table:: Transition State Search
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``ts``
      -  CLI for transition state calculation for ORCA
   -  -  ``modred``
      -  CLI for running ORCA modred jobs
   -  -  ``irc``
      -  CLI for running ORCA IRC jobs
   -  -  ``scan``
      -  CLI for running ORCA scan jobs

.. list-table:: Multiscale QM/MM Calculations
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``qmmm``
      -  CLI for multiscale QM/MM calculations using ORCA

.. list-table:: Direct Input File Execution
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``inp``
      -  Run an ORCA input job as it is

************
 Next Steps
************

For specific calculation types, see the detailed tutorials:

-  Submit Structure Optimization and Singlet Point Jobs
-  Submit Transition State Search Jobs
-  Submit QM/MM Multiscale Jobs
-  Submit Direct Input File ORCA Jobs
